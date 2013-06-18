%% Examples
%
%   cb = maxwell(grid, epsilon, J); % Simple example.
%   
%   cb = maxwell(grid, epsilon, J, 'option', val); % Option pairs allowed.
%
% Available option pairs
% * mu
% * E0
% * max_iters
% * err_thresh
% * vis_progress
%
    
function [server_url, name] = maxwell_upload(grid, epsilon, J, varargin)
    server_url = 'http://raven1.stanford.edu:8008/';
 
    % Parse input and option parameters.
    [omega, s_prim, s_dual, mu, epsilon, E0, J, max_iters, err_thresh] = ...
        my_parse_inputs(grid, epsilon, J, varargin{:});

    % Generate a random (and hopefully unique) ID.
    id = [datestr(now, 'HHMMSSFFF'), '-', num2str(randi([1e6 1e7-1]))];

    % Choose a prefix for the filename. 
    name = ['maxwell-', id, '.'];
    prefix = [tempdir, name];
    url = [server_url, name];

    % Write the grid file.
    gridfile = [prefix, 'grid']; 
    hdf5write(gridfile, 'omega_r', real(omega), 'WriteMode', 'overwrite');
    hdf5write(gridfile, 'omega_i', imag(omega), 'WriteMode', 'append');
    xyz = 'xyz';
    for k = 1 : length(xyz)
        hdf5write(gridfile, ['sp_', xyz(k), 'r'], ...
                    real(s_prim{k}), 'WriteMode', 'append');
        hdf5write(gridfile, ['sp_', xyz(k), 'i'], ...
                    imag(s_prim{k}), 'WriteMode', 'append');
        hdf5write(gridfile, ['sd_', xyz(k), 'r'], ...
                    real(s_dual{k}), 'WriteMode', 'append');
        hdf5write(gridfile, ['sd_', xyz(k), 'i'], ...
                    imag(s_dual{k}), 'WriteMode', 'append');
    end
    hdf5write(gridfile, 'max_iters', int64(max_iters), 'WriteMode', 'append');
    hdf5write(gridfile, 'err_thresh', double(err_thresh), 'WriteMode', 'append');

    % Write the other files (if needed).
    my_write(prefix, 'e', epsilon);
    my_write(prefix, 'J', J);
    my_write(prefix, 'm', mu);
    my_write(prefix, 'E', E0);

    
    % Upload files.
    files = dir([prefix, '*']);
    my_upload({files(:).name}, tempdir, server_url);

    % Upload a (empty) request file.
    request_file = fopen([prefix, 'request'], 'w');
    fprintf(request_file, 'all files uploaded');
    my_upload({strrep([prefix, 'request'], '/tmp/', '')}, tempdir, server_url);
    
    % Delete files.
    for f = 'eJmE'
        for g = 'xyz'
            for h = 'ri'
                delete([prefix, f, '_', g, h])
            end
        end
    end
    delete([prefix, 'grid'])
    delete([prefix, 'request'])
end

function my_write(prefix, name, field)
    xyz = 'xyz';
    for k = 1 : 3
        file_prefix = [prefix, name, '_', xyz(k)];
        my_write_data([file_prefix, 'r'], 'data', real(field{k}));
        my_write_data([file_prefix, 'i'], 'data', imag(field{k}));
    end
end
 
function my_write_data(filename, dset_name, data)
    if ndims(data) ~= 3 
        return % No need to write, default values still present.
    end
    N = 3; % 3D arrays assumed.
    data = permute((data), [ndims(data):-1:1]); % Make data row-major.
    dims = fliplr(size(data)); % Size of the array.
    chunk_dims = [1 dims(2:3)]; % Should heavily affect compression.

    % Create file.
    file = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

    % Create the dataspace.
    space = H5S.create_simple(N, dims, []);

    % Set dataspace properties.
    dcpl = H5P.create('H5P_DATASET_CREATE');
    H5P.set_deflate(dcpl, 1); % Deflation level: 0 (none) to 9 (most).
    H5P.set_chunk(dcpl, chunk_dims);

    % Create dataset and write to file.
    dset = H5D.create(file, dset_name, 'H5T_IEEE_F32BE', space, dcpl); 
    H5D.write(dset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', single(real(data)));

    % Close resources.
    H5P.close(dcpl);
    H5D.close(dset);
    H5S.close(space);
    H5F.close(file) % Close file, flushing to storage.
end

