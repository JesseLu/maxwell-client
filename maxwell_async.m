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
    
function [callback] = maxwell(grid, epsilon, J, varargin)
 
    % Parse input and option parameters.
    [omega, s_prim, s_dual, mu, epsilon, E0, J, max_iters, err_thresh] = ...
        my_parse_inputs(grid, epsilon, J, varargin{:});

    % Generate a random (and hopefully unique) ID.
    id = [datestr(now, 'HHMMSSFFF'), '-', num2str(round(1e6*rand(1)))];

    % Choose a prefix for the filename. 
    prefix = ['maxwell-', id, '.'];
    if ~isempty(tempdir) % Means temporary directory is accessible.
        prefix = [tempdir, prefix];
    end
    prefix



        %
        % Create input file.
        %

    %% Make sure hdf5 compression is available.
    if ~H5Z.filter_avail('H5Z_FILTER_DEFLATE') || ...
        ~H5ML.get_constant_value('H5Z_FILTER_CONFIG_ENCODE_ENABLED') || ...
        ~H5ML.get_constant_value('H5Z_FILTER_CONFIG_DECODE_ENABLED') || ...
        ~H5Z.get_filter_info('H5Z_FILTER_DEFLATE')
        error('HDF5 gzip filter not available!') 
    end
    
    xyz = 'xyz'
    gridfile = [prefix, 'grid'];
    hdf5write(gridfile, 'username', 'user', 'WriteMode', 'overwrite');
    hdf5write(gridfile, 'password', 'pwd', 'WriteMode', 'append');
    hdf5write(gridfile, 'omega_r', double(real(omega)), 'WriteMode', 'append');
    hdf5write(gridfile, 'omega_i', double(imag(omega)), 'WriteMode', 'append');
    hdf5write(gridfile, 'max_iters', int64(max_iters), 'WriteMode', 'append');
    hdf5write(gridfile, 'err_thresh', double(err_thresh), 'WriteMode', 'append');

