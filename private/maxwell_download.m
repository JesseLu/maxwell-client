%% maxwell_download
% Downloads result files from simulation server.

function [E, H, err, state, s] = maxwell_download(server_url, name)

    E = [];
    H = [];
    err = [];

    url = [server_url, name];
    endings = {'request', 'status', 'log', 'finished'};
    for k = 1 : length(endings)
        [s{k}, status(k)] = urlread([url, endings{k}]);
    end

    status_as_str = sprintf('%d%d%d%d', status); 
    switch status_as_str
        case '1000'
            state = 'queued';
        case '0010'
            state = 'loading';
        case '0110'
            state = 'executing';
        case '0111'
            state = 'finished';
        otherwise
            state = sprintf('unknown (%s)', status_as_str);
    end

    if strcmp(state, 'executing') | strcmp (state, 'finished')
        % Get the error from the status file.
        err = str2num(s{2});
    end

    if strcmp(state, 'finished')
        % Download all files.
        [a, b, c] = ndgrid('EH', 'xyz', 'ri');
        for k = 1 : numel(a)
            files{k} = [name, a(k), '_', b(k), c(k)];
        end
        my_download(files, tempdir, server_url);

        % Load files.
        for k = 2 : 2 : numel(files)
            F{k/2} = double(h5read([tempdir, files{k-1}], '/data')) + ...
                1i * double(h5read([tempdir, files{k}], '/data'));
            F{k/2} = permute(F{k/2}, [ndims(F{k/2}):-1:1]); % Convert to column-major.
        end
        E = F(1:3);
        H = F(4:6);

        % Delete files.
        for k = 1 : numel(files)
            delete([tempdir, files{k}]);
        end
    
    end

