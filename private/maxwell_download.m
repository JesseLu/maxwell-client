function [E, H, err, state, s] = maxwell_download(server_url, name)

    E = [];
    H = [];
    err = [];

    url = [server_url, name]
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
    
    end

