function my_update_documentation()
% Update documentation for Maxwell.
    
    maxwell_dir = strrep(mfilename('fullpath'), ['private/', mfilename], '');
    files = dir(maxwell_dir);

    % Find example_*.m and maxwell_*.m files
    maxwell_files = {};
    example_files = {};
    for k = 1 : length(files)
        if strfind(files(k).name, 'maxwell_') == 1 & ...
            files(k).name(end-1:end) == '.m'
            maxwell_files{end+1} = files(k).name;
        elseif strfind(files(k).name, 'example') == 1 & ...
            files(k).name(end-1:end) == '.m'
            example_files{end+1} = files(k).name;
        end
    end

    % Generate maxwell_help file.
    try 
        delete([maxwell_dir, 'maxwell_help.m']);
    end

    return

    % Publish documentation.
    doc_files = [maxwell_files, example_files, {'maxwell_help.m'}];
    for k = 1 : length(doc_files)
        publish([maxwell_dir, doc_files{k}], 'evalCode', false);
    end


