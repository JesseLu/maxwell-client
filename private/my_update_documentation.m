function my_update_documentation()
% Update documentation for Maxwell.
    
    maxwell_dir = strrep(mfilename('fullpath'), ['private/', mfilename], '');
    files = dir(maxwell_dir);

        %
        % Find example_*.m and maxwell_*.m files
        %

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
    doc_files = [maxwell_files, example_files];


        %
        % Generate maxwell_help file.
        %

    % Get by-line (second line) from files.
    for k = 1 : length(doc_files)
        f = fopen([maxwell_dir, doc_files{k}], 'r');
        title = fgetl(f); % Usually not needed.
        byline{k} = fgetl(f);
        try
            byline{k}(1:2) = []; % Delete first 2 characters (should be '% ').
        end 
        fclose(f);
    end

    % Get template.
    f = fopen('my_help_template.m', 'r');
    my_help = char(fread(f))';
    fclose(f);

    % Make one-liners
    for k = 1 : length(byline)
        oneliners{k} = sprintf('%%\n%% * <%s %s> - %s\n', ...
                                strrep(doc_files{k}, '.m', '.html'), ...
                                strrep(doc_files{k}, '.m', ''), ...
                                byline{k});
    end

    my_help = strrep(my_help, '<<<functions>>>', ...
                        [oneliners{1:length(maxwell_files)}]);
    my_help = strrep(my_help, '<<<examples>>>', ...
                        [oneliners{length(maxwell_files)+1:end}]);

    f = fopen([maxwell_dir, 'maxwell_help.m'], 'w');
    fwrite(f, my_help);
    fclose(f);


        %
        % Publish documentation.
        %

    doc_files = [doc_files, {'maxwell_help.m'}];
    for k = 1 : length(doc_files)
        publish([maxwell_dir, doc_files{k}], 'evalCode', false);
    end


