function s3_download(filenames, local_dir)
    my_disp = @(s) my_display_status(s, 'text');
    for k = 1 : length(filenames)
        [inputStream{k}, file{k}] = my_open_connection(filenames{k}, local_dir);
    end

    my_disp('Receiving...');

    my_stream_send(inputStream, file, 'received', my_disp);

    for k = 1 : length(filenames) % Close the files and connections.
        inputStream{k}.close()
        file{k}.close()
    end

end


function [inputStream, file] = my_open_connection(filename, local_dir)
    url = ['https://s3.amazonaws.com/maxwell-out/', filename];
    urlConnection = my_urlreadwrite(url);
    urlConnection.connect();
    inputStream = urlConnection.getInputStream();
    file = java.io.FileOutputStream([local_dir, filesep, filename]);
end
                 

