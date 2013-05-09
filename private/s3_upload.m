%% s3_upload
% Upload files to Maxwell's simulation input S3 bucket.

function s3_upload(filenames, local_dir)
    my_disp = @(s) my_display_status(s, 'text');
    url = 'http://localhost:8000';

    for k = 1 : length(filenames)
        % Open connections and send headers.
        [infile{k}, outputStream{k}, footer{k}] = ...
            my_open_connection(filenames{k}, local_dir);
    end

    % Stream the files.
    my_stream_send(infile, outputStream, 'sent', my_disp);

    for k = 1 : length(filenames)
        % Send the footers.
        outputStream{k}.write(footer{k}.getBytes(), 0, footer{k}.length); 

        % Close off the messages.
        infile{k}.close(); 
        outputStream{k}.flush();
    end
end

%% Open a connection (POST).
function [infile, outputStream, footer] = my_open_connection(filename, local_dir)
    url = 'https://s3.amazonaws.com/maxwell-in';
    params = {  'key', filename, ...
                'acl', 'public-read', ...
                'content-type', 'application/octet-stream', ...
                'AWSAccessKeyId', 'AKIAI6BMVTBFEAMEMHQQ', ...
                'policy', 'ewogICJleHBpcmF0aW9uIjogIjIwMjAtMDEtMDFUMTI6MDA6MDAuMDAwWiIsCiAgImNvbmRpdGlvbnMiOiBbCiAgICB7ImJ1Y2tldCI6ICJtYXh3ZWxsLWluIiB9LAogICAgeyJhY2wiOiAicHVibGljLXJlYWQiIH0sCiAgICBbInN0YXJ0cy13aXRoIiwgIiRrZXkiLCAiIl0sCiAgICBbInN0YXJ0cy13aXRoIiwgIiRDb250ZW50LVR5cGUiLCAiYXBwbGljYXRpb24vb2N0ZXQtc3RyZWFtIl0sCiAgXQp9Cg==', ...
                'signature', 'QEPdrd7W6aOG8yO/tMOeyWWeNFY='};

    % Create a urlConnection.
    [urlConnection, errorid, errormsg] = my_urlreadwrite(url);
    if isempty(urlConnection)
        error(['Could not connect to url: ', url]);
    end

    urlConnection.setDoOutput(true); % Sets the request mode to POST.
    boundary = '*** maxwell_client boundary ***';
    urlConnection.setRequestProperty('Content-Type',...
        ['multipart/form-data; boundary=', boundary]);

    eol = [char(13),char(10)]; % End-of-line character.

    % Build the header, body and footer of the POST request.
    header = [];
    for k = 1 : 2 : length(params) % Form data for text parameters.
        header = [header, '--', boundary, eol, ...
                    'Content-Disposition: form-data; name="', params{k}, '"', ...
                    eol, eol, params{k+1}, eol];
    end
    header = java.lang.String([header, '--', boundary, eol, ...
                'Content-Disposition: form-data; name="file"; ', ...
                    'filename="', filename, '"', eol, ...
                'Content-Type: application/octet-stream', eol, eol]); 
                % Form data for binary data (the simulation file.
    file = java.io.File([local_dir, filesep, filename]);
    footer = java.lang.String([eol, '--', boundary, '--', eol]);

    % We used a streaming connection, crucial for large files.
    total_length = header.length + file.length() + footer.length;
    urlConnection.setFixedLengthStreamingMode(total_length);
    outputStream = java.io.DataOutputStream(urlConnection.getOutputStream);

    outputStream.write(header.getBytes(), 0, header.length); % Send the header.

    infile = java.io.FileInputStream(file);
end
