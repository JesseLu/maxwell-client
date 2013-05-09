    function S3_post
    my_disp = @(s) my_display_status(s, 'text');
url = 'http://localhost:8000';
url = 'http://s3.amazonaws.com/maxwell-in';
filename = 'small';
params = {'key', 'matlabtest', ...
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
                'Content-Disposition: form-data; name="file"; filename="matlabtest"', eol, ...
                'Content-Type: application/octet-stream', eol, eol]); 
                % Form data for binary data (the simulation file.
    file = java.io.File(filename);
    footer = java.lang.String([eol, '--', boundary, '--', eol]);

    % We used a streaming connection, crucial for large files.
    total_length = header.length + file.length() + footer.length;
    urlConnection.setFixedLengthStreamingMode(total_length);
    outputStream = java.io.DataOutputStream(urlConnection.getOutputStream);

    outputStream.write(header.getBytes(), 0, header.length); % Send the header.

    % fprintf('Sending...'); % Send the file.
    my_disp('Sending...'); % Send the file.
    infile = java.io.FileInputStream(file);
    stream_send({infile}, {outputStream}, 'sent', my_disp);

    outputStream.write(footer.getBytes(), 0, footer.length); % Send the footer.

    infile.close(); % Close off the message.
    outputStream.flush();
    end

function stream_send (in, out, action_name, display_fun)

    copier = MaxwellCopier; % Requires the Maxwell.jar library to be loaded.

    if length(in) ~= length(out)
        error('Unequal number of input and output streams.');
    end


    N = length(in);
    running = true;
    start_time = tic;
    status_time = start_time;
    prevlen = 0;

    while any(running)
        for k = 1 : N % Transfer some data.
            running = copier.copy(in{k}, out{k});
        end

        if toc(status_time) > 0.3 || all(~running) % Periodically give updates.
            megabytes = copier.total_bytes_transferred / 1e6;
            status_line = sprintf('%1.2f MB %s (%1.2f MB/s)', ...
                megabytes, action_name, megabytes/toc(start_time));
            display_fun(status_line);
%             fprintf([repmat('\b', 1, prevlen), status_line]); % Write-over.
%             prevlen = length(status_line);
            status_time = tic;
        end
    end
    % fprintf('\n');
end

function my_display_status(status_text, option)
% Used to display text on command line or title area of a plot.
    max_length = 60;
    switch option
        case 'text'
            if length(status_text) > max_length
                fprintf(status_text(1:max_length));
            else 
                fprintf([status_text, ...
                        repmat(' ', 1, max_length - length(status_text))]);
            end
            fprintf(repmat('\b', 1, max_length));
        case 'plot'
            title(status_text);
            drawnow
        otherwise
            error('Unrecognized option, must be either ''plot'' or ''text''');
    end
end
