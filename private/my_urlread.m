%% maxwell_urlread
% Used to GET or POST to a url.
% Reads out the response in a stream-it-on-one-line fashion.
function [status, response] = ...
    my_urlread(urlChar, method, params, verbosity, varargin)

    % This function requires Java.
    if ~usejava('jvm')
       error(message('MATLAB:urlread:NoJvm'));
    end

    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier;

    % Be sure the proxy settings are set.
    com.mathworks.mlwidgets.html.HTMLPrefs.setProxySettings

%     % Check number of inputs and outputs.
%     error(nargchk(1,3,nargin))
%     error(nargoutchk(0,2,nargout))
    if ~ischar(urlChar)
        error('MATLAB:urlread:InvalidInput','The first input, the URL, must be a character array.');
    end
    if (nargin > 1) && ~strcmpi(method,'get') && ~strcmpi(method,'post')
        error('MATLAB:urlread:InvalidInput','Second argument must be either "get" or "post".');
    end

    % Do we want to throw errors or catch them?
    if nargout == 2
        catchErrors = true;
    else
        catchErrors = false;
    end

    % Set default outputs.
    output = '';
    status = 0;

    % GET method.  Tack param/value to end of URL.
    if (nargin > 1) && strcmpi(method,'get')
        if mod(length(params),2) == 1
            error('MATLAB:urlread:InvalidInput','Invalid parameter/value pair arguments.');
        end
        for i=1:2:length(params)
            if (i == 1), separator = '?'; else separator = '&'; end
            param = char(java.net.URLEncoder.encode(params{i}));
            value = char(java.net.URLEncoder.encode(params{i+1}));
            urlChar = [urlChar separator param '=' value];
        end
    end

    % Create a urlConnection.
    [urlConnection,errorid,errormsg] = my_urlreadwrite(urlChar, varargin{:});
    if isempty(urlConnection)
        if catchErrors, return
        else error(errorid,errormsg);
        end
    end

    % POST method.  Write param/values to server.
    if (nargin > 1) && strcmpi(method,'post')
%        try
            urlConnection.setDoOutput(true);
            urlConnection.setRequestProperty( ...
                'Content-Type','application/x-www-form-urlencoded');
            printStream = java.io.PrintStream(urlConnection.getOutputStream);
            for i=1:2:length(params)
                if (i > 1), printStream.print('&'); end
                param = char(java.net.URLEncoder.encode(params{i}));
                value = char(java.net.URLEncoder.encode(params{i+1}));
                printStream.print([param '=' value]);
            end
            printStream.close;
%         catch
%             if catchErrors, return
%             else error('MATLAB:urlread:ConnectionFailed','Could not POST to URL.');
%             end
%         end
    end

    % Read the data from the connection.
        
    % Receive the stream.
    inputStream = urlConnection.getInputStream();
    reader = java.io.BufferedReader(java.io.InputStreamReader(inputStream));

    response = '';
    if verbosity == 1 % Stream character-by-character.
        c = reader.read();
        while c ~= -1
            fprintf('%s', char(c))
            response = [response, char(c)];
            c = reader.read();
        end

    elseif verbosity == 0 % Stream line-by-line with overwriting.
        prevlen = 0;
        c = reader.readLine();
        while ~(isempty(c) && isnumeric(c)) % Null in Matlab: empty double array.
            warning off; c = sprintf(char(c)); warning on;
            fprintf('%s', sprintf(repmat('\b', 1, prevlen)));
            warning off; fprintf('%s', c); warning on;
            prevlen = length(c) - ...
                        2*length(strfind(c, sprintf('\b')));
            response = [response, c, sprintf('\n')];

            c = reader.readLine();
        end
        fprintf('\n');
    end
    reader.close();
    inputStream.close();
    status = 1;
end

