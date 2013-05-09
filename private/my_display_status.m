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
