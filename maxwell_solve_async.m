%% maxwell_solve_async
% Asynchronous solve for matlab
function [cb, vis_progress] = maxwell_solve_async(grid, epsilon, J, varargin)

    my_axis = gca;
    start_time = tic;
    [server_url, prefix, vis_progress] = maxwell_upload(grid, epsilon, J, varargin{:});

    % Persistent variables for the callback function.
    p_is_done = false;
    p_E = [];
    p_H = [];
    p_err = [];
    p_state = [];
    
    function [is_done, E, H, err] = maxwell_callback()
        if ~p_is_done
            [p_E, p_H, p_err, p_state, s] = maxwell_download(server_url, prefix);
            if strcmp(p_state, 'finished')
                p_is_done = true;
            end
        end

        is_done = p_is_done;
        E = p_E;
        H = p_H;
        err = p_err;
        state = p_state;

        % Show the progress.
        if isempty(err)
            progress_text = sprintf('[%s] err: ----, iter: 0, seconds: %1.1f', ...
                                    state, toc(start_time));
        else
            if is_done
                progress_text = sprintf('[%s] err: %e, iter: %d', ...
                                        state, err(end), length(err));
            else
                progress_text = sprintf('[%s] err: %e, iter: %d, seconds: %1.1f', ...
                                        state, err(end), length(err), toc(start_time));
            end
        end

        if strcmp(vis_progress, 'text') | strcmp(vis_progress, 'both')
            % Normalized text prints constant length of 80.
            norm_p_text = [progress_text, ...
                    repmat(' ', 1, 80 - length(progress_text))];
            fprintf(norm_p_text);
        end

        if strcmp(vis_progress, 'plot') | strcmp(vis_progress, 'both')
            axes(my_axis);
            title(progress_text);
            xlabel('iterations');
            ylabel('error');
            if isempty(err)
                semilogy(1);
            else 
                semilogy(err, 'b.-');
                hold on
                semilogy(length(err), err(end), 'rx');
                hold off
            end
        end

    end

    cb = @maxwell_callback;
end


