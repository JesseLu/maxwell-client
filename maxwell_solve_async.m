%% maxwell_solve_async
% Initiate an electromagnetic solve and immediately return a callback function.

%%% Syntax
% * |cb_fun = maxwell_solve_async(grid, eps, J)|
%   uploads a simulation and returns a callback function |cb_fun|
%   which is used to query progress and return the results.
%   The callback function's syntax is |[done, E, H] = cb_fun()|
%   where |done| is a boolean variable that is set to true
%   only when |E| and |H| are the final solution fields.
%
% * |cb_fun = maxwell_solve_async(grid, [eps mu], J)|
%   is the same as above except that it allows |mu ~= 1|.
%
% The optional input parameters for |maxwell_solve| are also valid here.

%%% Description
% |maxwell_solve_async| is an asynchronous version of |maxwell_solve|
% in that it exits as soon as the simulation has been uploaded to 
% the server, but before it is complete.
% For this reason, instead of returning the solution fields,
% a callback function is returned which is used to query the
% progress of the simulation and to obtain the solution fields.
% The callback function can be used in this simple way:
%
%   while ~cb_fun(); end; % Wait until simulation finishes.
%   [~, E, H] = cb_fun(); % Get solution fields.
%
% The benefit of such a function is that it frees the local Matlab process
% to perform additional work.
%

function [cb, vis_progress] = maxwell_solve_async(grid, eps_mu, J, varargin)

        %
        % Get current axis (for plotting) and start time (for timing).
        %

    my_axis = gca;
    start_time = tic;


        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1);
    end
    my_validate_field(eps, grid.shape, 'eps', mfilename);
    my_validate_field(mu, grid.shape, 'mu', mfilename);

    my_validate_field(J, grid.shape, 'J', mfilename);

    % Optional parameter-value pairs.
    options = my_parse_options(struct(  'E0', {my_default_field(grid.shape, 0)}, ...
                                        'vis_progress', 'plot', ...
                                        'max_iters', 1e5, ...
                                        'err_thresh', 1e-6), ...
                                varargin);
    my_validate_field(options.E0, grid.shape, 'E0', mfilename);
    validateattributes(options.vis_progress, {'char'}, ...
                {}, mfilename, 'vis_progress');
    validateattributes(options.max_iters, {'numeric'}, ...
        {'integer', 'positive', 'scalar', 'real'}, mfilename, 'max_iters');
    validateattributes(options.err_thresh, {'double'}, ...
        {'positive', 'scalar', '<', 1}, mfilename, 'err_thresh');


        %
        % Check if the simulation is 2D (can be solved locally).
        %
   
    if any(grid.shape == 1)
        % Compute E-field.
        [A, ~, b] = maxwell_axb(grid, [eps mu], options.E0, J);
        x = A \ b;
        N = prod(grid.shape);
        for k = 1 : 3
            E{k} = reshape(x((k-1)*N+1:k*N), grid.shape);
        end
        err = norm(A*x-b)/norm(b);
        
        % Compute H-field.
        H = my_E2H(grid, mu, E);

        % Return solution.
        vis_progress = 'none';
        cb = @() my_simple_callback(true, E, H, err);
        return
    end


        %
        % Upload simulation.
        %

    [server_url, prefix, vis_progress] = maxwell_upload(grid, eps, mu, J, ...
                                    options.E0, options.max_iters, ...
                                    options.err_thresh, options.vis_progress);


        %
        % Set up callback function.
        %

    % Persistent variables for the callback function.
    p_is_done = false;
    p_E = [];
    p_H = [];
    p_err = [];
    p_state = [];
    
    function [is_done, E, H, err] = maxwell_callback()
    % Queries server to inform user of the state of the simulation.
        if ~p_is_done % Not done, keep trying.
            [p_E, p_err, p_state, s] = maxwell_download(server_url, prefix);
            if strcmp(p_state, 'finished')
                p_is_done = true;
                p_H = my_E2H(grid, mu, p_E);
            end
        end

        % Update static variables.
        is_done = p_is_done;
        E = p_E;
        H = p_H;
        err = p_err;
        state = p_state;

        % Show the progress.
        if isempty(err) % Simulation not yet started.
            progress_text = sprintf('[%s] err: ----, iter: 0, seconds: %1.1f', ...
                                    state, toc(start_time));
        else 
            if is_done % Simulation complete.
                progress_text = sprintf('[%s] err: %e, iter: %d', ...
                                        state, err(end), length(err));
            else % Simulation in progress.
                progress_text = sprintf('[%s] err: %e, iter: %d, seconds: %1.1f', ...
                                        state, err(end), length(err), toc(start_time));
            end
        end

        if strcmp(vis_progress, 'text') | strcmp(vis_progress, 'both')
            % Normalized text progress output prints constant length of 60.
            norm_p_text = [progress_text, ...
                    repmat(' ', 1, 60 - length(progress_text))];
            fprintf(norm_p_text);
        end

        if strcmp(vis_progress, 'plot') | strcmp(vis_progress, 'both')
            % Plot the progress in log-scale.
            axes(my_axis);
            if isempty(err)
                semilogy(1, 'bx');
            else 
                semilogy(err, 'b.-');
            end
            title(progress_text);
            xlabel('iterations');
            ylabel('error');

            % Add a dotted line showing the error threshold.
            hold on
            a = axis;
            semilogy(a(1:2), options.err_thresh * [1 1], 'k--');
            axis([a(1:2), options.err_thresh/10, a(4)]);
            hold off
        end

    end

    cb = @maxwell_callback;
end


function [varargout] = my_simple_callback(varargin)
    varargout = varargin;
end

