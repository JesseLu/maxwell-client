% Parse inputs for the maxwell_upload function.
function [omega, s_prim, s_dual, mu, epsilon, E0, J, max_iters, err_thresh, vis_progress] = ...
        parse_inputs(grid, epsilon, J, varargin);

        %
        % Process option pairs.
        %

    % Here are the default options
    maxwell_options = struct(   'mu', {{1, 1, 1}}, ...
                                'E0', {{0, 0, 0}}, ...
                                'max_iters', 5e4, ...
                                'err_thresh', 1e-6', ...
                                'vis_progress', 'plot'); % Options: none, text, plot, both.

    for k = 2 : 2 : length(varargin)
        maxwell_options = setfield(maxwell_options, varargin{k-1}, varargin{k});
    end
      
    mu = maxwell_options.mu;
    E0 = maxwell_options.E0;
    max_iters = maxwell_options.max_iters;
    err_thresh = maxwell_options.err_thresh;
    vis_progress = maxwell_options.vis_progress;
    
    omega = grid.omega;
    s_prim = grid.s_prim;
    s_dual = grid.s_dual;


        %
        % Verify inputs.
        %

    % Check omega.
    if numel(omega) ~= 1 
        error('OMEGA must be a scalar.');
    end

    
    shape = size(epsilon{1}); % Make sure all 3D fields have this shape.
    % Expand mu and E0 if needed.
    for k = 1 : 3
        if numel(mu{k}) == 1
            mu{k} = mu{k} * ones(shape);
        end
        if numel(E0{k}) == 1
            E0{k} = E0{k} * ones(shape);
        end
    end

    % Check shapes of mu, epsilon, J.
    % Specifically, make sure all component fields have shape xx-yy-zz.
    if any([numel(mu), numel(epsilon), numel(E0), numel(J)] ~= 3)
        error('All 3D vector fields (EPSILON, J, MU, E0) must have three cell elements.');
    end

    fields = [mu, epsilon, E0, J];
    for k = 1 : length(fields)
       if any(size(fields{k}) ~= shape)
            error('All 3D vector fields (EPSILON, J, MU, E0) must have the same shape.');
        end
    end


    % Check shapes of s_prim and s_dual.
    % Specifically, each array of each must have length xx, yy, and zz respectively.
    if any([numel(s_prim), numel(s_dual)] ~= 3)
        error('S_PRIM, and S_DUAL must each have three cell elements.');
    end
    for k = 1 : 3
        s_prim{k} = s_prim{k}(:);
        s_dual{k} = s_dual{k}(:);
        if (length(s_prim{k}) ~= shape(k)) || (length(s_dual{k}) ~= shape(k))
            error('The lengths of S_PRIM and S_DUAL vectors must be xx, yy, and zz, in that order.')
        end
    end

    % Make sure the value for max_iters is valid.
    if mod(max_iters,1) ~= 0 || max_iters <= 0 || ~isreal(max_iters) || max_iters > 1e8
        error('MAX_ITERS must be a positive integer less than 1e8.');
    end

    % Make sure the value for err_thresh is valid.
    if ~isfloat(err_thresh) || err_thresh <= 0 || err_thresh >= 1 || ~isreal(err_thresh)
        error('ERR_THRESH must be a positive number between 0 and 1.');
    end

    % Make sure the value for view_progress is valid.
    if all(~strcmp(vis_progress, {'plot', 'text', 'none', 'both'}))
        error('VIS_PROGRESS must be either ''plot'', ''text'', or ''none''.');
    end


