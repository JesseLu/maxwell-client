%% maxopt_case_2wbeam
% Sets up a frequency-doubling cavity optimization.

function [fun, x0] = maxopt_case_2wbeam(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    
    beam_height = 220;
    beam_width = 440;
    hole_pos = 400 * [0.5:7.5];
    hole_xlen = 120 * ones(size(hole_pos));
    hole_ylen = 200 * ones(size(hole_pos));

    hole_params = [hole_pos; hole_xlen; hole_ylen];
    x0 = [beam_height; beam_width; hole_params(:)];

    wvlen = [1550, 775];
    eps_val = [3.2^2, 3.5^2];


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, ~, E, H, grid, eps] = solve_structure(varargin{:});
    end

    flt = options.flatten;
    switch type
        case 'get_fields'
            fun = @(x) get_fields(wvlen, eps_val, x, flt, false);
        case 'fval'
            fun = @(x) solve_structure(wvlen, eps_val, x, flt, false);
        case 'grad_f'
            fun = @(x) solve_structure(wvlen, eps_val, x, flt, true);
        otherwise
            error('Invalid type.');
    end
end

function [Fval, grad_F, E, H, grid, eps] = ...
                solve_structure(wvlen, eps_val, varargin)
% Simulate all structures.
    for k = 1 : length(wvlen)
        subplot(length(wvlen), 1, k);
        [fval(k), grad_f{k}, E{k}, H{k}, grid{k}, eps{k}] = ...
                    solve_one_structure(wvlen(k), eps_val(k), varargin{:});
    end

    % Find the worst one.
    [Fval, ind] = max(fval);
    grad_F = grad_f{ind};

    % Pretty print.
    fprintf('fvals: ');
    for k = 1 : length(wvlen)
        if k == ind
            fprintf('[%e] ', fval(k));
        else
            fprintf('%e ', fval(k));
        end
    end
    fprintf('\n');
end

function [fval, grad_f, E, H, grid, eps] = ...
                solve_one_structure(wvlen, eps_val, params, flatten, calc_grad)
% Simulate a nanobeam structure.


    [grid, eps, J] = make_structure(wvlen, eps_val, params, flatten);

        %
        % Use central point source as the excitation.
        %

    [x, y, z] = maxwell_pos2ind(grid, 'Ey', [0 0 0]);
    J{2}(x, y, z) = 1;
    J{2}(x-1, y, z) = 1;

        
        %
        % Solve.
        %

    [E, H] = maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); 


        %
        % Measure power reflected back to the center (figure of merit).
        %

    function [fval, grad_E] = fitness(E)
    % Calculates figure of merit and its derivative.
        % Figure of merit.
        E_meas = [E{2}(x, y, z); E{2}(x-1, y, z)];
        fval = -sum(abs(E_meas)); % This is the figure of merit.

        % Field gradient.
        grad_E = my_default_field(grid.shape, 0); 
        grad_E{2}(x, y, z) = -E{2}(x, y, z) / abs(E{2}(x, y, z));
        grad_E{2}(x-1, y, z) = -E{2}(x-1, y, z) / abs(E{2}(x-1, y, z));
    end
        
    [fval, grad_E] = fitness(E);


        % 
        % Calculate structural gradient.
        %

    if ~calc_grad % Skip if not needed.
        grad_f = nan;
        return
    end

    function [eps] = make_eps(params)
    % Function handle for creating the structure.
        [~, eps] = make_structure(wvlen, eps_val, params, flatten);
    end

    % Calculate the structural gradient.
    grad_f = maxopt_field_gradient(grid, E, @fitness, params, @make_eps, ...
                'solver_fun', @(eps) maxwell_solve(grid, eps, J), ...
                'check_gradients', false);
end


function [grid, eps, J] = make_structure(wvlen, eps_val, params, flatten)

    height = params(1);
    width = params(2);
    hole_xpos = abs(params(3:3:end)) + 60; % Leave a 120 nm gap in the center.
    hole_xlen = abs(params(4:3:end));
    hole_ylen = abs(params(5:3:end));


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2*pi / wvlen;
    delta = 20;
    x = -4000 : delta : 4000;
    y = -1000 : delta : 1000;
    z = -1000 : delta : 1000;

    if flatten
        z = 0;
    end

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z);


        %
        % Setup the structure.
        %

    % Structure constants.
    air_eps = 1;

    my_box = @(pos, siz) maxwell_box_smooth(pos, siz, 'smooth_dist', delta);
    
    % Draw the beam.
    eps = maxwell_shape(grid, eps, eps_val, ...
                        my_box([0 0 0], [1e9 width height]));

    % Draw rectangular holes.
    for k = 1 : length(hole_xpos)
        for l = [-1, 1]
            hole_pos = [l*hole_xpos(k) 0 0];
            hole_size = [hole_xlen(k) hole_ylen(k) 2*height];
            hole = my_box(hole_pos, hole_size);
            eps = maxwell_shape(grid, eps, air_eps, hole); 
        end
    end

    % maxwell_view(grid, eps, [], 'y', [nan nan 0]);
end
