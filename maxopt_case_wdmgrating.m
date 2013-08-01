%% maxopt_case_wdmgrating
% Used to optimize a wavelength-splitting grating coupler.


function [fun, x0] = maxopt_case_wdmgrating(type, varargin)

        %
        % Parse inputs.
        %

    validateattributes(type, {'char'}, {'vector'}, 'type', mfilename);

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Return recommended starting parameters.
        %

    x0 = zeros(100, 2);
    wvlens = [1300 1500];


        %
        % Return appropriate function handle.
        %
        
    function [E, H, grid, eps] = get_fields(varargin)
        [~, E, H, grid, eps] = solve_structure(varargin{:});
    end

    switch type
        case 'get_fields'
            fun = @(x) get_fields(wvlens, x, options.flatten);
        case 'fval'
            fun = @(x) solve_structure(wvlens, x, options.flatten);
        otherwise
            error('Invalid type.');
    end
end


function [fval, E, H, grid, eps] = solve_structure(wvlens, params, flatten)

    N = length(wvlens);

    if flatten
        my_plane = [nan 0 nan];
    else
        my_plane = [nan nan 0];
    end

    % Initiate solves.
    for k = 1 : N
        subplot(1, N, k);
        [cb{k}, grid{k}, eps{k}, E_out{k}] = ...
            start_structure_solve(wvlens(k), params, flatten);
    end

    % Wait for solves to finish.
    for k = 1 : N
        while ~cb{k}(); end
        [~, E{k}, H{k}] = cb{k}();

        % Visualize.
        subplot(1, N, k);
        maxwell_view(grid{k}, eps{k}, E{k}, 'y', my_plane);
    end

%     done = false * ones(N, 1);
%     while ~all(done)
%         for k = 1 : N
%             [done(k), E{k}, H{k}] = cb{k}();
%         end
%     end

    % General fitness function.
    [vec, unvec] = my_vec(grid{1}.shape);
    function [fval, grad_E] = fitness(E, E_ref, P_in)
        overlap = vec(E_ref)' * vec(E);
        fval = (-1/P_in) * abs(overlap)^2;
        grad_E = (-1/P_in) * overlap * vec(E_ref);
    end
    
    % Compute fitness functions.
    for k = 1 : N
        fitness_fun{k} = @(E) fitness(E, E_out{k}{k}, 1);
        fval(k) = fitness_fun{k}(E{k})
    end
end
    
    

function [cb, grid, eps, E_out, H_out] = ...
                start_structure_solve(wvlen, params, flatten)
% Simulate the structure for the specific parameters.

        %
        % Get the structural parameters.
        %

    n = numel(params)/2;
    x_shift = params(1:n);
    y_shift = params(n+1:end);
    r_shift = zeros(n, 1); % No radius shift.


%     % With radius shifts.
%     x_shift = params(1:round(n/3));
%     y_shift = params(round(n/3)+1:round(2*n/3));
%     r_shift = params(round(2*n/3)+1:end);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 900 nm.
    delta = 40;
    omega = 2*pi/wvlen;

    my_size = 2500;
    x = -my_size : delta : my_size+500;
    y = -my_size : delta : my_size;
    z = -1000 : delta : 1000;

    if flatten
        y = -500;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z); 


        %
        % Setup the structure.
        %

    % Material permittivities.
    si = 13;
    ox = 2.25;
    air = 1;

    % Structural constants.
    wg_pos = [-500 500];
    wg_width = 220;
    height = 440;
    radius = 100;
    a = 400;

    % Draw the silicon slab. 
    eps = maxwell_shape(grid, eps, si, ...
                        maxwell_box([0 0 0], [4e3 4e3 height]));

    % Draw the silicon waveguides. 
    eps = maxwell_shape(grid, eps, si, ...
                maxwell_box([max(x) wg_pos(1) 0], [2*max(x) wg_width height]));
    eps = maxwell_shape(grid, eps, si, ...
                maxwell_box([max(x) wg_pos(2) 0], [2*max(x) wg_width height]));

    % Draw the holes.
    k = 1;
    for i = -4.5 : 4.5
        for j = -4.5 : 4.5
            pos = a * [i, j] + [x_shift(k), y_shift(k)];
            r = abs(radius + r_shift(k));
            my_cyl = maxwell_cyl_smooth([pos 0], r, 2*height, ...
                                        'smooth_dist', delta/2);
            eps = maxwell_shape(grid, eps, air, my_cyl);
            k = k + 1;
        end
    end

    % Draw oxide slab.
    eps = maxwell_shape(grid, eps, ox, ...
                maxwell_box([0 0 min(z)-height/2+delta], ...
                            [inf inf 2*abs(min(z))]));

%     subplot 121; maxwell_view(grid, eps, [], 'y', [nan wg_pos(1) nan]);
%     subplot 122; maxwell_view(grid, eps, [], 'y', [nan nan 0]);

    if flatten
        my_plane = [nan 0 nan];
    else
        my_plane = [nan nan 0];
    end
    maxwell_view(grid, eps, [], 'y', my_plane);

        %
        % Excite with a Gaussian wave.
        %

    J = maxwell_gaussian(grid, eps, [0 0 500], [2*my_size 2*my_size -inf], ...
                        'y', 500, 2000);

        %
        % Initiate solve.
        %
    
    cb = maxwell_solve_async(grid, eps, J);


        %
        % Calculate output modes.
        %

    for k = 1 : 2
        [J{k}, E_out{k}, H_out{k}] = ...
                maxwell_wgmode(grid, eps, [2100 wg_pos(k) 0], [+inf 1e3 1e3]);
    end
    
%     subplot 121; 
%     maxwell_view(grid, eps, E, 'y', [nan -500 nan], 'field_phase', nan); 
%     subplot 122; 
%     maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); 

end

