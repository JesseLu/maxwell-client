%% example2_metalmode.m
% Mode of a metallic resonator

function [omega, E, H, grid, eps] = example2_metalmode(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 25, ...
                                        'flatten', false, ...
                                        'hires_delta', [3 3 6], ...
                                        'view_only', false, ...
                                        'sim_only', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    omega = 2*pi/900;
    extra = 400;
    % extra = 0;
    x = -200-extra : options.delta : 200+extra;
    y = -200-extra : options.delta : 200+extra;
    z = -300-extra : options.delta : 100+extra;
    if options.flatten
        x = 0;
    end

    hires_option = {[0 0 -100], [120 120 220], [options.hires_delta]};

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z, ...
                                        'hires_box', hires_option, ...
                                        'growth_rate', 1.3);


    % Structure constants.
    radius = 50;
    height = 200;
    smooth_dist = 1e-5;

    gaas = 3.5^2;
    silver = -40-4i;
    air = 1;
    my_inf = 1e4;


        % 
        % Draw the GaAs and the silver shapes.
        %

    fprintf('Building shapes [%dx%dx%d grid]...\n', grid.shape);

    eps = maxwell_shape(grid, eps, silver, ...
                        maxwell_smooth_box([0 0 -height/2], [my_inf, my_inf, height], ...
                                            'smooth_dist', smooth_dist), ...
                        'upsample_ratio', 1, 'f_avg', @my_pass);
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_smooth_cyl([0 0 -height], radius, 2*height, ...
                                            'smooth_dist', smooth_dist), ...
                        'upsample_ratio', 1, 'f_avg', @my_pass);
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_smooth_box([0 0 z(1)], [my_inf my_inf 2*(-z(1)-height)], ...
                                            'smooth_dist', smooth_dist), ...
                        'upsample_ratio', 1, 'f_avg', @my_pass);

    if options.view_only
        maxwell_view(grid, eps, [], 'y', [0 nan nan]);
        omega = nan;
        E = nan;
        H = nan;
        return
    end



        %
        % Get excitation.
        %

    c = round(grid.shape/2);
    if options.flatten % Use plane wave to solve 2D mode.
        J{2}(:, :, end-12) = 1;
        J{2}(c(1), c(2), c(3)) = 1;
        % J = maxwell_wgmode(grid, eps, [0 0 100], [1e4 1e4 -inf], 'mode_number', 2, 'view', true);
        fprintf('=== 2D solve ===\n');
    else
        % To get the current excitation for 3D, use the 2D mode.
        % We do this via a recursive call.
        [~, E] = example2_metalmode('flatten', true); 
        J{2}(:,round(grid.shape(2)/2), :) = E{2};
        fprintf('=== 3D solve ===\n');
    end


        %
        % Solve.
        %

    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J); % Use this solution as an initial guess.
    maxwell_view(grid, eps, E, 'y', [0 nan nan]);

    if options.sim_only
        omega = grid.omega;
        return
    end

    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E, 'eig_err_thresh', 1e-10); % Use this solution as an initial guess.
end

function [res] = my_pass(z, ~)
    res = z;
end

function [res] = staircase_fun(z, dir)
    res = round(sum(z(:))/numel(z)); % Staircase the following calls.
end
