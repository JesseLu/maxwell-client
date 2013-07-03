%% example2_metalmode.m
% Mode of a metallic resonator

function [omega, E, H, grid, eps] = example2_metalmode(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'delta', 25, ...
                                        'flatten', false, ...
                                        'hires_delta', [1 1 2], ...
                                        'view_only', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    omega = 2*pi/940;
    x = -200 : options.delta : 200;
    y = -200 : options.delta : 200;
    z = -300 : options.delta : 100;
    if options.flatten
        x = 0;
    end

    hires_option = {[0 0 -100], [100 100 200]-1, [options.hires_delta]};

    [grid, eps, ~, J] = maxwell_grid(omega, x, y, z, 'hires_box', hires_option);

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

    fprintf('Building shapes...\n');

    eps = maxwell_shape(grid, eps, silver, ...
                        maxwell_smooth_box([0 0 -height/2], [my_inf, my_inf, height], ...
                                            'smooth_dist', smooth_dist), ...
                        'upsample_ratio', 1, 'f_avg', @my_pass);
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_smooth_cyl([0 0 -height], radius, 2*height, ...
                                            'smooth_dist', smooth_dist), ...
                        'upsample_ratio', 1, 'f_avg', @my_pass);
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_smooth_box([0 0 -2*height], [my_inf my_inf 2*height], ...
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
        % Point excitation in center for initial excitation.
        %

    c = round(grid.shape/2); % Center.
    % J{1}(:, :, end-12) = 1;
    J{2}(c(1), c(2), c(3)) = 1;
    % J{1}(c(1)+[-5:5], c(2)+[-5:5], c(3)) = 1;


        %
        % Solve.
        %

    function my_vis()
    end
    omega = grid.omega;
    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J, 'vis_progress', 'both'); % Use this solution as an initial guess.
    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E, 'eig_err_thresh', 1e-14); % Use this solution as an initial guess.
    fprintf('\n');
end

function [res] = my_pass(z, ~)
    res = z;
end

function [res] = staircase_fun(z, dir)
    res = round(sum(z(:))/numel(z)); % Staircase the following calls.
end
