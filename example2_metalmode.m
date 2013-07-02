%% example2_metalmode.m
% Mode of a metallic resonator

function [omega, E, H, grid, eps] = example2_metalmode(varargin)
        %
        % Create grid.
        %

    delta = 10;
    % [grid, eps, J] = maxwell_grid(2*pi/1200, -150:delta:150, -150:delta:150, -100);
    [grid, eps, J] = maxwell_grid(2*pi/940, -150:delta:150, -150:delta:150, -250:delta:50);

    if isempty(varargin)
        % Structure constants.
        radius = 100;
        height = 200;
        smooth_dist = 10;

        gaas = 3.5^2;
        silver = -40+4i;
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
        maxwell_view(grid, eps, [], 'y', [nan nan -100]);
%         omega = nan;
%         E = nan;
%         H = nan;
%         return

    else
        eps = varargin{1};
    end


        %
        % Point excitation in center for initial excitation.
        %

    c = round(grid.shape/2); % Center.
    J{1}(c(1)+[-10:10], c(2)+[-10:10], c(3)) = 1;


        %
        % Solve.
        %

    function my_vis()
    end
    omega = grid.omega;
    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J, 'vis_progress', 'both'); % Use this solution as an initial guess.
    % [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E); % Use this solution as an initial guess.
end

function [res] = my_pass(z, ~)
    res = z;
end

function [res] = staircase_fun(z, dir)
    res = round(sum(z(:))/numel(z)); % Staircase the following calls.
end
