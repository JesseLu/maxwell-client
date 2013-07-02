%% example2_metalmode.m
% Mode of a metallic resonator

function [omega, E, H, grid] = example2_metalmode()
        %
        % Create grid.
        %

    delta = 1;
    [grid, eps, J] = maxwell_grid(2*pi/940, -150:delta:150, -150:delta:150, -100);
    % [grid, eps, J] = maxwell_grid(2*pi/940, -150:150, -150:150, -250:50);

    % Structure constants.
    radius = 100;
    height = 200;

    gaas = 3.5^2;
    silver = -40-4i;


        % 
        % Draw the GaAs and the silver shapes.
        %

    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_cyl([0 0 -height], radius, 2*height));
    % Staircase the following calls.
    eps = maxwell_shape(grid, eps, silver, ...
                        maxwell_box([0 0 -height/2], [inf, inf, height]));
    eps = maxwell_shape(grid, eps, gaas, ...
                        maxwell_cyl([0 0 -height/2], radius, height));
    maxwell_view(grid, eps, [], 'y', [nan nan -100]);


        %
        % Point excitation in center for initial excitation.
        %

    c = round(grid.shape/2); % Center.
    J{1}(c(1), c(2), c(3)) = 1;


        %
        % Solve.
        %

    function my_vis()
    end
    [E, H] =  maxwell_solve(grid, eps, J); % Use this solution as an initial guess.
    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E); % Use this solution as an initial guess.
end

