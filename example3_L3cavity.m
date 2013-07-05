% Simulate L3 photonic crystal cavity.

function [omega, E, H, grid, eps] = example3_L3cavity(varargin)

    [omega, E, H, grid, eps] = deal(nan);

    omega = 0.078; % Our guess.
    eps = getfield(load('l3.mat'), 'eps');
    dims = size(eps{1});

    for k = 1 : 3
        xyz_pos{k} = [0:dims(k)] - round(dims(k)/2);
    end

    [grid, ~, ~, J] = maxwell_grid(omega, xyz_pos{:});
    c = round(dims/2);
    J{2}(c(1), c(2), c(3)) = 1;

    [E, H] =  maxwell_solve(grid, eps, J, 'err_thresh', 1e-6, 'vis_progress', 'both'); % Use this solution as an initial guess.



