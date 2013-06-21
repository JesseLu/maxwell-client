%% Purpose
% The goal of this library is to provide a open-source toolbox
% for solving the linear electromagnetic wave equation.

%% Features
% * Open-source 
% * Simple
% * Extensible.
% * Fast (gpu solver in cloud)

%% Functionality
% * Dynamix grid.
% * Isotropic materials (no anisotropy for now).
% * Epsilon and mu both covered.
% * E and H fields.
% * Extensible subgrid smoothing.
% * Only J (no M).


% Most general use-case

[grid, J, eps] = maxwell_grid(0.3, -50:50, -80:80, -20:20);
[grid, J, eps, mu] = maxwell_grid(0.3, -50:50, -80:80, -20:20); % If we want mu.
[grid, J, eps, mu] = maxwell_grid(0.3, -50:50, -80:80, -20:20, options); % If we want mu.


eps = maxwell_epsilon(grid, eps, 10, my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size), options);

J = maxwell_wgmode(grid, [eps, mu], 'x+', [40 40 20], [0 50 30], 1);
[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], 'x+', [40 40 20], [0 50 30], 1);

E = maxwell_solve(grid, eps, J, 'vis_progress', 'both');
[E, H] = maxwell_solve(grid, [eps, mu], J, 'vis_progress', 'both');

maxwell_flux(grid, E, 'x+', [40 40 20], [0 50 30])
maxwell_flux(grid, E, E1)

maxwell_totalview(grid, [eps, mu], [E, H]);
