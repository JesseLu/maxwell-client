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

[grid, J, eps] = maxwell_grid(1550, -50:50, -80:80, -20:20);
[grid, J, eps, mu] = maxwell_grid(1550, -50:50, -80:80, -20:20); % If we want mu.
[grid, J, eps, mu] = maxwell_grid(1550, -50:50, -80:80, -20:20, options); 
... = maxwell_grid(1550, -50:50, -80:80, -20:20, 'hires_box', [dx dy dz], 'hires_box', delta);

eps = maxwell_epsilon(grid, eps, 10, my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size), options);
% Options include upsample_ratio, functions for averaging and eps/mu-modification.

[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], [+inf, 50, 50], [40 40 20]);
[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], [+inf, 50, 50], [40 40 20], 1); % Optional mode number.

E = maxwell_solve(grid, eps, J, 'vis_progress', 'both');
[E, H] = maxwell_solve(grid, [eps, mu], J, 'vis_progress', 'both');
 
P = maxwell_flux(grid, E, [20 -inf 3], [40 40 20])
P = maxwell_flux(grid, E, [E1, H1]); % Special case where [E1, H1] is output from maxwell_wgmode  

maxwell_view(grid, E{2}, 'y', 50);
% Various options for absolute value, movie, phase, with structure, transperancy, etc...
