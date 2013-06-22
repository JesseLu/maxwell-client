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

[grid, eps, mu, J] = maxwell_grid(1550, -50:50, -80:80, -20:20, ...
                'nopml', , 'num_pml_cells', , 'hires_box', {start end delta});
... = maxwell_grid(1550, -50:50, -80:80, -20:20, 'hires_box', [dx dy dz], 'hires_box', delta);

eps = maxwell_epsilon(grid, eps, 10, my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size));
[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], my_box(box_pos, box_size), options);
% Options include upsample_ratio, functions for averaging and eps/mu-modification.

[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], [40 40 20], [+inf, 50, 50]);
[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], [40 40 20], [+inf, 50, 50], 1); % Optional mode number.

E = maxwell_solve(grid, eps, J, 'vis_progress', 'both');
[E, H] = maxwell_solve(grid, [eps, mu], J, 'vis_progress', 'both');
 
P = maxwell_flux(grid, E, [40 40 20], [20 -inf 3])
P = maxwell_flux(grid, E, [E1, H1]); % Special case where [E1, H1] is output from maxwell_wgmode  

maxwell_view(grid, eps, E, 'y', [4 12 inf]);
% Various options for absolute value, movie, phase, with structure, transperancy, etc...
