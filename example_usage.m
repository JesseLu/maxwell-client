%% Purpose
% The goal of this library is to provide a open-source toolbox
% for solving the linear electromagnetic wave equation.

%% Features
% * Open-source 
% * Simple
% * Extensible.
% * Fast (gpu solver in cloud)

%% Functionality
% * Isotropic materials (no anisotropy for now).
% * Epsilon and mu both covered.
% * E and H fields.
% * Extensible subgrid smoothing.
% * Only J (no M).


% Most general use-case

% Initialize simulation grid.
[grid, eps, mu, J] = maxwell_grid(2*pi/1550, -50:50, -80:80, -20:20, ...
                                    'nopml', 'xyz', 'num_pml_cells', 10);

[grid, eps, mu, J] = maxwell_grid(0.08 - 0.003i, -50:50, -80:80, -20:20);


% Place object in grid.
eps = maxwell_epsilon(grid, eps, 10, maxwell_box(box_pos, box_size));

[eps, mu] = maxwell_epsilon(grid, [eps, mu], [10, 3], maxwell_box(box_pos, box_size), ...
                            'upsample_ratio', 6, 'f_avg', @mean, 'f_rep', @f);

% Construct excitation source.
[J, E1, H1] = maxwell_wgmode(grid, [eps, mu], [40 40 20], [+inf 50 50], ...
                            'mode_number', 1);

[J, E1, H1] = maxwell_gaussian(grid, [eps, mu], [0 0 20], [80 80 -inf], ...
                                'focus', [0 0 0], 'waist', 40);

% Solve.
E = maxwell_solve(grid, eps, J);

[E, H] = maxwell_solve(grid, [eps, mu], J, 'vis_progress', 'plot', ...
                            'E0', E0, 'max_iters', 1e5, 'err_thresh', 1e-6);

cb = maxwell_solve_async(grid, eps, J);
while ~cb(); end
[~, E, H] = cb();

[omega, E, H] = maxwell_solve_eigenmode(grid, eps, E0, ...
                                        'eig_iters', 10, 'eig_err_thresh', 1e-6, ...
                                        'max_iters', 1e5, 'err_thresh', 1e-6);
 
% Check the simulation.
[A, x, b] = maxwell_axb(grid, eps, E, J);
err = norm(A*x-b) / norm(b);

[A, x, b] = maxwell_axb(grid, [eps mu], [E H], J);
err = norm(A*x-b) / norm(b);

% Calculate output powers.
P = maxwell_flux(grid, [E H], [40 40 20], [20 -inf 3])

P = maxwell_flux(grid, [E H], [E1, H1]); % Special case where [E1, H1] is output from maxwell_wgmode  

% Visualize.
maxwell_view(grid, eps, E, 'y', [nan nan 50]);

maxwell_view(grid, eps, [], 'y', [nan nan 50]);

maxwell_view(grid, [], E, 'y', [nan nan 50]);

