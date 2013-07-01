%% example1_ring.m
% Waveguide-coupled ring resonator.

    %
    % Create grid.
    %

% Make a grid for a wavelength of 1550 nm.
% The grid size if 6 x 6 x 3 microns with a resolution of 50 nm.
[grid, eps] = maxwell_grid(2*pi/1.55, -3:0.025:3, -3:0.025:3, 0); % Use this for 2D.
% [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);

% Structure constants.
height = 0.5;
ring_radii = [1.5 1.0];
wg_width = 0.5;
wg_offset = 0;
si_eps = 13;
air_eps = 1;


    % 
    % Setup the waveguide (which will be used to couple to the ring).
    %

% Draw coupling waveguide.
eps = maxwell_shape(grid, eps, si_eps, ...
                    maxwell_box([0 wg_offset 0], [inf, wg_width, height]));
maxwell_view(grid, eps, [], 'z', [nan nan 0]); % Visualize the waveguide.

J = maxwell_wgmode(grid, eps, [-1 wg_offset 0], [+inf 2 2], 'mode_number', 2, 'pause_and_view', true);
[E, H] =  maxwell_solve(grid, eps, J);
subplot 111
maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', nan); % Visualize the excited waveguide.


    %
    % Setup the ring.
    %

% Draw ring.
eps = maxwell_shape(grid, eps, si_eps, ...
                    maxwell_cyl([0 0 0], ring_radii(1), height));
eps = maxwell_shape(grid, eps, air_eps, ...
                    maxwell_cyl([0 0 0], ring_radii(2), height));

