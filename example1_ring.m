%% example1_ring.m
% Waveguide-coupled ring resonator.

% Make a grid for a wavelength of 1550 nm.
% The grid size if 6 x 6 x 3 microns with a resolution of 50 nm.
[grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);

% Structure constants.
height = 0.5;
ring_radius = [1.5 1.0];
wg_width = 0.5;
si_eps = 13;
air_eps = 1;

% eps = maxwell_shape(grid, eps, si_eps, ...
%                     maxwell_cyl([0 0 0], ring_radii(1), height));
% eps = maxwell_shape(grid, eps, air_eps, ...
%                     maxwell_cyl([0 0 0], ring_radii(2), height));
eps = maxwell_shape(grid, eps, si_eps, ...
                    maxwell_box([0 0 0], [inf, wg_width, height]));
