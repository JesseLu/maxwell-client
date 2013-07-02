%% example1_ring.m
% Waveguide-coupled ring resonator.

        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    % The grid size if 6 x 6 x 3 microns with a resolution of 50 nm.
    % [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, 0); % Use this for 2D.
    [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);

    % Structure constants.
    height = 0.2;
    ring_radii = [1.3 1.0];
    wg_width = 0.4;
    wg_offset = -1.55;
    si_eps = 13;
    air_eps = 1;


        % 
        % Setup the waveguide (which will be used to couple to the ring).
        %

    % Draw coupling waveguide.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 wg_offset 0], [inf, wg_width, height]));



        %
        % Setup the ring.
        %

    % Draw ring.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_cyl([0 0 0], ring_radii(1), height));
    eps = maxwell_shape(grid, eps, air_eps, ...
                        maxwell_cyl([0 0 0], ring_radii(2), height));


    J = maxwell_wgmode(grid, eps, [-2 wg_offset 0], [+inf 2 2], 'mode_number', 1);

    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); % Visualize the excited waveguide.
