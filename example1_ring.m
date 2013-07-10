%% example1_ring.m
% Waveguide-coupled ring resonator.


function [omega, E, H, grid, eps] = example1_ring(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct('sim_only', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    if options.flatten
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, 0); % Use this for 2D.
        m = 2;
    else
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -3:0.05:3, -1.5:0.05:1.5);
        m = 1;
    end


        %
        % Setup the ring.
        %

    % Structure constants.
    height = 0.2;
    ring_radii = [1.4 1.0];
    si_eps = 13;
    air_eps = 1;

    % Draw ring.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_cyl([0 0 0], ring_radii(1), height));
    eps = maxwell_shape(grid, eps, air_eps, ...
                        maxwell_cyl([0 0 0], ring_radii(2), height));


        %
        % Solve for initial excitation.
        %

    % Excitation for the fundamental mode (of the ring's waveguide).
    J = maxwell_wgmode(grid, eps, [0 mean(ring_radii) 0], [+inf 2 2], 'mode_number', m);

    fprintf('Initial excitation -- ');
    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); % Visualize the excited waveguide.

    if options.sim_only
        omega = grid.omega;
        return
    end


        % 
        % Solve for the eigenmode.
        %
    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E); % Use this solution as an initial guess.
