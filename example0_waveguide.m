%% example0_waveguide.m
% Waveguide example.


function [E, H, grid, eps] = example0_waveguide(varargin)


        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'flatten', false), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2 * pi/1.55;
    x = -2 : 0.05 : 2;
    y = -1 : 0.05 : 1;
    z = -1 : 0.05 : 1;
    if options.flatten
        z = 0;
        mode_num = 1;
    else
        mode_num = 1;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z);
                            % 'hires_box', {[0 -.4 0], [1 .4 1], [.05 .05 .05]});


        %
        % Setup the waveguide.
        %

    % Structure constants.
    wg_height = 0.2;
    wg_width = 0.4;
    si_eps = 13;

    % Draw waveguide.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_smooth_box([0 0 0], [1e9 wg_width wg_height], ...
                                            'smooth_dist', 0.02));

        %
        % Solve for initial excitation.
        %

    [J, ~, ~, beta]  = maxwell_wgmode(grid, eps, [-1 0 0], [+inf 3 3], 'mode_number', mode_num);
    beta

    fprintf('Initial excitation -- ');
    [E, H] =  maxwell_solve(grid, eps, J);
    [A, x, b] = maxwell_axb(grid, eps, E, J);
    norm(A*x - b) / norm(b)
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); % Visualize the excited waveguide.


        % 
        % Measure the output powers.
        %

    % Solve wgmode for filtering out the mode we care about.
    [~, E1, H1] = maxwell_wgmode(grid, eps, [0 0 0], [+inf 3 3], 'mode_number', 1);

    P0 = maxwell_flux(grid, [E H], [0 0 0], [+inf 100 100]);
    P1 = maxwell_flux(grid, [E H], [E1 H1]);
    fprintf('Output powers at x = 0,\n');
    fprintf('Total power: %1.5f\n', P0);
    fprintf('Power in mode: %1.5f\n', P1);

    params = {'field_phase', nan, 'show_grid', true};
    maxwell_view(grid, eps, E, 'y', [nan nan 0], params{:});
