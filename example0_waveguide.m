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
    if options.flatten
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -1:0.05:1, 0); % Use this for 2D.
        mode_num = 1;
    else
        [grid, eps] = maxwell_grid(2*pi/1.55, -3:0.05:3, -1:0.05:1, -1:0.05:1);
        mode_num = 1;
    end



        %
        % Setup the waveguide.
        %

    % Structure constants.
    wg_height = 0.2;
    wg_width = 0.4;
    si_eps = 6;

    % Draw waveguide.
    eps = maxwell_shape(grid, eps, si_eps, ...
                        maxwell_box([0 0 0], [inf wg_width wg_height]));


        %
        % Solve for initial excitation.
        %

    % Excitation for the fundamental mode (of the ring's waveguide).
    J = maxwell_wgmode(grid, eps, [-2 0 0], [+inf 1 1], 'mode_number', mode_num);

    fprintf('Initial excitation -- ');
    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan nan 0], 'field_phase', 0); % Visualize the excited waveguide.


        % 
        % Measure the output powers.
        %

    xm = 0:0.1:2;
    for k = 1 : length(xm)
        P(k) = maxwell_flux(grid, [E H], [xm(k) 0 0], [+inf 100 100]);
    end

    P = P(:);
    plot([real(P), imag(P), abs(P)], '.-');

