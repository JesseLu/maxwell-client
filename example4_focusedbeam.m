function [E, H, grid, eps] = example4_focusedbeam(varargin)

        %
        % Parse inputs.
        %

    options = my_parse_options(struct(  'flatten', false, ...
                                        'flen', 0), ...
                                varargin, mfilename);


        %
        % Create grid.
        %

    % Make a grid for a wavelength of 1550 nm.
    omega = 2 * pi/1.55;
    x = -5 : 0.1 : 5;
    y = -5 : 0.1 : 5;
    z = -3 : 0.1 : 3;
    if options.flatten
        y = 0;
        mode_num = 1;
    else
        mode_num = 1;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z);
                            % 'hires_box', {[0 -.4 0], [1 .4 1], [.05 .05 .05]});

    J = maxwell_gaussian(grid, eps, [0 0 2], [8 8 -inf], 'x', options.flen, 0.9);
   
%     mode_fun = gaussian([0 0 0], 0.2, 1);
%     J = maxwell_fsmode(grid, eps, [0 0 2], [8 8 -inf], mode_fun, 'focal_length', options.flen);

    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'x', [nan 0 nan], 'field_phase', nan); % Visualize the excited waveguide.
