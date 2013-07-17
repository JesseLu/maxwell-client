function [E, H, grid, eps] = example4_focusedbeam(type, varargin)

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

    [grid, eps] = maxwell_grid(omega, x, y, z, ...
                            'hires_box', {[0 0 0], [1 .4 1], [.02 .02 .02]});
    switch type
        case 'gaussian'
            J = maxwell_gaussian(grid, eps, [0 0 2], [8 8 -inf], 'y', options.flen, 0.9);
        case 'donut'
            mode_fun = zdonut([0 0 0], 0.8);
            J = maxwell_fsmode(grid, eps, [0 0 2], [8 8 -inf], mode_fun, 'focal_length', options.flen);
        otherwise
            error('Type must either ''gaussian'' or ''donut''.');
    end
   

    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [nan 0 nan], 'field_phase', nan); % Visualize the excited waveguide.
end


function [fun] = zdonut(center, width)

    function [E] = mode_fun(w, x, y, z)
        r = sqrt(   (x - center(1)).^2 + ...
                    (y - center(2)).^2) + 1e-10;
        E = (w == 1) .* (y./r) .* sin(pi*r/width/2) .* (r < 2*width) + ...
            (w == 2) .* (x./r) .* sin(pi*r/width/2) .* (r < 2*width);
    end

    fun = @mode_fun;
end
        
