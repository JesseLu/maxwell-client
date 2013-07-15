function example4_focusedbeam(varargin)

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
    x = -3 : 0.1 : 3;
    y = -3 : 0.1 : 3;
    z = -4 : 0.1 : 1;
    if options.flatten
        x = 0;
        mode_num = 1;
    else
        mode_num = 1;
    end

    [grid, eps] = maxwell_grid(omega, x, y, z);
                            % 'hires_box', {[0 -.4 0], [1 .4 1], [.05 .05 .05]});

    mode_fun = gaussian([0 0 0], 2, 2);
    J = maxwell_fpmode(grid, eps, [0 0 0], [5 5 +inf], mode_fun);

    [E, H] =  maxwell_solve(grid, eps, J);
    maxwell_view(grid, eps, E, 'y', [0 nan nan], 'field_phase', 0); % Visualize the excited waveguide.
end

function [fun] = gaussian(center, fwhm, pol)
    sigma = fwhm / (2 * sqrt(2 * log(2)));

    function [E] = mode_fun(w, x, y, z)
        r = sqrt(   (x - center(1)).^2 + ...
                    (y - center(2)).^2 + ...
                    (z - center(3)).^2);
        E = (w == pol) .* ...
            1/(sigma*sqrt(2*pi)) * exp(-r.^2 / (2*sigma^2));
    end

    fun = @mode_fun;
    % fun = @(w, x, y, z) (w == pol) * ones(size(x));
end
        
