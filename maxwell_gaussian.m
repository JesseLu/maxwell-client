%% maxwell_gaussian
% Free-space Gaussian mode (0th-order).


function [J] = maxwell_gaussian(grid, eps_mu, plane_pos, plane_size, ...
                                polarization, focal_length, beam_diameter)


        %
        % Validate and parse inputs.
        %

    validateattributes(polarization, {'char'}, ...
                    {'scalar'}, mfilename, 'polarization');
    if ~any(polarization == 'xyz')
        error('Polarization must be either ''x'', ''y'', or ''z''.');
    end

    validateattributes(beam_diameter, {'numeric'}, ...
                    {'positive', 'finite', 'real'}, mfilename, 'beam_diameter');

    mode_fun = gaussian([0 0 0], beam_diameter, find(polarization == 'xyz'));
    J = maxwell_fsmode(grid, eps_mu, plane_pos, plane_size, mode_fun, ...
                        'focal_length', focal_length);
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
        
