
function [cyl_fun] = maxwell_smooth_cyl(center, radius, height, varargin)


        %
        % Validate and parse inputs.
        %

    validateattributes(center, {'double'}, {'numel', 3, 'nonnan', 'finite'}, ...
                        mfilename, 'center');

    validateattributes(radius, {'double'}, {'scalar', 'positive'}, ...
                        mfilename, 'radius');

    validateattributes(height, {'double'}, {'scalar', 'positive'}, ...
                        mfilename, 'height');

    % Optional parameters.
    options = my_parse_options(struct('smooth_dist', 1), ...
                                varargin, mfilename);
    validateattributes(options.smooth_dist, {'numeric'}, ...
        {'positive', 'scalar'}, mfilename, 'smooth_dist');


        %
        % Create function handle.
        %

    box_size = [2*radius, 2*radius, height] + 2*options.smooth_dist;
    bounding_box = {center - box_size/2, ...
                    center + box_size/2};
        
    function [out] = f(x, y, z)
        if nargin == 0 % Asking for bounding box.
            out = bounding_box;
        else
            dist_r = (radius - sqrt((x-center(1)).^2 + (y-center(2)).^2));
            dist_z = (height/2 - abs(z-center(3)));
            out = my_val_clamp(dist_r, options.smooth_dist) .* ...
                    my_val_clamp(dist_z, options.smooth_dist);
        end
    end

    cyl_fun = @f;
end


