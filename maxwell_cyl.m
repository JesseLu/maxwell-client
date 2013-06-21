function [cyl_fun] = maxwell_cyl(center, radius, height)

    box_size = [radius, radius, height];
    bounding_box = {center - box_size/2, ...
                    center + box_size/2};
        
    function [out] = f(x, y, z)
        if nargin == 0 % Asking for bounding box.
            out = bounding_box;
        else
            out =   (x.^2 + y.^2 < radius^2) & ...
                    (abs(z-center(3)) < height/2);
        end
    end

    cyl_fun = @f;
end
