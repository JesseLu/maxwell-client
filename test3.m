function [eps] = test(box_pos, box_size)
    [grid, eps, J] = maxwell_grid(0.3, -20:20, -20:20, -20:20);
    % eps = maxwell_shape(grid, eps, 10, my_box(box_pos, box_size));
    % eps = maxwell_shape(grid, eps, 10, maxwell_box(box_pos, box_size));
    eps = maxwell_shape(grid, eps, 2, maxwell_cyl(box_pos, 6, 6));
end


function [fun] = my_box(box_pos, box_size)

    function [out] = f(x, y, z)
        if nargin == 0 % Asking for bounding box.
            out{1} = box_pos - box_size/2;
            out{2} = box_pos + box_size/2;
            return
        end
        out =   (abs(x-box_pos(1)) < box_size(1)/2) & ...
                (abs(y-box_pos(2)) < box_size(2)/2) & ...
                (abs(z-box_pos(3)) < box_size(3)/2);
    end
    fun = @f;
end


