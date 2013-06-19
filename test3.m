function [E, H, err, grid, epsilon] = test(box_pos, box_size)
    [grid, epsilon, J] = maxwell_grid(0.3, -50:50, -80:80, -20:20);
        
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

    epsilon = maxwell_epsilon(grid, epsilon, @f, 10);
end
