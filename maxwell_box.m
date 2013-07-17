%% maxwell_box
% Box of constant epsilon/mu within the simulation grid.

%%% Syntax
%
% * |fun = maxwell_box(center, box_size)| 
%   returns function handle |fun| which describes a rectangular prism 
%   centered at |center| and of size |box_size|.

%%% Source code  
function [box_fun] = maxwell_box(center, box_size)

    validateattributes(center, {'double'}, ...
        {'vector', 'numel', 3}, mfilename, 'center');
    validateattributes(box_size, {'double'}, ...
        {'vector', 'positive', 'numel', 3}, mfilename, 'box_size');

    bounding_box = {center - box_size/2, ...
                    center + box_size/2};
        
    function [out] = f(x, y, z)
        if nargin == 0 % Asking for bounding box.
            out = bounding_box;
        else
            out =   (abs(x-center(1)) < box_size(1)/2) & ...
                    (abs(y-center(2)) < box_size(2)/2) & ...
                    (abs(z-center(3)) < box_size(3)/2);
        end
    end

    box_fun = @f;
end
