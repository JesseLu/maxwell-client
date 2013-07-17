%% maxwell_smooth_box
% Used to insert a "smoothed" box of constant epsilon/mu into the simulation grid.

%%% Syntax
%
% * |fun = maxwell_smooth_box(center, box_size)| 
%   returns function handle |fun| which describes a rectangular prism 
%   centered at |center| and of size |box_size|.
%


%%% Source code
function [box_fun] = maxwell_smooth_box(center, box_size, varargin)


        %
        % Validate and parse inputs.
        %

    validateattributes(center, {'double'}, {'numel', 3, 'nonnan', 'finite'}, ...
                        mfilename, 'center');

    validateattributes(box_size, {'double'}, {'numel', 3, 'positive', 'finite'}, ...
                        mfilename, 'box_size');

    % Optional parameters.
    options = my_parse_options(struct('smooth_dist', 1), ...
                                varargin, mfilename);
    validateattributes(options.smooth_dist, {'numeric'}, ...
        {'positive', 'scalar'}, mfilename, 'smooth_dist');


        %
        % Create function handle.
        %


    bounding_box = {center - box_size/2 - options.smooth_dist, ...
                    center + box_size/2 + options.smooth_dist};
        
    function [out] = f(x, y, z)
        if nargin == 0 % Asking for bounding box.
            out = bounding_box;
        else
            s = options.smooth_dist;
            out =   my_val_clamp(box_size(1)/2 - abs(x-center(1)), s) .* ...
                    my_val_clamp(box_size(2)/2 - abs(y-center(2)), s) .* ...
                    my_val_clamp(box_size(3)/2 - abs(z-center(3)), s);
        end
    end

    box_fun = @f;
end


