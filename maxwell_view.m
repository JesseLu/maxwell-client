%% maxwell_view
% Slice view of a field and/or structure.

%%% Syntax
%
% * |maxwell_view(grid, eps, E, 'x', [nan 50 nan])|
%   visualizes the |y = 50| slice of the x-components of |eps| and |E|.
%
% * |maxwell_view(grid, eps, [], 'x', [nan 50 nan])|
%   visualizes only the structure.
%
% * |maxwell_view(grid, [], E, 'x', [nan 50 nan])|
%   visualizes only the field.
%
% |mu| and |H| can be visualized by replacing either |eps| or |E|.

function maxwell_view(grid, mat, F, dir, slice_ind)

