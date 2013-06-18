function [grid, epsilon, J] = maxwell_grid(omega, x, y, z, varargin)

    % Check omega.
    if numel(omega) ~= 1 | isnan(omega) 
        error('OMEGA must be a scalar.');
    end

    grid.omega = omega;
    grid.origin = [x(1), y(1), z(1)];

    pos = {x(:), y(:), z(:)};
    for k = 1 : 3
        if length(x) == 1
            grid.s_prim{k} = Inf;
            grid.s_prim{k} = Inf;
        else
            grid.s_prim{k} = diff(pos{k});
            grid.s_dual{k} = my_dual(pos{k});
        end
    end

    dims = [length(grid.s_prim{1}), ...
            length(grid.s_prim{2}), ...
            length(grid.s_prim{3})];

    epsilon = {ones(dims), ones(dims), ones(dims)};
    J = {zeros(dims), zeros(dims), zeros(dims)};


function [s] = my_dual(w)
    w_avg = (w + [w(2:end); (w(end) + (w(2)-w(1)))]) ./ 2;
    s = diff(w_avg);
