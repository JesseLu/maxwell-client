%% stretched_coordinates
% Produces s-parameters needed to implement stretched-coordinate PML.

%% Description
% Produces the s-parameters needed to implement the stretched-coordinate 
% perfectly-matched layer (PML) boundary.
%
% Steven Johnson has a great reference on this.
%

function [s_prim, s_dual] = make_scpml(omega, origin, s_prim, s_dual, num_cells)

%% Input parameters
% * |omega| is the angular frequency of the simulation.
% * |t_pml| represents the depth, in grid points, of the pml in each direction.
%   It is also a 3-element vector. 
%   For no pml in a particular direction, set that element to 0.

%% Output parameters
% * |s_prim, s_dual| the primal and dual s-parameters for the grid.

%% Example
%
%   omega = 0.08;
%   dims = [80 40 20];
%   t_pml = [10 10 10];
%   [s_prim, s_dual] = stretched_coordinates(omega, dims, t_pml);

%% Source code

    % Position functions.
    w_p = origin + [0 cumsum(s_prim)];
    w_s = mean(w_p(1:2)) + [0 cumsum(s_dual)];

    % Define the borders.
    if isempty(find(w_p == 0))
        % 0 position does not occur on the primary grid,
        % act as if it occurs on the dual grid.
        border = [(w_s(1) - s_dual(end)), w_s(end)];

    else
        % 0-position occurs on the primary grid.
        border = [w_p(1), w_p(end)];
    end
    % Helper functions.
    pos = @(z) (z > 0) .* z; % Only take positive values.
    l = @(u, n, t) pos(t - u) + pos(u - (n - t)); % Distance to nearest pml boundary.

    % Compute the stretched-coordinate grid spacing values.
    for k = 1 : 3
        if t_pml(k) > 0 % PML requested in this direction.
            s_prim{k} = 1 - i * (4 / omega) * ...
                            (l(0:dims(k)-1, dims(k), t_pml(k)) / t_pml(k)).^4;
            s_dual{k} = 1 - i * (4 / omega) * ...
                            (l(0.5:dims(k)-0.5, dims(k), t_pml(k)) / t_pml(k)).^4;

        else % No PML requested in this direction 
            s_prim{k} = ones(1, dims(k));
            s_dual{k} = ones(1, dims(k));
        end
    end

