%% maxwell_grid
% Create a simulation domain.

%%% Syntax
%
% * |[grid, eps] = maxwell_grid(wavelength, x, y, z)|
%   produces a grid for a specific wavelength with
%   grid points located at the points specified at |x|, |y|, and |z|.
%   An initial variable for the permittivity |eps| is returned.
%   Perfectly matched layer (PML) absorbing boundaries
%   are automatically included in the grid.
%   Use the |'nopml'| option (below) to create grids without PMLs.
%
% * |[grid, eps, mu, J] = maxwell_grid(...)|
%   also returns initial variables for the permeability, |mu|,
%   and the excitation current, |J|.
%
% * |... = maxwell_grid(..., 'nopml', xyz, 'num_pml_cells', n_pml)
%   allows for grids without PMLs in one or many directions of the grid.
%   |xyz| can be any combination of |'x'|, |'y'|, or |'z'|,
%   in order to not include PMLs in the x-, y-, or z-directions.
%   Additionally, the |'num_pml_cells'| option can be used to 
%   change the number of grid points used for the PML layer
%   (defaults to 10). 
%

function [grid, epsilon, J] = maxwell_grid(omega, x, y, z, varargin)

    % Check omega.
    if numel(omega) ~= 1 | isnan(omega) 
        error('OMEGA must be a scalar.');
    end

    % Here are the default options
    grid_options = struct(  'nopml', '', ...
                            'num_pml_cells', 10); 

    for k = 2 : 2 : length(varargin)
        grid_options = setfield(grid_options, varargin{k-1}, varargin{k});
    end

    if numel(grid_options.num_pml_cells) == 1
        grid_options.num_pml_cells = grid_options.num_pml_cells * ones(1, 3);
    end
      
    grid.omega = omega;
    grid.origin = [x(1), y(1), z(1)];

    pos = {x(:), y(:), z(:)};
    xyz = 'xyz';
    for k = 1 : 3
        if length(x) == 1
            grid.s_prim{k} = Inf;
            grid.s_prim{k} = Inf;
        else
            grid.s_prim{k} = diff(pos{k});
            grid.s_dual{k} = my_dual(pos{k});
        end

        % Add pml.
        if ~any(grid_options.nopml == xyz(k))
            [grid.s_prim{k}, grid.s_dual{k}] = ...
                stretched_coordinates(grid.omega, grid.origin(k), ...
                grid.s_prim{k}, grid.s_dual{k}, grid_options.num_pml_cells(k));
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
