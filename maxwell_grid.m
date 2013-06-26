%% maxwell_grid
% Create a simulation domain.

%%% Syntax
%
% * |[grid, eps] = maxwell_grid(omega, x, y, z)|
%   produces a grid for a specific angular frequency with
%   grid points located at the points specified at |x|, |y|, and |z|.
%   An initial variable for the permittivity |eps| is returned.
%   Perfectly matched layer (PML) absorbing boundaries
%   are automatically included in the grid.
%   Use the |'nopml'| option (below) to create grids without PMLs.
%
% * |... = maxwell_grid(2*pi/wavelength, ...)|
%   allows the user to define the angular frequency using a wavelength
%   parameter instead of an angular frequency parameter.
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

%%% Description
% |maxwell_grid| returns the |grid| variable along with initial values 
% for certain field vectors.
% |grid| defines the position and spacing of the simulation domain,
% and has the following important properties:
%
% * Periodic, wrap-around boundaries. 
%   This is very important in the case where PML has been removed 
%   (via the |'nopml'| option).
%
% * The number of unique grid points is one less than the length 
%   of the position arguments.
%   In other words, the number of unique grid points in the x-direction
%   is |length(x) - 1|.
%   This is a natural result of wrap-around boundaries.
%   

function [grid, eps, mu, J] = maxwell_grid(omega, x, y, z, varargin)

        %
        % Validate and parse inputs.
        %

    my_simple_check = @(var, var_name) ...
        validateattributes(var, {'double'}, {'nonnan', 'finite', 'vector'}, ...
                            mfilename, var_name);

    my_simple_check(omega, 'omega');
    validateattributes(omega, {'double'}, {'scalar'}, mfilename, 'omega');

    my_simple_check(x, 'x');
    my_simple_check(y, 'y');
    my_simple_check(z, 'z');

    % Optional arguments
    options = my_parse_options(struct('nopml', '', 'num_pml_cells', 10), ...
                                varargin);

    validateattributes(options.nopml, {'char'}, {}, mfilename, 'nopml');

    if length(options.num_pml_cells) == 1
        options.num_pml_cells = options.num_pml_cells * ones(3, 1);
    end
    validateattributes(options.num_pml_cells, {'numeric'}, ...
                        {'positive', 'integer', 'vector'}, mfilename, 'num_pml_cells');

 
        %
        % Initialize the grid structure.
        %

    grid = struct(  'omega', omega, ...
                    'origin', [x(1), y(1), z(1)], ...
                    'shape', [length(x), length(y), length(z)] - 1);
                    

        %
        % Compute the s-parameters for the grid (spacing between grid points).
        %

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

        % Add pml if needed.
        if ~any(options.nopml == xyz(k))
            [grid.s_prim{k}, grid.s_dual{k}] = ...
                stretched_coordinates(grid.omega, grid.origin(k), ...
                grid.s_prim{k}, grid.s_dual{k}, options.num_pml_cells(k));
        end
    end

    my_validate_grid(grid, mfilename); % Sanity check.

    
        %
        % Form the initial field vectors.
        %

    dims = grid.shape;
    eps = {ones(dims), ones(dims), ones(dims)};
    mu = {ones(dims), ones(dims), ones(dims)};
    J = {zeros(dims), zeros(dims), zeros(dims)};



function [s] = my_dual(w)
% Private function to compute s_dual.
    w_avg = (w + [w(2:end); (w(end) + (w(2)-w(1)))]) ./ 2;
    s = diff(w_avg);
