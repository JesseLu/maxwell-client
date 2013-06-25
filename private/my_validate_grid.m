function f = my_validate_grid(grid, name)
%     f = @val_grid;
% 
% function val_grid(grid)
    validateattributes(grid, {'struct'}, {'nonempty'}, name, 'grid', 1)

    if ~isfield(grid, 'omega')
        error(['MATLAB:' name, ':invalidType'], 'Expected input number 1, grid, to be a structure with field "omega".')
    else
        validateattributes(grid.omega, {'double'}, {'scalar'}, name, 'grid.omega', 1)
    end
