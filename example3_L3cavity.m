% Simulate L3 photonic crystal cavity.

function [omega, E, H, grid, eps] = example3_L3cavity(varargin)

        %
        % Get optional parameters.
        %

    options = my_parse_options(struct(  'flatten', false, ...
                                        'sim_only', false), ...
                                varargin, mfilename);

    
        %
        % Load the structure and create grid for it.
        %

    omega = 0.078;
    eps = getfield(load('l3.mat'), 'eps');
    dims = size(eps{1});

    if options.flatten % Make 2D, if needed.
        for k = 1 : 3
            eps{k} = eps{k}(:,:,round(dims(3)/2));
        end
        dims(3) = 1;
        omega = 0.063;
    end
        

    for k = 1 : 3
        xyz_pos{k} = [0:dims(k)] - round(dims(k)/2);
    end
    [grid, ~, ~, J] = maxwell_grid(omega, xyz_pos{:});

    
        %
        % Simulate with point source for initial guess field.
        %

    c = round(dims/2);
    J{2}(c(1), c(2), c(3)) = 1;

    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J); % Use this solution as an initial guess.

    if options.sim_only
        omega = grid.omega;
        return
    end


        %
        % Find the eigenmode, using previous result as initial guess field.
        %
    
    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E); 

