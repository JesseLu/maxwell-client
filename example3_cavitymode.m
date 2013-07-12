% Simulate L3 photonic crystal cavity.

function [omega, E, H, grid, eps] = example3_cavitymode(cavity_type, varargin)

        %
        % Get optional parameters.
        %

    options = my_parse_options(struct(  'flatten', false, ...
                                        'omega_guess', [], ...
                                        'central_Jy', [], ...
                                        'sim_only', false), ...
                                varargin, mfilename);

    
        %
        % Load the structure and create grid for it.
        %

    switch cavity_type
        case 'L3'
            filename = 'l3.mat';
            omega_guess = struct('D2', 0.063, 'D3', 0.078);
        case 'beam'
            filename = 'beam.mat';
            omega_guess = struct('D2', 0.062, 'D3', 0.080);
        case 'beam2w'
            filename = 'beam.mat';
            omega_guess = struct('D2', 0.118, 'D3', 0.166);
        otherwise
            error('cavity_type must either be ''L3'' or ''beam''.');
    end

    eps = getfield(load(filename), 'eps');
    omega = omega_guess.D3; % Guess frequency for 3D.
    dims = size(eps{1});

    if options.flatten % Make 2D, if needed.
        for k = 1 : 3
            eps{k} = eps{k}(:,:,round(dims(3)/2));
        end
        omega = omega_guess.D2; % Guess frequency for 2D.
        dims(3) = 1;
    end

    if ~isempty(options.omega_guess)
        omega = options.omega_guess;
    end
        

    for k = 1 : 3
        xyz_pos{k} = [0:dims(k)] - round(dims(k)/2);
    end
    [grid, ~, ~, J] = maxwell_grid(omega, xyz_pos{:});

    
        %
        % Obtain initial field by doing a simulation!
        %

    c = round(dims/2);
    if options.flatten 
        if strcmp(cavity_type, 'beam2w')
            [~, E] = example3_cavitymode('beam', 'flatten', true);
            J{2} = abs(E{2});
        else
            % Use point source to excite 2D mode.
            J{2}(c(1), c(2), c(3)) = 1;
        end
        fprintf('=== 2D solve ===\n');
    else
        % To get the current excitation for 3D, use the 2D mode.
        % We do this via a recursive call.
        [~, E] = example3_cavitymode(cavity_type, 'flatten', true); 
        J{2}(:,:,c(3)) = E{2};
        fprintf('=== 3D solve ===\n');
    end

    if ~isempty(options.central_Jy)
        J{2}(:,:,c(3)) = options.central_Jy;
    end

    fprintf('Solving for initial field... ');
    [E, H] =  maxwell_solve(grid, eps, J); % Use this solution as an initial guess.

    maxwell_view(grid, eps, E, 'y', [nan nan 0]);

    if options.sim_only
        omega = grid.omega;
        return
    end


        %
        % Find the eigenmode, using previous result as initial guess field.
        %
    
    [omega, E, H] =  maxwell_solve_eigenmode(grid, eps, E);

