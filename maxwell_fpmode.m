%% maxwell_fpmode
% Excite general free-space modes.


function [J] = maxwell_fpmode(grid, eps_mu, plane_pos, plane_size, mode_fun, ...
                                varargin)


        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    if isempty(mu)
        mu = my_default_field(grid.shape, 1); 
    end
    my_validate_field(eps, grid.shape, 'eps', mfilename);
    my_validate_field(mu, grid.shape, 'mu', mfilename);

    validateattributes(plane_pos, {'numeric'}, ...
                {'nonnan', 'finite', 'numel', 3}, mfilename, 'plane_pos');
    validateattributes(plane_size, {'numeric'}, ...
                {'nonnan', 'numel', 3}, mfilename, 'plane_size');
    if length(find(isinf(plane_size))) ~= 1
        error('plane_size must have exactly 1 element equal to either +inf or -inf.');
    end

    % Optional arguments.
    options = my_parse_options(struct(  'focal_length', 0, ...
                                        'view', false), ...
                                varargin, mfilename);
    validateattributes(options.focal_length, {'numeric'}, ...
                        {'real', 'nonnan'}, mfilename, 'focal_length');
    validateattributes(options.view, {'logical'}, ...
                        {'binary'}, mfilename, 'view');


        %
        % Find plane (sub-grid) on which to put the free-space mode.
        %

    % Determine desired direction of propagation.
    [p0, p1, prop_dir, prop_in_pos_dir] = ...
                                    my_find_plane(grid, plane_pos, plane_size);

    % Get position information.
    pos = my_s2pos(grid);

    % Cut out the bounded plane.
    sub_shape = p1 - p0 + 1;
    for k = 1 : 3
        for l = 1 : 3
            pos{l}{k} = pos{l}{k}(p0(k):p1(k));
        end
        e{k} = eps{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        m{k} = mu{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end

    % Get uniform values of epsilon and mu.
    eps = [e{1}(:); e{2}(:); e{3}(:)];
    mu =  [m{1}(:); m{2}(:); m{3}(:)];
    eps_val = mean(eps);
    mu_val = mean(mu);

    if any(eps ~= eps_val) || any(mu ~= mu_val)
        error('Material parameters must not vary over excitation plane.');
    end 
    
    % Step size in propagation direction.
    prop_step = real(grid.s_dual{prop_dir}(p0(prop_dir))); 


        %
        % Get the mode shape.
        %

    for k = 1 : 3
        [x, y, z] = ndgrid(pos{k}{1}, pos{k}{2}, pos{k}{3});
        E{k} = mode_fun(k, x, y, z);
    end

        
        %
        % Adjust to correct focal length.
        %

    E{k} = my_propagate_beam(real(grid.omega), prop_dir, ...
                                -options.focal_length, E, pos);


        %
        % Form current source.
        %

    J = my_default_field(grid.shape, 0);
    for k = 1 : 3
        J{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3)) = E{k};
    end
