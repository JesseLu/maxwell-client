
function [eps, mu] = maxwell_smooth_shape(grid, eps_mu, f)

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);
    my_validate_field(eps, grid.shape, 'eps', mfilename);

    if isempty(mu)
        validateattributes(val, {'double'}, {'scalar', 'nonnan', 'finite'}, ...
                           mfilename, 'eps_val'); 
    else
        my_validate_field(mu, grid.shape, 'mu', mfilename);
        validateattributes(val, {'double'}, {'numel', 2, 'nonnan', 'finite'}, ...
                           mfilename, '[eps_val mu_val]'); 
    end

    % Test if we can get a bounding box.
    % TODO: Check bounding box has non-zero (positive) volume.
    try 
        bnd_box = f();
        validateattributes(bnd_box, {'cell'}, {'numel', 2});
        validateattributes(bnd_box{1}, {'numeric'}, {'numel', 3});
        validateattributes(bnd_box{2}, {'numeric'}, {'numel', 3});
    catch
        bnd_box = {[-Inf -Inf -Inf], [Inf Inf Inf]};
    end

        %
        % Update material fields.
        %

    [eps_grid_pos, mu_grid_pos] = my_s2pos(grid); % Get grid positions.

    % Update components of epsilon.
    params = {bnd_box, f};
    for k = 1 : 3
        eps{k} = my_update(eps_grid_pos{k}, eps{k}, val(1), params{:}); 
        if ~isempty(mu)
            mu{k} = my_update(mu_grid_pos{k}, mu{k}, val(2), params{:});
        end
    end


function [mat] = my_update(pos, mat, val, box, f)
% Updates a single component of a field.

    % Determine box on which to evaluate f.
    for k = 1 : 3
        if box{1}(k) <= pos{k}(1)
            ind = 1;
        else
            ind = max(find(pos{k} <= box{1}(k)));
        end
        s{1}(k) = ind;

        if box{2}(k) >= pos{k}(end)
            ind = length(pos{k});
        else
            ind = min(find(pos{k} >= box{2}(k)));
        end
        s{2}(k) = ind;
        if s{1}(k) == s{2}(k) % 2D special case.
            s{2}(k) = s{2}(k) + 1;
        end
    end

    % Produce upsampled grid.
    for k = 1 : 3
        c{k} = s{1}(k) : s{2}(k);
    end

    % Obtain upsampled values.
    [x, y, z] = ndgrid(c{1}, c{2}, c{3});

    if multipoint
        inside_shape = f(x(:), y(:), z(:));
    else
        for k = 1 : numel(x)
            inside_shape(k) = f(x(k), y(k), z(k));
        end
    end
    inside_shape = reshape(inside_shape, size(x));

    % Downsample results by averaging.
    for i = 1 : (s{2}(1) - s{1}(1))
        for j = 1 : (s{2}(2) - s{1}(2))
            for k = 1 : (s{2}(3) - s{1}(3))
                fill_fraction(i, j, k) = f_avg(inside_shape(...
                                    (i-1)*up_ratio+[1:up_ratio], ...
                                    (j-1)*up_ratio+[1:up_ratio], ...
                                    (k-1)*up_ratio+[1:up_ratio]), dir);
            end
        end
    end

    % Make changes to the material.
    m = mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1);
    mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1) = ...
        reshape(f_rep(val, fill_fraction(:), m(:), dir), size(m));


