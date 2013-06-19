% Add shapes to epsilon data structure.
function [epsilon] = maxwell_epsilon(grid, epsilon, f, eps_val, varargin)

    % Here are the default options
    options = struct(  'upsample_ratio', 6);
    for k = 2 : 2 : length(varargin)
        options = setfield(options, varargin{k-1}, varargin{k});
    end

    % Test if we can get a bounding box.
    % TODO: Check bounding box has non-zero (positive) volume.
    try 
        bnd_box = f();
    catch
        bnd_box = {[-Inf -Inf -Inf], [Inf Inf Inf]};
    end

    % Test if we can give multiple points to f.
    try 
        out = f([0 1], [2 3], [4 5]);
        multipoint = true;
    catch
        multipoint = false;
    end

    % Build up grid info.
    origin_prim = grid.origin(:);
    origin_dual = grid.origin(:) + [real(grid.s_prim{1}(1))/2; ...
                                    real(grid.s_prim{2}(1))/2; ...
                                    real(grid.s_prim{3}(1))/2];
    for k = 1 : 3
        pos_prim{k} = origin_prim(k) + [0; cumsum(real(grid.s_prim{k}))];
        pos_dual{k} = origin_dual(k) + [0; cumsum(real(grid.s_dual{k}))];
    end

    for k = 1 : 3
        for l = 1 : 3
            eps_grid_pos{k}{l} = pos_dual{l};
        end
        eps_grid_pos{k}{k} = pos_prim{k};
    end

    % Update components of epsilon.
    for k = 1 : 3
        epsilon{k} = my_update(eps_grid_pos{k}, epsilon{k}, ...
                                bnd_box, f, eps_val, ...
                                options.upsample_ratio, multipoint);
    end


function [epsilon] = my_update(pos, epsilon, box, f, val, up_ratio, multipoint)
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
    end

    % Produce upsampled grid.
    for k = 1 : 3
        c{k} = zeros(up_ratio * (s{2}(k) - s{1}(k)), 1);
        cnt = 0;
        for l = s{1}(k) : s{2}(k)-1
            c{k}(cnt*up_ratio+[1:up_ratio]) = pos{k}(l) + ...
                ([0.5 : up_ratio] ./ (up_ratio * (pos{k}(l+1)-pos{k}(l))));
            cnt = cnt + 1;
        end
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

    % Downsample results by averaging (TODO: other option later).
    ds_fun = @(z) mean(z(:));
    for i = 1 : (s{2}(1) - s{1}(1))
        for j = 1 : (s{2}(2) - s{1}(2))
            for k = 1 : (s{2}(3) - s{1}(3))
                ds_eps(i, j, k) = ds_fun(inside_shape(...
                                    (i-1)*up_ratio+[1:up_ratio], ...
                                    (j-1)*up_ratio+[1:up_ratio], ...
                                    (k-1)*up_ratio+[1:up_ratio]));
            end
        end
    end

    % Make changes to epsilon.
    epsilon(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1) = ds_eps * val + ...
        (1 - ds_eps) .* epsilon(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1);
