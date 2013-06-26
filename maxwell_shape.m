%% maxwell_shape
% Generic function to add constant material shapes to the simulation domain.

%%% Syntax
%
% * |eps = maxwell_shape(grid, eps, eps_val, shape_fun)|
%   modifies the shape structure |eps| by inserting |eps_val| 
%   in the volume described by |shape_fun|.
%
% * |[eps, mu] = maxwell_shape(grid, [eps, mu], [eps_val mu_val], shape_fun)|
%   modifies both |eps| and |mu|.
%
% * |... = maxwell_shape(..., 'upsample_ratio', ratio)|
%   allows for an upsampling ratio of |ratio|, defaults to |ratio = 6|.
%
% * |... = maxwell_shape(..., 'f_avg', f_avg, 'f_rep', f_rep)|
%   allows for custom functions which determine the averaging function
%   for a grid point (|f_avg|) and how values of |eps| (and |mu|) 
%   are replaced (|f_rep|).

%%% Description
% |maxwell_shape| modifies
function [eps, mu] = maxwell_shape(grid, eps_mu, val, f, varargin)


        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    [eps, mu] = my_split(eps_mu, grid.shape, {'eps', 'mu'}, mfilename);

    if isempty(mu)
        validateattributes(val, {'double'}, {'scalar', 'nonnan', 'finite'}, ...
                           mfilename, 'eps_val'); 
    else
        validateattributes(val, {'double'}, {'numel', 2, 'nonnan', 'finite'}, ...
                           mfilename, '[eps_val mu_val]'); 
    end

    validateattributes(f, {'function_handle'}, {}, mfilename, 'f');

    % Optional arguments.
    options = my_parse_options(struct(  'upsample_ratio', 6, ...
                                        'f_avg', @(z) mean(z(:)),  ...
                                        'f_rep', @(val, ff, m) ff * val + ...
                                                                (1-ff) .* m), ...
                                varargin);
    validateattributes(options.upsample_ratio, {'numeric'}, ...
        {'positive', 'integer', 'scalar'}, mfilename, 'upsample_ratio');
    validateattributes(options.f_avg, {'function_handle'}, {}, ...
                        mfilename, 'f_avg');
    validateattributes(options.f_rep, {'function_handle'}, {}, ...
                        mfilename, 'f_rep');

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

    % Test if we can give multiple points to f.
    try 
        out = f([0 1], [2 3], [4 5]);
        multipoint = true;
    catch
        multipoint = false;
    end


        %
        % Prepare necessary grid information.
        %

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
            mu_grid_pos{k}{l} = pos_prim{l};
        end
        eps_grid_pos{k}{k} = pos_prim{k};
        mu_grid_pos{k}{k} = pos_dual{k};
    end


        %
        % Update material fields.
        %

    % Update components of epsilon.
    params = {bnd_box, f, options.upsample_ratio, multipoint, ...
                options.f_avg, options.f_rep};
    for k = 1 : 3
        eps{k} = my_update(eps_grid_pos{k}, eps{k}, val(1), params{:}); 
        if ~isempty(mu)
            mu{k} = my_update(mu_grid_pos{k}, mu{k}, mu_val, params{:});
        end
    end


function [mat] = my_update(pos, mat, val, box, f, up_ratio, multipoint, ...
                            f_avg, f_rep)
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

    % Downsample results by averaging.
    for i = 1 : (s{2}(1) - s{1}(1))
        for j = 1 : (s{2}(2) - s{1}(2))
            for k = 1 : (s{2}(3) - s{1}(3))
                fill_fraction(i, j, k) = f_avg(inside_shape(...
                                    (i-1)*up_ratio+[1:up_ratio], ...
                                    (j-1)*up_ratio+[1:up_ratio], ...
                                    (k-1)*up_ratio+[1:up_ratio]));
            end
        end
    end

    % Make changes to epsilon.
    m = mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1);
    mat(s{1}(1):s{2}(1)-1, s{1}(2):s{2}(2)-1, s{1}(3):s{2}(3)-1) = ...
        reshape(f_rep(val, fill_fraction(:), m(:)), size(m));
