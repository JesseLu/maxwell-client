function [E] = my_propagate_beam(omega, prop_dir, prop_dist, E, pos)

        %
        % Upsample pos (if not uniformly spaced).
        %
    
    % Get minimum grid spacing.
    min_delta = inf;
    for k = 1 : 3
        for l = 1 : 3
            min_delta = min([diff(pos{k}{l}(:)); min_delta]);
        end
    end
    
    % Get upsampled positions.
    for k = 1 : 3
        for l = 1 : 3
            n = ceil((pos{k}{l}(end) - pos{k}{l}(1)) / min_delta);
            p{k}{l} = pos{k}{l}(1) + min_delta * [0:n];
        end
    end

    % Get upsampled values.
    for k = 1 : 3
        for l = 1 : 3
            p0{l} = [pos{k}{l}(:); 1e9];
            E0_shape(l) = length(p0{l});
        end
        E0 = zeros(E0_shape);
        E0(1:end-1, 1:end-1, 1:end-1) = E{k};
        size(E0)

        [x0, y0, z0] = ndgrid(p0{1}, p0{2}, p0{3});
        [x1, y1, z1] = ndgrid(p{k}{1}, p{k}{2}, p{k}{3});
        e{k} = interp3(x0, y0, z0, E0, x1, y1, z1);
        % e{k} = interp3(p0{1}, p0{2}, p0{3}, E0, p{k}{1}, p{k}{2}, p{k}{3});
    end

        
        %
        % Convert to plane-wave basis.
        %


        %
        % Propagate some distance.
        %


        %
        % Convert back to real-space basis.
        %


        % 
        % Downsample back to original positions (if not uniformly spaced).
        %
