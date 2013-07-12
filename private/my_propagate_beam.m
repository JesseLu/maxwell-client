function [E] = my_propagate_beam(omega_eff, prop_dir, prop_dist, E, pos)
 

        %
        % Determine plane wave modes to use.
        %

    % Resolution is determined by extent in each spatial direction.
    for i = 1 : 3
        extent(i) = pos{1}{i}(end) - pos{1}{i}(1);
    end
    dk = 2*pi / max(extent);

    % Tabulate all plane wave modes.
    n = ceil(omega_eff/dk);
    for i = 1 : 3
        if i ~= prop_dir
            k_range{i} = dk * [-n:n];
        end
    end
    k_range{prop_dir} = 0;
    [k{1}, k{2}, k{3}] = ndgrid(k_range{1}, k_range{2}, k_range{3});
    k{prop_dir} = sqrt(omega_eff^2 - (k{1}.^2 + k{2}.^2 + k{3}.^2 - k{prop_dir}.^2));
    k = [k{1}(:), k{2}(:), k{3}(:)];

    % Eliminate those which are outside the light-cone
    k(find(imag(k(:,prop_dir))), :) = [];

    % k_tot = sqrt(k(:,1).^2 + k(:,2).^2 + k(:,3).^2) % Debug.

    % Obtain the polarizations.
    p{1} = [-k(:,2), k(:,1), zeros(size(k, 1), 1)];
    p{2} = [-k(:,1).*sin(k(:,3)), -k(:,2).*sin(k(:,3)), cos(k(:,3))];
    p{3} = k;

    % Normalize the magnitudes (of the polarizations).
    for i = 1 : length(p) 
        p{i} = p{i} ./ (sqrt(sum(p{i}.^2, 2)) * ones(1, 3));
        ind = find(isnan(p{i}(:,1)))
        if ~isempty(ind)
            if length(ind) > 1
                error('Too many nan values.');
            end
            p{i}(ind,:) = zeros(1, 3);
            p{i}(ind, i) = 1;
        end
    end



    % Check orthogonatlity.
    size(p{1})
    any(([sum(p{1}.*p{2}, 2); sum(p{2}.*p{3}, 2); sum(p{1}.*p{3}, 2)]) ~= 0)
    [sum(p{1}.*p{2}, 2); sum(p{2}.*p{3}, 2); sum(p{1}.*p{3}, 2)]'

        %
        % Convert to plane-wave basis.
        %


        %
        % Propagate some distance.
        %


        %
        % Convert back to real-space basis.
        %
