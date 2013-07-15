function [E1] = my_propagate_beam(omega_eff, prop_dir, prop_dist, E0, pos)
 

        %
        % Determine plane wave modes to use.
        %

    % Resolution is determined by extent in each spatial direction.
    for i = 1 : 3
        extent(i) = pos{1}{i}(end) - pos{1}{i}(1);
    end
    dk = pi / max(extent);

    % Tabulate all plane wave modes.
    n = ceil(omega_eff/dk);
    for i = 1 : 3
        if i ~= prop_dir
            k_range{i} = dk * [0:1];
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
    p{2} = [-k(:,1).*k(:,3), -k(:,2).*k(:,3), (k(:,1).^2+k(:,2).^2)];
    p{3} = k;

    % Normalize the magnitudes (of the polarizations).
    for i = 1 : length(p) 
        p{i} = p{i} ./ (sqrt(sum(p{i}.^2, 2)) * ones(1, 3));
        ind = find(isnan(p{i}(:,1)));
        if ~isempty(ind)
            if length(ind) > 1
                error('Too many nan values.');
            end
            p{i}(ind,:) = zeros(1, 3);
            p{i}(ind, i) = 1;
        end
    end

%     % Check orthogonatlity.
%     size(p{1})
%     any(([sum(p{1}.*p{2}, 2); sum(p{2}.*p{3}, 2); sum(p{1}.*p{3}, 2)]) ~= 0)
%     a = [sum(p{1}.*p{2}, 2), sum(p{2}.*p{3}, 2), sum(p{1}.*p{3}, 2)]
%     max(a(:))


        %
        % Convert to plane-wave basis.
        %


    for i = 1 : size(k, 1)   
        m = zeros(1, 3);
        for w = 1 : 3
            [x, y, z] = ndgrid(pos{w}{1}, pos{w}{2}, pos{w}{3});
            m(w) = m(w) + mean(E0{w}(:) .* exp( 1i * k(i,1) * x(:) + ...
                                                1i * k(i,2) * y(:) + ...
                                                1i * k(i,3) * z(:)));
        end

        for j = 1 : 2
            mag{j}(i) = sum(m .* p{j}(i,:));
        end
    end
    k
    [p{1}, mag{1}(:)]
    [p{2}, mag{2}(:)]



        %
        % Propagate some distance.
        %


        %
        % Convert back to real-space basis.
        %

    for w = 1 : 3
        E1{w} = 0 * E0{w};
        for i = 1 : size(k, 1)   
            [x, y, z] = ndgrid(pos{w}{1}, pos{w}{2}, pos{w}{3});
            for j = 1 : 2
                E1{w} = E1{w} + mag{j}(i) * p{j}(i,w) .* ...
                                    exp(-1i * k(i,1) * x + ...
                                        -1i * k(i,2) * y + ...
                                        -1i * k(i,3) * z);
            end
        end
        E1{w} = E1{w};
    end

    for w = 1 : 3
        norm(E0{w}(:) - E1{w}(:)) ./ norm(E0{w}(:))
    end
    plot(real([E0{2}(:), E1{2}(:)]), '.-');
    pause

