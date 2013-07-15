function [E1] = my_propagate_beam(omega_eff, prop_dir, prop_dist, E0, pos)
 

        %
        % Upsample positions.
        %

    % Upsampling determined by smallest interval.
    min_delta = inf;
    for i = 1 : 3
        min_delta = min([min_delta; diff(pos{1}{i}(:))]);
    end

    % Determine the upsampled spacing.
    for i = 1 : 3
        n(i) = 2 * round((pos{1}{i}(end) - pos{1}{i}(1)) / min_delta);
        delta(i) = (pos{1}{i}(end) - pos{1}{i}(1)) / n(i);
    end

    % Obtain upsampled positions.
    for i = 1 : 3
        for j = 1 : 3
            if length(pos{i}{j}) > 1
                pos_up{i}{j} = pos{i}{j}(1) : delta(j) : pos{i}{j}(end);
            else
                pos_up{i}{j} = pos{i}{j};
            end
        end
        up_shape(i) = length(pos_up{1}{i}); % Upsampling shape.
    end

    % Obtain upsampled fields.
    E0_up = my_interp(pos, E0, pos_up);


        %
        % Determine plane wave modes to use.
        %

    for i = 1 : 3
        if n(i) == 0
            k_range{i} = 0;
        else
            k_range{i} = fftshift(2*pi / up_shape(i) * [0:up_shape(i)-1]);
        end
    end

    % Determine k's in propagation direction.
    k_range{prop_dir} = 0;
    [k{1}, k{2}, k{3}] = ndgrid(k_range{1}, k_range{2}, k_range{3});
    k{prop_dir} = sqrt(omega_eff^2 - (k{1}.^2 + k{2}.^2 + k{3}.^2 - k{prop_dir}.^2));
    k = [k{1}(:), k{2}(:), k{3}(:)];

%     % Eliminate k-vectors which are outside the light-cone
%     k(find(imag(k(:,prop_dir))), :) = [];

    % k_tot = sqrt(k(:,1).^2 + k(:,2).^2 + k(:,3).^2) % Debug.

    % Obtain the polarizations.
    p{1} = [-k(:,2), k(:,1), zeros(size(k, 1), 1)];
    p{2} = [-k(:,1).*conj(k(:,3)), -k(:,2).*conj(k(:,3)), (k(:,1).^2+k(:,2).^2)];
    p{3} = k;

    % Normalize the magnitudes (of the polarizations).
    for i = 1 : length(p) 
        p{i} = p{i} ./ (sqrt(sum(abs(p{i}).^2, 2)) * ones(1, 3));
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
%     a = [sum(conj(p{1}).*p{2}, 2), sum(conj(p{2}).*p{3}, 2), sum(conj(p{1}).*p{3}, 2)]
%     max(a(:))
%     return


        %
        % Convert to plane-wave basis.
        %

    % Convert to spatial frequency basis.
    for i = 1 : 3
        temp = fftshift(fftn(E0_up{i}));
        E0_f(:,i) = temp(:);
    end

    % Convert to plane wave modes.
    for j = 1 : 3
        for i = 1 : size(k, 1)
            mag{j}(i) =  conj(p{j}(i,:)) * E0_f(i,:).';
        end
    end


        %
        % Propagate some distance.
        %


        %
        % Convert back to real-space basis.
        %

    for i = 1 : length(mag{j})
        E1_f(i,:) = mag{1}(i) * p{1}(i,:) + mag{2}(i) * p{2}(i,:);
        E1_f(i,:) = mag{1}(i) * p{1}(i,:) + mag{2}(i) * p{2}(i,:) + mag{3}(i) * p{3}(i,:);
    end
    norm(E1_f(:) - E0_f(:))

    for i = 1 : 3
        E1_up{i} = ifft(ifftshift(reshape(E1_f(:,i), up_shape)));
    end


        %
        % Downsample.
        %

    E1 = my_interp(pos_up, E1_up, pos);
    plot(abs([E0{2}(:), E1{2}(:)]), '.-')
    pause




function [E_up] = my_interp(pos, E, pos_up)
    for i = 1 : 3
        up_shape(i) = length(pos_up{1}{i});
    end

    % Obtain upsampled field values.
    for i = 1 : 3
        [u0{1}, u0{2}, u0{3}] = ndgrid(pos{i}{1}, pos{i}{2}, pos{i}{3});
        [u1{1}, u1{2}, u1{3}] = ndgrid(pos_up{i}{1}, pos_up{i}{2}, pos_up{i}{3});

        my_shape = [size(E{i}), ones(1, 3-ndims(E{i}))];
        u0(find(my_shape == 1)) = [];
        u1(find(my_shape == 1)) = [];

        switch length(u0)
            case 1
                interp_fun = @interp1;
            case 2
                interp_fun = @interp2;
            case 3
                interp_fun = @interp3;
        end

        E_up{i} = interp_fun(u0{:}, E{i}, u1{:});
        E_up{i} = reshape(E_up{i}, up_shape);
    end



