function [H] = my_E2H(grid, mu, E)
% Calculate H from E.

    % Inline function that takes (forward) derivatives of E.
    [sdx, sdy, sdz] = ndgrid(grid.s_dual{1}, grid.s_dual{2}, grid.s_dual{3});
    function [f] = d(f, dir)
        switch dir
            case 'x'
                f = (f([2:end, 1],:,:) - f) ./ sdx;
            case 'y'
                f = (f(:,[2:end, 1],:) - f) ./ sdy;
            case 'z'
                f = (f(:,:,[2:end, 1]) - f) ./ sdz;
            otherwise
                error('dir must be ''x'', ''y'', or ''z''.');
        end
    end

    % Calculate curl, assuming mu = 1.
    H{1} = (d(E{3}, 'y') - d(E{2}, 'z')) ./ (-1i * grid.omega);
    H{2} = (d(E{1}, 'z') - d(E{3}, 'x')) ./ (-1i * grid.omega);
    H{3} = (d(E{2}, 'x') - d(E{1}, 'y')) ./ (-1i * grid.omega);

    % If mu does not equal to 1, make a simple correction.
    if ~isempty(mu)
        for k = 1 : 3
            H{k} = H{k} ./ mu{k};
        end
    end

end
