function [E_grid_pos, H_grid_pos] = my_s2pos(grid)
% Convert s-parameters to positions.
    origin_prim = grid.origin(:);
    origin_dual = grid.origin(:) + [real(grid.s_prim{1}(1))/2; ...
                                    real(grid.s_prim{2}(1))/2; ...
                                    real(grid.s_prim{3}(1))/2];
    for k = 1 : 3
        pos_prim{k} = origin_prim(k) + [0; cumsum(real(grid.s_prim{k}))];
        pos_dual{k} = origin_dual(k) + [0; cumsum(real(grid.s_dual{k}))];
    end

    % Build up grid info.
    for k = 1 : 3
        for l = 1 : 3
            E_grid_pos{k}{l} = pos_dual{l};
            H_grid_pos{k}{l} = pos_prim{l};
        end
        E_grid_pos{k}{k} = pos_prim{k};
        H_grid_pos{k}{k} = pos_dual{k};
    end



