function [E, H, err] = test1
    dims = [100 160 40];
    for k = 1 : 3
        s{k} = ones(dims(k), 1);
        epsilon{k} = ones(dims);
        J{k} = zeros(dims);
    end
    J{2}(50, 80, 20) = 1;

    grid = struct('omega', 0.08, 's_prim', {s}, 's_dual', {s});

    modify_javapath();
    % [cb] = maxwell_solve_async(grid, epsilon, J, 'none', 'both');
    [E, H, err] = maxwell_solve(grid, epsilon, J, 'vis_progress', 'both');



end
