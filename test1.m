function [E, H, err, grid, eps] = test1
%     dims = [100 160 40];
%     for k = 1 : 3
%         s{k} = ones(dims(k), 1);
%         epsilon{k} = ones(dims);
%         J{k} = zeros(dims);
%     end
%     J{2}(50, 80, 20) = 1;
% 
%     grid = struct('omega', 0.08, 's_prim', {s}, 's_dual', {s});
%     [grid, eps, ~, J] = maxwell_grid(0.3, -100:100, -100:100, -20:20);
%     J{2}(100,100,20) = 1;
% 
    [grid, eps, ~, J] = maxwell_grid(0.3, -100:100, -100:100, 0);
    J{2}(100,100,1) = 1;

    [E, H, err] = maxwell_solve(grid, eps, J, 'vis_progress', 'both');

    [A, x, b] = maxwell_axb(grid, eps, E, J);
    E_err = norm(A*x-b) / norm(b);
    [A, x, b] = maxwell_axb(grid, eps, [E H], J);
    EH_err = norm(A*x-b) / norm(b);
    fprintf('Error: %e (E), %e (EH)\n', E_err, EH_err);

end
