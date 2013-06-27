function test2(grid, epsilon, E)
    maxwell_solve_eigenmode(grid, epsilon, E, 'eig_max_iters', 3, 'vis_progress', 'both');
