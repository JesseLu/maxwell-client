%% maxwell_solve_eigenmode
% Find the eigenmode based on an initial guess for E.

%%% Syntax
%
% * |[omega, E, H] = maxwell_solve_eigenmode(grid, eps, E0)|
%   returns the frequency (|omega|), and fields (|E| and |H|)
%   of the eigenmode "nearest" to the field |E0|.
%   |maxwell_solve_eigenmode| works by calling 
%   an underlying |maxwell_solve| via the Rayleigh Quotient Iteration algorithm.
%
% * |... = maxwell_solve_eigenmode(grid, [eps mu], E0)|
%   does the same except for |mu ~= 1|.
%
% * |... = maxwell_solve_eigenmode(..., 'vis_progress', vis_opt)|
%   controls the progress visualization where |vis_opt| can be
%   |none|, |plot|, |text|, or |both|. Defaults to |plot|.
%
% * |... = maxwell_solve_eigenmode(..., 'eig_max_iters', eig_n, ...
%                                       'eig_err_thresh', eig_err)|
%   sets the termination conditions for the eigenmode algorithm 
%   (Rayleigh quotient iteration).
%   Defaults to |eig_n = 10|, and |eig_err = 1e-6|.
%
% * |... = maxwell_solve_eigenmode(..., 'max_iters', n, 'err_thresh', err)|
%   sets the termination conditions for the underlying calls to |maxwell_solve|,

function [omega, E, H] = maxwell_solve_eigenmode(grid, epsilon, E0, varargin) 


%% Form matrices and function handles 
% We now form the necessary linear algebra components and function hanles
% to solve the system using |eigenmode_solver|.
%
% We actually used a modified electromagnetic wave equation where $F = \sqrt{\epsilon} E$,
%
% $$ \frac{1}{\sqrt{\epsilon}}\nabla\times\mu^{-1}\nabla\times\frac{1}{\sqrt{\epsilon}}F - \omega^2 F = 0$$
%

    % Use the parse function to get select parameters.
    % Notice the E0 takes the place of J.
    [omega, ~, ~, ~, epsilon, ~, E0] = my_parse_inputs(grid, epsilon, E0, varargin{:});

    % Get ingredient matrices and vectors.
    [A, v] = maxwell_axb(grid, epsilon, E0, varargin{:}, 'E0', E0); 
    e = [epsilon{1}(:); epsilon{2}(:); epsilon{3}(:)];

    % Helper functions.
    dims = size(epsilon{1});
    n = prod(dims);
    unvec = @(z) {reshape(z(1:n), dims), reshape(z(n+1:2*n), dims), reshape(z(2*n+1:3*n), dims)};
    vec = @(z) [z{1}(:); z{2}(:); z{3}(:)]; 

    % Compose function handles.
    mult_A = @(x) e.^-0.5 .* (A * (e.^-0.5 .* x));
    mult_A_dag = @(x) (e.^-0.5 .* (A.' * (e.^-0.5 .* conj(x)))).';

    function [x] = solve_A_shifted(lambda, b)
    % Solves for the F-field.
        grid.omega = sqrt(lambda);
        J = unvec(-i * omega * b);
        E = maxwell_solve(grid, epsilon, J, varargin{:});
        x = sqrt(e) .* vec(E);
    end

    function my_vis(lambda, b)
        fprintf('omega: %e\n', sqrt(lambda));
    end

    [lambda, v] = my_solve_eigenmode(mult_A, @solve_A_shifted, @my_vis, v, 10, 1e-6);

    % Get back E and H.

end
