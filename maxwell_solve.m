%% maxwell_solve
% Solves the electromagnetic wave equation.

%%% Syntax
% 
% * |[E H] = maxwell_solve(grid, eps, J)|
%   returns the E- and H-fields of the solution to the electromagnetic
%   wave equation.
%
% * |[E H] = maxwell_solve(grid, [eps mu], J)|
%   same as above except for |mu ~= 1|.
%%
% * |... = maxwell_solve(..., 'vis_progress', vis_opt)|
%   controls the progress visualization where |vis_opt| can be
%   |none|, |plot|, |text|, or |both|. Defaults to |plot|.
%
% * |... = maxwell_solve(..., 'E0', E0)|
%   allows one to control the initial value of E.
%   Defaults to |E0 = 0|.
%
% * |... = maxwell_solve(..., 'max_iters', n, 'err_thresh', err)|
%   sets the termination conditions for the solve.
%   Defaults to |n = 5e4| and |err = 1e-6|.
%
% 
function [E, H, err] = maxwell_solve(grid, epsilon, J, varargin)
    [cb, vis_progress] = maxwell_solve_async(grid, epsilon, J, varargin{:});
    text_output = strcmp(vis_progress, 'text') | strcmp(vis_progress, 'both');
    while ~cb()
        if text_output
            fprintf(repmat('\b', 1, 80));
        end
    end
    if text_output; fprintf(repmat('\b', 1, 80)); end
    [~, E, H, err] = cb();
    if text_output; fprintf('\n'); end
    
