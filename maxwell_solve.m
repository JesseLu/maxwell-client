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
    
