
function [cb] = maxwell_solve_async(grid, epsilon, J)

    [server_url, prefix] = maxwell_upload(grid, epsilon, J);
    cb = @() maxwell_download(server_url, prefix);
