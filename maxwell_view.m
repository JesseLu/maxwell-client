%% maxwell_view
% Slice view of a field and/or structure.

%%% Syntax
%
% * |maxwell_view(grid, eps, E, 'x', [nan 50 nan])|
%   visualizes the |y = 50| slice of the x-components of |eps| and |E|.
%
% * |maxwell_view(grid, eps, [], 'x', [nan 50 nan])|
%   visualizes only the structure.
%
% * |maxwell_view(grid, [], E, 'x', [nan 50 nan])|
%   visualizes only the field.
%
% * |maxwell_view(grid, mu, H, 'x', [nan 50 nan], 'grid_type', 'dual');
%   allows for |mu| and |H| to be visualized in their correct locations,
%   accounting for the proper shift in the Yee cell.
%   |'grid_type'| defaults to |'prim'|, which assumes we are visualizing
%   E or eps.
%
% * |maxwell_view(..., 'field_phase', phase)|
%   allows one to adjust the phase of the field, 
%   where |phase = 0| is the real component,
%   |phase = pi/2| is the imaginary component,
%   |phase = nan| denotes the absolute value of the field,
%   and |phase = inf| will play a movie.
%   Defaults to 0.
%
% * |maxwell_view(..., 'clims', clims)|
%   allows the user to determine the minimum and maximum values of 
%   the colorbar.
%
% * |maxwell_view(..., 'reverse_structure', true)|
%   allows the user to reverse the grayscale/shading used to draw the structure.


function maxwell_view(grid, mat, F, dir, slice_ind, varargin)

        %
        % Validate and parse inputs.
        %

    my_validate_grid(grid, mfilename);

    if isempty(mat) && isempty(F)
        error('Either eps (mu) or E (H) must be non-empty.');
    end

    if ~isempty(mat)
        my_validate_field(mat, grid.shape, 'eps (or mu)', mfilename);
    end

    if ~isempty(F)
        my_validate_field(F, grid.shape, 'E (or H)', mfilename);
    end

    validateattributes(dir, {'char'}, {'scalar'}, 'dir', mfilename);
    if ~any(dir == 'xyz')
        error('Input dir must be either ''x'', ''y'', or ''z''.');
    end

    validateattributes(slice_ind, {'numeric'}, {'numel', 3}, ...
                        'slice_ind', mfilename);
    if sum(isnan(slice_ind)) ~= 2
        error('Exactly two elements of slice_ind must be nan.');
    end

    % Optional arguments.
    options = my_parse_options(struct(  'grid_type', 'prim', ...
                                        'clims', [], ...
                                        'field_phase', 0, ...
                                        'reverse_structure', false), ...
                                varargin, mfilename);

    validateattributes(options.grid_type, {'char'}, {}, 'grid_type', mfilename);
    if ~strcmp(options.grid_type, 'prim') && ~strcmp(options.grid_type, 'dual')
        error('grid_type option must be either ''prim'' or ''dual''.');
    end

    validateattributes(options.field_phase, {'numeric'}, {'scalar'}, ...
                        'field_phase', mfilename);

    validateattributes(options.reverse_structure, {'logical'}, {'binary'}, ...
                        'reverse_structure', mfilename);

        %
        % Determine slice.
        %

    slice_comp = find(~isnan(slice_ind)); 
    field_comp = find(dir == 'xyz');

    [eps_grid_pos, mu_grid_pos] = my_s2pos(grid); % Get grid positions.
    
    % Determine the position grid for the specific field and field component
    % that we wish to visualize.
    if strcmp(options.grid_type, 'prim')
        pos = eps_grid_pos{field_comp};
    else
        pos = mu_grid_pos{field_comp};
    end

    % Find index which is closest to the slice position.
    [~, ind] = min(abs(pos{slice_comp} - slice_ind(slice_comp)));
    if ind == length(pos{slice_comp}) % Special case related to extra pos.
        ind = 1
    end

    % Extract out the relevant slice.
    if ~isempty(mat)
        m_data = squeeze(real(my_slice(mat{field_comp}, slice_comp, ind)));
    end

    if ~isempty(F)
        F_data = squeeze(my_slice(F{field_comp}, slice_comp, ind));
    end

    % Determine grids which will be used for plotting the slice.
    perp_comp = find(isnan(slice_ind));
    xyz = 'xyz';
    xlabel = xyz(perp_comp(1));
    ylabel = xyz(perp_comp(2));
    x = pos{perp_comp(1)}(1:end-1); % Clip off extra position at end.
    y = pos{perp_comp(2)}(1:end-1);
    

        %
        % Plot.
        %

    if ~isempty(F) % Just field.

        % Prepare data.
        if isnan(options.field_phase)
            data = abs(F_data);
            cmax = max(data(:));
        elseif isinf(options.field_phase)
            num_frames = 24;
            for k = 1 : num_frames
                data{k} = real(exp(i*2*pi*(k-1)/num_frames) .* F_data);
            end
            cmax = max(abs(F_data(:)))/2;
        else
            data = real(exp(i*options.field_phase) .* F_data);
            cmax = max(abs(data(:)));
        end

        % Prepare alpha data.
        if isempty(mat) 
            alpha_data = [];
        else
            alpha_data = m_data - min(m_data(:));
            if options.reverse_structure
                alpha_data = 1 - alpha_data;
            end
            alpha_data = -alpha_data ./ max(alpha_data(:));
        end
        
        % Allow for user-override of cmax.
        if ~isempty(options.clims)
            clims = options.clims;
        else 
            clims = cmax * [-1 1];
        end

        % Plot.
        if ~isinf(options.field_phase)
            my_plot(x, y, data, alpha_data, {xlabel, ylabel}, clims);
            colormap('jet');
        else
            cnt = 0;
            while true
                frame_start = tic;
                my_plot(x, y, data{mod(cnt, num_frames)+1}, alpha_data, ...
                        {xlabel, ylabel}, clims);
                colormap('jet');
                pause_time = toc(frame_start) - 2 / num_frames;
                if pause_time > 0
                    pause(pause_time);
                end
                cnt = cnt + 1;
                drawnow;
            end
        end

    elseif ~isempty(mat) && isempty(F) % Just material.
        my_plot(x, y, m_data, [], {xlabel, ylabel}, options.clims);
        cmap = colormap('gray');
        if ~options.reverse_structure
            colormap(cmap(end:-1:1,:)); % Reversed grayscale colormap.
        end
    else
        error('Could not figure out what to plot.'); % Should never get here.
    end

function [data] = my_slice(arr, comp, ind)
% Take a slice out of arr.
    switch comp
        case 1
            data = arr(ind,:,:);
        case 2
            data = arr(:,ind,:);
        case 3
            data = arr(:,:,ind);
        otherwise
            error('comp needs to be either 1, 2, or 3.');
    end


function my_plot(x, y, data, alpha_data, labels, clims)
% Plot that allows for variable grid spacing.

    pcolor(x, y, data.');

    % Add transparency effect if needed.
    if ~isempty(alpha_data) && ~any(isnan(alpha_data(:)))
        alpha(alpha_data(2:end, 2:end).'); % 2:end is because of pcolor function.
        set(gca, 'Color', 'black', 'Alim', [-4 0]);
    end

   shading('flat');
   xlabel(labels{1});
   ylabel(labels{2});
   colorbar();
   axis equal tight;
   set(gca, 'YDir', 'normal');

   if ~isempty(clims)
       caxis(clims);
   end
   drawnow;


