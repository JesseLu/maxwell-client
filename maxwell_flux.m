%% maxwell_flux
% Computes the Poynting vector through a finite plane of the simulation.
% TODO: Not working yet!

%%% Syntax
%
% * |P = maxwell_flux(grid, [E H], plane_pos, plane_size)|
%   returns the integral of the Poynting vector over a plane
%   in the simulation domain.
%
% * |P = maxwell_flux(grid, [E H], [E1 H1])|
%   uses |[E1 H1]|, which is assumed to be the output from 
%   the |maxwell_wgmode| function (or something equivalent),
%   in order to define the plane and also filter |[E H]|.
%   This is useful for computing the amount of power in
%   only a single mode.
%   Note that this may fail if large evanescent fields are present.

function [P] = maxwell_flux(grid, E_H, varargin)

        %
        % Validate and parse inputs.
        %

    
    my_validate_grid(grid, mfilename);

    [E, H] = my_split(E_H, grid.shape, {'E', 'H'}, mfilename);
    my_validate_field(E, grid.shape, 'E', mfilename);
    my_validate_field(H, grid.shape, 'H', mfilename);

    switch length(varargin)
        case 1
            E1_H1 = varargin{1};
            [E1, H1] = my_split(E1_H1, grid.shape, {'E1', 'H1'}, mfilename);
            my_validate_field(E1, grid.shape, 'E1', mfilename);
            my_validate_field(H1, grid.shape, 'H1', mfilename);

            [plane_pos, plane_size] = deal([], []);
        case 2
            plane_pos = varargin{1};
            plane_size = varargin{2};
            validateattributes(plane_pos, {'numeric'}, ...
                    {'nonnan', 'finite', 'numel', 3}, mfilename, 'plane_pos');
            validateattributes(plane_size, {'numeric'}, ...
                    {'nonnan', 'numel', 3}, mfilename, 'plane_size');

            [E1, H1] = deal([], []);
        otherwise
            error('Invalid number of input parameters.');
    end


        %
        % Filter, if needed.
        %

    [vec, unvec] = my_vec(grid.shape);

    function [F] = my_project(F1, F2)
    % Project z1 onto (normalized) z2.
        F = unvec(vec(F1) * dot(vec(F1), vec(F2)) / norm(vec(F2));
    end

    % Filter.
    if ~isempty(E1)
        E = my_project(E, E1);
        H = my_project(H, H1);
    end


        %
        % Cut out slice, if needed.
        %

    if ~isempty(plane_pos)
        [p0, p1, prop_dir, prop_pos] = my_find_plane(grid, plane_pos, plane_size);

        % Cut out the plane.
        for k = 1 : 3
            sp{k} = grid.s_prim{k}(p0(k):p1(k));
            sd{k} = grid.s_dual{k}(p0(k):p1(k));
            E{k} = E{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
            H{k} = H{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
        end
    end


    % TODO: Fix this, need to find prop_dir for filtered case.
    % Always need to ignore the field in the direction of propagation.
    a_dir = mod(prop_dir, 3) + 1;
    b_dir = mod(prop_dir+1, 3) + 1;
    Ea = E{a_dir};
    Ha = H{a_dir};
    Eb = E{b_dir};
    Hb = H{b_dir};


        %
        % Calculate the power.
        %

    P = dot(Ea(:), Ha(:)) + dot(Eb(:), Hb(:));
