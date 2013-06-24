%% maxwell_flux
% Computes the Poynting vector through a finite plane of the simulation.

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
