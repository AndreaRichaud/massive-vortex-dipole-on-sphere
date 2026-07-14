function [M_excl, rho_ref] = Compute_excluded_mass_sphere(psi_a, m_a, N_a, geom)
%COMPUTE_EXCLUDED_MASS_SPHERE
% Compute the excluded mass associated with vortex cores on a sphere.
%
% The excluded mass is defined as the surface integral of the positive part
% of (rho_ref - rho), multiplied by the atomic mass m_a:
%
%   M_excl = m_a * \int max(rho_ref - rho(theta,phi), 0) dS
%
% where
%   rho(theta,phi) = N_a * |psi_a(theta,phi)|^2
%
% and rho_ref is chosen as the WEIGHTED mean of the CENTRAL 50% of the
% SORTED rho values. The central window is selected by number of elements,
% while the average inside that window is weighted with geom.W.
%
% INPUTS
%   psi_a   : wavefunction of component a, normalized to 1 on the sphere
%   m_a     : atomic mass of component a
%   N_a     : atom number of component a
%   geom    : spherical geometry struct, with integration weights geom.W
%
% OUTPUTS
%   M_excl  : excluded mass [kg]
%   rho_ref : reference density [1/m^2]
%
% NOTE
%   Negative contributions rho_ref - rho < 0 are ignored.

    % Surface density of component a [1/m^2]
    rho = N_a * abs(psi_a).^2;

    % Flatten arrays
    rho_vec = rho(:);
    w_vec   = geom.W(:);

    % Keep only finite entries with positive weights
    valid = isfinite(rho_vec) & isfinite(w_vec) & (w_vec > 0);
    rho_vec = rho_vec(valid);
    w_vec   = w_vec(valid);

    if isempty(rho_vec)
        error('Compute_excluded_mass_sphere:InvalidInput', ...
              'No valid grid points found in rho or geom.W.');
    end

    % Reference density: weighted mean over the central 50% of sorted rho
    rho_ref = weighted_central_window_mean(rho_vec, w_vec, 0.5);

    % Positive density deficit only
    rho_deficit = max(rho_ref - rho, 0);

    % Excluded mass [kg]
    M_excl = m_a * sum(rho_deficit(:) .* geom.W(:));

    % Remove tiny imaginary parts from numerical noise
    M_excl = real(M_excl);
    rho_ref = real(rho_ref);

end


function x_mean = weighted_central_window_mean(x, w, frac_keep)
%WEIGHTED_CENTRAL_WINDOW_MEAN
% Sort x in ascending order, keep the central fraction frac_keep of the
% total number of elements, and compute the weighted mean over that subset.
%
% INPUTS
%   x         : values
%   w         : weights
%   frac_keep : fraction of elements to keep (e.g. 0.5 for central 50%)
%
% OUTPUT
%   x_mean    : weighted mean over the selected central window

    x = x(:);
    w = w(:);

    if numel(x) ~= numel(w)
        error('weighted_central_window_mean:SizeMismatch', ...
              'x and w must have the same number of elements.');
    end

    if isempty(x)
        error('weighted_central_window_mean:EmptyInput', ...
              'Input arrays must not be empty.');
    end

    if frac_keep <= 0 || frac_keep > 1
        error('weighted_central_window_mean:InvalidFraction', ...
              'frac_keep must be in the interval (0,1].');
    end

    % Sort values and corresponding weights
    [x_sorted, idx] = sort(x, 'ascend');
    w_sorted = w(idx);

    n = numel(x_sorted);
    n_keep = max(1, round(frac_keep * n));

    % Central window indices
    i_start = floor((n - n_keep)/2) + 1;
    i_end   = i_start + n_keep - 1;

    x_sel = x_sorted(i_start:i_end);
    w_sel = w_sorted(i_start:i_end);

    w_tot = sum(w_sel);
    if w_tot <= 0
        error('weighted_central_window_mean:ZeroTotalWeight', ...
              'Sum of weights in selected window must be positive.');
    end

    x_mean = sum(x_sel .* w_sel) / w_tot;
end