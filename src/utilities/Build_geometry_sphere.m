function geom = Build_geometry_sphere(R_sphere, L_band, method)

%BUILD_GEOMETRY_SPHERE Construct the SSHT sampling grid and spectral metadata.
%
%   geom = BUILD_GEOMETRY_SPHERE(R_sphere, L_band, method) builds the
%   spherical sampling grid, integration weights, and Laplace--Beltrami
%   eigenvalues required by the solver.
%
%   Inputs
%   ------
%   R_sphere
%       Radius of the spherical shell [m].
%   L_band
%       Spherical-harmonic band limit used by SSHT.
%   method
%       SSHT sampling theorem. Default: 'MW'.
%
%   Output
%   ------
%   geom
%       Structure containing angular grids, quadrature weights, spectral
%       indices, and geometry parameters.
%
%   External dependency: SSHT MATLAB library.

    if nargin < 3
        method = 'MW';
    end

    [thetas, phis, n, ntheta, nphi] = ssht_sampling(L_band, 'Method', method, 'Grid', true);

    dtheta = pi / ntheta;
    dphi   = 2*pi / nphi;
    W = R_sphere^2 * sin(thetas) * dtheta * dphi;

    n_lm = L_band^2;
    ell_vec  = zeros(n_lm,1);
    m_vec    = zeros(n_lm,1);
    lap_eigs = zeros(n_lm,1);

    ind = 1;
    for ell = 0:(L_band-1)
        for m = -ell:ell
            ell_vec(ind)  = ell;
            m_vec(ind)    = m;
            lap_eigs(ind) = ell*(ell+1)/R_sphere^2;
            ind = ind + 1;
        end
    end

    geom = struct();
    geom.R_sphere = R_sphere;
    geom.L_band   = L_band;
    geom.method   = method;

    geom.n       = n;
    geom.ntheta  = ntheta;
    geom.nphi    = nphi;

    geom.theta   = thetas;
    geom.phi     = phis;
    geom.W       = W;

    geom.n_lm     = n_lm;
    geom.ell_vec  = ell_vec;
    geom.m_vec    = m_vec;
    geom.lap_eigs = lap_eigs;
end
