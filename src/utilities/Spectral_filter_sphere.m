function psi_filt = Spectral_filter_sphere(psi, geom, p)
%SPECTRAL_FILTER_SPHERE Apply an exponential spherical-harmonic filter.
%
%   psi_filt = SPECTRAL_FILTER_SPHERE(psi, geom, p) transforms psi to
%   spherical-harmonic space, multiplies degree ell by
%   exp[-(ell/(L-1))^p], and transforms the filtered field back to the
%   SSHT sampling grid.
%
%   Inputs
%   ------
%   psi
%       Scalar field sampled on the SSHT grid.
%   geom
%       Geometry structure returned by Build_geometry_sphere.
%   p
%       Positive filter exponent.
%
%   Output
%   ------
%   psi_filt
%       Filtered field on the original grid.

    reality_flag = isreal(psi);

    flm = ssht_forward(psi, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality_flag);

    sigma = zeros(size(flm));

    ind = 1;
    for ell = 0:(geom.L_band-1)
        filt = exp(-(ell/(geom.L_band-1))^p);
        for m = -ell:ell
            sigma(ind) = filt;
            ind = ind + 1;
        end
    end

    flm = sigma .* flm;

    psi_filt = ssht_inverse(flm, geom.L_band, ...
        'Method', geom.method, ...
        'Reality', reality_flag);
end