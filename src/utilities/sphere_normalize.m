function psi = sphere_normalize(psi, W)
%SPHERE_NORMALIZE Normalize a wavefunction on the spherical grid.
%
%   psi = SPHERE_NORMALIZE(psi, W) rescales psi so that
%   sum(abs(psi).^2 .* W) = 1, where W contains the surface-integration
%   weights associated with the SSHT grid.

    norm_psi = sqrt(sum(abs(psi(:)).^2 .* W(:)));
    psi = psi / norm_psi;

end