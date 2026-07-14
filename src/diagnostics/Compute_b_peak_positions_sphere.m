function [theta_N, phi_N, theta_S, phi_S] = Compute_b_peak_positions_sphere(psi_b, geom)
%COMPUTE_B_PEAK_POSITIONS_SPHERE
% Compute the center-of-mass position of the two disconnected peaks of psi_b,
% one in the Northern hemisphere and one in the Southern hemisphere.
%
% INPUTS
%   psi_b : wavefunction of component b
%   geom  : geometry struct, containing at least:
%           geom.theta, geom.phi, geom.W
%
% OUTPUTS
%   theta_N, phi_N : spherical coordinates of the Northern peak
%   theta_S, phi_S : spherical coordinates of the Southern peak
%
% NOTES
%   The positions are computed as weighted centers of mass on the unit sphere,
%   using weights rho = |psi_b|^2 * dA, separately on the two hemispheres.
%

    TH = geom.theta;
    PH = geom.phi;
    W  = geom.W;

    rho = abs(psi_b).^2;

    % Unit vectors on the sphere
    X = sin(TH) .* cos(PH);
    Y = sin(TH) .* sin(PH);
    Z = cos(TH);

    % Masks for the two hemispheres
    mask_N = (TH <= pi/2);
    mask_S = (TH >  pi/2);

    % ---------- Northern hemisphere ----------
    weight_N = rho .* W .* mask_N;
    mass_N = sum(weight_N(:));

    if mass_N > 0
        XN = sum(weight_N(:) .* X(:));
        YN = sum(weight_N(:) .* Y(:));
        ZN = sum(weight_N(:) .* Z(:));

        rN = sqrt(XN^2 + YN^2 + ZN^2);

        if rN > 0
            XN = XN / rN;
            YN = YN / rN;
            ZN = ZN / rN;

            theta_N = acos(max(min(ZN,1),-1));
            phi_N   = mod(atan2(YN, XN), 2*pi);
        else
            theta_N = NaN;
            phi_N   = NaN;
        end
    else
        theta_N = NaN;
        phi_N   = NaN;
    end

    % ---------- Southern hemisphere ----------
    weight_S = rho .* W .* mask_S;
    mass_S = sum(weight_S(:));

    if mass_S > 0
        XS = sum(weight_S(:) .* X(:));
        YS = sum(weight_S(:) .* Y(:));
        ZS = sum(weight_S(:) .* Z(:));

        rS = sqrt(XS^2 + YS^2 + ZS^2);

        if rS > 0
            XS = XS / rS;
            YS = YS / rS;
            ZS = ZS / rS;

            theta_S = acos(max(min(ZS,1),-1));
            phi_S   = mod(atan2(YS, XS), 2*pi);
        else
            theta_S = NaN;
            phi_S   = NaN;
        end
    else
        theta_S = NaN;
        phi_S   = NaN;
    end

end