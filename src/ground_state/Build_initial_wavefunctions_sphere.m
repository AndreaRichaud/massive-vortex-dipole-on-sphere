function [psi_a, psi_b] = Build_initial_wavefunctions_sphere( ...
    geom, ...
    theta_1, phi_1, q_1, ...
    theta_2, phi_2, q_2, ...
    sigma_psi_b, xi_core_a)

%BUILD_INITIAL_WAVEFUNCTIONS_SPHERE
% Initial states for a two-component GPE on a sphere.
%
% Component a:
%   seeded with a vortex dipole at (theta_1,phi_1) and (theta_2,phi_2)
%
% Component b:
%   seeded as a sum of two Gaussians localized at the two vortex cores
%
% INPUTS
%   geom        : geometry struct from Build_geometry_sphere
%   theta_1     : polar angle of vortex 1 [rad]
%   phi_1       : azimuthal angle of vortex 1 [rad]
%   q_1         : charge of vortex 1 (typically +1)
%   theta_2     : polar angle of vortex 2 [rad]
%   phi_2       : azimuthal angle of vortex 2 [rad]
%   q_2         : charge of vortex 2 (typically -1)
%   sigma_psi_b : linear width [m] of component-b seed around each core
%   xi_core_a   : linear core size [m] used to deplete density of psi_a
%
% OUTPUTS
%   psi_a       : initial state for component a, normalized to 1
%   psi_b       : initial state for component b, normalized to 1
%
% NOTES
%   This is a practical imaginary-time seed, not an exact stationary solution.
%   It is especially natural for a neutral dipole q_1 + q_2 = 0.
%

    % -----------------------------
    % Basic checks
    % -----------------------------
    validateattributes(theta_1, {'numeric'}, {'scalar','real','>=',0,'<=',pi});
    validateattributes(theta_2, {'numeric'}, {'scalar','real','>=',0,'<=',pi});
    validateattributes(phi_1,   {'numeric'}, {'scalar','real'});
    validateattributes(phi_2,   {'numeric'}, {'scalar','real'});
    validateattributes(q_1,     {'numeric'}, {'scalar','integer'});
    validateattributes(q_2,     {'numeric'}, {'scalar','integer'});
    validateattributes(sigma_psi_b, {'numeric'}, {'scalar','real','positive'});
    validateattributes(xi_core_a,   {'numeric'}, {'scalar','real','positive'});

    % This implementation is meant for a vortex dipole
    if ~(q_1 == +1 && q_2 == -1)
        error('Build_initial_wavefunctions_sphere:InvalidCharges', ...
            'This implementation assumes a dipole with q_1=+1 and q_2=-1.');
    end

    % Wrap azimuths
    phi_1 = mod(phi_1, 2*pi);
    phi_2 = mod(phi_2, 2*pi);

    TH = geom.theta;
    PH = geom.phi;
    R  = geom.R_sphere;

    % -----------------------------
    % Unit vector field on the sphere
    % -----------------------------
    rx = sin(TH).*cos(PH);
    ry = sin(TH).*sin(PH);
    rz = cos(TH);

    % Vortex-center unit vectors
    n1 = [sin(theta_1)*cos(phi_1), sin(theta_1)*sin(phi_1), cos(theta_1)];
    n2 = [sin(theta_2)*cos(phi_2), sin(theta_2)*sin(phi_2), cos(theta_2)];

    % -----------------------------
    % Angular geodesic distances from the two cores
    % -----------------------------
    cos_gamma_1 = rx*n1(1) + ry*n1(2) + rz*n1(3);
    cos_gamma_1 = min(max(cos_gamma_1, -1), 1);
    gamma_1 = acos(cos_gamma_1);

    cos_gamma_2 = rx*n2(1) + ry*n2(2) + rz*n2(3);
    cos_gamma_2 = min(max(cos_gamma_2, -1), 1);
    gamma_2 = acos(cos_gamma_2);

    % -----------------------------
    % Phase seed for component a
    % -----------------------------
    % Stereographic complex coordinate:
    %   zeta = tan(theta/2) * exp(i phi)
    %
    % Dipole phase:
    %   S = arg(zeta - zeta1) - arg(zeta - zeta2)
    %
    % where:
    %   zeta1 = tan(theta_1/2) * exp(i phi_1)
    %   zeta2 = tan(theta_2/2) * exp(i phi_2)

    zeta  = tan(TH/2)      .* exp(1i * PH);
    zeta1 = tan(theta_1/2) .* exp(1i * phi_1);
    zeta2 = tan(theta_2/2) .* exp(1i * phi_2);

    phase_a = angle(zeta - zeta1) - angle(zeta - zeta2);
    phase_a = angle(exp(1i * phase_a));   % wrap to (-pi, pi]

    % -----------------------------
    % Amplitude depletion for component a
    % -----------------------------
    % Convert linear core size to angular core size
    xi_ang = xi_core_a / R;

    amp_a = tanh(gamma_1 / xi_ang) .* tanh(gamma_2 / xi_ang);

    psi_a = amp_a .* exp(1i * phase_a);

    % -----------------------------
    % Initial seed for component b
    % -----------------------------
    % Sum of two geodesic Gaussians centered at the two vortex cores.
    sigma_b_ang = sigma_psi_b / R;

    psi_b = exp( -gamma_1.^2 / (2*sigma_b_ang^2) ) ...
          + exp( -gamma_2.^2 / (2*sigma_b_ang^2) );

    % -----------------------------
    % Normalize both wavefunctions to 1
    % -----------------------------
    psi_a = sphere_normalize(psi_a, geom.W);
    psi_b = sphere_normalize(psi_b, geom.W);

end