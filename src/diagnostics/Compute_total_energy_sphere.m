function [Ene] = Compute_total_energy_sphere( ...
    psi_a, psi_b, ...
    m_a, m_b, ...
    g_a, g_b, g_ab, ...
    N_a, N_b, ...
    L_z, ...
    Mat_V_a, Mat_V_b, ...
    geom, hbar)

%COMPUTE_TOTAL_ENERGY_SPHERE
% Compute the total energy of two coupled GPEs on a sphere.
%
% INPUTS
%   psi_a, psi_b : wavefunctions normalized to 1 on the sphere
%   m_a, m_b     : masses
%   g_a, g_b     : intra-species 3D couplings
%   g_ab         : inter-species 3D coupling
%   N_a, N_b     : atom numbers
%   L_z          : effective radial thickness used for 3D -> 2D reduction
%   Mat_V_a      : external potential acting on component a
%   Mat_V_b      : external potential acting on component b
%   geom         : spherical geometry struct
%   hbar         : reduced Planck constant
%
% OUTPUT
%   Ene         : total energy [J]

    % Effective 2D couplings
    g_a_2D  = g_a  / L_z;
    g_b_2D  = g_b  / L_z;
    g_ab_2D = g_ab / L_z;

    % Laplace-Beltrami terms
    Laplacian_a = spherical_laplacian(psi_a, geom, 'Reality', isreal(psi_a));
    Laplacian_b = spherical_laplacian(psi_b, geom, 'Reality', isreal(psi_b));

    % Kinetic energy
    Kinetic_a = -N_a * hbar^2/(2*m_a) * sum(conj(psi_a(:)) .* Laplacian_a(:) .* geom.W(:));
    Kinetic_b = -N_b * hbar^2/(2*m_b) * sum(conj(psi_b(:)) .* Laplacian_b(:) .* geom.W(:));

    % Potential energy
    Pot_a = N_a * sum(Mat_V_a(:) .* abs(psi_a(:)).^2 .* geom.W(:));
    Pot_b = N_b * sum(Mat_V_b(:) .* abs(psi_b(:)).^2 .* geom.W(:));

    % Intra-species interaction energy
    Intra_a = g_a_2D * N_a^2 / 2 * sum(abs(psi_a(:)).^4 .* geom.W(:));
    Intra_b = g_b_2D * N_b^2 / 2 * sum(abs(psi_b(:)).^4 .* geom.W(:));

    % Inter-species interaction energy
    Inter_ab = g_ab_2D * N_a * N_b * sum(abs(psi_a(:)).^2 .* abs(psi_b(:)).^2 .* geom.W(:));

    Ene = real(Kinetic_a + Kinetic_b + Pot_a + Pot_b + Intra_a + Intra_b + Inter_ab);

end