function [Hpsi_a, Hpsi_b] = Imaginary_time_iteration_sphere( ...
    psi_a, psi_b, ...
    m_a, m_b, ...
    g_a, g_b, g_ab, ...
    N_a, N_b, ...
    L_z, ...
    Mat_V_a, Mat_V_b, ...
    geom, hbar)

%IMAGINARY_TIME_ITERATION_SPHERE
% Compute H psi for two coupled GPEs on a sphere.
%
% INPUTS
%   psi_a, psi_b : wavefunctions of components a and b
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
% OUTPUTS
%   Hpsi_a, Hpsi_b : Hamiltonian applied to psi_a and psi_b

    % Effective 2D couplings
    g_a_2D  = g_a  / L_z;
    g_b_2D  = g_b  / L_z;
    g_ab_2D = g_ab / L_z;

    % Laplace-Beltrami terms
    Laplacian_a = spherical_laplacian(psi_a, geom, 'Reality', isreal(psi_a));
    Laplacian_b = spherical_laplacian(psi_b, geom, 'Reality', isreal(psi_b));

    % Nonlinear intra-species terms
    Intra_a = g_a_2D * N_a * abs(psi_a).^2 .* psi_a;
    Intra_b = g_b_2D * N_b * abs(psi_b).^2 .* psi_b;

    % Nonlinear inter-species terms
    Inter_ab = g_ab_2D * N_b * abs(psi_b).^2 .* psi_a;
    Inter_ba = g_ab_2D * N_a * abs(psi_a).^2 .* psi_b;

    % External potentials
    Pot_a = Mat_V_a .* psi_a;
    Pot_b = Mat_V_b .* psi_b;

    % Hamiltonian applied to the wavefunctions
    Hpsi_a = -(hbar^2/(2*m_a)) * Laplacian_a + Intra_a + Inter_ab + Pot_a;
    Hpsi_b = -(hbar^2/(2*m_b)) * Laplacian_b + Intra_b + Inter_ba + Pot_b;

end