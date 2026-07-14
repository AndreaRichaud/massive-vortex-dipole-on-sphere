function [Hpsi_a, Hpsi_b] = Real_time_iteration_sphere( ...
    psi_a, psi_b, ...
    m_a, m_b, ...
    g_a, g_b, g_ab, ...
    N_a, N_b, ...
    L_z, ...
    Mat_V_a, Mat_V_b, ...
    geom, hbar)

%REAL_TIME_ITERATION_SPHERE Apply the coupled Gross--Pitaevskii Hamiltonian.
%
%   [Hpsi_a,Hpsi_b] = REAL_TIME_ITERATION_SPHERE(...) evaluates the
%   Hamiltonian action for the two condensate components on the spherical
%   grid. The function does not advance time; it is used by the real-time
%   integration scheme.
%
%   Wavefunctions are normalized to unity; atom numbers enter explicitly
%   in the nonlinear terms. Couplings supplied to this function are 3D
%   couplings and are reduced using the effective radial thickness L_z.

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

    % Hamiltonian action
    Hpsi_a = -(hbar^2/(2*m_a)) * Laplacian_a + Intra_a + Inter_ab + Pot_a;
    Hpsi_b = -(hbar^2/(2*m_b)) * Laplacian_b + Intra_b + Inter_ba + Pot_b;

end
