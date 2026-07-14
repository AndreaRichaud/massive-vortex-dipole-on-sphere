function [Lz_tot, Lz_a, Lz_b] = Total_angular_momentum_L_z_sphere( ...
    psi_a, psi_b, N_a, N_b, geom, hbar)
%TOTAL_ANGULAR_MOMENTUM_L_Z_SPHERE
% Compute the total angular momentum along z for a two-component condensate
% on the sphere.
%
% INPUTS
%   psi_a, psi_b : wavefunctions normalized to 1
%   N_a, N_b     : atom numbers
%   geom         : geometry struct
%   hbar         : reduced Planck constant
%
% OUTPUTS
%   Lz_tot : total angular momentum along z
%   Lz_a   : contribution from component a
%   Lz_b   : contribution from component b

    Lz_psi_a = L_z_operator_sphere(psi_a, geom, hbar);
    Lz_psi_b = L_z_operator_sphere(psi_b, geom, hbar);

    Lz_a = N_a * real(sum(conj(psi_a(:)) .* Lz_psi_a(:) .* geom.W(:)));
    Lz_b = N_b * real(sum(conj(psi_b(:)) .* Lz_psi_b(:) .* geom.W(:)));

    Lz_tot = Lz_a + Lz_b;
end