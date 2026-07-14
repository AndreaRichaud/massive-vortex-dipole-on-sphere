function [E_c, E_i, E_hd, E_qp, out] = Compute_Ec_Ei_sphere(psi, m, N, geom, hbar, varargin)
%COMPUTE_EC_EI_SPHERE
% Compressible / incompressible hydrodynamic kinetic energy on a sphere
% for a single condensate component.
%
% INPUTS
%   psi   : wavefunction normalized to 1 on the sphere
%   m     : mass
%   N     : atom number of this component
%   geom  : geometry struct, must contain at least:
%             geom.theta
%             geom.phi
%             geom.W
%             geom.R_sphere
%             geom.L_band
%             geom.method
%             geom.lap_eigs
%   hbar  : reduced Planck constant
%
% OPTIONAL NAME-VALUE PAIRS
%   'RhoFloorRel' : relative density floor used in u = j/sqrt(rho)
%                   default = 1e-10
%
% OUTPUTS
%   E_c   : compressible hydrodynamic kinetic energy [J]
%   E_i   : incompressible hydrodynamic kinetic energy [J]
%   E_hd  : direct hydrodynamic energy from u [J]
%   E_qp  : quantum pressure energy [J]
%   out   : diagnostics struct
%
% NOTES
%   Definitions:
%       rho = |psi|^2
%       j = (hbar/m) Im(psi^* grad_S psi)
%       u = sqrt(rho) v = j / sqrt(rho)
%
%   Hodge decomposition on the sphere:
%       u = grad_S Phi + n x grad_S Chi
%
%   with
%       Delta_S Phi = div_S u
%       Delta_S Chi = - curl_S(u)
%
%   Energies:
%       E_c  = (mN/2) ∫ |u_c|^2 dA
%       E_i  = (mN/2) ∫ |u_i|^2 dA
%       E_hd = (mN/2) ∫ |u|^2 dA
%       E_qp = N hbar^2/(2m) ∫ |grad_S sqrt(rho)|^2 dA
%

p = inputParser;
addParameter(p, 'RhoFloorRel', 1e-10, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, varargin{:});
rhoFloorRel = p.Results.RhoFloorRel;

theta = geom.theta;
W     = geom.W;
R     = geom.R_sphere;

if ~isequal(size(psi), size(theta), size(W))
    error('Compute_Ec_Ei_sphere:SizeMismatch', ...
        'psi, geom.theta, and geom.W must have the same size.');
end

sinT = sin(theta);
tolPole = 1e-12;
sinT_safe = sinT;
sinT_safe(abs(sinT_safe) < tolPole) = NaN;  % to detect poles cleanly

dtheta = infer_dtheta_from_geom(geom);

% ------------------------------------------------------------
% Density
% ------------------------------------------------------------
rho = abs(psi).^2;
rho_floor = rhoFloorRel * max(rho(:));
sqrt_rho_safe = sqrt(max(rho, rho_floor));

% ------------------------------------------------------------
% Derivatives of psi
%   d/dphi via L_z operator
%   d/dtheta via finite differences
% ------------------------------------------------------------
dpsi_dtheta = dtheta_fd(psi, dtheta);
dpsi_dphi   = dphi_spectral(psi, geom, hbar);

% Surface gradient components in orthonormal basis (e_theta, e_phi)
gradpsi_theta = dpsi_dtheta / R;
gradpsi_phi   = zeros(size(psi), 'like', psi);

mask = ~isnan(sinT_safe);
gradpsi_phi(mask) = dpsi_dphi(mask) ./ (R * sinT_safe(mask));
gradpsi_phi(~mask) = 0;

% ------------------------------------------------------------
% Current j = (hbar/m) Im(psi^* grad psi)
% ------------------------------------------------------------
j_theta = (hbar/m) * imag(conj(psi) .* gradpsi_theta);
j_phi   = (hbar/m) * imag(conj(psi) .* gradpsi_phi);

% ------------------------------------------------------------
% u = j / sqrt(rho)
% ------------------------------------------------------------
u_theta = j_theta ./ sqrt_rho_safe;
u_phi   = j_phi   ./ sqrt_rho_safe;

u_theta(~isfinite(u_theta)) = 0;
u_phi(~isfinite(u_phi))     = 0;

% ------------------------------------------------------------
% div_S u = 1/(R sinθ) [ dθ(sinθ uθ) + dφ(uφ) ]
% curl_S u = 1/(R sinθ) [ dθ(sinθ uφ) - dφ(uθ) ]
% ------------------------------------------------------------
A = sinT .* u_theta;
B = sinT .* u_phi;

dA_dtheta = dtheta_fd(A, dtheta);
dB_dtheta = dtheta_fd(B, dtheta);

duphi_dphi   = dphi_spectral(u_phi,   geom, hbar);
dutheta_dphi = dphi_spectral(u_theta, geom, hbar);

div_u = zeros(size(psi));
omega = zeros(size(psi));

div_u(mask) = (dA_dtheta(mask) + duphi_dphi(mask)) ./ (R * sinT_safe(mask));
omega(mask) = (dB_dtheta(mask) - dutheta_dphi(mask)) ./ (R * sinT_safe(mask));

% Regularize pole rows by copying nearest interior row
[div_u, omega] = regularize_poles(div_u, omega, sinT);

% Enforce real RHS for Poisson
div_u = real(div_u);
omega = real(omega);

% Remove numerical mean to ensure solvability
div_u = div_u - sum(div_u(:) .* W(:)) / sum(W(:));
omega = omega - sum(omega(:) .* W(:)) / sum(W(:));

% ------------------------------------------------------------
% Poisson solves on the sphere
% ------------------------------------------------------------
Phi = solve_poisson_sphere_spectral(div_u,  geom);
Chi = solve_poisson_sphere_spectral(-omega, geom);

% ------------------------------------------------------------
% Reconstruct u_c = grad_S Phi
% ------------------------------------------------------------
dPhi_dtheta = dtheta_fd(Phi, dtheta);
dPhi_dphi   = dphi_spectral(Phi, geom, hbar);

u_c_theta = dPhi_dtheta / R;
u_c_phi   = zeros(size(Phi));
u_c_phi(mask) = dPhi_dphi(mask) ./ (R * sinT_safe(mask));
u_c_phi(~mask) = 0;

% ------------------------------------------------------------
% Reconstruct u_i = n x grad_S Chi
%   (u_i)_theta = - (1/(R sinθ)) dφ Chi
%   (u_i)_phi   =   (1/R) dθ Chi
% ------------------------------------------------------------
dChi_dtheta = dtheta_fd(Chi, dtheta);
dChi_dphi   = dphi_spectral(Chi, geom, hbar);

u_i_theta = zeros(size(Chi));
u_i_theta(mask) = - dChi_dphi(mask) ./ (R * sinT_safe(mask));
u_i_theta(~mask) = 0;

u_i_phi = dChi_dtheta / R;

% ------------------------------------------------------------
% Quantum pressure
% ------------------------------------------------------------
sqrt_rho = sqrt(rho);
dsqrt_dtheta = dtheta_fd(sqrt_rho, dtheta);
dsqrt_dphi   = dphi_spectral(sqrt_rho, geom, hbar);

grad_s_theta = dsqrt_dtheta / R;
grad_s_phi   = zeros(size(sqrt_rho));
grad_s_phi(mask) = dsqrt_dphi(mask) ./ (R * sinT_safe(mask));
grad_s_phi(~mask) = 0;

% ------------------------------------------------------------
% Energies
% ------------------------------------------------------------
dens_c  = abs(u_c_theta).^2 + abs(u_c_phi).^2;
dens_i  = abs(u_i_theta).^2 + abs(u_i_phi).^2;
dens_hd = abs(u_theta).^2   + abs(u_phi).^2;
dens_qp = abs(grad_s_theta).^2 + abs(grad_s_phi).^2;

E_c  = 0.5 * m * N * sum(dens_c(:)  .* W(:));
E_i  = 0.5 * m * N * sum(dens_i(:)  .* W(:));
E_hd = 0.5 * m * N * sum(dens_hd(:) .* W(:));
E_qp = N * hbar^2/(2*m) * sum(dens_qp(:) .* W(:));

E_c  = real(E_c);
E_i  = real(E_i);
E_hd = real(E_hd);
E_qp = real(E_qp);

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
out = struct();
out.rho       = rho;
out.j_theta   = j_theta;
out.j_phi     = j_phi;
out.u_theta   = u_theta;
out.u_phi     = u_phi;
out.div_u     = div_u;
out.omega     = omega;
out.Phi       = Phi;
out.Chi       = Chi;
out.u_c_theta = u_c_theta;
out.u_c_phi   = u_c_phi;
out.u_i_theta = u_i_theta;
out.u_i_phi   = u_i_phi;

out.E_c       = E_c;
out.E_i       = E_i;
out.E_hd      = E_hd;
out.E_qp      = E_qp;
out.E_hd_rec  = E_c + E_i;
out.rel_rec_err = abs(out.E_hd - out.E_hd_rec) / abs(out.E_hd);
end


% =====================================================================
% Helpers
% =====================================================================

function dphi_f = dphi_spectral(f, geom, hbar)
% d/dphi f = (i/hbar) L_z f
dphi_f = 1i / hbar * L_z_operator_sphere(f, geom, hbar);
end

function df = dtheta_fd(f, dtheta)
% Centered finite differences in theta, one-sided at boundaries.
df = zeros(size(f), 'like', f);
df(2:end-1,:) = (f(3:end,:) - f(1:end-2,:)) / (2*dtheta);
df(1,:)       = (-3*f(1,:) + 4*f(2,:) - f(3,:)) / (2*dtheta);
df(end,:)     = ( 3*f(end,:) - 4*f(end-1,:) + f(end-2,:)) / (2*dtheta);
end

function dtheta = infer_dtheta_from_geom(geom)
th = geom.theta(:,1);
dth = diff(th);
dtheta = median(dth(:));
end

function f = solve_poisson_sphere_spectral(rhs, geom)
% Solve Delta_S f = rhs, setting the l=0 mode to zero.
rhs_lm = ssht_forward(rhs, geom.L_band, ...
    'Method', geom.method, ...
    'Reality', true);

lap = geom.lap_eigs(:);   % + l(l+1)/R^2
f_lm = zeros(size(rhs_lm));

mask = lap > 0;
f_lm(mask) = - rhs_lm(mask) ./ lap(mask);
f_lm(~mask) = 0;

f = ssht_inverse(f_lm, geom.L_band, ...
    'Method', geom.method, ...
    'Reality', true);
end

function [div_u, omega] = regularize_poles(div_u, omega, sinT)
% Copy nearest interior row into polar rows if present.
pole_rows = find(all(abs(sinT) < 1e-12, 2));
if isempty(pole_rows)
    return
end

nTh = size(div_u,1);

for r = pole_rows(:).'
    if r == 1 && nTh >= 2
        div_u(r,:) = div_u(2,:);
        omega(r,:) = omega(2,:);
    elseif r == nTh && nTh >= 2
        div_u(r,:) = div_u(end-1,:);
        omega(r,:) = omega(end-1,:);
    end
end
end
