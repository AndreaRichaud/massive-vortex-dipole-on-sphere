function [eigVals, M] = eig_linear_matrix(Mb, ma, na, R, hbar, theta)
% EIG_LINEAR_MATRIX
% Builds the 8x8 linearized matrix M and returns its eigenvalues.
%
% Inputs:
%   Mb    - mass of component b / effective vortex mass parameter
%   ma    - mass of component a atom
%   na    - 2D density of component a
%   R     - sphere radius
%   hbar  - reduced Planck constant
%   theta - polar angle (radians), with 0 < theta < pi/2
%
% Outputs:
%   eigVals - eigenvalues of M
%   M       - 8x8 matrix
%
% Example:
%   [eigVals, M] = eig_linear_matrix(1.0, 1.0, 0.1, 2.0, 1.0, pi/6);

% Basic input checks
if nargin ~= 6
    error('Function requires exactly 6 inputs: Mb, ma, na, R, hbar, theta.');
end

if ~isscalar(Mb) || ~isscalar(ma) || ~isscalar(na) || ...
        ~isscalar(R) || ~isscalar(hbar) || ~isscalar(theta)
    error('All inputs must be scalars.');
end

if Mb <= 0 || ma <= 0 || na <= 0 || R <= 0 || hbar <= 0
    error('Mb, ma, na, R, and hbar must be positive.');
end

if theta <= 0 || theta >= pi/2
    warning('theta is usually expected to satisfy 0 < theta < pi/2.');
end

M = [

0, 0, 0, 0, 2/(Mb*R^2), 0, 0, 0;

(2*sqrt(2*pi)*sqrt(na*(-Mb+2*ma*na*pi*R^2)/ma)*hbar/(Mb*R))*(1/sin(theta)), ...
0, 0, 0, 0, (2/(Mb*R^2))*(1/sin(theta)^2), 0, 0;

0, 0, 0, 0, 0, 0, 2/(Mb*R^2), 0;

0, 0, -(2*sqrt(2*pi)*sqrt(na*(-Mb+2*ma*na*pi*R^2)/ma)*hbar/(Mb*R))*(1/sin(theta)), ...
0, 0, 0, 0, (2/(Mb*R^2))*(1/sin(theta)^2);

-(na*pi*hbar^2*(1/cos(theta)^2)* ...
(-5*Mb + 12*ma*na*pi*R^2 - 2*(Mb-2*ma*na*pi*R^2)*cos(2*theta) ...
- 4*sqrt(2*pi)*R*sqrt(ma*na*(-Mb+2*ma*na*pi*R^2))*sin(theta)^2)) ...
/(2*ma*Mb), ...
0, -(na*pi*hbar^2*(1/cos(theta)^2))/(2*ma), 0, ...
0, -(2*sqrt(2*pi)*sqrt(na*(-Mb+2*ma*na*pi*R^2)/ma)*hbar/(Mb*R))*(1/sin(theta)), ...
0, 0;

0, -(na*pi*hbar^2*tan(theta)^2)/(2*ma), ...
0, (na*pi*hbar^2*tan(theta)^2)/(2*ma), ...
0, 0, 0, 0;

-(na*pi*hbar^2*(1/cos(theta)^2))/(2*ma), ...
0, ...
-(na*pi*hbar^2*(1/cos(theta)^2)* ...
(-5*Mb + 12*ma*na*pi*R^2 - 2*(Mb-2*ma*na*pi*R^2)*cos(2*theta) ...
- 4*sqrt(2*pi)*R*sqrt(ma*na*(-Mb+2*ma*na*pi*R^2))*sin(theta)^2)) ...
/(2*ma*Mb), ...
0, 0, 0, 0, ...
(2*sqrt(2*pi)*sqrt(na*(-Mb+2*ma*na*pi*R^2)/ma)*hbar/(Mb*R))*(1/sin(theta));

0, (na*pi*hbar^2*tan(theta)^2)/(2*ma), ...
0, -(na*pi*hbar^2*tan(theta)^2)/(2*ma), ...
0, 0, 0, 0

];

% Eigenvalues
eigVals = eig(M);
end