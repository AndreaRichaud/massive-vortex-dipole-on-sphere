function [Mat_V_a, Mat_V_b] = Build_pinning_potential_sphere( ...
    geom, ...
    theta_1, phi_1, theta_2, phi_2, ...
    Pinning_flag, V0_pin_a, sigma_pin)

%BUILD_PINNING_POTENTIAL_SPHERE
% Builds a pinning potential on the sphere for component a only.
%
% INPUTS
%   geom         : geometry struct from Build_geometry_sphere
%   theta_1,phi_1: angular position of vortex 1
%   theta_2,phi_2: angular position of vortex 2
%   Pinning_flag : 0 or 1
%   V0_pin_a     : pinning strength for component a [J]
%   sigma_pin    : angular width of each Gaussian [rad]
%
% OUTPUTS
%   Mat_V_a      : pinning potential acting on component a [J]
%   Mat_V_b      : potential acting on component b [J]
%
% NOTES
%   The Gaussian is defined using the geodesic angular distance gamma on the
%   sphere:
%
%       V ~ exp( - gamma^2 / (2 sigma_pin^2) )
%
%   This is preferable to using naive Euclidean distances in (theta,phi).

validateattributes(theta_1, {'numeric'}, {'scalar','real','>=',0,'<=',pi});
validateattributes(theta_2, {'numeric'}, {'scalar','real','>=',0,'<=',pi});
validateattributes(phi_1,   {'numeric'}, {'scalar','real'});
validateattributes(phi_2,   {'numeric'}, {'scalar','real'});
validateattributes(Pinning_flag, {'numeric','logical'}, {'scalar'});
validateattributes(V0_pin_a, {'numeric'}, {'scalar','real'});
validateattributes(sigma_pin, {'numeric'}, {'scalar','real','positive'});

Mat_V_a = zeros(size(geom.theta));
Mat_V_b = zeros(size(geom.theta));

if Pinning_flag == 0
    return
end

% Wrap phi in [0,2*pi) only for clarity
phi_1 = mod(phi_1, 2*pi);
phi_2 = mod(phi_2, 2*pi);

TH = geom.theta;
PH = geom.phi;

% Geodesic angular distance from vortex 1
cos_gamma_1 = sin(TH).*sin(theta_1).*cos(PH - phi_1) + cos(TH).*cos(theta_1);
cos_gamma_1 = min(max(cos_gamma_1, -1), 1);
gamma_1 = acos(cos_gamma_1);

% Geodesic angular distance from vortex 2
cos_gamma_2 = sin(TH).*sin(theta_2).*cos(PH - phi_2) + cos(TH).*cos(theta_2);
cos_gamma_2 = min(max(cos_gamma_2, -1), 1);
gamma_2 = acos(cos_gamma_2);

% Repulsive pinning barriers for component a
Mat_V_a = V0_pin_a * exp( -gamma_1.^2 / (2*sigma_pin^2) ) ...
    + V0_pin_a * exp( -gamma_2.^2 / (2*sigma_pin^2) );

% Component b does not feel the pinning
Mat_V_b = zeros(size(TH));

% figure
%     subplot(1,2,1)
% ssht_plot_sphere(Mat_V_a, geom.L_band, ...
%     'Method', geom.method, ...
%     'Type', 'colour', ...
%     'ColourBar', true, ...
%     'Lighting', false);
% title('Pinning potential V_a on the sphere');
% subplot(1,2,2)
% ssht_plot_sphere(Mat_V_b, geom.L_band, ...
%     'Method', geom.method, ...
%     'Type', 'colour', ...
%     'ColourBar', true, ...
%     'Lighting', false);
% title('Pinning potential V_b on the sphere');


end


