%% Analyze the initial-position sweep
% Load the saved real-time trajectories and diagnostics for different
% initial polar angles and generate the comparison plots used in the analysis.
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));
output_root = fullfile(repository_root, 'output');

base_output_folder = fullfile(output_root, 'Sweep_theta_initial');

% Original values used in the sweep
theta_1_values_all = linspace(0.1, pi/2 - 0.1, 10);

% Exclude the last one because that simulation failed
theta_0_values = theta_1_values_all(1:end-1);

nCases = numel(theta_0_values);

% Layout for 9 panels
nRows = 3;
nCols = 3;

%% =========================================================
%  1) Vortex trajectories on the sphere
% ==========================================================
fig1 = figure('Position', [100 100 1400 1200], 'Color', 'w');
tl1 = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(nCases,1);

for k = 1:nCases

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    ax(k) = nexttile;

    theta_pi_str = sprintf('%.3f', theta_0/pi);

    if ~exist(mat_file, 'file')
        axis(ax(k), 'off')
        title(ax(k), sprintf('$\\theta_0 = %s\\,\\pi$\\newline file not found', theta_pi_str), ...
            'Interpreter', 'latex', 'FontSize', 14)
        continue
    end

    S = load(mat_file, ...
        'vec_position_b_North', ...
        'vec_position_b_South');

    theta_N = S.vec_position_b_North(:,1);
    phi_N   = S.vec_position_b_North(:,2);

    theta_S = S.vec_position_b_South(:,1);
    phi_S   = S.vec_position_b_South(:,2);

    % Convert spherical -> Cartesian on unit sphere
    xN = sin(theta_N).*cos(phi_N);
    yN = sin(theta_N).*sin(phi_N);
    zN = cos(theta_N);

    xS = sin(theta_S).*cos(phi_S);
    yS = sin(theta_S).*sin(phi_S);
    zS = cos(theta_S);

    axes(ax(k)); %#ok<LAXES>

    % Reference sphere
    [XS,YS,ZS] = sphere(80);
    surf(XS,YS,ZS, ...
        'FaceAlpha', 0.08, ...
        'EdgeColor', 'none', ...
        'FaceColor', [0.8 0.8 0.8]);
    hold on

    % Trajectories
    plot3(xN, yN, zN, 'r', 'LineWidth', 1.8)
    plot3(xS, yS, zS, 'b', 'LineWidth', 1.8)

    % Initial points
    plot3(xN(1), yN(1), zN(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
    plot3(xS(1), yS(1), zS(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5)

    % Final points
    plot3(xN(end), yN(end), zN(end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 5)
    plot3(xS(end), yS(end), zS(end), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 5)

    axis equal
    axis vis3d
    axis off
    view(35,25)

    title(sprintf('$\\theta_0 = %s\\,\\pi$', theta_pi_str), ...
        'Interpreter', 'latex', 'FontSize', 15)

    hold off
end

sgtitle(tl1, 'Vortex trajectories on the sphere for different initial $\theta_0$', ...
    'Interpreter', 'latex', 'FontSize', 20)

% Link camera properties across all valid axes
valid_ax = ax(isgraphics(ax));

hlink = linkprop(valid_ax, ...
    {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});

% Store link handle so it is not destroyed
setappdata(fig1, 'CameraLink', hlink);

% Enable interactive rotation
rotate3d(fig1, 'on');


%% =========================================================
%  2) Final relative variations of E and Lz
% ==========================================================
fig2 = figure('Position', [100 100 1200 1000], 'Color', 'w');
tl2 = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:nCases

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    nexttile
    axis off

    theta_pi_str = sprintf('%.3f', theta_0/pi);

    if ~exist(mat_file, 'file')
        text(0.5, 0.5, sprintf('$\\theta_0 = %s\\,\\pi$\\newline file not found', theta_pi_str), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',14)
        continue
    end

    S = load(mat_file, 'vec_Ene', 'vec_L_z');

    vec_E = S.vec_Ene(:);
    vec_L = S.vec_L_z(:);

    if isempty(vec_E) || isempty(vec_L)
        text(0.5, 0.5, sprintf('$\\theta_0 = %s\\,\\pi$\\newline empty data', theta_pi_str), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',14)
        continue
    end

    E0 = vec_E(1);
    Ef = vec_E(end);

    L0 = vec_L(1);
    Lf = vec_L(end);

    dE_percent = 100 * (Ef - E0) / E0;
    dL_percent = 100 * (Lf - L0) / L0;

    str = sprintf([...
        '$\\theta_0 = %s\\,\\pi$ \n\n' ...
        '$\\Delta E = %.3f\\%%$ \n' ...
        '$\\Delta L_z = %.3f\\%%$'], ...
        theta_pi_str, dE_percent, dL_percent);

    text(0.5, 0.5, str, ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',16)
end

sgtitle(tl2, 'Final relative variations of $E$ and $L_z$', ...
    'Interpreter','latex', 'FontSize',20)


%% =========================================================
%  3) Plot E_PV(t) for all theta_0
% ==========================================================
fig3 = figure('Position', [100 100 1000 700], 'Color', 'w');
hold on
grid on
box on

cmap = lines(nCases);

for k = 1:nCases

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        fprintf('File not found for theta_0 = %.5f\n', theta_0);
        continue
    end

    S = load(mat_file, ...
        'vec_t', ...
        'vec_position_b_North', ...
        'vec_position_b_South', ...
        'm_a','m_b','N_a','N_b','hbar','geom');

    t = S.vec_t(:);

    thetaN = S.vec_position_b_North(:,1);
    phiN   = S.vec_position_b_North(:,2);

    thetaS = S.vec_position_b_South(:,1);
    phiS   = S.vec_position_b_South(:,2);

    % unwrap azimuthal angles
    phiN_u = unwrap(phiN);
    phiS_u = unwrap(phiS);

    % time derivatives
    dthetaN = gradient(thetaN, t);
    dphiN   = gradient(phiN_u, t);

    dthetaS = gradient(thetaS, t);
    dphiS   = gradient(phiS_u, t);

    % parameters
    R    = S.geom.R_sphere;
    Mb   = S.m_b * S.N_b;
    hbar = S.hbar;
    m_a  = S.m_a;
    N_a  = S.N_a;

    % kinetic energy
    T_N = (Mb/4) * R^2 .* (dthetaN.^2 + sin(thetaN).^2 .* dphiN.^2);
    T_S = (Mb/4) * R^2 .* (dthetaS.^2 + sin(thetaS).^2 .* dphiS.^2);
    T   = T_N + T_S;

    % chi (interaction term)
    DeltaPhi = phiN_u - phiS_u;

    arg = 2 ...
        - 2*cos(thetaN).*cos(thetaS) ...
        - 2*sin(thetaN).*sin(thetaS).*cos(DeltaPhi);

    arg = max(arg, eps);
    chi = 0.5 * log(arg);

    % dipolar energy
    na = N_a / (4*pi*R^2);
    prefE = 2 * hbar^2 * pi * na / m_a;

    xi = R / 50;
    E_const = prefE * log(R/xi);

    E_dip = E_const + prefE * chi;

    % total PV energy
    E_PV = T + E_dip;

    plot(t, E_PV, 'LineWidth', 1.7, ...
        'Color', cmap(k,:), ...
        'DisplayName', sprintf('$\\theta_0 = %.3f\\,\\pi$', theta_0/pi));
end

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$E_{\mathrm{PV}}(t)$ [J]', 'Interpreter', 'latex', 'FontSize', 16)
title('Point-vortex energy reconstructed from trajectories', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'best')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)


%% =========================================================
%  4) Angular distance between vortices vs time
% ==========================================================
fig4 = figure('Position', [100 100 1000 700], 'Color', 'w');
hold on
grid on
box on

cmap = lines(nCases);

for k = 1:nCases

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        fprintf('File not found for theta_0 = %.5f\n', theta_0);
        continue
    end

    S = load(mat_file, ...
        'vec_t', ...
        'vec_position_b_North', ...
        'vec_position_b_South');

    t = S.vec_t(:);

    thetaN = S.vec_position_b_North(:,1);
    phiN   = S.vec_position_b_North(:,2);

    thetaS = S.vec_position_b_South(:,1);
    phiS   = S.vec_position_b_South(:,2);

    cos_gamma = cos(thetaN).*cos(thetaS) + ...
        sin(thetaN).*sin(thetaS).*cos(phiN - phiS);

    % numerical safety
    cos_gamma = min(max(cos_gamma, -1), 1);

    gamma = acos(cos_gamma);

    plot(t, gamma, 'LineWidth', 1.7, ...
        'Color', cmap(k,:), ...
        'DisplayName', sprintf('$\\theta_0 = %.3f\\,\\pi$', theta_0/pi));
end

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\gamma(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 16)
title('Angular distance between the two vortices', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend('Interpreter', 'latex', 'FontSize', 11, 'Location', 'best')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)


%% =========================================================
%  5) All trajectories on a single spherical surface
%      + analytical prediction in simplified regime
% ==========================================================
fig5 = figure('Position', [100 100 950 800], 'Color', 'w');
ax5 = axes(fig5);
hold(ax5, 'on')

% Reference sphere
[XS,YS,ZS] = sphere(100);
hsurf = surf(ax5, XS, YS, ZS, ...
    'FaceAlpha', 0.07, ...
    'EdgeColor', 'none', ...
    'FaceColor', [0.8 0.8 0.8], ...
    'HandleVisibility', 'off');
% extra safety: keep sphere out of legend
set(get(get(hsurf,'Annotation'),'LegendInformation'), ...
    'IconDisplayStyle', 'off');

cmap = lines(nCases);

for k = 1:nCases

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        fprintf('File not found for theta_0 = %.5f\n', theta_0);
        continue
    end

    S = load(mat_file, ...
        'vec_t', ...
        'vec_position_b_North', ...
        'vec_position_b_South', ...
        'm_a', 'N_a', 'm_b', 'N_b', 'hbar', 'geom');

    t = S.vec_t(:);

    theta_N = S.vec_position_b_North(:,1);
    phi_N   = S.vec_position_b_North(:,2);

    theta_S = S.vec_position_b_South(:,1);
    phi_S   = S.vec_position_b_South(:,2);

    % Numerical trajectories
    xN = sin(theta_N).*cos(phi_N);
    yN = sin(theta_N).*sin(phi_N);
    zN = cos(theta_N);

    xS = sin(theta_S).*cos(phi_S);
    yS = sin(theta_S).*sin(phi_S);
    zS = cos(theta_S);

    thisColor = cmap(k,:);

    % lighter color for analytical curves
    analColor = 0.55*thisColor + 0.45*[1 1 1];

    % Numerical: same color for both trajectories of the same simulation
    plot3(ax5, xN, yN, zN, '-', 'LineWidth', 2.0, 'Color', thisColor, ...
        'DisplayName', sprintf('$\\theta_0 = %.3f\\,\\pi$', theta_0/pi));

    plot3(ax5, xS, yS, zS, '-', 'LineWidth', 2.0, 'Color', thisColor, ...
        'HandleVisibility', 'off');

    % Initial points
    plot3(ax5, xN(1), yN(1), zN(1), 'o', ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', thisColor, ...
        'MarkerSize', 4, 'HandleVisibility', 'off');

    plot3(ax5, xS(1), yS(1), zS(1), 'o', ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', thisColor, ...
        'MarkerSize', 4, 'HandleVisibility', 'off');

    % Final points
    plot3(ax5, xN(end), yN(end), zN(end), 's', ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', thisColor, ...
        'MarkerSize', 4, 'HandleVisibility', 'off');

    plot3(ax5, xS(end), yS(end), zS(end), 's', ...
        'MarkerFaceColor', thisColor, 'MarkerEdgeColor', thisColor, ...
        'MarkerSize', 4, 'HandleVisibility', 'off');

    %% ----- Analytical simplified model -----
    R    = S.geom.R_sphere;
    hbar = S.hbar;
    m_a  = S.m_a;
    N_a  = S.N_a;
    m_b  = S.m_b;
    N_b  = S.N_b;

    n_a = N_a / (4*pi*R^2);
    M_b = m_b * N_b;

    % Use initial polar angles and azimuths
    thetaN0 = theta_N(1);
    phiN0   = phi_N(1);

    thetaS0 = theta_S(1);
    phiS0   = phi_S(1);

    % Radicand of the analytical expression
    radicand = 2*pi*n_a * (2*m_a*n_a*pi*R^2 - M_b) / m_a;

    if radicand >= 0 && abs(cos(thetaN0)) > eps && abs(cos(thetaS0)) > eps

        pref = hbar/M_b/R;
        bracket = (2*n_a*pi*R - sqrt(radicand));

        % North analytic trajectory
        phidotN = pref * bracket / cos(thetaN0);

        phiN_an   = phiN0 + phidotN * t;
        thetaN_an = thetaN0 * ones(size(t));

        xN_an = sin(thetaN_an).*cos(phiN_an);
        yN_an = sin(thetaN_an).*sin(phiN_an);
        zN_an = cos(thetaN_an);

        % South analytic trajectory
        % sign opposite to the one given in the original expression
        phidotS = - pref * bracket / cos(thetaS0);

        phiS_an   = phiS0 + phidotS * t;
        thetaS_an = thetaS0 * ones(size(t));

        xS_an = sin(thetaS_an).*cos(phiS_an);
        yS_an = sin(thetaS_an).*sin(phiS_an);
        zS_an = cos(thetaS_an);

        % Plot analytical prediction: dashed and lighter
        plot3(ax5, xN_an, yN_an, zN_an, '--', ...
            'LineWidth', 1.2, 'Color', analColor, 'HandleVisibility', 'off');

        plot3(ax5, xS_an, yS_an, zS_an, '--', ...
            'LineWidth', 1.2, 'Color', analColor, 'HandleVisibility', 'off');

    else
        fprintf(['Analytical model not plotted for theta_0 = %.5f ' ...
                 '(radicand < 0 or cos(theta)=0).\n'], theta_0);
    end
end

axis(ax5, 'equal')
axis(ax5, 'vis3d')
axis(ax5, 'off')
view(ax5, 35, 25)
rotate3d(fig5, 'on')

title(ax5, 'All vortex trajectories on the sphere', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend(ax5, 'Interpreter', 'latex', 'FontSize', 11, 'Location', 'eastoutside')
set(ax5, 'FontSize', 14)

%% =========================================================
%  6) Sinusoidal fit of gamma(t) for each theta_0
%      + extraction of omega(theta_0)
%      + comparison with the analytical model
% ==========================================================

omega_fit_theta   = nan(size(theta_0_values));
A_fit_theta       = nan(size(theta_0_values));
gamma0_fit_theta  = nan(size(theta_0_values));
phi_fit_theta     = nan(size(theta_0_values));

fig6 = figure('Position', [100 100 1400 1200], 'Color', 'w');
tl6 = tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(theta_0_values)

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    nexttile
    hold on
    grid on
    box on

    if ~exist(mat_file, 'file')
        title(sprintf('$\\theta_0 = %.3f\\,\\pi$\\newline file not found', theta_0/pi), ...
            'Interpreter', 'latex', 'FontSize', 14)
        axis off
        continue
    end

    S = load(mat_file, ...
        'vec_t', ...
        'vec_position_b_North', ...
        'vec_position_b_South');

    t = S.vec_t(:);

    thetaN = S.vec_position_b_North(:,1);
    phiN   = S.vec_position_b_North(:,2);

    thetaS = S.vec_position_b_South(:,1);
    phiS   = S.vec_position_b_South(:,2);

    if isempty(t) || numel(t) < 5 || ...
       any(~isfinite(t)) || any(~isfinite(thetaN)) || any(~isfinite(phiN)) || ...
       any(~isfinite(thetaS)) || any(~isfinite(phiS))
        title(sprintf('$\\theta_0 = %.3f\\,\\pi$\\newline invalid data', theta_0/pi), ...
            'Interpreter', 'latex', 'FontSize', 14)
        axis off
        continue
    end

    % Angular distance gamma(t)
    cos_gamma = cos(thetaN).*cos(thetaS) + ...
                sin(thetaN).*sin(thetaS).*cos(phiN - phiS);
    cos_gamma = min(max(cos_gamma, -1), 1);
    gamma = acos(cos_gamma);

    % Remove repeated times
    [t, ia] = unique(t, 'stable');
    gamma = gamma(ia);
    thetaN = thetaN(ia);

    % Initial guesses
    gamma0_guess = mean(gamma);
    A_guess = 0.5 * (max(gamma) - min(gamma));

    % Frequency guess from FFT
    dt = mean(diff(t));
    y = gamma - mean(gamma);
    n = numel(y);

    if n < 8 || all(abs(y) < 1e-12)
        gamma0_fit_theta(k) = gamma0_guess;
        A_fit_theta(k) = 0;
        omega_fit_theta(k) = 0;
        phi_fit_theta(k) = 0;

        plot(t, gamma, 'k', 'LineWidth', 1.5)
        title(sprintf('$\\theta_0 = %.3f\\,\\pi$\\newline quasi-constant $\\gamma(t)$', theta_0/pi), ...
            'Interpreter', 'latex', 'FontSize', 13)
        xlabel('$t$ [s]', 'Interpreter', 'latex')
        ylabel('$\gamma$ [rad]', 'Interpreter', 'latex')
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11)
        hold off
        continue
    end

    Y = fft(y);
    freqs = (0:n-1)'/(n*dt);   % Hz
    half_idx = 2:floor(n/2);   % skip zero-frequency component

    [~, idx_max] = max(abs(Y(half_idx)));
    f_guess = freqs(half_idx(idx_max));
    omega_guess = 2*pi*f_guess;

    if ~isfinite(omega_guess) || omega_guess <= 0
        T_total = t(end) - t(1);
        omega_guess = 2*pi / max(T_total, eps);
    end

    phi_guess = 0;

    % Model: p = [gamma0, A, omega, phi]
    model = @(p,tt) p(1) + p(2)*sin(p(3)*tt + p(4));

    % Objective function
    objfun = @(p) sum((gamma - model(p,t)).^2);

    % Initial parameter vector
    p0 = [gamma0_guess, A_guess, omega_guess, phi_guess];

    % Fit with fminsearch
    opts = optimset('Display','off', 'MaxFunEvals', 1e5, 'MaxIter', 1e5);
    pfit = fminsearch(objfun, p0, opts);

    % Standardize sign: A >= 0
    if pfit(2) < 0
        pfit(2) = -pfit(2);
        pfit(4) = pfit(4) + pi;
    end

    % Standardize frequency: omega >= 0
    if pfit(3) < 0
        pfit(3) = -pfit(3);
        pfit(4) = -pfit(4);
    end

    gamma0_fit_theta(k) = pfit(1);
    A_fit_theta(k)      = pfit(2);
    omega_fit_theta(k)  = pfit(3);
    phi_fit_theta(k)    = pfit(4);

    % Plot data + fit
    tt_fine = linspace(t(1), t(end), 1000);
    gamma_fit_curve = model(pfit, tt_fine);

    plot(t, gamma, 'ko', 'MarkerSize', 3, 'DisplayName', 'data')
    plot(tt_fine, gamma_fit_curve, 'r-', 'LineWidth', 1.8, 'DisplayName', 'fit')

    title(sprintf('$\\theta_0 = %.3f\\,\\pi$\\newline $\\omega = %.4g\\ \\mathrm{rad/s}$', ...
        theta_0/pi, omega_fit_theta(k)), ...
        'Interpreter', 'latex', 'FontSize', 13)

    xlabel('$t$ [s]', 'Interpreter', 'latex')
    ylabel('$\gamma$ [rad]', 'Interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11)

    hold off
end

sgtitle(tl6, 'Sinusoidal fits of $\gamma(t)$ for different initial $\theta_0$', ...
    'Interpreter', 'latex', 'FontSize', 20);


%% =========================================================
%  7) omega(theta_0): numerical fit vs analytical prediction
% ==========================================================

% Use parameters from the first valid file
sampleFound = false;
for k = 1:nCases
    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if exist(mat_file, 'file')
        S0 = load(mat_file, 'm_a', 'm_b', 'N_a', 'N_b', 'hbar', 'geom');
        sampleFound = true;
        break
    end
end

if ~sampleFound
    error('No valid .mat file found in Sweep_theta_initial.');
end

m_a  = S0.m_a;
m_b  = S0.m_b;
N_a  = S0.N_a;
N_b  = S0.N_b;
hbar = S0.hbar;
R    = S0.geom.R_sphere;

na = N_a / (4*pi*R^2);
Mb = m_b * N_b;
tolEig = 1e-8;

omega_analytic_theta = nan(size(theta_0_values));
modeIndex_best_theta = nan(size(theta_0_values));
theta_eq_values      = nan(size(theta_0_values));

omega_candidates_all_theta = cell(size(theta_0_values));

for k = 1:numel(theta_0_values)

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        continue
    end

    S = load(mat_file, ...
        'vec_position_b_North', ...
        'vec_t');

    thetaN = S.vec_position_b_North(:,1);
    t = S.vec_t(:);

    % remove repeated times consistently
    [~, ia] = unique(t, 'stable');
    thetaN = thetaN(ia);

    % better estimate of equilibrium angle
    theta_eq_values(k) = mean(thetaN);

    try
        [eigVals, ~] = eig_linear_matrix(Mb, m_a, na, R, hbar, theta_eq_values(k));

        omega_candidates = imag(eigVals);
        omega_candidates = omega_candidates(omega_candidates > tolEig);
        omega_candidates = sort(omega_candidates, 'ascend');

        if ~isempty(omega_candidates)
            omega_candidates = unique(round(omega_candidates, 12));
        end

        omega_candidates_all_theta{k} = omega_candidates;

        if isfinite(omega_fit_theta(k)) && ~isempty(omega_candidates)
            [~, idxBest] = min(abs(omega_candidates - omega_fit_theta(k)));
            omega_analytic_theta(k) = omega_candidates(idxBest);
            modeIndex_best_theta(k) = idxBest;
        end

    catch
        omega_analytic_theta(k) = NaN;
        modeIndex_best_theta(k) = NaN;
        omega_candidates_all_theta{k} = [];
    end
end

%% =========================================================
%  7) Extended linearized-theory analysis vs theta_0
%      - compute M_ex(theta_0)
%      - build ordered spectra Im(lambda)
%      - MANUAL branch selection
%      - plot omega(theta_0) with double analytical prediction
%      - plot M_ex/M_a vs theta_0
%      - plot all Im(lambda) branches vs theta_0
%% ==========================================================

% Use parameters from the first valid file
sampleFound = false;
for k = 1:nCases
    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if exist(mat_file, 'file')
        S0 = load(mat_file, 'm_a', 'm_b', 'N_a', 'N_b', 'hbar', 'geom');
        sampleFound = true;
        break
    end
end

if ~sampleFound
    error('No valid .mat file found in Sweep_theta_initial.');
end

m_a  = S0.m_a;
m_b  = S0.m_b;
N_a  = S0.N_a;
N_b  = S0.N_b;
hbar = S0.hbar;
R    = S0.geom.R_sphere;

tolEig = 1e-8;

% ---------------------------------------------------------
% Containers
% ---------------------------------------------------------
M_ex_theta      = nan(size(theta_0_values));
M_a_theta       = nan(size(theta_0_values));
theta_eq_values = nan(size(theta_0_values));

imagEig_cell_theta_Mb         = cell(size(theta_0_values));
imagEig_cell_theta_Mb_plusMex = cell(size(theta_0_values));

omega_analytic_branch_theta_Mb         = nan(size(theta_0_values));
omega_analytic_branch_theta_Mb_plusMex = nan(size(theta_0_values));

% ---------------------------------------------------------
% Loop over theta_0:
%   - compute M_ex(theta_0)
%   - build full ordered spectra Im(lambda)
% ---------------------------------------------------------
for k = 1:numel(theta_0_values)

    theta_0 = theta_0_values(k);
    output_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        continue
    end

    S = load(mat_file, ...
        'psi_a', ...
        'vec_position_b_North', ...
        'vec_t', ...
        'm_a', 'm_b', 'N_a', 'N_b', 'hbar', 'geom');

    % Basic parameters
    m_a_case  = S.m_a;
    m_b_case  = S.m_b;
    N_a_case  = S.N_a;
    N_b_case  = S.N_b;
    hbar_case = S.hbar;
    R_case    = S.geom.R_sphere;

    M_b_case = m_b_case * N_b_case;
    M_a_case = m_a_case * N_a_case;
    n_a_case = N_a_case / (4*pi*R_case^2);

    M_a_theta(k) = M_a_case;

    % Better estimate of equilibrium angle
    t_case = S.vec_t(:);
    thetaN_case = S.vec_position_b_North(:,1);

    [~, ia] = unique(t_case, 'stable');
    thetaN_case = thetaN_case(ia);

    theta_eq_values(k) = mean(thetaN_case);

    % Compute excluded mass M_ex(theta_0)
    psi_a_case = S.psi_a;
    geom_case  = S.geom;

    if ~isequal(size(psi_a_case), size(geom_case.W))
        sz  = size(psi_a_case);
        wsz = size(geom_case.W);

        if numel(sz) == 3 && isequal(sz(1:2), wsz)
            psi_a_case = psi_a_case(:,:,end);
        elseif numel(sz) == 4 && isequal(sz(1:2), wsz)
            psi_a_case = psi_a_case(:,:,end,end);
        else
            error(['psi_a has incompatible size in %s. ' ...
                   'Size(psi_a) = %s, size(geom.W) = %s'], ...
                   mat_file, mat2str(size(psi_a_case)), mat2str(size(geom_case.W)));
        end
    end

    [M_ex_theta(k), ~] = Compute_excluded_mass_sphere(psi_a_case, m_a_case, N_a_case, geom_case);

    % Full ordered spectrum for Mb
    try
        [eigVals_Mb, ~] = eig_linear_matrix(M_b_case, m_a_case, n_a_case, R_case, hbar_case, theta_eq_values(k));
        imagEig_cell_theta_Mb{k} = sort(imag(eigVals_Mb(:)), 'ascend');
    catch ME
        warning('eig_linear_matrix failed for Mb at theta_0 = %.5f: %s', theta_0, ME.message);
        imagEig_cell_theta_Mb{k} = [];
    end

    % Full ordered spectrum for Mb + M_ex
    try
        if ~isfinite(M_ex_theta(k))
            error('Excluded mass is not finite at theta_0 = %.5f.', theta_0);
        end

        [eigVals_Mb_plusMex, ~] = eig_linear_matrix(M_b_case + M_ex_theta(k), ...
            m_a_case, n_a_case, R_case, hbar_case, theta_eq_values(k));

        imagEig_cell_theta_Mb_plusMex{k} = sort(imag(eigVals_Mb_plusMex(:)), 'ascend');
    catch ME
        warning('eig_linear_matrix failed for Mb+M_ex at theta_0 = %.5f: %s', theta_0, ME.message);
        imagEig_cell_theta_Mb_plusMex{k} = [];
    end
end

%% =========================================================
%  Manual branch selection
%  These indices refer to the FULL ordered spectrum:
%     imagVals = sort(imag(eigVals(:)), 'ascend')
%% =========================================================
branchIndex_theta_Mb         = 7;   % <-- set manually
branchIndex_theta_Mb_plusMex = 7;   % <-- set manually

% Extract the chosen branch for all theta_0
for k = 1:numel(theta_0_values)

    imagVals_Mb = imagEig_cell_theta_Mb{k};
    if ~isempty(imagVals_Mb) && numel(imagVals_Mb) >= branchIndex_theta_Mb
        omega_analytic_branch_theta_Mb(k) = imagVals_Mb(branchIndex_theta_Mb);
    end

    imagVals_MbMex = imagEig_cell_theta_Mb_plusMex{k};
    if ~isempty(imagVals_MbMex) && numel(imagVals_MbMex) >= branchIndex_theta_Mb_plusMex
        omega_analytic_branch_theta_Mb_plusMex(k) = imagVals_MbMex(branchIndex_theta_Mb_plusMex);
    end
end

fprintf('\n=========================================================\n');
fprintf('Manual branch selection vs theta_0\n');
fprintf('---------------------------------------------------------\n');
fprintf('Mb        : branch index = %d\n', branchIndex_theta_Mb);
fprintf('Mb + M_ex : branch index = %d\n', branchIndex_theta_Mb_plusMex);
fprintf('=========================================================\n\n');

%% =========================================================
%  Plot omega vs initial theta_0 with double analytical prediction
%% ==========================================================
fig7 = figure('Position', [100 100 1000 700], 'Color', 'w');
hold on
grid on
box on

xvals = theta_0_values / pi;

valid_num       = isfinite(omega_fit_theta);
valid_an_Mb     = isfinite(omega_analytic_branch_theta_Mb);
valid_an_MbMex  = isfinite(omega_analytic_branch_theta_Mb_plusMex);

plot(xvals(valid_num), omega_fit_theta(valid_num), 'o-', ...
    'LineWidth', 1.8, 'MarkerSize', 7, ...
    'DisplayName', 'Numerical fit')

plot(xvals(valid_an_Mb), omega_analytic_branch_theta_Mb(valid_an_Mb), '--', ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('Analytical branch %d ($M_b$)', branchIndex_theta_Mb))

plot(xvals(valid_an_MbMex), omega_analytic_branch_theta_Mb_plusMex(valid_an_MbMex), '-.', ...
    'LineWidth', 2.0, ...
    'DisplayName', sprintf('Analytical branch %d ($M_b + M_{\\mathrm{ex}}$)', branchIndex_theta_Mb_plusMex))

xlabel('$\theta_0/\pi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\omega$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('Angular frequency extracted from sinusoidal fits of $\gamma(t)$', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend('Interpreter', 'latex', 'FontSize', 13, 'Location', 'best')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xlim([min(xvals)-0.01, max(xvals)+0.01])
ylim([0,inf])

%% =========================================================
%  Plot M_ex / M_a vs theta_0
%% ==========================================================
fig8 = figure('Position', [120 120 900 650], 'Color', 'w');
hold on
grid on
box on

valid_mex = isfinite(M_ex_theta) & isfinite(M_a_theta) & (M_a_theta > 0);

plot(xvals(valid_mex), (M_ex_theta(valid_mex) ./ M_a_theta(valid_mex)), 'o-', ...
    'LineWidth', 1.8, 'MarkerSize', 7)

xlabel('$\theta_0/\pi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$M_{\mathrm{ex}}/M_a$', 'Interpreter', 'latex', 'FontSize', 16)
title('Excluded mass fraction vs initial $\theta_0$', ...
    'Interpreter', 'latex', 'FontSize', 18)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xlim([min(xvals)-0.01, max(xvals)+0.01])

%% =========================================================
%  Convert full spectra from cells to matrices for plotting
%% ==========================================================
nEig_theta_Mb = 0;
for k = 1:numel(theta_0_values)
    nEig_theta_Mb = max(nEig_theta_Mb, numel(imagEig_cell_theta_Mb{k}));
end

nEig_theta_Mb_plusMex = 0;
for k = 1:numel(theta_0_values)
    nEig_theta_Mb_plusMex = max(nEig_theta_Mb_plusMex, numel(imagEig_cell_theta_Mb_plusMex{k}));
end

imagEig_all_theta_Mb = nan(nEig_theta_Mb, numel(theta_0_values));
for k = 1:numel(theta_0_values)
    vals = imagEig_cell_theta_Mb{k};
    if ~isempty(vals)
        imagEig_all_theta_Mb(1:numel(vals), k) = vals(:);
    end
end

imagEig_all_theta_Mb_plusMex = nan(nEig_theta_Mb_plusMex, numel(theta_0_values));
for k = 1:numel(theta_0_values)
    vals = imagEig_cell_theta_Mb_plusMex{k};
    if ~isempty(vals)
        imagEig_all_theta_Mb_plusMex(1:numel(vals), k) = vals(:);
    end
end

%% =========================================================
%  Plot all Im(lambda) branches vs theta_0
%  for both Mb and Mb + M_ex
%% ==========================================================
fig9 = figure('Position', [100 100 1250 520], 'Color', 'w');
tl9 = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Left panel: Im(lambda) for Mb
nexttile
hold on
grid on
box on

for j = 1:size(imagEig_all_theta_Mb,1)
    plot(xvals, imagEig_all_theta_Mb(j,:), '-', 'LineWidth', 1.2);
end

if branchIndex_theta_Mb <= size(imagEig_all_theta_Mb,1)
    plot(xvals, imagEig_all_theta_Mb(branchIndex_theta_Mb,:), 'k--', 'LineWidth', 2.5);
end

xlabel('$\theta_0/\pi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${\rm Im}(\lambda)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('$M_b$', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
xlim([min(xvals), max(xvals)])

% Right panel: Im(lambda) for Mb + M_ex
nexttile
hold on
grid on
box on

for j = 1:size(imagEig_all_theta_Mb_plusMex,1)
    plot(xvals, imagEig_all_theta_Mb_plusMex(j,:), '-', 'LineWidth', 1.2);
end

if branchIndex_theta_Mb_plusMex <= size(imagEig_all_theta_Mb_plusMex,1)
    plot(xvals, imagEig_all_theta_Mb_plusMex(branchIndex_theta_Mb_plusMex,:), 'k--', 'LineWidth', 2.5);
end

xlabel('$\theta_0/\pi$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${\rm Im}(\lambda)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('$M_b + M_{\mathrm{ex}}$', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
xlim([min(xvals), max(xvals)])

title(tl9, 'Imaginary parts of all eigenvalues from linearized theory', ...
    'Interpreter', 'latex', 'FontSize', 18)
