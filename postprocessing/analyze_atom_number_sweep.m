%% Analyze the core-filling atom-number sweep
% Load the saved real-time trajectories and diagnostics for different N_b
% values and generate the comparison plots used in the analysis.
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));
output_root = fullfile(repository_root, 'output');

base_output_folder = fullfile(output_root, 'Sweep_Nb');

% Values to display
N_b_values = 2000:1000:10000;

% Create figure
fig = figure('Position', [100 100 1400 1200], 'Color', 'w');
tl = tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(numel(N_b_values),1);

for k = 1:numel(N_b_values)

    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    ax(k) = nexttile;

    if ~exist(mat_file, 'file')
        axis(ax(k), 'off')
        title(ax(k), sprintf('$N_b = %d$\\newline file not found', N_b), ...
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

    title(sprintf('$N_b = %d$', N_b), ...
        'Interpreter', 'latex', 'FontSize', 15)

    hold off
end

sgtitle(tl, 'Vortex trajectories on the sphere for different $N_b$', ...
    'Interpreter', 'latex', 'FontSize', 20)

% Link camera properties across all valid axes
valid_ax = ax(isgraphics(ax));

hlink = linkprop(valid_ax, ...
    {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});

% Store link handle so it is not destroyed
setappdata(fig, 'CameraLink', hlink);

% Enable interactive rotation
rotate3d(fig, 'on');


fig = figure('Position', [100 100 1200 1000], 'Color', 'w');
tl = tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(N_b_values)

    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    nexttile
    axis off

    if ~exist(mat_file, 'file')
        text(0.5, 0.5, sprintf('$N_b = %d$\\newline file not found', N_b), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','center', ...
            'FontSize',14)
        continue
    end

    S = load(mat_file, 'vec_Ene', 'vec_L_z');

    vec_E = S.vec_Ene(:);
    vec_L = S.vec_L_z(:);

    if isempty(vec_E) || isempty(vec_L)
        text(0.5, 0.5, sprintf('$N_b = %d$\\newline empty data', N_b), ...
            'Interpreter','latex', ...
            'HorizontalAlignment','center', ...
            'FontSize',14)
        continue
    end

    E0 = vec_E(1);
    Ef = vec_E(end);

    L0 = vec_L(1);
    Lf = vec_L(end);

    dE_percent = 100 * (Ef - E0) / E0;

    
    dL_percent = 100 * (Lf - L0) / L0;
  

    % Create formatted text
    str = sprintf([...
        '$N_b = %d$ \n\n' ...
        '$\\Delta E = %.3f\\%%$ \n' ...
        '$\\Delta L_z = %.3f\\%%$'], ...
        N_b, dE_percent, dL_percent);

    text(0.5, 0.5, str, ...
        'Interpreter','latex', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', ...
        'FontSize',16)

end

sgtitle(tl, 'Final relative variations of $E$ and $L_z$', ...
    'Interpreter','latex', 'FontSize',20)


%% Plot E_PV(t) for all N_b

fig = figure('Position', [100 100 900 650], 'Color', 'w');
hold on
grid on
box on

cmap = lines(numel(N_b_values));

for k = 1:3:numel(N_b_values)

    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        fprintf('File not found for N_b = %d\n', N_b);
        continue
    end

    S = load(mat_file, ...
        'vec_t', ...
        'vec_position_b_North', ...
        'vec_position_b_South', ...
        'm_a','m_b','N_a','hbar','geom');

    t = S.vec_t(:);

    thetaN = S.vec_position_b_North(:,1);
    phiN   = S.vec_position_b_North(:,2);

    thetaS = S.vec_position_b_South(:,1);
    phiS   = S.vec_position_b_South(:,2);

    % --- unwrap ---
    phiN_u = unwrap(phiN);
    phiS_u = unwrap(phiS);

    % --- derivatives ---
    dthetaN = gradient(thetaN, t);
    dphiN   = gradient(phiN_u, t);

    dthetaS = gradient(thetaS, t);
    dphiS   = gradient(phiS_u, t);

    % --- parameters ---
    R    = S.geom.R_sphere;
    Mb   = S.m_b * N_b;
    hbar = S.hbar;
    m_a  = S.m_a;
    N_a  = S.N_a;

    % --- kinetic energy ---
    T_N = (Mb/4) * R^2 .* (dthetaN.^2 + sin(thetaN).^2 .* dphiN.^2);
    T_S = (Mb/4) * R^2 .* (dthetaS.^2 + sin(thetaS).^2 .* dphiS.^2);
    T   = T_N + T_S;

    % --- chi (interaction) ---
    DeltaPhi = phiN_u - phiS_u;

    arg = 2 ...
        - 2*cos(thetaN).*cos(thetaS) ...
        - 2*sin(thetaN).*sin(thetaS).*cos(DeltaPhi);

    arg = max(arg, eps);
    chi = 0.5 * log(arg);

    % --- dipolar energy ---
    na = N_a / (4*pi*R^2);
    prefE = 2 * hbar^2 * pi * na / m_a;

    xi = R / 50;
    E_const = prefE * log(R/xi);

    E_dip = E_const + prefE * chi;

    % --- total PV energy ---
    E_PV = T + E_dip;

    plot(t, E_PV, 'LineWidth', 1.8, ...
        'Color', cmap(k,:), ...
        'DisplayName', sprintf('$N_b = %d$', N_b));
end

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$E_{\mathrm{PV}}(t)$ [J]', 'Interpreter', 'latex', 'FontSize', 16)
title('Point-vortex energy reconstructed from trajectories', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)

%% Plot angular distance between the two vortices vs time

fig = figure('Position', [100 100 900 650], 'Color', 'w');
hold on
grid on
box on

cmap = lines(numel(N_b_values));

for k = 1:numel(N_b_values)

    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~exist(mat_file, 'file')
        fprintf('File not found for N_b = %d\n', N_b);
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

    if isempty(t) || isempty(thetaN) || isempty(thetaS) || ...
       any(~isfinite(t)) || any(~isfinite(thetaN)) || any(~isfinite(phiN)) || ...
       any(~isfinite(thetaS)) || any(~isfinite(phiS))
        fprintf('Invalid data for N_b = %d\n', N_b);
        continue
    end

    % Angular distance on the sphere
    cos_gamma = cos(thetaN).*cos(thetaS) + ...
                sin(thetaN).*sin(thetaS).*cos(phiN - phiS);

    % Numerical safety
    cos_gamma = min(max(cos_gamma, -1), 1);

    gamma = acos(cos_gamma);

    plot(t, gamma, 'LineWidth', 1.8, ...
        'Color', cmap(k,:), ...
        'DisplayName', sprintf('$N_b = %d$', N_b));
end

xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\gamma(t)$ [rad]', 'Interpreter', 'latex', 'FontSize', 16)
title('Angular distance between the two vortices', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)



%% Fit sinusoidally gamma(t) for each N_b and extract omega(N_b)

omega_fit   = nan(size(N_b_values));
A_fit       = nan(size(N_b_values));
gamma0_fit  = nan(size(N_b_values));
phi_fit     = nan(size(N_b_values));

fig_fit = figure('Position', [100 100 1400 1200], 'Color', 'w');
tl = tiledlayout(3,3, 'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:numel(N_b_values)

    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    nexttile
    hold on
    grid on
    box on

    if ~exist(mat_file, 'file')
        title(sprintf('$N_b = %d$\\newline file not found', N_b), ...
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
        title(sprintf('$N_b = %d$\\newline invalid data', N_b), ...
            'Interpreter', 'latex', 'FontSize', 14)
        axis off
        continue
    end

    % Angular distance gamma(t)
    cos_gamma = cos(thetaN).*cos(thetaS) + ...
                sin(thetaN).*sin(thetaS).*cos(phiN - phiS);
    cos_gamma = min(max(cos_gamma, -1), 1);
    gamma = acos(cos_gamma);

    % Remove possible repeated times
    [t, ia] = unique(t, 'stable');
    gamma = gamma(ia);

    % Initial guesses
    gamma0_guess = mean(gamma);
    A_guess = 0.5 * (max(gamma) - min(gamma));

    % Frequency guess from FFT
    dt = mean(diff(t));
    y = gamma - mean(gamma);
    n = numel(y);

    if n < 8 || all(abs(y) < 1e-12)
        % Degenerate case: almost constant signal
        gamma0_fit(k) = gamma0_guess;
        A_fit(k) = 0;
        omega_fit(k) = 0;
        phi_fit(k) = 0;

        plot(t, gamma, 'k', 'LineWidth', 1.5)
        title(sprintf('$N_b = %d$\\newline quasi-constant $\\gamma(t)$', N_b), ...
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

    gamma0_fit(k) = pfit(1);
    A_fit(k)      = pfit(2);
    omega_fit(k)  = pfit(3);
    phi_fit(k)    = pfit(4);

    % Plot data + fit
    tt_fine = linspace(t(1), t(end), 1000);
    gamma_fit_curve = model(pfit, tt_fine);

    plot(t, gamma, 'ko', 'MarkerSize', 3, 'DisplayName', 'data')
    plot(tt_fine, gamma_fit_curve, 'r-', 'LineWidth', 1.8, 'DisplayName', 'fit')

    title(sprintf('$N_b = %d$\\newline $\\omega = %.4g\\ \\mathrm{rad/s}$', ...
        N_b, omega_fit(k)), ...
        'Interpreter', 'latex', 'FontSize', 13)

    xlabel('$t$ [s]', 'Interpreter', 'latex')
    ylabel('$\gamma$ [rad]', 'Interpreter', 'latex')
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11)

    hold off
end

sgtitle(tl, 'Sinusoidal fits of $\gamma(t)$', ...
    'Interpreter', 'latex', 'FontSize', 18);


%% =========================================================
% Analytical prediction from linearized theory:
% for each N_b, choose the analytical mode closest to omega_fit
%
% We compute two analytical curves:
%   1) using Mb
%   2) using Mb + M_ex, where M_ex is the excluded mass of psi_a
%% =========================================================

reference_file = fullfile(base_output_folder, 'N_b_2000', ...
    'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

if ~isfile(reference_file)
    error('Reference file not found: %s', reference_file);
end

S_ref = load(reference_file, ...
    'm_a', 'm_b', 'N_a', 'hbar', 'geom', 'vec_position_b_North');

m_a    = S_ref.m_a;
m_b    = S_ref.m_b;
N_a    = S_ref.N_a;
hbar   = S_ref.hbar;
R      = S_ref.geom.R_sphere;
thetaN = S_ref.vec_position_b_North(:,1);

na     = N_a / (4*pi*R^2);   % 2D density of component a
theta0 = mean(thetaN, 'omitnan'); % representative polar angle
tolEig = 1e-8;               % tolerance for zero / numerical modes

% =========================================================
% Compute excluded mass M_ex for each N_b
% =========================================================
M_ex_values = nan(size(N_b_values));

for k = 1:numel(N_b_values)

    Nb = N_b_values(k);

    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', Nb));
    mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_with_Ec_Ei.mat');

    if ~isfile(mat_file)
        warning('File not found: %s', mat_file);
        continue
    end

    S = load(mat_file);

    % Required variables
    if ~isfield(S, 'psi_a') || ~isfield(S, 'N_a') || ~isfield(S, 'm_a') || ~isfield(S, 'geom')
        warning('Missing one or more required variables in %s', mat_file);
        continue
    end

    psi_a_case = S.psi_a;
    N_a_case   = S.N_a;
    m_a_case   = S.m_a;
    geom_case  = S.geom;

    % If psi_a is time dependent, take final state
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

    [M_ex_values(k), ~] = Compute_excluded_mass_sphere(psi_a_case, m_a_case, N_a_case, geom_case);

end


%% =========================================================
% Analytical curves built from the SAME ordered spectrum used in Fig. 8
%
% Strategy:
%   1) for each N_b, sort all Im(lambda) in ascending order
%   2) choose a reference branch index from the first N_b where omega_fit
%      is available, selecting the positive Im(lambda) closest to omega_fit
%   3) keep that SAME index for all N_b
%
% This ensures that Fig. 6 follows one branch of the ordered spectrum,
% consistently with Fig. 8.
%% =========================================================

% ---------------------------------------------------------
% Containers for full ordered spectra
% ---------------------------------------------------------
imagEig_cell_Mb         = cell(size(N_b_values));
imagEig_cell_Mb_plusMex = cell(size(N_b_values));

omega_analytic_branch_Mb         = nan(size(N_b_values));
omega_analytic_branch_Mb_plusMex = nan(size(N_b_values));

branchIndex_Mb         = NaN;   % index in FULL sorted spectrum
branchIndex_Mb_plusMex = NaN;   % index in FULL sorted spectrum

k_ref_Mb         = NaN;
k_ref_Mb_plusMex = NaN;

% ---------------------------------------------------------
% Build ordered spectra for all N_b
% ---------------------------------------------------------
for k = 1:numel(N_b_values)

    Nb   = N_b_values(k);
    Mb   = m_b * Nb;
    M_ex = M_ex_values(k);

    % --------------------------------------------
    % Case 1: Mb
    % --------------------------------------------
    try
        [eigVals_Mb, ~] = eig_linear_matrix(Mb, m_a, na, R, hbar, theta0);
        imagEig_cell_Mb{k} = sort(imag(eigVals_Mb(:)), 'ascend');
    catch ME
        warning('eig_linear_matrix failed for Mb at N_b = %d: %s', Nb, ME.message);
        imagEig_cell_Mb{k} = [];
    end

    % --------------------------------------------
    % Case 2: Mb + M_ex
    % --------------------------------------------
    try
        if ~isfinite(M_ex)
            error('Excluded mass is not finite for N_b = %g.', Nb);
        end

        [eigVals_Mb_plusMex, ~] = eig_linear_matrix(Mb + M_ex, m_a, na, R, hbar, theta0);
        imagEig_cell_Mb_plusMex{k} = sort(imag(eigVals_Mb_plusMex(:)), 'ascend');
    catch ME
        warning('eig_linear_matrix failed for Mb+M_ex at N_b = %d: %s', Nb, ME.message);
        imagEig_cell_Mb_plusMex{k} = [];
    end
end

%% =========================================================
% Choose fixed branch index in the FULL ordered spectrum
% using the first valid omega_fit as reference
%% =========================================================

% ---------------------------------------------------------
% Mb
% ---------------------------------------------------------
for k = 1:numel(N_b_values)

    if ~isfinite(omega_fit(k)) || isempty(imagEig_cell_Mb{k})
        continue
    end

    imagVals = imagEig_cell_Mb{k};

    % Restrict to positive oscillatory modes only, but keep their
    % indices in the FULL sorted spectrum
    idxPos = find(imagVals > tolEig);

    if isempty(idxPos)
        continue
    end

    [~, loc] = min(abs(imagVals(idxPos) - omega_fit(k)));

    k_ref_Mb = k;
    branchIndex_Mb = idxPos(loc);   % FULL-spectrum index
    break
end

% ---------------------------------------------------------
% Mb + M_ex
% ---------------------------------------------------------
for k = 1:numel(N_b_values)

    if ~isfinite(omega_fit(k)) || isempty(imagEig_cell_Mb_plusMex{k})
        continue
    end

    imagVals = imagEig_cell_Mb_plusMex{k};

    idxPos = find(imagVals > tolEig);

    if isempty(idxPos)
        continue
    end

    [~, loc] = min(abs(imagVals(idxPos) - omega_fit(k)));

    k_ref_Mb_plusMex = k;
    branchIndex_Mb_plusMex = idxPos(loc);   % FULL-spectrum index
    break
end

%% =========================================================
% Extract that SAME branch for all N_b
%% =========================================================

for k = 1:numel(N_b_values)

    % ----- Mb -----
    if isfinite(branchIndex_Mb) && ~isempty(imagEig_cell_Mb{k})
        imagVals = imagEig_cell_Mb{k};
        if numel(imagVals) >= branchIndex_Mb
            omega_analytic_branch_Mb(k) = imagVals(branchIndex_Mb);
        end
    end

    % ----- Mb + M_ex -----
    if isfinite(branchIndex_Mb_plusMex) && ~isempty(imagEig_cell_Mb_plusMex{k})
        imagVals = imagEig_cell_Mb_plusMex{k};
        if numel(imagVals) >= branchIndex_Mb_plusMex
            omega_analytic_branch_Mb_plusMex(k) = imagVals(branchIndex_Mb_plusMex);
        end
    end
end

%% =========================================================
% Print selected branch information
%% =========================================================
fprintf('\n=========================================================\n');
fprintf('Fixed branch selection from FULL ordered spectrum\n');
fprintf('---------------------------------------------------------\n');

if isfinite(branchIndex_Mb)
    fprintf('Mb        : branch index = %d   selected at N_b = %d\n', ...
        branchIndex_Mb, N_b_values(k_ref_Mb));
else
    fprintf('Mb        : no valid branch found\n');
end

if isfinite(branchIndex_Mb_plusMex)
    fprintf('Mb + M_ex : branch index = %d   selected at N_b = %d\n', ...
        branchIndex_Mb_plusMex, N_b_values(k_ref_Mb_plusMex));
else
    fprintf('Mb + M_ex : no valid branch found\n');
end

fprintf('=========================================================\n\n');

%% =========================================================
% Convert cell spectra to matrices for Fig. 8-style plotting
%% =========================================================
nEig_Mb = 0;
for k = 1:numel(N_b_values)
    nEig_Mb = max(nEig_Mb, numel(imagEig_cell_Mb{k}));
end

nEig_Mb_plusMex = 0;
for k = 1:numel(N_b_values)
    nEig_Mb_plusMex = max(nEig_Mb_plusMex, numel(imagEig_cell_Mb_plusMex{k}));
end

imagEig_all_Mb = nan(nEig_Mb, numel(N_b_values));
for k = 1:numel(N_b_values)
    vals = imagEig_cell_Mb{k};
    if ~isempty(vals)
        imagEig_all_Mb(1:numel(vals), k) = vals(:);
    end
end

imagEig_all_Mb_plusMex = nan(nEig_Mb_plusMex, numel(N_b_values));
for k = 1:numel(N_b_values)
    vals = imagEig_cell_Mb_plusMex{k};
    if ~isempty(vals)
        imagEig_all_Mb_plusMex(1:numel(vals), k) = vals(:);
    end
end

%% =========================================================
% Plot imaginary parts of all eigenvalues vs N_b
%% =========================================================
fig_spec = figure('Position', [100 100 1200 500], 'Color', 'w');
tl = tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---------------------------------------------------------
% Left panel: Im(lambda) for Mb
% ---------------------------------------------------------
nexttile
hold on
grid on
box on

for j = 1:size(imagEig_all_Mb,1)
    plot(N_b_values, imagEig_all_Mb(j,:), '-', 'LineWidth', 1.2);
end

xlabel('$N_b$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${\rm Im}(\lambda)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('$M_b$', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
xlim([min(N_b_values), max(N_b_values)])

% ---------------------------------------------------------
% Right panel: Im(lambda) for Mb + M_ex
% ---------------------------------------------------------
nexttile
hold on
grid on
box on

for j = 1:size(imagEig_all_Mb_plusMex,1)
    plot(N_b_values, imagEig_all_Mb_plusMex(j,:), '-', 'LineWidth', 1.2);
end

xlabel('$N_b$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('${\rm Im}(\lambda)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('$M_b + M_{\mathrm{ex}}$', 'Interpreter', 'latex', 'FontSize', 17)
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13)
xlim([min(N_b_values), max(N_b_values)])

title(tl, 'Imaginary parts of all eigenvalues from linearized theory', ...
    'Interpreter', 'latex', 'FontSize', 18)

%% =========================================================
% Rebuild Fig. 6 by manually selecting the ordered-spectrum branch index
%
% Required variables already in workspace:
%   N_b_values
%   omega_fit
%   imagEig_all_Mb
%   imagEig_all_Mb_plusMex
%% =========================================================

% ---------------------------------------------------------
% USER CHOICE: branch indices to plot
% ---------------------------------------------------------
branchIndex_Mb         = 7;   % <-- choose here
branchIndex_Mb_plusMex = 7;   % <-- choose here

% ---------------------------------------------------------
% Safety checks
% ---------------------------------------------------------
if ~exist('imagEig_all_Mb', 'var') || ~exist('imagEig_all_Mb_plusMex', 'var')
    error(['imagEig_all_Mb and/or imagEig_all_Mb_plusMex not found in workspace. ' ...
           'Run first the block that computes the full ordered spectra.']);
end

if ~exist('omega_fit', 'var') || ~exist('N_b_values', 'var')
    error('omega_fit and/or N_b_values not found in workspace.');
end

if branchIndex_Mb < 1 || branchIndex_Mb > size(imagEig_all_Mb,1)
    error('branchIndex_Mb = %d is out of range. Valid range: 1 ... %d', ...
          branchIndex_Mb, size(imagEig_all_Mb,1));
end

if branchIndex_Mb_plusMex < 1 || branchIndex_Mb_plusMex > size(imagEig_all_Mb_plusMex,1)
    error('branchIndex_Mb_plusMex = %d is out of range. Valid range: 1 ... %d', ...
          branchIndex_Mb_plusMex, size(imagEig_all_Mb_plusMex,1));
end

% ---------------------------------------------------------
% Extract chosen branches
% ---------------------------------------------------------
omega_analytic_branch_Mb         = imagEig_all_Mb(branchIndex_Mb, :);
omega_analytic_branch_Mb_plusMex = imagEig_all_Mb_plusMex(branchIndex_Mb_plusMex, :);

% Optional: ignore non-positive frequencies if you only want oscillatory modes
% omega_analytic_branch_Mb(omega_analytic_branch_Mb <= 0) = NaN;
% omega_analytic_branch_Mb_plusMex(omega_analytic_branch_Mb_plusMex <= 0) = NaN;

% ---------------------------------------------------------
% Valid entries
% ---------------------------------------------------------
valid_num             = isfinite(omega_fit);
valid_an_branch_Mb    = isfinite(omega_analytic_branch_Mb);
valid_an_branch_MbMex = isfinite(omega_analytic_branch_Mb_plusMex);

% ---------------------------------------------------------
% Plot Fig. 6-like figure
% ---------------------------------------------------------
fig_omega_manual = figure('Position', [100 100 900 650], 'Color', 'w');
hold on
grid on
box on

plot(N_b_values(valid_num), omega_fit(valid_num), 'o-', ...
    'LineWidth', 1.8, 'MarkerSize', 7)

plot(N_b_values(valid_an_branch_Mb), omega_analytic_branch_Mb(valid_an_branch_Mb), '--', ...
    'LineWidth', 2.0)

plot(N_b_values(valid_an_branch_MbMex), omega_analytic_branch_Mb_plusMex(valid_an_branch_MbMex), '-.', ...
    'LineWidth', 2.0)

xlabel('$N_b$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$\omega$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 16)
title('Angular frequency extracted from sinusoidal fits of $\gamma(t)$', ...
    'Interpreter', 'latex', 'FontSize', 18)

legend({ ...
    'Numerical fit', ...
    sprintf('Ordered-spectrum branch %d ($M_b$)', branchIndex_Mb), ...
    sprintf('Ordered-spectrum branch %d ($M_b + M_{\\mathrm{ex}}$)', branchIndex_Mb_plusMex)}, ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xlim([min(N_b_values), max(N_b_values)])

% ---------------------------------------------------------
% Print selected values for quick inspection
% ---------------------------------------------------------
fprintf('\n=========================================================\n');
fprintf('Manual branch selection for Fig. 6\n');
fprintf('---------------------------------------------------------\n');
fprintf('Mb branch index         = %d\n', branchIndex_Mb);
fprintf('Mb + M_ex branch index  = %d\n', branchIndex_Mb_plusMex);
fprintf('=========================================================\n\n');

%% =========================================================
% Plot M_ex / M_a vs N_b
%% =========================================================

% Total mass of component a
M_a = N_a * m_a;

% Valid entries
valid_mex = isfinite(M_ex_values);

fig_mex = figure('Position', [100 100 900 650], 'Color', 'w');
hold on
grid on
box on

plot(N_b_values(valid_mex), M_ex_values(valid_mex) / M_a, 'o-', ...
    'LineWidth', 1.8, 'MarkerSize', 7)

xlabel('$N_b$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$M_{\mathrm{ex}}/M_a$', 'Interpreter', 'latex', 'FontSize', 16)
title('Excluded mass fraction vs $N_b$', ...
    'Interpreter', 'latex', 'FontSize', 18)

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14)
xlim([min(N_b_values), max(N_b_values)])
