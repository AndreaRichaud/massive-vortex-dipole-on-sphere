function Run_real_time_dynamics_residual_gravity( ...
    N_b, g_fraction, input_mat_file, output_folder)

%RUN_REAL_TIME_DYNAMICS_RESIDUAL_GRAVITY Evolve the dipole under residual gravity.
%
%   RUN_REAL_TIME_DYNAMICS_RESIDUAL_GRAVITY(N_b, g_fraction,
%   input_mat_file, output_folder) loads a converged two-component state,
%   adds the gravitational potentials
%
%       V_g,a(theta) = m_a a_res R cos(theta),
%       V_g,b(theta) = m_b a_res R cos(theta),
%
%   with a_res = g_fraction * g_earth, and performs the real-time evolution.
%
%   Inputs
%   ------
%   N_b
%       Atom number of the core-filling component b.
%   g_fraction
%       Residual acceleration in units of terrestrial gravity.
%   input_mat_file
%       MAT file containing the converged initial state.
%   output_folder
%       Directory used for results, frames, figures, and the execution log.
%
%   Main file written
%   -----------------
%   Real_time_dynamics_on_sphere_residual_g_with_Ec_Ei.mat
%
%   External dependency: SSHT MATLAB library.

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

temp_real_folder = fullfile(output_folder, 'Temp_real');
if ~exist(temp_real_folder, 'dir')
    mkdir(temp_real_folder);
end

log_file = fullfile(output_folder, 'log_movie_vortexon_residual_g.txt');
if exist(log_file, 'file')
    delete(log_file);
end

diary(log_file);
diary on;

t_start = tic;

try
    fprintf('============================================================\n');
    fprintf('Run_real_time_dynamics_residual_gravity started\n');
    fprintf('N_b = %d\n', N_b);
    fprintf('g_fraction = %.6e\n', g_fraction);
    fprintf('Output folder: %s\n', output_folder);
    fprintf('Start time: %s\n', datestr(now));
    fprintf('============================================================\n');

    %--------------------------------------------------------------%
    % Path setup
    %--------------------------------------------------------------%
    this_file_folder = fileparts(mfilename('fullpath'));
    addpath(genpath(this_file_folder));

    %--------------------------------------------------------------%
    % Load input MAT
    %--------------------------------------------------------------%
    if ~exist(input_mat_file, 'file')
        error('Input file not found: %s', input_mat_file);
    end

    S = load(input_mat_file);

    % unpack variables used below
    hbar       = S.hbar;
    m_a        = S.m_a;
    m_b        = S.m_b;
    N_a        = S.N_a;
    L_z        = S.L_z;
    g_a        = S.g_a;
    g_b        = S.g_b;
    g_ab       = S.g_ab;
    geom       = S.geom;
    dt         = S.dt;
    psi_a      = S.psi_a;
    psi_b      = S.psi_b;
    theta_1    = S.theta_1;
    phi_1      = S.phi_1;
    theta_2    = S.theta_2;
    phi_2      = S.phi_2;
    V0_pin_a   = S.V0_pin_a;
    sigma_pin  = S.sigma_pin;

    % Sanity check: ensure loaded N_b matches requested N_b
    if isfield(S, 'N_b') && S.N_b ~= N_b
        warning('Loaded N_b (%d) differs from function input N_b (%d). Using input N_b.', S.N_b, N_b);
    end

    fprintf('Ground-state file loaded successfully.\n');

    %--------------------------------------------------------------%
    % No pinning during real-time evolution
    %--------------------------------------------------------------%
    Pinning_flag = 0;

    [Mat_V_a, Mat_V_b] = Build_pinning_potential_sphere( ...
        geom, ...
        theta_1, phi_1, theta_2, phi_2, ...
        Pinning_flag, V0_pin_a, sigma_pin);

   %--------------------------------------------------------------%
% Residual-gravity potential acting on both components
%--------------------------------------------------------------%
validateattributes(g_fraction, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, ...
    mfilename, 'g_fraction');

g_earth = 9.80665;                  % standard gravity [m/s^2]
a_res   = g_fraction * g_earth;     % residual acceleration [m/s^2]
R       = geom.R_sphere;            % sphere radius [m]

% geom.theta is the polar-angle grid.
V_g_a = m_a * a_res * R .* cos(geom.theta);   % [J]
V_g_b = m_b * a_res * R .* cos(geom.theta);   % [J]

if ~isequal(size(V_g_a), size(psi_a))
    error(['V_g_a and psi_a have incompatible sizes: ' ...
        'size(V_g_a) = [%s], size(psi_a) = [%s].'], ...
        num2str(size(V_g_a)), num2str(size(psi_a)));
end

if ~isequal(size(V_g_b), size(psi_b))
    error(['V_g_b and psi_b have incompatible sizes: ' ...
        'size(V_g_b) = [%s], size(psi_b) = [%s].'], ...
        num2str(size(V_g_b)), num2str(size(psi_b)));
end

Mat_V_a = Mat_V_a + V_g_a;
Mat_V_b = Mat_V_b + V_g_b;

fprintf('Residual gravity added to components a and b.\n');
fprintf('g_earth = %.8f m/s^2\n', g_earth);
fprintf('a_res = %.8e m/s^2 = %.3e g_earth\n', ...
    a_res, g_fraction);
fprintf('max|V_g_a| = %.8e J\n', max(abs(V_g_a(:))));
fprintf('max|V_g_b| = %.8e J\n', max(abs(V_g_b(:))));

    %--------------------------------------------------------------%
    % Real-time evolution parameters
    %--------------------------------------------------------------%
    c_p = 1;

    n_iterations     = 1e5;
    sample_frequency = 1000;

    [vec_Ene, vec_L_z, vec_position_b_North, vec_position_b_South, ...
        vec_t, vec_Ec, vec_Ei, vec_Ehd, vec_Eqp] = ...
        Initialize_output_arrays_sphere(n_iterations, sample_frequency);

    fprintf('Real-time evolution setup completed.\n');
    fprintf('n_iterations = %d\n', n_iterations);
    fprintf('sample_frequency = %d\n', sample_frequency);

    %--------------------------------------------------------------%
    % Real-time evolution loop
    %--------------------------------------------------------------%
    for i = 1:n_iterations

        [psi_a, psi_b] = Real_time_evolve_sphere( ...
            psi_a, psi_b, ...
            m_a, m_b, ...
            g_a, g_b, g_ab, ...
            N_a, N_b, ...
            L_z, ...
            Mat_V_a, Mat_V_b, ...
            geom, hbar, dt);

        if i == 1 || mod(i, sample_frequency) == 0

            elapsed_time = toc(t_start);

            vec_t(c_p) = dt*i;

            vec_Ene(c_p) = Compute_total_energy_sphere( ...
                psi_a, psi_b, ...
                m_a, m_b, ...
                g_a, g_b, g_ab, ...
                N_a, N_b, ...
                L_z, ...
                Mat_V_a, Mat_V_b, ...
                geom, hbar);

            vec_L_z(c_p) = Total_angular_momentum_L_z_sphere( ...
                psi_a, psi_b, N_a, N_b, geom, hbar);

            % Relative percentage variations with respect to the first sampled value
            if c_p == 1
                dE_percent  = 0;
                dLz_percent = 0;
            else
                E0  = vec_Ene(1);
                Lz0 = vec_L_z(1);

                dE_percent = 100 * (vec_Ene(c_p) - E0) / E0;

                dLz_percent = 100 * (vec_L_z(c_p) - Lz0) / Lz0;
              
            end

            fprintf(['Avanzamento = %.2f %% in %.0f s | ' ...
                'dE = %+10.3e %% | dL_z = %+10.3e %%\n'], ...
                100*i/n_iterations, elapsed_time, dE_percent, dLz_percent);

            [theta_N, phi_N, theta_S, phi_S] = ...
                Compute_b_peak_positions_sphere(psi_b, geom);

            vec_position_b_North(c_p,1) = theta_N;
            vec_position_b_North(c_p,2) = phi_N;
            vec_position_b_South(c_p,1) = theta_S;
            vec_position_b_South(c_p,2) = phi_S;

            [vec_Ec(c_p), vec_Ei(c_p), vec_Ehd(c_p), vec_Eqp(c_p), out_hydro] = ...
                Compute_Ec_Ei_sphere(psi_a, m_a, N_a, geom, hbar); %#ok<NASGU>

            fig = Plot_frame_sphere(psi_a, psi_b, geom.L_band, geom.method, ...
                'Visible', 'off', ...
                'FigureTitle', sprintf('$t = %.3f\\ \\mathrm{s}$', dt*i));

            fig_name = fullfile(temp_real_folder, sprintf('%05d.jpg', c_p));
            saveas(fig, fig_name);
            close(fig);

            c_p = c_p + 1;
        end
    end

    %--------------------------------------------------------------%
    % Trim arrays
    %--------------------------------------------------------------%
    c_p = c_p - 1;

    vec_Ene              = vec_Ene(1:c_p);
    vec_L_z              = vec_L_z(1:c_p);
    vec_position_b_North = vec_position_b_North(1:c_p,:);
    vec_position_b_South = vec_position_b_South(1:c_p,:);
    vec_t                = vec_t(1:c_p);
    vec_Ec               = vec_Ec(1:c_p);
    vec_Ei               = vec_Ei(1:c_p);
    vec_Ehd              = vec_Ehd(1:c_p);
    vec_Eqp              = vec_Eqp(1:c_p);

    % force column vectors where convenient
    vec_Ene = vec_Ene(:);
    vec_L_z = vec_L_z(:);
    vec_t   = vec_t(:);
    vec_Ec  = vec_Ec(:);
    vec_Ei  = vec_Ei(:);
    vec_Ehd = vec_Ehd(:);
    vec_Eqp = vec_Eqp(:);

    %--------------------------------------------------------------%
    % Save real-time dynamics MAT
    %--------------------------------------------------------------%
    output_mat_file = fullfile(output_folder, 'Real_time_dynamics_on_sphere_residual_g_with_Ec_Ei.mat');

    save(output_mat_file, ...
        'N_b', ...
        'g_fraction', 'g_earth', 'a_res', 'V_g_a', 'V_g_b', ...
        'input_mat_file', ...
        'n_iterations', 'sample_frequency', ...
        'Pinning_flag', ...
        'psi_a', 'psi_b', ...
        'Mat_V_a', 'Mat_V_b', ...
        'vec_Ene', 'vec_L_z', ...
        'vec_position_b_North', 'vec_position_b_South', ...
        'vec_t', 'vec_Ec', 'vec_Ei', 'vec_Ehd', 'vec_Eqp', ...
        'm_a', 'm_b', 'N_a', 'L_z', ...
        'g_a', 'g_b', 'g_ab', ...
        'geom', 'hbar', 'dt', ...
        'theta_1', 'phi_1', 'theta_2', 'phi_2', ...
        'V0_pin_a', 'sigma_pin');

    fprintf('Saved MAT file: %s\n', output_mat_file);

    %% Plot vortex trajectories on the sphere
    R = geom.R_sphere;

    theta_N = vec_position_b_North(:,1);
    phi_N   = vec_position_b_North(:,2);

    theta_S = vec_position_b_South(:,1);
    phi_S   = vec_position_b_South(:,2);

    xN = sin(theta_N).*cos(phi_N);
    yN = sin(theta_N).*sin(phi_N);
    zN = cos(theta_N);

    xS = sin(theta_S).*cos(phi_S);
    yS = sin(theta_S).*sin(phi_S);
    zS = cos(theta_S);

    fig_traj = figure('Visible','off','Position',[100 100 900 800]);

    [XS,YS,ZS] = sphere(120);
    surf(XS,YS,ZS, ...
        'FaceAlpha',0.08, ...
        'EdgeColor','none', ...
        'FaceColor',[0.8 0.8 0.8]);

    hold on;
    plot3(xN,yN,zN,'r','LineWidth',2);
    plot3(xS,yS,zS,'b','LineWidth',2);

    plot3(xN(1),yN(1),zN(1),'ro','MarkerFaceColor','r','MarkerSize',7);
    plot3(xS(1),yS(1),zS(1),'bo','MarkerFaceColor','b','MarkerSize',7);

    plot3(xN(end),yN(end),zN(end),'rs','MarkerFaceColor','r','MarkerSize',7);
    plot3(xS(end),yS(end),zS(end),'bs','MarkerFaceColor','b','MarkerSize',7);

    axis equal;
    grid on;

    xlabel('$x/R$','Interpreter','latex','FontSize',16);
    ylabel('$y/R$','Interpreter','latex','FontSize',16);
    zlabel('$z/R$','Interpreter','latex','FontSize',16);

    title('Vortex trajectories on the sphere','Interpreter','latex','FontSize',18);

    set(gca,'TickLabelInterpreter','latex','FontSize',14);
    view(35,25);

    legend({'sphere','North vortex','South vortex'}, ...
        'Interpreter','latex','FontSize',14, ...
        'Location','best');

    hold off;

    saveas(fig_traj, fullfile(output_folder, 'vortex_trajectories.fig'));
    saveas(fig_traj, fullfile(output_folder, 'vortex_trajectories.png'));
    close(fig_traj);

    %% Plot temporal evolution of E(t), L_z(t), and Ec/Ei
    fig_diag1 = figure('Visible','off','Position',[100 100 900 900]);

    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(vec_t, vec_Ene, 'k', 'LineWidth', 1.8);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$E(t)$ [J]','Interpreter','latex','FontSize',16);
    title('$E(t)$','Interpreter','latex','FontSize',18);
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    nexttile;
    plot(vec_t, vec_L_z, 'k', 'LineWidth', 1.8);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$L_z(t)$ [J\,s]','Interpreter','latex','FontSize',16);
    title('$L_z(t)$','Interpreter','latex','FontSize',18);
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    nexttile;
    plot(vec_t, vec_Ec, 'r', 'LineWidth', 1.8); hold on;
    plot(vec_t, vec_Ei, 'b', 'LineWidth', 1.8);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$E_{\mathrm{c}},\,E_{\mathrm{i}}$ [J]','Interpreter','latex','FontSize',16);
    title('Compressible vs incompressible','Interpreter','latex','FontSize',18);
    legend({'$E_{\mathrm{c}}$','$E_{\mathrm{i}}$'}, ...
        'Interpreter','latex','FontSize',14,'Location','best');
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    sgtitle('Real-time diagnostics','Interpreter','latex','FontSize',18);

    saveas(fig_diag1, fullfile(output_folder, 'real_time_diagnostics.fig'));
    saveas(fig_diag1, fullfile(output_folder, 'real_time_diagnostics.png'));
    close(fig_diag1);

    %% Reconstruct point-vortex-model L_z(t) and E(t) from numerical trajectories
    R    = geom.R_sphere;
    Mb   = m_b * N_b;
    qN   = +1;
    qS   = -1;
    g    = - N_a * hbar / 2;
    na   = N_a / (4*pi*R^2);

    xi = R/50;

    thetaN = vec_position_b_North(:,1);
    phiN   = vec_position_b_North(:,2);

    thetaS = vec_position_b_South(:,1);
    phiS   = vec_position_b_South(:,2);

    t = vec_t;

    phiN_u = unwrap(phiN);
    phiS_u = unwrap(phiS);

    dthetaN = gradient(thetaN, t);
    dphiN   = gradient(phiN_u, t);

    dthetaS = gradient(thetaS, t);
    dphiS   = gradient(phiS_u, t);

    Lz_N_PV = (Mb/2) * R^2 .* sin(thetaN).^2 .* dphiN + g*qN .* (1 - cos(thetaN));
    Lz_S_PV = (Mb/2) * R^2 .* sin(thetaS).^2 .* dphiS + g*qS .* (1 - cos(thetaS));
    vec_Lz_PV = Lz_N_PV + Lz_S_PV;

    T_N_PV = (Mb/4) * R^2 .* (dthetaN.^2 + sin(thetaN).^2 .* dphiN.^2);
    T_S_PV = (Mb/4) * R^2 .* (dthetaS.^2 + sin(thetaS).^2 .* dphiS.^2);
    vec_T_PV = T_N_PV + T_S_PV;

    DeltaPhi = phiN_u - phiS_u;
    arg = 2 ...
        - 2*cos(thetaN).*cos(thetaS) ...
        - 2*sin(thetaN).*sin(thetaS).*cos(DeltaPhi);

    arg = max(arg, eps);

    chi_dip = 0.5 * log(arg);

    prefE = 2 * hbar^2 * pi * na / m_a;

    E_const  = prefE * log(R/xi);
    vec_Edip = E_const + prefE * chi_dip;

    vec_E_PV = vec_T_PV + vec_Edip;

    % force column shape
    vec_Lz_PV = vec_Lz_PV(:);
    vec_E_PV  = vec_E_PV(:);

    %% Comparison plots
    fig_diag2 = figure('Visible','off','Position',[100 100 1000 800]);
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    plot(t, vec_Ene, 'k', 'LineWidth', 1.8); hold on;
    plot(t, vec_E_PV, '--r', 'LineWidth', 1.6);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$E$ [J]','Interpreter','latex','FontSize',16);
    title('Energy','Interpreter','latex','FontSize',18);
    legend({'numerical','point-vortex from trajectory'}, ...
        'Interpreter','latex','FontSize',13,'Location','best');
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    nexttile;
    plot(t, vec_L_z, 'k', 'LineWidth', 1.8); hold on;
    plot(t, vec_Lz_PV, '--b', 'LineWidth', 1.6);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$L_z$ [J\,s]','Interpreter','latex','FontSize',16);
    title('Angular momentum','Interpreter','latex','FontSize',18);
    legend({'numerical','point-vortex from trajectory'}, ...
        'Interpreter','latex','FontSize',13,'Location','best');
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    nexttile;
    plot(t, vec_Ene - vec_E_PV, 'k', 'LineWidth', 1.6);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$\Delta E$ [J]','Interpreter','latex','FontSize',16);
    title('$E_{\mathrm{num}}-E_{\mathrm{PV}}$','Interpreter','latex','FontSize',18);
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    nexttile;
    plot(t, vec_L_z - vec_Lz_PV, 'k', 'LineWidth', 1.6);
    grid on;
    xlabel('$t$ [s]','Interpreter','latex','FontSize',16);
    ylabel('$\Delta L_z$ [J\,s]','Interpreter','latex','FontSize',16);
    title('$L_{z,\mathrm{num}}-L_{z,\mathrm{PV}}$','Interpreter','latex','FontSize',18);
    set(gca,'TickLabelInterpreter','latex','FontSize',14);

    sgtitle('Diagnostics vs reconstructed point-vortex model', ...
        'Interpreter','latex','FontSize',18);

    saveas(fig_diag2, fullfile(output_folder, 'pv_comparison_diagnostics.fig'));
    saveas(fig_diag2, fullfile(output_folder, 'pv_comparison_diagnostics.png'));
    close(fig_diag2);

    fprintf('Run_real_time_dynamics_residual_gravity completed successfully in %.2f s.\n', toc(t_start));

    diary off;

catch ME
    fprintf('\nERROR in Run_real_time_dynamics_residual_gravity\n');
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    diary off;
    rethrow(ME);
end
end
