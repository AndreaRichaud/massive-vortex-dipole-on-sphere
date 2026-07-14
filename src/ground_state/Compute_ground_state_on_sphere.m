function Compute_ground_state_on_sphere(N_b, output_folder)

%COMPUTE_GROUND_STATE_ON_SPHERE Compute a pinned stationary vortex-dipole state.
%
%   COMPUTE_GROUND_STATE_ON_SPHERE(N_b, output_folder) obtains an
%   imaginary-time stationary state of two coupled Gross-Pitaevskii
%   equations on a spherical shell. Component b fills the two vortex cores
%   of component a.
%
%   Inputs
%   ------
%   N_b
%       Atom number of the core-filling component b.
%   output_folder
%       Directory used for simulation data, diagnostic figures, frames,
%       and the execution log.
%
%   Files written
%   -------------
%   Vortexon_dipole_on_sphere.mat
%       Converged wavefunctions, geometry, and physical parameters.
%   log_compute_ground_state.txt
%       Text log of the imaginary-time calculation.
%   Temp_img/
%       Intermediate diagnostic frames.
%
%   External dependency: SSHT MATLAB library.

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    temp_img_folder = fullfile(output_folder, 'Temp_img');
    if ~exist(temp_img_folder, 'dir')
        mkdir(temp_img_folder);
    end

    log_file = fullfile(output_folder, 'log_compute_ground_state.txt');
    if exist(log_file, 'file')
        delete(log_file);
    end

    diary(log_file);
    diary on;

    t_start = tic;

    try
        fprintf('============================================================\n');
        fprintf('Compute_ground_state_on_sphere started\n');
        fprintf('N_b = %d\n', N_b);
        fprintf('Output folder: %s\n', output_folder);
        fprintf('Start time: %s\n', datestr(now));
        fprintf('============================================================\n');

        %--------------------------------------------------------------%
        % Path setup
        %--------------------------------------------------------------%
        this_file_folder = fileparts(mfilename('fullpath'));
        addpath(genpath(this_file_folder));

        %--------------------------------------------------------------%
        % Physical constants
        %--------------------------------------------------------------%
        hbar = 1.0545718e-34;      % Planck's constant / 2pi [J s]
        amu  = 1.66054e-27;        % atomic mass unit [kg]
        a0   = 5.291777e-11;       % Bohr radius [m]

        % Species masses
        m_a = 23 * amu;            % Sodium-23
        m_b = 39 * amu;            % Potassium-39

        % Atom numbers
        N_a = 500000;

        % Geometry
        R_sphere = 20e-6;          % [m]
        L_z      = 2e-6 / sqrt(2*pi);

        % Interactions
        a_aa = 52.0 * a0;
        a_bb =  7.6 * a0;

        g_a  = 4*pi*hbar^2*a_aa/m_a;
        g_b  = 4*pi*hbar^2*a_bb/m_b;
        g_ab = 3*sqrt(g_a*g_b);

        % Numerical parameters
        n_iterations = 1e4;
        dt           = 1e-6;
        L_band       = 224;

        % Build geometry
        geom = Build_geometry_sphere(R_sphere, L_band, 'MW');

        fprintf('Geometry built successfully.\n');
        fprintf('L_band = %d\n', L_band);
        fprintf('n_iterations = %d\n', n_iterations);
        fprintf('dt = %.3e\n', dt);

        %--------------------------------------------------------------%
        % Initial conditions
        %--------------------------------------------------------------%
        q_1     = +1;
        theta_1 = pi/2 - 0.4;
        phi_1   = 3*pi/4;

        q_2     = -1;
        theta_2 = pi/2 + 0.4;
        phi_2   = 3*pi/4;

        % Pinning potential
        Pinning_flag = 1;
        V0_pin_a     = 5e-31;   % [J]
        sigma_pin    = 0.04;    % [rad]

        [Mat_V_a, Mat_V_b] = Build_pinning_potential_sphere( ...
            geom, ...
            theta_1, phi_1, theta_2, phi_2, ...
            Pinning_flag, V0_pin_a, sigma_pin);

        % Initial wavefunctions
        sigma_psi_b = 1e-6;     % [m]
        xi_core_a   = 0.8e-6;   % [m]

        [psi_a, psi_b] = Build_initial_wavefunctions_sphere( ...
            geom, ...
            theta_1, phi_1, q_1, ...
            theta_2, phi_2, q_2, ...
            sigma_psi_b, xi_core_a);

        fprintf('Initial conditions generated.\n');

        %--------------------------------------------------------------%
        % Diagnostics allocation
        %--------------------------------------------------------------%
        sample_frequency = 1000;
        n_samples        = ceil(n_iterations / sample_frequency) + 1;

        vec_mu_a = zeros(n_samples,1);
        vec_mu_b = zeros(n_samples,1);
        vec_Ene  = zeros(n_samples,1);
        vec_t    = zeros(n_samples,1);

        i_sampling = 1;

        %--------------------------------------------------------------%
        % Imaginary-time evolution
        %--------------------------------------------------------------%
        for i = 1:n_iterations

            [psi_a, psi_b, Hpsi_a, Hpsi_b] = Imaginary_time_evolve_sphere( ...
                psi_a, psi_b, ...
                m_a, m_b, ...
                g_a, g_b, g_ab, ...
                N_a, N_b, ...
                L_z, ...
                Mat_V_a, Mat_V_b, ...
                geom, hbar, dt);

            if i == 1 || mod(i, sample_frequency) == 0

                elapsed_time = toc(t_start);
                fprintf('Avanzamento = %.2f %% in %.0f s\n', ...
                    100*i/n_iterations, elapsed_time);

                vec_mu_a(i_sampling) = real(sum(conj(psi_a(:)) .* Hpsi_a(:) .* geom.W(:)));
                vec_mu_b(i_sampling) = real(sum(conj(psi_b(:)) .* Hpsi_b(:) .* geom.W(:)));

                vec_Ene(i_sampling) = Compute_total_energy_sphere( ...
                    psi_a, psi_b, ...
                    m_a, m_b, ...
                    g_a, g_b, g_ab, ...
                    N_a, N_b, ...
                    L_z, ...
                    Mat_V_a, Mat_V_b, ...
                    geom, hbar);

                vec_t(i_sampling) = i * dt;

                fig = Plot_frame_sphere(psi_a, psi_b, geom.L_band, geom.method, ...
                    'Visible', 'off', ...
                    'FigureTitle', sprintf('$t = %.3f\\ \\mathrm{s}$', i*dt));

                fig_name = fullfile(temp_img_folder, sprintf('%05d.jpg', i_sampling));
                saveas(fig, fig_name);
                close(fig);

                i_sampling = i_sampling + 1;
            end
        end

        %--------------------------------------------------------------%
        % Trim unused preallocated entries
        %--------------------------------------------------------------%
        vec_mu_a = vec_mu_a(1:i_sampling-1);
        vec_mu_b = vec_mu_b(1:i_sampling-1);
        vec_Ene  = vec_Ene(1:i_sampling-1);
        vec_t    = vec_t(1:i_sampling-1);

        %--------------------------------------------------------------%
        % Final state figure
        %--------------------------------------------------------------%
        fig_final = Plot_frame_sphere(psi_a, psi_b, geom.L_band, geom.method, ...
            'Visible', 'off', ...
            'FigureTitle', sprintf('Final state at $t = %.3e\\ \\mathrm{s}$', n_iterations*dt));

        saveas(fig_final, fullfile(output_folder, 'final_state.fig'));
        saveas(fig_final, fullfile(output_folder, 'final_state.png'));
        close(fig_final);

        %--------------------------------------------------------------%
        % Diagnostics figure
        %--------------------------------------------------------------%
        fig_diag = figure('Visible', 'off', 'Position', [100, 100, 1200, 900]);

        tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        plot(vec_t, vec_mu_a, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$\mu_a$ [J]', 'Interpreter', 'latex', 'FontSize', 15);
        title('$\mu_a(t)$', 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

        nexttile;
        plot(vec_t, vec_mu_b, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$\mu_b$ [J]', 'Interpreter', 'latex', 'FontSize', 15);
        title('$\mu_b(t)$', 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

        nexttile;
        plot(vec_t, vec_Ene, 'LineWidth', 1.8);
        grid on;
        xlabel('$t$ [s]', 'Interpreter', 'latex', 'FontSize', 15);
        ylabel('$E$ [J]', 'Interpreter', 'latex', 'FontSize', 15);
        title('$E(t)$', 'Interpreter', 'latex', 'FontSize', 16);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

        sgtitle('Imaginary-time diagnostics', 'Interpreter', 'latex', 'FontSize', 18);

        saveas(fig_diag, fullfile(output_folder, 'diagnostics.fig'));
        saveas(fig_diag, fullfile(output_folder, 'diagnostics.png'));
        close(fig_diag);

        %--------------------------------------------------------------%
        % Save workspace output
        %--------------------------------------------------------------%
        output_mat_file = fullfile(output_folder, 'Vortexon_dipole_on_sphere.mat');

        save(output_mat_file, ...
            'hbar', 'amu', 'a0', ...
            'm_a', 'm_b', ...
            'N_a', 'N_b', ...
            'R_sphere', 'L_z', ...
            'a_aa', 'a_bb', ...
            'g_a', 'g_b', 'g_ab', ...
            'n_iterations', 'dt', 'L_band', ...
            'geom', ...
            'q_1', 'theta_1', 'phi_1', ...
            'q_2', 'theta_2', 'phi_2', ...
            'Pinning_flag', 'V0_pin_a', 'sigma_pin', ...
            'Mat_V_a', 'Mat_V_b', ...
            'sigma_psi_b', 'xi_core_a', ...
            'psi_a', 'psi_b', ...
            'vec_mu_a', 'vec_mu_b', 'vec_Ene', 'vec_t', ...
            'sample_frequency');

        fprintf('Saved MAT file: %s\n', output_mat_file);
        fprintf('Compute_ground_state_on_sphere completed successfully in %.2f s.\n', toc(t_start));

        diary off;

    catch ME
        fprintf('\nERROR in Compute_ground_state_on_sphere\n');
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        diary off;
        rethrow(ME);
    end
end
