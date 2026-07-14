%% Sweep over residual gravitational acceleration
% This script starts from the converged N_b = 6000 state produced by the
% core-filling atom-number sweep and evolves it for several fractions of
% terrestrial gravity. Results are stored under output/Sweep_g_residual/.
%
% Run run_sweep_core_filling_atom_number.m first, or otherwise provide the
% required N_b = 6000 ground-state file in the expected output directory.
%
% Requirements: SSHT MATLAB library and Parallel Computing Toolbox.
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));
Check_ssht_dependency;
output_root = fullfile(repository_root, 'output');
if ~exist(output_root, 'dir')
    mkdir(output_root);
end

N_b = 6000;
g_residual_fractions = [1e-7, 1e-6, 1e-5, 1e-4];

base_output_folder = fullfile(output_root, 'Sweep_g_residual');
input_mat_file = fullfile(output_root, 'Sweep_Nb', 'N_b_6000', ...
    'Vortexon_dipole_on_sphere.mat');

if ~exist(input_mat_file, 'file')
    error(['Previously converged state not found:\n  %s\n' ...
        'Run the core-filling atom-number sweep first or place the state in output/Sweep_Nb/N_b_6000.'], ...
        input_mat_file);
end

if ~exist(base_output_folder, 'dir')
    mkdir(base_output_folder);
end

parfor k = 1:numel(g_residual_fractions)
    g_fraction = g_residual_fractions(k);

    % Folder labels such as g_residual_1e-07.
    case_label = sprintf('g_residual_%1.0e', g_fraction);
    output_folder = fullfile(base_output_folder, case_label);

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    Run_real_time_dynamics_residual_gravity( ...
        N_b, g_fraction, input_mat_file, output_folder);
end
