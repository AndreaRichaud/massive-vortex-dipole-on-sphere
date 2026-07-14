%% Sweep over the atom number of the core-filling component
% For every N_b value, this script computes a pinned imaginary-time state
% and then performs the unpinned real-time evolution of the massive vortex
% dipole. Results are stored under output/Sweep_Nb/N_b_<value>/.
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

N_b_values = 1000:1000:10000;
base_output_folder = fullfile(output_root, 'Sweep_Nb');

if ~exist(base_output_folder, 'dir')
    mkdir(base_output_folder);
end

parfor k = 1:numel(N_b_values)
    N_b = N_b_values(k);
    output_folder = fullfile(base_output_folder, sprintf('N_b_%d', N_b));

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    Compute_ground_state_on_sphere(N_b, output_folder);
    Run_real_time_dynamics_sphere(N_b, output_folder);
end
