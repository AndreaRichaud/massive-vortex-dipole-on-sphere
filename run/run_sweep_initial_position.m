%% Sweep over the initial vortex position
% At fixed N_b, this script computes and evolves massive vortex-dipole
% states for several initial polar angles. Results are stored under
% output/Sweep_theta_initial/theta1_<value>/.
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

N_b = 3000;

% 10 values of theta_1 between 0.1 and pi/2 - 0.1
theta_1_values = linspace(0.1, pi/2 - 0.1, 10);

% Base output folder
base_output_folder = fullfile(output_root, 'Sweep_theta_initial');

if ~exist(base_output_folder, 'dir')
    mkdir(base_output_folder);
end

parfor k = 1:numel(theta_1_values)

    theta_1 = theta_1_values(k);

    % Folder name
    output_folder = fullfile(base_output_folder, ...
        sprintf('theta1_%0.5f', theta_1));

    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    Compute_ground_state_on_sphere_given_theta(N_b, theta_1, output_folder);
    Run_real_time_dynamics_sphere(N_b, output_folder);
end
