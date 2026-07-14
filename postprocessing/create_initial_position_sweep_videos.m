%% Create videos for the initial-position sweep
% Convert the JPEG frames stored in each Temp_real directory into MPEG-4
% videos. Simulation data are not modified.
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

target_duration_sec = 10;   % desired video duration
video_quality = 100;        % for MPEG-4 writer

for k = 1:numel(theta_1_values_all)

    theta_0 = theta_1_values_all(k);

    case_folder = fullfile(base_output_folder, sprintf('theta1_%0.5f', theta_0));
    temp_real_folder = fullfile(case_folder, 'Temp_real');

    if ~exist(temp_real_folder, 'dir')
        fprintf('Skipping theta_0 = %.5f: Temp_real not found\n', theta_0);
        continue
    end

    % Read jpg frames
    frame_files = dir(fullfile(temp_real_folder, '*.jpg'));

    if isempty(frame_files)
        fprintf('Skipping theta_0 = %.5f: no jpg frames found\n', theta_0);
        continue
    end

    % Sort by filename
    [~, idx] = sort({frame_files.name});
    frame_files = frame_files(idx);

    n_frames = numel(frame_files);

    % Choose frame rate to get about target_duration_sec
    fps = max(1, round(n_frames / target_duration_sec));

    output_video_file = fullfile(case_folder, sprintf('theta1_%0.5f.mp4', theta_0));

    fprintf('Creating video for theta_0 = %.5f\n', theta_0);
    fprintf('  Frames: %d\n', n_frames);
    fprintf('  FPS: %d\n', fps);
    fprintf('  Output: %s\n', output_video_file);

    v = VideoWriter(output_video_file, 'MPEG-4');
    v.FrameRate = fps;
    v.Quality = video_quality;
    open(v);

    for j = 1:n_frames
        img_path = fullfile(temp_real_folder, frame_files(j).name);
        img = imread(img_path);
        writeVideo(v, img);
    end

    close(v);

    actual_duration = n_frames / fps;
    fprintf('  Done. Video duration: %.2f s\n\n', actual_duration);
end
