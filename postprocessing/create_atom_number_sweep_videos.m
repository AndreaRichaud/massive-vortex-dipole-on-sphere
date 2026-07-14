%% Create videos for the core-filling atom-number sweep
% Convert the JPEG frames stored in each Temp_real directory into MPEG-4
% videos. Simulation data are not modified.
clear
close all
clc

script_directory = fileparts(mfilename('fullpath'));
repository_root = fileparts(script_directory);
addpath(genpath(fullfile(repository_root, 'src')));
output_root = fullfile(repository_root, 'output');

base_folder = fullfile(output_root, 'Sweep_Nb');

% Find all N_b_* folders
dir_list = dir(fullfile(base_folder, 'N_b_*'));
dir_list = dir_list([dir_list.isdir]);

target_duration_sec = 10;   % desired video duration
video_quality = 100;        % for MPEG-4 writer

for k = 1:numel(dir_list)

    case_folder = fullfile(base_folder, dir_list(k).name);
    temp_real_folder = fullfile(case_folder, 'Temp_real');

    if ~exist(temp_real_folder, 'dir')
        fprintf('Skipping %s: Temp_real not found\n', dir_list(k).name);
        continue
    end

    % Read jpg frames
    frame_files = dir(fullfile(temp_real_folder, '*.jpg'));

    if isempty(frame_files)
        fprintf('Skipping %s: no jpg frames found\n', dir_list(k).name);
        continue
    end

    % Sort by filename
    [~, idx] = sort({frame_files.name});
    frame_files = frame_files(idx);

    n_frames = numel(frame_files);

    % Choose frame rate to get about target_duration_sec
    fps = max(1, round(n_frames / target_duration_sec));

    output_video_file = fullfile(case_folder, sprintf('%s.mp4', dir_list(k).name));

    fprintf('Creating video for %s\n', dir_list(k).name);
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
