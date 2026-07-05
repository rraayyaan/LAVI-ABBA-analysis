% Removes 50Hz artifacts from LFP epochs
% Run this immediately after the python epoching script
% The script assumes FieldTrip toolbox is on the MATLAB path

clear; close all; clc;

ft_defaults;  % initialize FieldTrip

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices';

% Setup paths
animal_id = 'r19';
condition = '1_habituation';
area = 'A1';
brain_state = 'NREM';

epoch_folder = fullfile(DATApath, [animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['4-minute epochs (' area ' ' brain_state ')']);

% Sampling frequency
fs = 1000;

% Define artifact frequencies
artF = [50, 100, 150, 200, 250]; % 50 Hz plus harmonics 
rband = [0.5 0.5 0.5 1 1]; % band width to be replaced
nband = [0.5 0.5 1 1 2]; % neighbouring frequency band

% List all matrix files
mat_files = dir(fullfile(epoch_folder, sprintf('%s_%s_matrix_*.mat', animal_id, condition)));
fprintf('Found %d epoch matrices to process\n', length(mat_files));

% Process each epoch matrix
for i = 1:length(mat_files)
    fname = fullfile(epoch_folder, mat_files(i).name);
    fprintf('Processing %s...\n', mat_files(i).name);
    
    % Load the epoch matrix
    data = load(fname);
    matrix_orig = data.matrix;
    
    % Prepare for FieldTrip format - single trial with all channels
    MEG = struct();
    MEG.trial{1} = matrix_orig; % data matrix (channels x time)
    MEG.time{1} = (0:size(matrix_orig,2)-1)/fs; % time vector
    MEG.fsample = fs; % sampling frequency
    MEG.label = cellstr(strcat('chan', string(1:size(matrix_orig,1))'));  % channel labels
    
    % Apply DFT filter to remove 50Hz artifacts using spectral interpolation
    filt = ft_preproc_dftfilter(MEG.trial{1}, fs, artF, ...
        'dftreplace', 'neighbour_fft', ...
        'dftbandwidth', rband, ...
        'dftneighbourwidth', nband);
    
    % Update data structure with filtered signal
    MEG2 = MEG;
    MEG2.trial{1} = filt;
    
    % Save filtered version (original + filtered)
    output_fname = strrep(mat_files(i).name, '.mat', '_filtered.mat');
    save(fullfile(epoch_folder, output_fname), 'MEG2', 'matrix_orig', 'fs');
    fprintf('  Saved filtered data: %s\n', output_fname);
end

fprintf('\n=== Processing complete ===\n');
fprintf('All epochs now have 50Hz + harmonics removed via spectral interpolation.\n');
fprintf('Filtered files end with "_filtered.mat" and contain MEG2 structure.\n');
