% Calculates LAVI values and 100 pink noise simulations (null reference) per channel across the frequencies of interest

clear; close all; clc;

% Add LAVI toolbox to your MATLAB path
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)'; % Location of the LAVI toolbox
addpath(LAVIpath); % Add the LAVI toolbox to the MATLAB path

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% LAVI parameters
foi         = logspace(log10(1), log10(40), 96);    % frequencies of interest
fsample     = 1000;                                 % sampling frequency
lag         = 1.5;                                  % lag between the signal and its copy (in cycles, default = 1.5)
width       = 5;                                    % wavelet width (in cycles, default = 5)
pink_reps   = 100;                                  % number of simulations created per channel. Default = 20.

% Number of files (epochs)
numFiles    = 6;

% Data parameters
animal_id = 'r14';
condition = '1_habituation';
area = 'A1';
brain_state = 'AW';

epoch_folder = fullfile(DATApath, [animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['4-minute epochs (' area ' ' brain_state ')']);

% Checks number of files in epoch folder
epoch_files = dir(fullfile(epoch_folder, sprintf('%s_%s_matrix_*_filtered.mat', animal_id, condition)));
if length(epoch_files) ~= numFiles
    error('Number of files in epoch folder (%d) does not match number of requested files (%d)', length(epoch_files), numFiles);
end

% Preallocate cell arrays to hold results for all epochs
all_LAVI = cell(1, numFiles);
all_PINK = cell(1, numFiles);

% Loops through each file
for fileIndex = 0:numFiles-1
    % Constructs the .mat path for each file
    dataFileName = fullfile(epoch_folder, sprintf('%s_%s_matrix_%d_filtered.mat', animal_id, condition, fileIndex));

    % Loads the .mat file
    loaded_struct = load(dataFileName);

    % Accesses the filtered data from FieldTrip structure
    myData = loaded_struct.MEG2.trial{1}; % Filtered data (channels x timepoints)
    disp(['Processing file: ', dataFileName]);

    % Checks dimensions and type
    disp(size(myData)); % displays the size of the data matrix
    disp(class(myData)); % displays the class/type of the data

    durs = size(myData, 2) / fsample;  % duration (in seconds)
    choi = 1:size(myData, 1);          % chooses channels of interest

    % Calculates LAVI of the data
    cfg = [];
    cfg.foi     = foi;
    cfg.fs      = fsample;
    cfg.lag     = lag;
    cfg.width   = width;
    cfg.verbose = 1;

    dat = myData(choi, :);
    if any(isnan(dat(:)))
        warning('Data contains NaNs. Calculating TFR by convolution in the time domain');
    end
    [LAVI, cfg] = Prepare_LAVI(cfg, dat);

    % Generates pink noise matching the data and calculates its LAVI values
    cfg = [];
    cfg.Pink_reps = pink_reps;
    cfg.durs      = durs;
    cfg.foi       = foi;
    PINK = computePinkLAVI(cfg, dat(choi, :)); % Dimord: rep_freq_chan

    % Stores in cell arrays
    all_LAVI{fileIndex+1} = LAVI;
    all_PINK{fileIndex+1} = PINK;

    % Saves results (LAVI and PINK)
    outputDir = fullfile(DATApath, [animal_id ' matrices'], ...
        [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
        ['LAVI processing (' area ' ' brain_state ')']);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    outputFile = fullfile(outputDir, sprintf('matrix%d.mat', fileIndex));

    % Initialises matrices
    numChannels = size(LAVI, 1);
    LAVI_matrix = zeros(numChannels, size(LAVI, 2));
    PINK_matrix = zeros(cfg.Pink_reps, size(LAVI, 2), numChannels);

    % Fills the matrices with LAVI and PINK data
    for chi = 1:numChannels
        LAVI_matrix(chi, :) = LAVI(chi, :);
        PINK_matrix(:, :, chi) = PINK(:, :, chi);
    end

    % Creates a structure to hold the results
    results.PINK = PINK_matrix;
    results.LAVI = LAVI_matrix;

    % Saves the results structure to a single .mat file
    save(outputFile, 'results');
    disp(['LAVI and PINK results saved successfully in: ', outputFile]);
end

% Combines results into single arrays across epochs
numChannels = size(all_LAVI{1}, 1);
LAVI_matrix = zeros(numChannels, size(all_LAVI{1}, 2), numFiles);
PINK_matrix = zeros(pink_reps, size(all_LAVI{1}, 2), numChannels, numFiles);

for i = 1:numFiles
    LAVI_matrix(:, :, i) = all_LAVI{i};
    PINK_matrix(:, :, :, i) = all_PINK{i};
end

% Defines a single output file path
outputDir_combined = fullfile(DATApath, [animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['LAVI processing (' area ' ' brain_state ')']);

if ~exist(outputDir_combined, 'dir')
    mkdir(outputDir_combined);
end

outputFile_combined = fullfile(outputDir_combined, ...
    sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

% Saves directly
save(outputFile_combined, 'LAVI_matrix', 'PINK_matrix');
disp(['All LAVI and PINK results saved successfully in: ', outputFile_combined]);
