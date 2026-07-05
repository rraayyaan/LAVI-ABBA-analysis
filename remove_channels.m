% Removes a hardcoded channel from copy of LAVI/PINK .mat file
% Channel to remove is set manually (no thresholding)

clear; close all; clc;

% Define path
DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Define parameters
animal_id = 'r16';
condition = '1_habituation';
area = 'HPC';
brain_state = 'REM';

% Channel to remove
removechannel = 3;

% Create a copy - do not overwrite original
overwrite_original = false;

% Construct paths
lavi_folder = fullfile([animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['LAVI processing (' area ' ' brain_state ')']);

inputFileName = fullfile(DATApath, lavi_folder, ...
    sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

if overwrite_original
    outputFileName = inputFileName;
else
    outputFileName = fullfile(DATApath, lavi_folder, ...
        sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', animal_id, condition, area, brain_state));
end

% Loads file
fprintf('Loading: %s\n', inputFileName);
S = load(inputFileName);

fprintf('File contains:\n');
disp(fieldnames(S));

LAVI_matrix = S.LAVI_matrix; % chan x freq x epoch
PINK_matrix = S.PINK_matrix; % rep x freq x chan x epoch

fprintf('Original LAVI_matrix size: %s\n', mat2str(size(LAVI_matrix)));
fprintf('Original PINK_matrix size: %s\n', mat2str(size(PINK_matrix)));

numChans = size(LAVI_matrix, 1);

if removechannel < 1 || removechannel > numChans
    error('removechannel (%d) is out of range 1-%d', removechannel, numChans);
end

fprintf('Removing channel: %d of %d\n', removechannel, numChans);

% Builds keep mask for channels
keepChan = true(1, numChans);
keepChan(removechannel) = false;

% Removes the channel from both matrices
LAVI_matrix = LAVI_matrix(keepChan, :, :); % channel dimension = 1
PINK_matrix = PINK_matrix(:, :, keepChan, :); % channel dimension = 3

fprintf('Cut LAVI_matrix size: %s\n', mat2str(size(LAVI_matrix)));
fprintf('Cut PINK_matrix size: %s\n', mat2str(size(PINK_matrix)));

% Saves new file with the same variable names so downstream code works unchanged
save(outputFileName, 'LAVI_matrix', 'PINK_matrix');
fprintf('Saved copy: %s\n', outputFileName);