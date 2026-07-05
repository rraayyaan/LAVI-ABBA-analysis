% Generates a single plot showing ABBA detected bands (without LAVI profiles) for all epochs of a specified channel
% Each row shows the ABBA bands from a single epoch
% One row plots the bands detected using the averaged (across all epochs) LAVI profile

clear; close all; clc;

% Add LAVI toolbox
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
addpath(LAVIpath);

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Define dataset
animal_id = 'r16';
condition = '1_habituation';
area = 'BLA';
brain_state = 'REM';

foi = logspace(log10(1), log10(40), 96); % frequencies of interest

% Set dataset path
lavi_folder = fullfile([animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['LAVI processing (' area ' ' brain_state ')']);

% Loads non-averaged data
dataFileName = fullfile(DATApath, lavi_folder, ...
    sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

loaded_struct = load(dataFileName);

LAVI_matrix = loaded_struct.LAVI_matrix; % chan x freq x epochs
PINK_matrix = loaded_struct.PINK_matrix; % rep x freq x chan x epochs

% SELECT CHANNEL TO PLOT
channel_idx = 1; % change to specify the channel

num_epochs = size(LAVI_matrix, 3);

% alpha_range is the frequency range in which we expact to find alpha
% The band of the peak in this band (alpha) will be assigned the index 0
% Bands with lower frequency than alpha band will be assigned with negative indices
% Bands with higher frequency than alpha band will be assigned with positive indices
% alpha_range is only used after band detection for relative numbering of identified bands
alpha_range = [6 8];

% Preallocates array to store the significance vectors from each epoch
sigVect_all_epochs = cell(1, num_epochs);

for ep = 1:num_epochs
    
    % Extracts the LAVI values for the current epoch and channel
    LAVI_ep = LAVI_matrix(channel_idx, :, ep); % size: 1 (channel) x freq
    
    % Extracts the pink noise surrogates for the current epoch and channel
    PINK_ep = PINK_matrix(:, :, channel_idx, ep); % size: rep x freq (channel, epoch)

    % Reshapes PINK_ep from (rep x freq) to (channel x freq x rep) so it matches the format ABBA expects
    pink_ep = permute(PINK_ep, [3, 2, 1]); % size: 1 (channel) x freq x rep
    
    % Computes min and max across the number of pink-noise simulations (rep)
    sig_lim_ep = cat(3, min(pink_ep, [], 3), max(pink_ep, [], 3)); % size: 1 (channel) x freq x 2 (min, max)
    
    % Runs ABBA to get significance vector for this epoch
    [~, ~, sigVect_ep] = ABBA(LAVI_ep, foi, alpha_range, sig_lim_ep, 0);

    % Stores the significance vector for this epoch
    sigVect_all_epochs{ep} = sigVect_ep{1}; % size: 1 (channel) x freq
end

% Averages LAVI across epochs
% Uses ABBA for the epoch-averaged data
LAVI_avg = mean(LAVI_matrix, 3); % size: channel x freq

% Averages pink-noise surrogates across epochs
PINK_avg = mean(PINK_matrix, 4); % size: rep x freq x channel

% Extracts the averaged LAVI values for the current channel
LAVI_avg_ch = LAVI_avg(channel_idx, :); % size: 1 (channel) x freq

% Extracts the averaged pink-noise simulations for the current channel
PINK_avg_ch = PINK_avg(:, :, channel_idx); % size: rep x freq

% Reshapes PINK_avg_ch from (rep x freq) to (channel x freq x rep) so it matches the format ABBA expects
pink_avg = permute(PINK_avg_ch, [3, 2, 1]); % size: 1 (channel) x freq x rep

% Builds the significance envelope - min and max across pink noise simulations (rep) at each frequency
sig_lim_avg = cat(3, min(pink_avg, [], 3), max(pink_avg, [], 3)); % size: 1 (channel) x freq x 2 (min, max)

% Runs ABBA to get the significance vector for the epoch-averaged data 
[~, ~, sigVect_avg_cell] = ABBA(LAVI_avg_ch, foi, alpha_range, sig_lim_avg, 0);

% Stores the significance vector for the epoch-averaged data
sigVect_avg = sigVect_avg_cell{1}; % size: 1 (channel) x freq

% Plots bands for all epochs and average
% Colours
cachol = [83 174 244] / 255; % blue (transient)
yarok = [108 192 12] / 255; % green (sustained)

y_epoch = 1:num_epochs;
y_avg   = num_epochs + 0.5;

figure(1); clf; hold on
set(gcf, 'position', [500, 100, 900, 120 * (num_epochs + 1)]);

% Plot bands for each epoch
for ep = 1:num_epochs
    sigVect_ep = sigVect_all_epochs{ep};
    
    posind = sigVect_ep > 0;
    negind = sigVect_ep < 0;
    
    scatter(foi(posind), ones(1, sum(posind)) * y_epoch(ep), 40, yarok, 'filled', 's');
    scatter(foi(negind), ones(1, sum(negind)) * y_epoch(ep), 40, cachol, 'filled', 's');
end

% Plot bands for average (at the bottom)
posind_avg = sigVect_avg > 0;
negind_avg = sigVect_avg < 0;

scatter(foi(posind_avg), ones(1, sum(posind_avg)) * y_avg, 40, yarok, 'filled', 's', ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);
scatter(foi(negind_avg), ones(1, sum(negind_avg)) * y_avg, 40, cachol, 'filled', 's', ...
    'MarkerFaceAlpha', 0.7, 'MarkerEdgeAlpha', 0.7);

% Axes
set(gca, 'xscale', 'log', 'xtick', [2:2:10, 20:10:foi(end)]);
xlim([1, 40]);
yticks([1:num_epochs, y_avg]);
yticklabels([arrayfun(@(x) sprintf('Epoch %d', x), 1:num_epochs, 'UniformOutput', false), 'Average']);
ylabel('Epoch');
xlabel('Frequency (Hz)');
title(sprintf('%s-%s-%s-%s | Channel %d | Significant bands (all epochs and average)', ...
    animal_id, condition, upper(area), upper(brain_state), channel_idx));

grid on;

h_pos = scatter(foi(1), 0, 40, yarok, 'filled', 's');
h_neg = scatter(foi(1), 0, 40, cachol, 'filled', 's');
legend([h_pos, h_neg], {'Sustained)', 'Transient)'}, ...
    'Location', 'northeastoutside');
