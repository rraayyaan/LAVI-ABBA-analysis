% Generates ABBA plots for all epochs of one specified channel
% ABBA results for each epoch are shown as a separate subplot

clear; close all; clc;

% Add LAVI toolbox
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
addpath(LAVIpath);

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Define dataset
animal_id = 'r14';
condition = '1_habituation';
area = 'A1';
brain_state = 'AW';

foi = logspace(log10(1), log10(40), 96); % frequencies of interest

% SELECT CHANNEL TO PLOT
channel_to_plot = 3;

% Load per-epoch data
dataFileName = fullfile(DATApath, [animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['LAVI processing (' area ' ' brain_state ') THRESH BOTH'], ...
    sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', animal_id, condition, area, brain_state));
loaded_struct = load(dataFileName);
LAVI_matrix = loaded_struct.LAVI_matrix; % chan x freq x epoch
PINK_matrix = loaded_struct.PINK_matrix; % rep x freq x chan x epoch
num_epochs = size(LAVI_matrix, 3);

% Colours
cachol = [83 174 244] / 255; % blue
varod  = [230 89 106] / 255; % pink
yarok  = [108 192 12] / 255; % green
figure; clf;
cols = 3;
rows = ceil(num_epochs / cols);
set(gcf, 'position', [500 300 600 150 * rows]);
for ep = 1:num_epochs
    
    % Extracts single channel and single epoch
    LAVI_ep = squeeze(LAVI_matrix(channel_to_plot, :, ep)); % 1 x freq
    PINK_ep = squeeze(PINK_matrix(:, :, channel_to_plot, ep)); % rep x freq
    
    % Prepares for ABBA
    pink = permute(PINK_ep, [3 2 1]); % chan x freq x rep
    sig_lim = cat(3, min(pink,[],3), max(pink,[],3));
    
    % alpha_range: the frequency range in which we expact to find alpha
    % The band of the peak in this band (alpha) will be assigned the index 0
    % Bands with lower frequency than alpha band will be assigned with negative indices
    % Bands with higher frequency than alpha band will be assigned with positive indices
    % alpha_range is only used after band detection for relative numbering of identified bands
    alpha_range = [6 8];
    
    % Runs ABBA to get significance vector for this epoch
    [~,~,sigVect] = ABBA(LAVI_ep, foi, alpha_range, sig_lim, 0);
    
    subplot(rows, cols, ep); hold on
    
    % Plots pink noise
    plot(foi, PINK_ep', 'color', varod);
    
    % Plots LAVI profiles
    plot(foi, LAVI_ep, 'k', 'linewidth', 1.5);
    
    % Significant points
    posind = sigVect{1} > 0;
    negind = sigVect{1} < 0;
    
    scatter(foi(posind), ones(1,sum(posind))*0.1, 35, yarok, 'fill', 's');
    scatter(foi(negind), ones(1,sum(negind))*0.1, 35, cachol, 'fill', 's');
    
    set(gca,'xscale','log','xtick',[2:2:10,20:10:foi(end)])
    ylim([0.1 0.8])
    xlim([1 40])
    
    title(sprintf('Epoch %d (Ch %d)', ep, channel_to_plot));
    xlabel('Frequency (Hz)');
    if mod(ep, cols) == 1
        ylabel('LAVI');
    end
end
