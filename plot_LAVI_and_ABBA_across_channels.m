% Plots LAVI and ABBA results for all channels of a specified animal for a single specified epoch
% Generates a single plot with the LAVI profiles of all channels for one specified epoch
% Generates ABBA plots for all channels for one specified epoch - ABBA results for each epoch are shown as a separate subplot
% Use these plots to identify unviable channels for exclusion

clear; close all; clc;

% Add LAVI toolbox
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
addpath(LAVIpath);

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Define dataset
animal_id = 'r16';
condition = '1_habituation';
area = 'HPC';
brain_state = 'REM';

foi = logspace(log10(1), log10(40), 96); % frequencies of interest

% Loads per-epoch data
lavi_folder = fullfile([animal_id ' matrices'], ...
    [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
    ['LAVI processing (' area ' ' brain_state ')']);

% Loads non-averaged data
dataFileName = fullfile(DATApath, lavi_folder, sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

% Loads the file
loaded_struct = load(dataFileName);

% Gets the field names
fieldNames = fieldnames(loaded_struct);
LAVI_matrix = loaded_struct.LAVI_matrix; % channel x freq x epochs
PINK_matrix = loaded_struct.PINK_matrix; % rep x freq x channel x epochs

% SELECT EPOCH TO PLOT
epoch_idx = 1;

LAVI_epoch = squeeze(LAVI_matrix(:, :, epoch_idx)); % channel x freq
PINK_epoch = squeeze(PINK_matrix(:, :, :, epoch_idx)); % rep x freq x channel

% Plots LAVI for the specified epoch
figure(646); clf; hold on
set(gcf,'position',[680 400 800 600]);
plot(foi, LAVI_epoch);
set(gca,'xscale','log','xtick',[2:2:10,20:10:foi(end)]);
xlim([1, 40]) % match frequencies of interest
ylim([0,1]) % match LAVI value limits

% Creates default labels for the legend based on the number of channels
numChannels = size(LAVI_epoch, 1);
defaultLabels = arrayfun(@(x) sprintf('Channel %d', x), 1:numChannels, 'UniformOutput', false);
legend(defaultLabels, 'Location', 'northeastoutside');

% Axes
xlabel('Frequency (Hz)');
ylabel('LAVI value');
title(sprintf('%s-%s-%s-%s | Epoch %d', animal_id, condition, upper(area), upper(brain_state), epoch_idx));

% Calculates LAVI of the pink-noise simulations (null reference)
 pink = permute(PINK_epoch,[3,2,1]);
 sig_lim = cat(3,min(pink,[],3), max(pink,[],3));
 alpha_range = [6 8];
 [borders2,~,sigVect2] = ABBA(LAVI_epoch, foi, alpha_range, sig_lim, 0);

% Saves significance vector for this epoch to the same file
% Uses epoch-specific field name to avoid overwriting other epochs
sigVect2_fieldname = sprintf('sigVect2_epoch%d', epoch_idx);
S.(sigVect2_fieldname) = sigVect2;
save(dataFileName, '-struct', 'S', '-append');

% Plots LAVI and ABBA for the selected epoch
% Each channel is plotted in a separate subplot
    % Black: LAVI profile of the selected epoch
    % Pink: LAVI profile of pink simulations for the selected epoch
    % Green: significant sustained bands
    % Blue: significant transient bands

figure(456); clf; hold on
cols = 3;
rows = ceil(size(LAVI_epoch, 1) / cols);
set(gcf, 'position', [580 480 560 120 * rows]);
cachol = [83 174 244] / 255; % blue
varod = [230 89 106] / 255; % pink
yarok = [108 192 12] / 255; % green

for chi = 1:size(LAVI_epoch,1)
    subplot(rows,cols,chi); hold on
    if exist('pink')==1 && ~isempty(PINK_epoch) % plot pink simulations
        plot(foi, squeeze(PINK_epoch(:,:,chi)),'color',ones(1,3)*0.5,'color',varod);
    end

    plot(foi, LAVI_epoch(chi,:),'k','linewidth',1.5);
    ylim([0.1 0.8])
    posind = sigVect2{chi}>0;
    negind = sigVect2{chi}<0;
    scatter(foi(posind), ones(1,sum(posind))*0.1,35,yarok,'fill','s');
    scatter(foi(negind), ones(1,sum(negind))*0.1,35,cachol,'fill','s');
    
    set(gca,'xscale','log','xtick',[2:2:10,20:10:foi(end)])
    if mod(chi,cols)==1, ylabel('LAVI'); end
    xlabel('Frequency (Hz)')
    title(sprintf('Channel %u',chi))
end

sgtitle(sprintf('%s-%s-%s-%s | Epoch %d', animal_id, condition, upper(area), upper(brain_state), epoch_idx));
