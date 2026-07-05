% Plots sustained bands, transient bands, and no-bands from grand-averaged data.

clear; close all; clc;

% Paths
DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
addpath(LAVIpath);

% Parameters
animal_ids = {'r14','r16','r19','r20'};
condition  = '1_habituation';
area       = 'PFC'; % A1, BLA, HPC, PFC
states     = {'AW', 'QW', 'NREM', 'REM'};

num_states = numel(states);
foi = logspace(log10(1), log10(40), 96); % frequencies of interest - log-spaced 1-40 Hz
num_frequencies = numel(foi);
alpha_range = [6 8]; % alpha band - used only for naming bands after band detection

% Colours
colors = {[0.8, 0.6, 0.0], ...   % AW
          [0.0, 0.6, 0.0], ...   % QW
          [1.0, 0.0, 0.0], ...   % NREM
          [0.0, 0.0, 1.0]};      % REM

% Storage:
% Each cell holds a [1 x num_frequencies] vector (channel-averaged) for one state/animal combination
% Cell arrays (rather than plain matrices) let entries stay empty when data is missing for a given state/animal
sustained_pct = cell(num_states, numel(animal_ids));
transient_pct = cell(num_states, numel(animal_ids));
noband_pct    = cell(num_states, numel(animal_ids));

% For every state/animal combination:
% Loads data and classifies each frequency for every epoch/channel as sustained, transient, or no-band.
for state_idx = 1:num_states
    for animal_idx = 1:numel(animal_ids)
        animal_id = animal_ids{animal_idx};
        
        % Builds path to LAVI processing folder
        lavi_folder = fullfile([animal_id ' matrices'], ...
            [animal_id ' ' condition ' matrices'], area, [area ' ' states{state_idx}], ...
            ['LAVI processing (' area ' ' states{state_idx} ')']);
        file_cut = fullfile(DATApath, lavi_folder, ...
            sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', animal_id, condition, area, states{state_idx}));
        
        % Skip if the file doesn't exist
        if ~isfile(file_cut)
            warning('Missing CUT file: %s', file_cut);
            continue;
        end

        % Loads LAVI (rhythmicity profile) and PINK (null references) matrices
        loaded_file = load(file_cut);
        if ~isfield(loaded_file, 'LAVI_matrix') || ~isfield(loaded_file, 'PINK_matrix')
            warning('Missing required variables in file: %s', file_cut);
            continue;
        end

        LAVI_matrix = loaded_file.LAVI_matrix;   % [chan x freq x epoch]
        PINK_matrix = loaded_file.PINK_matrix;   % [rep x freq x chan x epoch]
        [num_chan_lavi, num_freq_lavi, num_epochs_lavi] = size(LAVI_matrix);
        [num_rep_pink, num_freq_pink, num_chan_pink, num_epochs_pink] = size(PINK_matrix);
        
        % Sanity check - LAVI and PINK matrices must agree on freq/chan/epoch dimensions
        if num_freq_lavi ~= num_freq_pink || num_chan_lavi ~= num_chan_pink || num_epochs_lavi ~= num_epochs_pink
            warning('Skipping file due to dimension mismatch: %s', file_cut);
            continue;
        end
        
        % Per-channel and per-frequency percentage of epochs classified as sustained/transient/no-band
        chan_sustained = nan(num_chan_lavi, num_freq_lavi);
        chan_transient = nan(num_chan_lavi, num_freq_lavi);
        chan_noband    = nan(num_chan_lavi, num_freq_lavi);
        
        for ch = 1:num_chan_lavi
        % Per-epoch boolean flag: 1 = band present at that frequency for this channel
            sustained_ep = nan(num_freq_lavi, num_epochs_lavi);
            transient_ep = nan(num_freq_lavi, num_epochs_lavi);
            noband_ep    = nan(num_freq_lavi, num_epochs_lavi);
            
            for ep = 1:num_epochs_lavi
                lavi_vec = squeeze(LAVI_matrix(ch, :, ep));
                if all(isnan(lavi_vec))
                    continue; % no valid data for this channel/epoch
                end

                % Builds surrogate significance limits (min/max across repetitions) for this channel/epoch
                % Then detects oscillatory bands using ABBA
                pink = squeeze(PINK_matrix(:, :, ch, ep));
                pink = permute(pink, [3 2 1]); % chan_freq_rep
                sig_lim = cat(3, min(pink,[],3), max(pink,[],3));
                [borders, ~, ~] = ABBA(lavi_vec, foi, alpha_range, sig_lim, 0);
                sust_mask = false(1, num_freq_lavi);
                tran_mask = false(1, num_freq_lavi);
                
                if ~isempty(borders{1})
                    bands = borders{1};
                
                    % Each row describes one detected band:
                    % Columns 1 and 2: start and end frequency index
                    % Column 9: direction
                        % 1 = sustained
                        % -1 = transient
                    % Column 11: significance
                    for b = 1:size(bands,1)
                        BegI = bands(b,1);
                        EndI = bands(b,2);
                        Dir  = bands(b,9);
                        Sig  = bands(b,11);
                        if Sig ~= 1
                            continue; % skip non-significant bands
                        end
                    
                        if Dir == 1
                            sust_mask(BegI:EndI) = true;
                        
                        elseif Dir == -1
                            tran_mask(BegI:EndI) = true;
                        end
                    end
                end

                % Frequencies not covered by any significant band are classified as no-band
                noband_mask = ~(sust_mask | tran_mask);
                sustained_ep(:, ep) = sust_mask;
                transient_ep(:, ep) = tran_mask;
                noband_ep(:, ep)    = noband_mask;
            end

            % Convert per-epoch flags into a percentage of epochs per frequency averaged over epochs
            chan_sustained(ch, :) = mean(sustained_ep, 2, 'omitnan')' * 100;
            chan_transient(ch, :) = mean(transient_ep, 2, 'omitnan')' * 100;
            chan_noband(ch, :)    = mean(noband_ep, 2, 'omitnan')' * 100;
        end

        % Averages across channels to get one value per frequency for this state/animal combination
        sustained_pct{state_idx, animal_idx} = mean(chan_sustained, 1, 'omitnan');
        transient_pct{state_idx, animal_idx} = mean(chan_transient, 1, 'omitnan');
        noband_pct{state_idx, animal_idx}    = mean(chan_noband, 1, 'omitnan');
    end
end

% Grand average across animals
% For each state stacks all per-frequency vectors for each animal and average
% Ignores animals with missing data
sustained_mean = nan(num_states, num_frequencies);
transient_mean = nan(num_states, num_frequencies);
noband_mean    = nan(num_states, num_frequencies);
for state_idx = 1:num_states
    tempS = nan(numel(animal_ids), num_frequencies);
    tempT = nan(numel(animal_ids), num_frequencies);
    tempN = nan(numel(animal_ids), num_frequencies);
    for animal_idx = 1:numel(animal_ids)
        if ~isempty(sustained_pct{state_idx, animal_idx})
            tempS(animal_idx, :) = sustained_pct{state_idx, animal_idx};
        end
        if ~isempty(transient_pct{state_idx, animal_idx})
            tempT(animal_idx, :) = transient_pct{state_idx, animal_idx};
        end
        if ~isempty(noband_pct{state_idx, animal_idx})
            tempN(animal_idx, :) = noband_pct{state_idx, animal_idx};
        end
    end

    sustained_mean(state_idx, :) = mean(tempS, 1, 'omitnan');
    transient_mean(state_idx, :) = mean(tempT, 1, 'omitnan');
    noband_mean(state_idx, :)    = mean(tempN, 1, 'omitnan');
end


% Figure 1: all states shown on one plot
% Sustained/transient/no-band each shown as a separate line
figure('Position', [100 100 900 650]);
hold on;
set(gca, 'XScale', 'log');
box on;
legendHandles = gobjects(num_states*3,1);
legendNames = cell(num_states*3,1);
idx = 1;

for state_idx = 1:num_states
    hS = semilogx(foi, sustained_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', colors{state_idx}, 'LineStyle', '-', ...
        'DisplayName', [upper(states{state_idx}) ' (Sustained)']);
    hT = semilogx(foi, transient_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', 'k', 'LineStyle', ':', ...
        'DisplayName', [upper(states{state_idx}) ' (Transient)']);
    hN = semilogx(foi, noband_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', ...
        'DisplayName', [upper(states{state_idx}) ' (No-band)']);
    legendHandles(idx) = hS;
    legendNames{idx} = get(hS, 'DisplayName');
    idx = idx + 1;
    legendHandles(idx) = hT;
    legendNames{idx} = get(hT, 'DisplayName');
    idx = idx + 1;
    legendHandles(idx) = hN;
    legendNames{idx} = get(hN, 'DisplayName');
    idx = idx + 1;
end

xlabel('Frequency (Hz)', 'FontSize', 14);
set(gca, 'XTick', [1 2 4 6 8 10 20 40], ...
    'XTickLabel', {'1', '2', '4', '6', '8', '10', '20', '40'}, 'FontSize', 12);
ylabel('Relative occurrence (%)', 'FontSize', 12);
title(sprintf('Sustained vs Transient vs No-band Across Brain States: %s', area), 'FontSize', 14);
ylim([0 105]);
xlim([0.8 42]);
grid on;
hold off;

% Figure 2: same data shown with one subplot per state
figure('Position', [100 100 1200 800]);
for state_idx = 1:num_states
    ax = subplot(2, 2, state_idx);
    hold(ax, 'on');
    set(ax, 'XScale', 'log');
    box(ax, 'on');
    semilogx(ax, foi, sustained_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', colors{state_idx}, 'LineStyle', '-', ...
        'DisplayName', 'Sustained');
    semilogx(ax, foi, transient_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', 'k', 'LineStyle', ':', ...
        'DisplayName', 'Transient');
    semilogx(ax, foi, noband_mean(state_idx, :), ...
        'LineWidth', 2, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', ...
        'DisplayName', 'No-band');
    xlabel(ax, 'Frequency (Hz)', 'FontSize', 12);
    set(ax, 'XTick', [1 2 4 6 8 10 20 40], ...
        'XTickLabel', {'1', '2', '4', '6', '8', '10', '20', '40'}, 'FontSize', 10);
    ylabel(ax, 'Relative occurrence (%)', 'FontSize', 12);
    title(ax, sprintf('%s', upper(states{state_idx})), 'FontSize', 14, 'FontWeight', 'bold');
    ylim(ax, [0 105]);
    xlim(ax, [0.8 42]);
    grid(ax, 'on');
    set(ax, 'FontSize', 10);
    % legend(ax, 'show', 'Location', 'best', 'FontSize', 10);
    hold(ax, 'off');
end

sgtitle(sprintf('Sustained vs Transient vs No-band Across Brain States: %s', area), ...
    'FontSize', 16, 'FontWeight', 'bold');
