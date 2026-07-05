% Computes and plots raw power averaged across all animals
% Generates one plot per brain region with all brain states

clear; close all; clc;


DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices';

% Define parameters
animals = {'r14','r16','r19','r20'};
condition = '1_habituation';
areas = {'A1', 'BLA', 'HPC', 'PFC'};
states = {'AW','QW','REM','NREM'};

fs = 1000; % sampling frequency
freq_range = [1 40]; % frequencies of interest

% PSD settings
window = 4 * fs; % Welch window length = 4-second segments - pwelch applies default Hamming window
noverlap = window/2; % 50% overlap between segments
nfft = 4096; % FFT length (next power of 2 >= window length)

% Builds the frequency vector once before loops (avoids recomputing it each iteration)
dummy_signal = zeros(1, window);
[~, fvec_full] = pwelch(dummy_signal, window, noverlap, nfft, fs); % native pwelch frequency bins
fvec_sel = logspace(log10(1), log10(40), 96)';   % 96 log-spaced frequencies at 1 and 40 Hz

% STORAGE
% Size: [n_animals x n_frequencies]
all_state_power = struct();
for s = 1:length(states)
    all_state_power.(states{s}) = struct();
    for area_idx = 1:length(areas)
        all_state_power.(states{s}).(areas{area_idx}) = [];
    end
end

% Colours for brain states
state_colours = [
    0.8, 0.6, 0.0;    % AW - dark yellow
    0.0, 0.6, 0.0;    % QW - green
    1.0, 0.0, 0.0;    % REM - red
    0.0, 0.0, 1.0     % NREM - blue
];

% Loop order: area - animal - brain state
% For each animal/area/state combination loads all 4-minute epoch files and computes a PSD per channel per epoch
% Then averages across epochs per channel, across channels per animal, and across animals.
% Result is one grand-averaged spectrum per state/area combination.
for area_idx = 1:length(areas)
    area = areas{area_idx};
    for a = 1:length(animals)
        animal_id = animals{a};
        fprintf('Processing animal: %s, area: %s\n', animal_id, area);
        for s = 1:length(states)
            brain_state = states{s};

            % Builds path to the folder containing the 4-minute epoch files
            epoch_folder = fullfile(DATApath, ...
                [animal_id ' matrices'], ...
                [animal_id ' ' condition ' matrices'], ...
                area, ...
                [area ' ' brain_state], ...
                ['4-minute epochs (' area ' ' brain_state ')']);
            if ~exist(epoch_folder, 'dir')
                warning('Missing folder: %s', epoch_folder);
                continue;
            end
            
            files = dir(fullfile(epoch_folder, '*_filtered.mat'));
            
            % Accumulates epochs per channel across all files
            % Holds one row per epoch found for a channel
            channel_epoch_psds = {};
            for f = 1:length(files)
                data = load(fullfile(epoch_folder, files(f).name));
                
                % Extracts data
                signal = data.MEG2.trial{1}; % channels x time
                n_channels = size(signal, 1);
                
                % Expands cell array to accommodate the number of channels
                if length(channel_epoch_psds) < n_channels
                    channel_epoch_psds{n_channels} = [];
                end

                for ch = 1:n_channels
                    % Computes Welch PSD for this channel/epoch
                    pxx_full = pwelch(signal(ch,:), window, noverlap, nfft, fs);
                
                    % Interpolates onto log-spaced target grid
                    pxx_interp = interp1(fvec_full, pxx_full, fvec_sel, 'pchip');
                
                    % Appends epoch as a new row for channel
                    channel_epoch_psds{ch} = [channel_epoch_psds{ch}; pxx_interp'];
                end
            end
            
            % Averages across all found epochs for each channel
            channel_means = [];
            for ch = 1:length(channel_epoch_psds)
                if isempty(channel_epoch_psds{ch})
                    continue;
                end
                
                % Averages across all epochs found for this channel
                channel_means = [channel_means; mean(channel_epoch_psds{ch}, 1)];
            end
            
            % Averages across all found channels for each animal
            if ~isempty(channel_means)
                mean_across_channels = mean(channel_means, 1);
                
                % Stores one row per animal
                all_state_power.(brain_state).(area) = ...
                    [all_state_power.(brain_state).(area); mean_across_channels];
            end
        end
    end
end

% Plots all brain states for each brain region
figure;
for area_idx = 1:length(areas)
    subplot(2,2,area_idx)
    
    ax = gca;
    set(ax, 'XScale', 'log');
    hold(ax, 'on');
    
    area = areas{area_idx};
    for s = 1:length(states)
        state_data = all_state_power.(states{s});
        area_data = state_data.(area);
        if isempty(area_data)
            continue;
        end

        % Computes mean across animals
        grand_mean = mean(area_data, 1);
        
        % Computes SEM across animals
        n_animals = size(area_data, 1); % animals x frequencies
        sem = std(area_data, 0, 1) / sqrt(n_animals);
        
        % Plots shaded SEM bands
        upper = grand_mean + sem;
        lower = grand_mean - sem;
        fill(ax, [fvec_sel; flipud(fvec_sel)], ...
            [upper(:); flipud(lower(:))], ...
            state_colours(s,:), ...
            'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, fvec_sel, grand_mean, 'Color', state_colours(s,:), 'LineWidth', 2, 'DisplayName', states{s});
    end

    xlim(ax, freq_range);
    set(ax, 'XTick', [1 2 4 6 8 10 20 40], 'XTickLabel', {'1','2','4','6','8','10','20','40'}, 'FontSize', 14);
    xlabel(ax, 'Frequency (Hz)', 'FontSize', 14);
    ylabel(ax, 'Power (V²/Hz)', 'FontSize', 14);
    title(ax, [area], 'FontSize', 14);
    legend(ax, 'show', 'Location', 'northeast');
    grid(ax, 'on');
end

sgtitle('Raw Power Spectra (1–40 Hz) - All Brain States');
