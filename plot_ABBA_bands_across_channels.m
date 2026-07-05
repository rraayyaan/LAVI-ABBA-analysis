% Generates a single plot showing ABBA detected bands (without LAVI profiles) for all channels from all animals
% Each row shows the averaged ABBA bands (across epochs) from a single channel
% Bands are computed per epoch then averaged within each channel across all its epochs

clear; close all; clc;

% Add LAVI toolbox
LAVIpath = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
addpath(LAVIpath);

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Define dataset
animal_ids = {'r14','r16','r19','r20'};
condition  = '1_habituation';
area       = 'PFC';
brain_state = 'NREM';

foi = logspace(log10(1), log10(40), 96); % frequencies of interest

% alpha_range is the frequency range in which we expact to find alpha
% The band of the peak in this band (alpha) will be assigned the index 0
% Bands with lower frequency than alpha band will be assigned with negative indices
% Bands with higher frequency than alpha band will be assigned with positive indices
% alpha_range is only used after band detection for relative numbering of identified bands
alpha_range = [6 8];

% STORAGE FOR PER-CHANNEL AVERAGED BANDS
% For each channel:
    % chan_sustained(ch, f) = percentage of epochs where frequency f has a sustained band
    % chan_transient(ch, f) = percentage of epochs where frequency f has a transient band
    % chan_noband(ch, f) = percentage of epochs where frequency f has no band
all_chan_sustained = {};
all_chan_transient = {};
all_chan_noband    = {};

% Loops across animals
for animal_idx = 1:numel(animal_ids)

    animal_id = animal_ids{animal_idx};

    lavi_folder = fullfile([animal_id ' matrices'], ...
        [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
        ['LAVI processing (' area ' ' brain_state ')']);

    dataFileName = fullfile(DATApath, lavi_folder, ...
        sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

    if ~isfile(dataFileName)
        warning('Missing file: %s', dataFileName);
        continue;
    end

    loaded = load(dataFileName);
    LAVI_matrix = loaded.LAVI_matrix; % channel x freq x epoch
    PINK_matrix = loaded.PINK_matrix; % rep x freq x channel x epoch

    [num_chan, num_freq, num_epochs] = size(LAVI_matrix);

    % Preallocates per-channel band occurrence across epochs
    chan_sustained = nan(num_chan, num_freq);
    chan_transient = nan(num_chan, num_freq);
    chan_noband    = nan(num_chan, num_freq);

    % Loops across channels
    for ch = 1:num_chan

        sustained_ep = nan(num_freq, num_epochs);
        transient_ep = nan(num_freq, num_epochs);
        noband_ep    = nan(num_freq, num_epochs);

        % Loops across epochs
        for ep = 1:num_epochs

            lavi_vec = squeeze(LAVI_matrix(ch, :, ep)); % 1 x freq
            if all(isnan(lavi_vec))
                continue;
            end

            pink = squeeze(PINK_matrix(:, :, ch, ep)); % rep x freq
            pink = permute(pink, [3 2 1]); % 1 x freq x rep

            sig_lim = cat(3, min(pink,[],3), max(pink,[],3));

            [borders, ~, ~] = ABBA(lavi_vec, foi, alpha_range, sig_lim, 0);

            % Flags each frequency as sustained or transient based on this epoch's detected bands
            sust_mask = false(1, num_freq);
            tran_mask = false(1, num_freq);

            if ~isempty(borders{1})
                bands = borders{1};
                for b = 1:size(bands,1)
                    BegI = bands(b,1);
                    EndI = bands(b,2);
                    Dir  = bands(b,9);
                    Sig  = bands(b,11);

                    if Sig ~= 1
                        continue;
                    end

                    if Dir == 1
                        sust_mask(BegI:EndI) = true;
                    elseif Dir == -1
                        tran_mask(BegI:EndI) = true;
                    end
                end
            end

            noband_mask = ~(sust_mask | tran_mask);

            sustained_ep(:, ep) = sust_mask';
            transient_ep(:, ep) = tran_mask';
            noband_ep(:, ep)    = noband_mask';
        end

        % Averages across epochs for this channel (percentage of epochs)
        chan_sustained(ch, :) = mean(sustained_ep, 2, 'omitnan')' * 100;
        chan_transient(ch, :) = mean(transient_ep, 2, 'omitnan')' * 100;
        chan_noband(ch, :)    = mean(noband_ep, 2, 'omitnan')' * 100;
    end

    % Stores per-channel results for this animal
    all_chan_sustained{end+1} = chan_sustained;
    all_chan_transient{end+1} = chan_transient;
    all_chan_noband{end+1}    = chan_noband;
end

% Collects all channels across all animals
all_sust = cat(1, all_chan_sustained{:}); % (total channels) x freq
all_tran = cat(1, all_chan_transient{:});
all_nob  = cat(1, all_chan_noband{:});

num_total_chan = size(all_sust, 1);

% Computes net band presence per channel: [sustained - transient]
    % positive = more often sustained
    % negative = more often transient
sig_prob_channels = all_sust - all_tran; % num_total_chan x freq

% Plots one row per channel
cachol = [83 174 244] / 255; % transient (blue)
yarok  = [108 192 12] / 255; % sustained (green)

y_channels = 1:num_total_chan;

figure; hold on
set(gcf, 'position', [500, 100, 900, 60 * num_total_chan]);

% Plots each channel
for ch = 1:num_total_chan
    sigVect = sig_prob_channels(ch, :);

    posind = sigVect > 0;
    negind = sigVect < 0;

    scatter(foi(posind), ones(1,sum(posind))*y_channels(ch), 40, yarok, 'filled', 's');
    scatter(foi(negind), ones(1,sum(negind))*y_channels(ch), 40, cachol, 'filled', 's');
end

set(gca, 'xscale', 'log', 'xtick', [2:2:10, 20:10:40]);
xlim([1 40]);

yticks(1:num_total_chan);
yticklabels(arrayfun(@(i) sprintf('Chan %d', i), 1:num_total_chan, 'UniformOutput', false));

ylabel('Channel');
xlabel('Frequency (Hz)');
title(sprintf('%s | %s-%s | All animals, all channels\n(ABBA per epoch → avg per channel)', ...
    condition, upper(area), upper(brain_state)));

grid on;

% Legend
h_pos = scatter(foi(1), 0, 40, yarok, 'filled', 's');
h_neg = scatter(foi(1), 0, 40, cachol, 'filled', 's');
legend([h_pos, h_neg], {'Sustained', 'Transient'}, ...
    'Location', 'northeastoutside');
