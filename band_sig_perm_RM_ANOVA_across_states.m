% Finds per-frequency significant differences in band presence across states
% Computes permutation repeated-measures ANOVA with FDR Correction for band distribution across brain states by shuffling state labels within each animal
% States are the repeated within-subject factor
% Obtains an observed F statistic, estimates permutation F distribution, computes permutation p-values, and then applies FDR correction.
% Post-hoc pairwise comparisons with Holm-Bonferroni correction for each significant frequency across animals
% Loops across all band-types for all areas and compares states (AW, QW, NREM, REM) within each area

clear; close all; clc;

rng(0); % Set random seed for reproducibility

% Parameters
areas    = {'A1','BLA','HPC','PFC'};
measures = {'sustained','transient','noband'};

% Set save path
save_path = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Sig band distributions HOLM ACROSS STATES';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Preallocate cache for data per area:
% Data is identical across measures within an area
% Data is loaded once per area and reused across the three measures instead of re-reading files from disk
sustained_all = cell(numel(areas), 1);
transient_all = cell(numel(areas), 1);
noband_all    = cell(numel(areas), 1);

% Initialise results CSV files
% combined_csv_file: post-hoc pairwise comparison results
% omnibus_csv_file: permutation based RM-ANOVA results
combined_csv_file = fullfile(save_path, 'posthoc_results_permutation_RM_ANOVA_measures_across_states.csv');
csv_headers = {'Frequency_Hz','Area','Measure','State1','State2', ...
               'N_Animals', ...
               'Mean_Diff_Pct','SE_Diff_Pct','CI_Lower_95_Pct','CI_Upper_95_Pct', ...
               'T_Statistic','DF', ...
               'Partial_Eta_Squared', ...
               'P_raw','P_holm','Significant'};

fid = fopen(combined_csv_file, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', csv_headers{:});
fclose(fid);

% Initialise omnibus ANOVA results CSV
omnibus_csv_file = fullfile(save_path, 'omnibus_anova_results_across_states.csv');
omnibus_headers  = {'Area','Measure','Frequency_Hz','F_Observed', ...
                    'P_permutation','P_FDR','Significant_FDR', 'Eta2p_Omnibus'};
fid = fopen(omnibus_csv_file, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s\n', omnibus_headers{:});
fclose(fid);

% Parameters
animal_ids      = {'r14','r16','r19','r20'};
states          = {'AW','QW','NREM','REM'};
num_states      = numel(states);
num_animals     = numel(animal_ids);
foi             = logspace(log10(1), log10(40), 96);
num_frequencies = numel(foi);
alpha           = 0.05;
n_perm          = 10000;

% Colours
colors = {[0.8, 0.6, 0.0], ...   % AW
          [0.0, 0.6, 0.0], ...   % QW
          [0.0, 0.0, 1.0], ...   % NREM
          [1.0, 0.0, 0.0]};      % REM

% All unique state combinations for post-hoc comparisons
pair_names_cell = {'AW vs QW','AW vs NREM','AW vs REM','QW vs NREM','QW vs REM','NREM vs REM'};
state_pairs = nchoosek(1:num_states, 2);
n_pairs     = size(state_pairs, 1);

% Outer loop: measures - figures are created per measure showing all areas
for measure_idx = 1:numel(measures)
    measure = measures{measure_idx};

    % Pre-allocate per-area storage for this measure's figures
    all_mean_by_state = cell(numel(areas), 1);   % [states × freqs] per area
    all_sem_by_state  = cell(numel(areas), 1);
    all_sig_matrix    = cell(numel(areas), 1);   % [n_pairs × freqs] per area
    all_data_matrix   = cell(numel(areas), 1);   % [animals × states × freqs]

    % Inner loop: areas
    for area_idx = 1:numel(areas)
        area = areas{area_idx};

        % Loads data:
        % On the second and third measure iterations retrieves the previously computed data instead of reloading files from disk
        if ~isempty(sustained_all{area_idx}) && ...
           ~isempty(transient_all{area_idx}) && ...
           ~isempty(noband_all{area_idx})
            fprintf('Using precomputed data for area: %s\n', area);
            sustained_pct = sustained_all{area_idx};
            transient_pct = transient_all{area_idx};
            noband_pct    = noband_all{area_idx};

        % Defines the data directory path and adds the LAVI toolbox folder to MATLAB search path on the first measure iteration for a new area
        else
            fprintf('Computing data for area: %s\n', area);
            DATApath  = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';
            LAVIpath  = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
            addpath(LAVIpath);
            condition   = '1_habituation';
            alpha_range = [6 8];
            sustained_pct = cell(num_states, num_animals);
            transient_pct = cell(num_states, num_animals);
            noband_pct    = cell(num_states, num_animals);

            % For every state/animal combination loads the LAVI/PINK matrices and classifies each frequency/epoch/channel as sustained, transient, or no-band
            for state_idx = 1:num_states
                for animal_idx = 1:num_animals
                    animal_id = animal_ids{animal_idx};
                    lavi_folder = fullfile([animal_id ' matrices'], ...
                        [animal_id ' ' condition ' matrices'], area, ...
                        [area ' ' states{state_idx}], ...
                        ['LAVI processing (' area ' ' states{state_idx} ')']);

                    file_cut = fullfile(DATApath, lavi_folder, ...
                        sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', ...
                        animal_id, condition, area, states{state_idx}));

                    if ~isfile(file_cut)
                        warning('Missing CUT file: %s', file_cut);
                        sustained_pct{state_idx, animal_idx} = nan(1, num_frequencies);
                        continue;
                    end

                    loaded_file = load(file_cut);
                    if ~isfield(loaded_file,'LAVI_matrix') || ~isfield(loaded_file,'PINK_matrix')
                        warning('Missing required variables in file: %s', file_cut);
                        sustained_pct{state_idx, animal_idx} = nan(1, num_frequencies);
                        continue;
                    end
                    LAVI_matrix = loaded_file.LAVI_matrix;
                    PINK_matrix = loaded_file.PINK_matrix;
                    [num_chan_lavi, ~, num_epochs_lavi] = size(LAVI_matrix);

                    chan_sustained = nan(num_chan_lavi, num_frequencies);
                    chan_transient = nan(num_chan_lavi, num_frequencies);
                    chan_noband    = nan(num_chan_lavi, num_frequencies);

                    for ch = 1:num_chan_lavi
                    % Per-epoch flags: 1 = band present at that frequency for this channel
                        sustained_ep = nan(num_frequencies, num_epochs_lavi);
                        transient_ep = nan(num_frequencies, num_epochs_lavi);
                        noband_ep    = nan(num_frequencies, num_epochs_lavi);

                        for ep = 1:num_epochs_lavi
                            lavi_vec = squeeze(LAVI_matrix(ch, :, ep));

                            if all(isnan(lavi_vec))
                                continue;
                            end

                            % Null reference significance limits from PINK simulations
                            pink    = squeeze(PINK_matrix(:, :, ch, ep));
                            pink    = permute(pink, [3 2 1]);
                            sig_lim = cat(3, min(pink,[],3), max(pink,[],3));

                            [borders, ~, ~] = ABBA(lavi_vec, foi, alpha_range, sig_lim, 0);

                            sust_mask = false(1, num_frequencies);
                            tran_mask = false(1, num_frequencies);

                            % For each row:
                            % Columns 1 and 2: start and end frequency index
                            % Column 9: direction
                                % 1 = sustained
                                % -1 = transient
                            % Column 11: significance
                            if ~isempty(borders{1})
                                bands = borders{1};
                                for b = 1:size(bands,1)
                                    BegI = bands(b,1);
                                    EndI = bands(b,2);
                                    Dir  = bands(b,9);
                                    Sig  = bands(b,11);
                                    if Sig ~= 1, continue; end
                                    if     Dir ==  1, sust_mask(BegI:EndI) = true;
                                    elseif Dir == -1, tran_mask(BegI:EndI) = true;
                                    end
                                end
                            end

                            % Frequencies not covered by any significant band are classified as no-band
                            noband_mask          = ~(sust_mask | tran_mask);
                            sustained_ep(:, ep)  = sust_mask;
                            transient_ep(:, ep)  = tran_mask;
                            noband_ep(:, ep)     = noband_mask;
                        end

                        % Averages across epochs and converts to percentage
                        chan_sustained(ch, :) = mean(sustained_ep, 2, 'omitnan')' * 100;
                        chan_transient(ch, :) = mean(transient_ep, 2, 'omitnan')' * 100;
                        chan_noband(ch, :)    = mean(noband_ep,    2, 'omitnan')' * 100;
                    end

                    % Average across channels
                    sustained_pct{state_idx, animal_idx} = mean(chan_sustained, 1, 'omitnan');
                    transient_pct{state_idx, animal_idx} = mean(chan_transient, 1, 'omitnan');
                    noband_pct{state_idx, animal_idx}    = mean(chan_noband,    1, 'omitnan');
                end
            end

            % Saves for reuse for the next measure (avoids reloading files)
            sustained_all{area_idx} = sustained_pct;
            transient_all{area_idx} = transient_pct;
            noband_all{area_idx}    = noband_pct;
        end

        % Measure selection
        switch measure
            case 'sustained', data = sustained_pct;
            case 'transient', data = transient_pct;
            case 'noband',    data = noband_pct;
            otherwise,        error('Unknown measure: %s', measure);
        end

        % Builds data matrix [animals × states × frequencies]
        data_matrix = nan(num_animals, num_states, num_frequencies);
        for state_idx = 1:num_states
            for animal_idx = 1:num_animals
                vec = data{state_idx, animal_idx};
                if ~isempty(vec) && ~all(isnan(vec))
                    data_matrix(animal_idx, state_idx, :) = vec;
                end
            end
        end
        all_data_matrix{area_idx} = data_matrix;

        % Averages and computes SEM across animals for plotting
        mean_by_state = nan(num_states, num_frequencies);
        sem_by_state  = nan(num_states, num_frequencies);
        for state_idx = 1:num_states
            temp = squeeze(data_matrix(:, state_idx, :));
            n_valid_animals = sum(~isnan(temp), 1);
            mean_by_state(state_idx, :) = mean(temp, 1, 'omitnan');
            sem_by_state(state_idx, :)  = std(temp, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);
        end
        all_mean_by_state{area_idx} = mean_by_state;
        all_sem_by_state{area_idx}  = sem_by_state;

        % Exports raw band occurrence percentages to CSV
        occurrence_csv_file = fullfile(save_path, 'raw_occurrence_percentages.csv');
        if area_idx == 1 && measure_idx == 1
            fid = fopen(occurrence_csv_file, 'w');
            fprintf(fid, 'Area,Measure,State,Animal,Frequency_Hz,Occurrence_Pct,Mean_Across_Animals,SEM_Across_Animals\n');
            fclose(fid);
        end

        fid = fopen(occurrence_csv_file, 'a');
        for state_idx = 1:num_states
            animal_matrix = nan(num_animals, num_frequencies);
            for animal_idx = 1:num_animals
                vec = data{state_idx, animal_idx};
                if ~isempty(vec) && ~all(isnan(vec))
                    animal_matrix(animal_idx, :) = vec;
                end
            end
            mean_across     = mean(animal_matrix, 1, 'omitnan');
            n_valid_animals = sum(~isnan(animal_matrix), 1);
            sem_across      = std(animal_matrix, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);

            for animal_idx = 1:num_animals
                vec = data{state_idx, animal_idx};
                if isempty(vec) || all(isnan(vec)), continue; end
                for f = 1:num_frequencies
                    fprintf(fid, '%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f\n', ...
                        area, measure, states{state_idx}, animal_ids{animal_idx}, ...
                        foi(f), vec(f), mean_across(f), sem_across(f));
                end
            end
        end
        fclose(fid);
        fprintf('Raw occurrence percentages appended for: %s | %s\n', area, measure);

        % PERMUTATION RM-ANOVA (with FDR)
        % For each frequency computes a one-way repeated-measures ANOVA across states (F statistic)
        % Then build a null distribution by repeatedly shuffling the state labels for each animal and recomputing the F statistic
        % The permutation p-value is how often the shuffled F statistic exceeds the observed F
        % FDR correction is then applied across frequencies
        fprintf('\nRunning permutation repeated-measures ANOVA...\n');
        fprintf('Area: %s | Measure: %s | Frequencies: %d | Permutations: %d\n', ...
            area, measure, num_frequencies, n_perm);

        p_vals    = nan(num_frequencies, 1);
        F_obs     = nan(num_frequencies, 1);
        eta2p_obs = nan(num_frequencies, 1);

        for f = 1:num_frequencies
            Y = squeeze(data_matrix(:, :, f)); % [animals × states]
            valid   = all(~isnan(Y), 2); % keeps only animals with data in every state
            Y       = Y(valid, :);
            n_valid = size(Y, 1);
            if n_valid < 3, continue; end % needs at least 3 animals

            % Observed repeated-measures ANOVA - state as within-subject factor
            grand_mean   = mean(Y, 'all');
            state_means  = mean(Y, 1);
            animal_means = mean(Y, 2);

            SS_total   = sum((Y - grand_mean).^2, 'all');
            SS_states  = n_valid * sum((state_means  - grand_mean).^2);
            SS_animals = num_states * sum((animal_means - grand_mean).^2);
            SS_error   = SS_total - SS_states - SS_animals;

            df_states = num_states - 1;
            df_error  = (n_valid - 1) * (num_states - 1);
            MS_states = SS_states / df_states;
            MS_error  = SS_error  / df_error;

            if MS_error == 0 || isnan(MS_error)
                F_obs(f)  = Inf;
                p_vals(f) = 0;
                continue;
            end

            F_obs(f)     = MS_states / MS_error;
            eta2p_obs(f) = SS_states / (SS_states + SS_error);

            % Permutation: shuffles state labels within each animal
            perm_F = zeros(n_perm, 1);
            for perm = 1:n_perm
                Y_perm = Y;
                for a = 1:n_valid
                    Y_perm(a, :) = Y(a, randperm(num_states));
                end
                gpmean = mean(Y_perm, 'all');
                sm     = mean(Y_perm, 1);
                am     = mean(Y_perm, 2);
                SSt  = sum((Y_perm - gpmean).^2, 'all');
                SSts = n_valid * sum((sm - gpmean).^2);
                SSam = num_states * sum((am - gpmean).^2);
                SSep = SSt - SSts - SSam;
                MSe  = SSep / df_error;
                if MSe > 0
                    perm_F(perm) = (SSts / df_states) / MSe;
                end
            end

            % Phipson & Smyth correction - avoids a p-value of exactly 0 from a finite number of permutations
            p_vals(f) = (sum(perm_F >= F_obs(f)) + 1) / (n_perm + 1);
            if mod(f, 20) == 0
                fprintf('  Frequency %d/%d complete (p = %.4f)\n', f, num_frequencies, p_vals(f));
            end
        end

        % FDR correction (Benjamini-Hochberg) - applied across all tested frequencies
        p_fdr     = mafdr(p_vals, 'BHFDR', true);
        sig_freqs = find(p_fdr < alpha);
        fprintf('\nPermutation RM-ANOVA + FDR correction complete.\n');
        fprintf('Significant frequencies (FDR < 0.05): %d / %d\n', ...
            length(sig_freqs), sum(~isnan(p_vals)));
        if ~isempty(sig_freqs)
            fprintf('Significant frequency range: %.1f – %.1f Hz\n', ...
                min(foi(sig_freqs)), max(foi(sig_freqs)));
        end

        % Writes omnibus ANOVA results to CSV
        fid = fopen(omnibus_csv_file, 'a');
        for f = 1:num_frequencies
            if isnan(p_vals(f)), continue; end
            fprintf(fid, '%s,%s,%.4f,%.4f,%.6f,%.6f,%d,%.4f\n', ...
                area, measure, foi(f), F_obs(f), p_vals(f), p_fdr(f), p_fdr(f) < alpha, eta2p_obs(f));
        end
        fclose(fid);

        % POST-HOC PAIRWISE COMPARISONS (Holm-Bonferroni)
        % Only runs at frequencies where the omnibus ANOVA was significant following FDR
        % For each frequency every state pair is compared with a paired t-test on the animal-level difference
        % Then Holm-Bonferroni corrects the six pairwise p-values for multiple comparisons
        fprintf('\nRunning post-hoc pairwise comparisons with Holm-Bonferroni correction...\n');
        fprintf('State pairs: %d | Holm-Bonferroni correction applied per frequency\n', n_pairs);

        p_posthoc_holm = nan(n_pairs, num_frequencies);
        t_crit = tinv(1 - alpha/2, num_animals - 1);

        if ~isempty(sig_freqs)
            for f_idx = sig_freqs'
                freq_hz = foi(f_idx);

                Y_posthoc = nan(num_animals, num_states);
                for a = 1:num_animals
                    for s = 1:num_states
                        vec = data{s, a};
                        if ~isempty(vec) && f_idx <= length(vec) && ~isnan(vec(f_idx))
                            Y_posthoc(a, s) = vec(f_idx);
                        end
                    end
                end

                valid     = all(~isnan(Y_posthoc), 2);
                Y_posthoc = Y_posthoc(valid, :);
                n_avail   = size(Y_posthoc, 1);
                if n_avail < 3, continue; end

                fprintf('\n--- Frequency: %.1f Hz (n = %d animals) ---\n', freq_hz, n_avail);

                p_raw      = nan(n_pairs, 1);
                mean_diffs = nan(n_pairs, 1);
                se_diffs   = nan(n_pairs, 1);
                t_stats    = nan(n_pairs, 1);
                eta2p      = nan(n_pairs, 1);

                % Paired t-test for every state pair
                for i = 1:n_pairs
                    s1   = state_pairs(i, 1);
                    s2   = state_pairs(i, 2);
                    diff = Y_posthoc(:, s1) - Y_posthoc(:, s2);
                    n_d           = numel(diff);
                    df_d          = n_d - 1;
                    mean_diffs(i) = mean(diff);
                    se_diffs(i)   = std(diff) / sqrt(n_d);
                    [~, p_raw(i), ~, stats] = ttest(diff);
                    t_stats(i) = stats.tstat;
                    eta2p(i)   = t_stats(i)^2 / (t_stats(i)^2 + df_d); % effect size from paired t
                end

                % Holm-Bonferroni correction
                [p_sorted, sort_idx] = sort(p_raw);
                holm_factors         = (n_pairs + 1 - (1:n_pairs)');
                p_holm_sorted        = min(cummax(p_sorted .* holm_factors), 1);
                p_holm               = nan(n_pairs, 1);
                p_holm(sort_idx)     = p_holm_sorted;
                p_posthoc_holm(:, f_idx) = p_holm;

                % Writes post-hoc results to CSV
                fid = fopen(combined_csv_file, 'a');
                for i = 1:n_pairs
                    s1          = state_pairs(i, 1);
                    s2          = state_pairs(i, 2);
                    ci_lower    = mean_diffs(i) - t_crit * se_diffs(i);
                    ci_upper    = mean_diffs(i) + t_crit * se_diffs(i);
                    p_adj       = p_holm(i);
                    is_sig      = (p_adj < alpha);

                    fprintf(fid, '%.4f,%s,%s,%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.6f,%.6f,%d\n', ...
                        freq_hz, area, measure, states{s1}, states{s2}, ...
                        n_avail, ...
                        mean_diffs(i), se_diffs(i), ci_lower, ci_upper, ...
                        t_stats(i), n_avail - 1, ...
                        eta2p(i), ...
                        p_raw(i), p_holm(i), is_sig);

                    if is_sig
                        fprintf('  %s vs %s: diff = %.2f%%, t(%d) = %.3f, eta2p = %.3f, p_holm = %.4f *\n', ...
                            states{s1}, states{s2}, mean_diffs(i), ...
                            n_avail - 1, t_stats(i), eta2p(i), p_holm(i));
                    end
                end
                fclose(fid);
            end
        else
            fprintf('No significant frequencies found — skipping post-hoc.\n');
        end

        fprintf('\nPost-hoc comparisons complete.\n');

        % Accumulates three-value significance matrix for this area heatmap:
        % NaN = frequency not tested post-hoc (omnibus not significant)
        % 0 = post-hoc tested but not significant
        % 1 = significant after Holm correction
        sig_matrix  = nan(n_pairs, num_frequencies);
        tested_mask = ~isnan(p_posthoc_holm);
        sig_matrix(tested_mask) = double(p_posthoc_holm(tested_mask) < alpha);
        all_sig_matrix{area_idx} = sig_matrix;

    end

    % Figure 1: 4-panel mean ± SEM line plots (one panel per area)
    % Shows each state's mean occurrence line plot with shaded SEM band
    fig_line = figure('Visible', 'off', 'Position', [50, 50, 1600, 1000]);

    for area_idx = 1:numel(areas)
        ax = subplot(2, 2, area_idx);
        hold(ax, 'on');
        set(ax, 'XScale', 'log', 'Box', 'on');

        mean_by_state = all_mean_by_state{area_idx};
        sem_by_state  = all_sem_by_state{area_idx};
        leg_h = gobjects(num_states, 1);

        for state_idx = 1:num_states
            x  = foi(:);
            y  = mean_by_state(state_idx, :)';
            e  = sem_by_state(state_idx, :)';
            ok = ~isnan(y) & ~isnan(e);

            % SEM shaded bands
            fill(ax, [x(ok); flipud(x(ok))], ...
                     [y(ok) - e(ok); flipud(y(ok) + e(ok))], ...
                 colors{state_idx}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');

            % Mean line
            leg_h(state_idx) = semilogx(ax, foi, mean_by_state(state_idx, :), ...
                'LineWidth', 2.2, 'Color', colors{state_idx}, ...
                'DisplayName', states{state_idx});
        end

        set(ax, 'XTick', [1 2 4 6 8 10 20 40], ...
                'XTickLabel', {'1','2','4','6','8','10','20','40'}, ...
                'XLim', [0.9 42], ...
                'FontSize', 15);
        xlabel(ax, 'Frequency (Hz)', 'FontSize', 19);
        if strcmp(measure, 'noband')
            ylab = 'Relative occurrence of no bands (%)';
        else
            ylab = 'Relative occurrence of bands (%)';
        end
        ylabel(ax, ylab, 'FontSize', 17);
        ylim(ax, [0 105]);
        title(ax, areas{area_idx}, 'FontSize', 20, 'FontWeight', 'bold');
        % legend(leg_h, 'Location', 'northeastoutside', 'FontSize', 10); % COMMENTED OUT
        grid(ax, 'on');
        hold(ax, 'off');
    end

    lineplot_filename = fullfile(save_path, sprintf('%s_grand_mean_across_states.png', measure));
    print(fig_line, '-dpng', '-r300', lineplot_filename);
    fprintf('Mean occurrence plot saved to: %s\n', lineplot_filename);
    close(fig_line);

    % Figure 2: 4-panel significance heatmaps (one panel per area)
    % Rows: state pairs
    % Columns: frequencies
    fig_heatmap = figure('Visible', 'off', 'Position', [50, 50, 1600, 1000]);

    for area_idx = 1:numel(areas)
        ax = subplot(2, 2, area_idx);

        sig_matrix  = all_sig_matrix{area_idx};
        plot_matrix = sig_matrix;
        plot_matrix(isnan(sig_matrix)) = -1;

        ytick_positions = 0.5:n_pairs-0.5;

        imagesc(ax, foi, 1:n_pairs, plot_matrix);
        set(ax, 'YDir', 'normal', 'XScale', 'log', 'TickLength', [0 0]);
        axis(ax, 'tight');

        % 3-color heatmap:
        colormap(ax, [0.75 0.75 0.75;   % grey  = untested (omnibus not significant)
                      1.00 1.00 1.00;   % white = tested but post-hoc not significant
                      0.00 0.00 0.00]); % black = post-hoc significant
        caxis(ax, [-1 1]);

        set(ax, 'XTick', [1 2 4 6 8 10 20 40], ...
                'XTickLabel', {'1','2','4','6','8','10','20','40'}, ...
                'YTick', ytick_positions, ...
                'YTickLabel', pair_names_cell, ...
                'XLim', [0.9 42], ...
                'YLim', [0.5, n_pairs + 0.5], ...
                'FontSize', 15);
        xlabel(ax, 'Frequency (Hz)', 'FontSize', 18);
        title(ax, areas{area_idx}, 'FontSize', 20, 'FontWeight', 'bold');
        grid(ax, 'on');
        box(ax, 'on');
    end

    cb_ax = axes(fig_heatmap, 'Position', [0.93 0.15 0.015 0.70], 'Visible', 'off');
    colormap(cb_ax, [0.75 0.75 0.75; 1 1 1; 0 0 0]);
    caxis(cb_ax, [-1 1]);
    cbar = colorbar(cb_ax, 'Location', 'eastoutside');
    set(cbar, 'Ticks',      [-1,                       0,            1], ...
              'TickLabels', {'Untested (RM-ANOVA NS)', 'Tested (NS)', 'Significant (p<0.05)'}, ...
              'FontSize',   11);

    heatmap_filename = fullfile(save_path, sprintf('%s_heatmap_across_states.png', measure));
    print(fig_heatmap, '-dpng', '-r300', heatmap_filename);
    fprintf('Heatmap saved to: %s\n', heatmap_filename);
    close(fig_heatmap);

end

fprintf('\n=== ALL COMBINATIONS COMPLETE ===\n');
fprintf('Total combinations processed: %d and %d = %d\n', ...
    numel(areas), numel(measures), numel(areas) * numel(measures));
fprintf('Combined post-hoc results saved to: %s\n', combined_csv_file);
fprintf('Omnibus ANOVA results saved to: %s\n', omnibus_csv_file);
