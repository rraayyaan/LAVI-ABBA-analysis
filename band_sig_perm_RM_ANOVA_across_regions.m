% Finds per-frequency significant differences in band presence across regions
% Computes permutation repeated-measures ANOVA with FDR Correction for band distribution across brain regions by shuffling region labels within each animal
% Regions are the repeated within-subject factor
% Obtains an observed F statistic, estimates permutation F distribution, computes permutation p-values, and then applies FDR correction.
% Post-hoc pairwise comparisons with Holm-Bonferroni correction for each significant frequency across animals
% Loops across all band-types for all states and compares regions (A1, BLA, HPC, PFC) within each state

clear; close all; clc;

rng(0); % Set random seed for reproducibility

% Parameters
states      = {'AW', 'QW', 'NREM', 'REM'};
band_types    = {'sustained', 'transient', 'noband'};
regions     = {'A1','BLA','HPC','PFC'};
animal_ids  = {'r14','r16','r19','r20'};

num_states      = numel(states);
num_regions     = numel(regions);
num_animals     = numel(animal_ids);
foi             = logspace(log10(1), log10(40), 96);
num_frequencies = numel(foi);
alpha           = 0.05;
n_perm          = 10000;

% Set save path
save_path = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Sig band distributions HOLM ACROSS REGIONS';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Colours
colors_regions = { ...
    [0.1250 0.5625 0.9375], ...   % A1
    [1.0000 0.4980 0.0549], ...   % BLA
    [0.1765 0.6275 0.1765], ...   % HPC
    [0.8392 0.1529 0.1608]};      % PFC

% All unique region combinations for post-hoc comparisons
pair_names_cell = {'A1 vs BLA', 'A1 vs HPC', 'A1 vs PFC', 'BLA vs HPC', 'BLA vs PFC', 'HPC vs PFC'};
region_pairs = nchoosek(1:num_regions, 2);
n_pairs      = size(region_pairs, 1);

% Preallocate cache for data per state:
% Data is identical across band-types within a state
% Data is loaded once per state and reused across the three band-types instead of re-reading files from disk
sustained_all = cell(num_states, 1);
transient_all = cell(num_states, 1);
noband_all    = cell(num_states, 1);

% Initialise results CSV files
% combined_csv_file: post-hoc pairwise comparison results
% omnibus_csv_file: permutation based RM-ANOVA results
combined_csv_file = fullfile(save_path, 'posthoc_results_permutation_RM_ANOVA_measures_across_regions.csv');
csv_headers = {'Frequency_Hz','State','Band-type','Region1','Region2', ...
               'N_Animals', ...
               'Mean_Diff_Pct','SE_Diff_Pct','CI_Lower_95_Pct','CI_Upper_95_Pct', ...
               'T_Statistic','DF', ...
               'Partial_Eta_Squared', ...
               'P_raw','P_holm','Significant'};
fid = fopen(combined_csv_file, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', csv_headers{:});
fclose(fid);

% Initialise omnibus ANOVA results CSV
omnibus_csv_file = fullfile(save_path, 'omnibus_anova_results_across_regions.csv');
omnibus_headers  = {'State','Band-type','Frequency_Hz','F_Observed', ...
                    'P_permutation','P_FDR','Significant_FDR','Partial_Eta_Squared'};
fid = fopen(omnibus_csv_file, 'w');
fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s\n', omnibus_headers{:});
fclose(fid);

fprintf('\n=== Starting loop across %d band-types and %d states ===\n', numel(band_types), numel(states));

% Outer loop: band types - figures are created per band type showing all states
for band_type_idx = 1:numel(band_types)
    band_type = band_types{band_type_idx};

    % Pre-allocate per-state storage for this band-type figures
    all_mean_by_region = cell(num_states, 1);   % [regions × freqs] per state
    all_sem_by_region  = cell(num_states, 1);
    all_sig_matrix     = cell(num_states, 1);   % [n_pairs × freqs] per state
    all_data_matrix    = cell(num_states, 1);   % [animals × regions × freqs]

    % Inner loop: states
    for state_idx = 1:num_states
        target_state = states{state_idx};

        % Loads data:
        % On the second and third band-type iterations retrieves the previously computed data instead of reloading files from disk
        if ~isempty(sustained_all{state_idx}) && ...
           ~isempty(transient_all{state_idx}) && ...
           ~isempty(noband_all{state_idx})
            fprintf('Using precomputed data for state: %s\n', target_state);
            sustained_pct = sustained_all{state_idx};
            transient_pct = transient_all{state_idx};
            noband_pct    = noband_all{state_idx};
        % Defines the data directory path and adds the LAVI toolbox folder to MATLAB search path on the first band-type iteration for a new state
        else

            fprintf('Computing data for state: %s\n', target_state);
            DATApath    = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';
            LAVIpath    = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)';
            addpath(LAVIpath);
            condition   = '1_habituation';
            alpha_range = [6 8];

            sustained_pct = cell(num_regions, num_animals);
            transient_pct = cell(num_regions, num_animals);
            noband_pct    = cell(num_regions, num_animals);

            % For every region/animal combination loads the LAVI/PINK matrices and classifies each frequency/epoch/channel as sustained, transient, or no-band
            for region_idx = 1:num_regions
                for animal_idx = 1:num_animals
                    animal_id   = animal_ids{animal_idx};
                    lavi_folder = fullfile([animal_id ' matrices'], ...
                        [animal_id ' ' condition ' matrices'], regions{region_idx}, ...
                        [regions{region_idx} ' ' target_state], ...
                        ['LAVI processing (' regions{region_idx} ' ' target_state ')']);
                    file_cut = fullfile(DATApath, lavi_folder, ...
                        sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', ...
                        animal_id, condition, regions{region_idx}, target_state));

                    if ~isfile(file_cut)
                        warning('Missing CUT file: %s', file_cut);
                        sustained_pct{region_idx, animal_idx} = nan(1, num_frequencies);
                        continue;
                    end

                    loaded_file = load(file_cut);
                    if ~isfield(loaded_file,'LAVI_matrix') || ~isfield(loaded_file,'PINK_matrix')
                        warning('Missing required variables in file: %s', file_cut);
                        sustained_pct{region_idx, animal_idx} = nan(1, num_frequencies);
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
                            noband_mask         = ~(sust_mask | tran_mask);
                            sustained_ep(:, ep) = sust_mask;
                            transient_ep(:, ep) = tran_mask;
                            noband_ep(:, ep)    = noband_mask;
                        end

                        % Averages across epochs and converts to percentage
                        chan_sustained(ch, :) = mean(sustained_ep, 2, 'omitnan')' * 100;
                        chan_transient(ch, :) = mean(transient_ep, 2, 'omitnan')' * 100;
                        chan_noband(ch, :)    = mean(noband_ep,    2, 'omitnan')' * 100;
                    end

                    % Average across channels
                    sustained_pct{region_idx, animal_idx} = mean(chan_sustained, 1, 'omitnan');
                    transient_pct{region_idx, animal_idx} = mean(chan_transient, 1, 'omitnan');
                    noband_pct{region_idx, animal_idx}    = mean(chan_noband,    1, 'omitnan');
                end
            end

            % Saves for reuse for the next band-type (avoids reloading files)
            sustained_all{state_idx} = sustained_pct;
            transient_all{state_idx} = transient_pct;
            noband_all{state_idx}    = noband_pct;
        end

        % Band-type selection
        switch band_type
            case 'sustained', data = sustained_pct;
            case 'transient', data = transient_pct;
            case 'noband',    data = noband_pct;
            otherwise,        error('Unknown band_type: %s', band_type);
        end

        % Builds data matrix [animals × regions × frequencies]
        data_matrix = nan(num_animals, num_regions, num_frequencies);
        for region_idx = 1:num_regions
            for animal_idx = 1:num_animals
                vec = data{region_idx, animal_idx};
                if ~isempty(vec) && ~all(isnan(vec))
                    data_matrix(animal_idx, region_idx, :) = vec;
                end
            end
        end
        
        all_data_matrix{state_idx} = data_matrix;

        % Averages and computes SEM across animals for plotting
        mean_by_region = nan(num_regions, num_frequencies);
        sem_by_region  = nan(num_regions, num_frequencies);
        for region_idx = 1:num_regions
            temp = squeeze(data_matrix(:, region_idx, :));
            n_valid_animals = sum(~isnan(temp), 1);
            mean_by_region(region_idx, :) = mean(temp, 1, 'omitnan');
            sem_by_region(region_idx, :)  = std(temp, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);
        end
        all_mean_by_region{state_idx} = mean_by_region;
        all_sem_by_region{state_idx}  = sem_by_region;

        % Exports raw band occurrence percentages to CSV
        occurrence_csv_file = fullfile(save_path, 'raw_occurrence_percentages.csv');
        if state_idx == 1 && band_type_idx == 1
            fid = fopen(occurrence_csv_file, 'w');
            fprintf(fid, 'State,Band-type,Region,Animal,Frequency_Hz,Occurrence_Pct,Mean_Across_Animals,SEM_Across_Animals\n');
            fclose(fid);
        end

        fid = fopen(occurrence_csv_file, 'a');
        for region_idx = 1:num_regions
            animal_matrix = nan(num_animals, num_frequencies);
            for animal_idx = 1:num_animals
                vec = data{region_idx, animal_idx};
                if ~isempty(vec) && ~all(isnan(vec))
                    animal_matrix(animal_idx, :) = vec;
                end
            end
            mean_across     = mean(animal_matrix, 1, 'omitnan');
            n_valid_animals = sum(~isnan(animal_matrix), 1);
            sem_across      = std(animal_matrix, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);

            for animal_idx = 1:num_animals
                vec = data{region_idx, animal_idx};
                if isempty(vec) || all(isnan(vec)), continue; end
                for f = 1:num_frequencies
                    fprintf(fid, '%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f\n', ...
                        target_state, band_type, regions{region_idx}, animal_ids{animal_idx}, ...
                        foi(f), vec(f), mean_across(f), sem_across(f));
                end
            end
        end
        fclose(fid);
        fprintf('Raw occurrence percentages appended for: %s | %s\n', target_state, band_type);

        % PERMUTATION RM-ANOVA (with FDR)
        % For each frequency computes a one-way repeated-measures ANOVA across regions (F statistic)
        % Then build a null distribution by repeatedly shuffling the region labels for each animal and recomputing the F statistic
        % The permutation p-value is how often the shuffled F statistic exceeds the observed F
        % FDR correction is then applied across frequencies
        fprintf('\nRunning permutation repeated-measures ANOVA\n');
        fprintf('State: %s | Band-type: %s | Frequencies: %d | Permutations: %d\n', ...
            target_state, band_type, num_frequencies, n_perm);

        p_vals = nan(num_frequencies, 1);
        F_obs = nan(num_frequencies, 1);
        eta2p_obs = nan(num_frequencies, 1);

        for f = 1:num_frequencies
            Y = squeeze(data_matrix(:, :, f)); % [animals × regions]
            valid   = all(~isnan(Y), 2); % keeps only animals with data in every region
            Y       = Y(valid, :);
            n_valid = size(Y, 1);
            if n_valid < 3, continue; end % needs at least 3 animals

            % Observed repeated-measures ANOVA - region as within-subject factor
            grand_mean   = mean(Y, 'all');
            region_means = mean(Y, 1);
            animal_means = mean(Y, 2);

            SS_total   = sum((Y - grand_mean).^2, 'all');
            SS_regions = n_valid * sum((region_means - grand_mean).^2);
            SS_animals = num_regions * sum((animal_means - grand_mean).^2);
            SS_error   = SS_total - SS_regions - SS_animals;

            df_regions = num_regions - 1;
            df_error   = (n_valid - 1) * (num_regions - 1);
            MS_regions = SS_regions / df_regions;
            MS_error   = SS_error   / df_error;

            if MS_error == 0 || isnan(MS_error)
                F_obs(f)  = Inf;
                p_vals(f) = 0;
                continue;
            end

            F_obs(f)     = MS_regions / MS_error;
            eta2p_obs(f) = SS_regions / (SS_regions + SS_error);

            % Permutation: shuffles region labels within each animal
            perm_F = zeros(n_perm, 1);
            for perm = 1:n_perm
                Y_perm = Y;
                for a = 1:n_valid
                    Y_perm(a, :) = Y(a, randperm(num_regions));
                end
                gpmean = mean(Y_perm, 'all');
                rm     = mean(Y_perm, 1);
                am     = mean(Y_perm, 2);
                SSt  = sum((Y_perm - gpmean).^2, 'all');
                SStr = n_valid * sum((rm - gpmean).^2);
                SSam = num_regions * sum((am - gpmean).^2);
                SSep = SSt - SStr - SSam;
                MSe  = SSep / df_error;
                if MSe > 0
                    perm_F(perm) = (SStr / df_regions) / MSe;
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
                target_state, band_type, foi(f), F_obs(f), p_vals(f), p_fdr(f), p_fdr(f) < alpha, eta2p_obs(f));
        end
        fclose(fid);

        % POST-HOC PAIRWISE COMPARISONS (Holm-Bonferroni)
        % Only runs at frequencies where the omnibus ANOVA was significant following FDR
        % For each frequency every region pair is compared with a paired t-test on the animal-level difference
        % Then Holm-Bonferroni corrects the six pairwise p-values for multiple comparisons
        fprintf('\nRunning post-hoc pairwise comparisons with Holm-Bonferroni correction...\n');
        fprintf('Region pairs: %d | Holm-Bonferroni correction applied per frequency\n', n_pairs);

        p_posthoc_holm = nan(n_pairs, num_frequencies);
        t_crit = tinv(1 - alpha/2, num_animals - 1);

        if ~isempty(sig_freqs)
            for f_idx = sig_freqs'
                freq_hz = foi(f_idx);

                Y_posthoc = nan(num_animals, num_regions);
                for a = 1:num_animals
                    for r = 1:num_regions
                        vec = data{r, a};
                        if ~isempty(vec) && f_idx <= length(vec) && ~isnan(vec(f_idx))
                            Y_posthoc(a, r) = vec(f_idx);
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

                % Paired t-test for every region pair
                for i = 1:n_pairs
                    r1   = region_pairs(i, 1);
                    r2   = region_pairs(i, 2);
                    diff = Y_posthoc(:, r1) - Y_posthoc(:, r2);
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
                    r1       = region_pairs(i, 1);
                    r2       = region_pairs(i, 2);
                    ci_lower = mean_diffs(i) - t_crit * se_diffs(i);
                    ci_upper = mean_diffs(i) + t_crit * se_diffs(i);
                    p_adj    = p_holm(i);
                    is_sig   = (p_adj < alpha);

                    fprintf(fid, '%.4f,%s,%s,%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.6f,%.6f,%d\n', ...
                        freq_hz, target_state, band_type, regions{r1}, regions{r2}, ...
                        n_avail, ...
                        mean_diffs(i), se_diffs(i), ci_lower, ci_upper, ...
                        t_stats(i), n_avail - 1, ...
                        eta2p(i), ...
                        p_raw(i), p_holm(i), is_sig);

                    if is_sig
                        fprintf('  %s vs %s: diff = %.2f%%, t(%d) = %.3f, eta2p = %.3f, p_holm = %.4f *\n', ...
                            regions{r1}, regions{r2}, mean_diffs(i), ...
                            n_avail - 1, t_stats(i), eta2p(i), p_holm(i));
                    end
                end
                fclose(fid);
            end
        else
            fprintf('No significant frequencies found — skipping post-hoc.\n');
        end

        fprintf('\nPost-hoc comparisons complete.\n');

        % Accumulates three-value significance matrix for this state heatmap:
        % NaN = frequency not tested post-hoc (omnibus not significant)
        % 0 = post-hoc tested but not significant
        % 1 = significant after Holm correction
        sig_matrix  = nan(n_pairs, num_frequencies);
        tested_mask = ~isnan(p_posthoc_holm);
        sig_matrix(tested_mask) = double(p_posthoc_holm(tested_mask) < alpha);
        all_sig_matrix{state_idx} = sig_matrix;
    end

    % Figure 1: 4-panel mean ± SEM line plots (one panel per state)
    % Shows each region's mean occurrence line plot with shaded SEM band
    fig_line = figure('Visible', 'off', 'Position', [50, 50, 1600, 1000]);

    for state_idx = 1:num_states
        ax = subplot(2, 2, state_idx);
        hold(ax, 'on');
        set(ax, 'XScale', 'log', 'Box', 'on');

        mean_by_region = all_mean_by_region{state_idx};
        sem_by_region  = all_sem_by_region{state_idx};
        leg_h = gobjects(num_regions, 1);

        for region_idx = 1:num_regions
            x  = foi(:);
            y  = mean_by_region(region_idx, :)';
            e  = sem_by_region(region_idx, :)';
            ok = ~isnan(y) & ~isnan(e);

            % SEM shaded bands
            fill(ax, [x(ok); flipud(x(ok))], ...
                     [y(ok) - e(ok); flipud(y(ok) + e(ok))], ...
                 colors_regions{region_idx}, 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');

            % Mean line
            leg_h(region_idx) = semilogx(ax, foi, mean_by_region(region_idx, :), ...
                'LineWidth', 2.2, 'Color', colors_regions{region_idx}, ...
                'DisplayName', regions{region_idx});
        end

        set(ax, 'XTick', [1 2 4 6 8 10 20 40], ...
                'XTickLabel', {'1','2','4','6','8','10','20','40'}, ...
                'XLim', [0.9 42], ...
                'FontSize', 15);
        xlabel(ax, 'Frequency (Hz)', 'FontSize', 19);
        if strcmp(band_type, 'noband')
            ylab = 'Relative occurrence of no bands (%)';
        else
            ylab = 'Relative occurrence of bands (%)';
        end
        ylabel(ax, ylab, 'FontSize', 17);
        ylim(ax, [0 105]);
        title(ax, states{state_idx}, 'FontSize', 20, 'FontWeight', 'bold');
        % legend(leg_h, 'Location', 'northeastoutside', 'FontSize', 10); % COMMENTED OUT
        grid(ax, 'on');
        hold(ax, 'off');
    end

    lineplot_filename = fullfile(save_path, sprintf('%s_grand_mean_across_regions.png', band_type));
    print(fig_line, '-dpng', '-r300', lineplot_filename);
    fprintf('Mean occurrence plot saved to: %s\n', lineplot_filename);
    close(fig_line);

    % Figure 2: 4-panel significance heatmaps (one panel per state)
    % Rows: region pairs
    % Columns: frequencies
    fig_heatmap = figure('Visible', 'off', 'Position', [50, 50, 1600, 1000]);

    for state_idx = 1:num_states
        ax = subplot(2, 2, state_idx);

        sig_matrix  = all_sig_matrix{state_idx};
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
        title(ax, states{state_idx}, 'FontSize', 20, 'FontWeight', 'bold');
        grid(ax, 'on');
        box(ax, 'on');
    end

    cb_ax = axes(fig_heatmap, 'Position', [0.93 0.15 0.015 0.70], 'Visible', 'off');
    colormap(cb_ax, [0.75 0.75 0.75; 1 1 1; 0 0 0]);
    caxis(cb_ax, [-1 1]);
    cbar = colorbar(cb_ax, 'Location', 'eastoutside');
    set(cbar, 'Ticks',      [-1,                        0,             1], ...
              'TickLabels', {'Untested (RM-ANOVA NS)', 'Tested (NS)', 'Significant (p<0.05)'}, ...
              'FontSize',   11);

    heatmap_filename = fullfile(save_path, sprintf('%s_heatmap_across_regions.png', band_type));
    print(fig_heatmap, '-dpng', '-r300', heatmap_filename);
    fprintf('Heatmap saved to: %s\n', heatmap_filename);
    close(fig_heatmap);

end

fprintf('\n=== ALL COMBINATIONS COMPLETE ===\n');
fprintf('Total combinations processed: %d and %d = %d\n', ...
    numel(band_types), numel(states), numel(band_types) * numel(states));
fprintf('Combined post-hoc results saved to: %s\n', combined_csv_file);
fprintf('Omnibus ANOVA results saved to: %s\n', omnibus_csv_file);