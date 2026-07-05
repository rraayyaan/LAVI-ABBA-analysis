% Finds per-frequency significant differences in rhythmicity across regions
% Computes grand-averaged LAVI profiles across regions per state
% Computes permutation repeated-measures ANOVA with FDR Correction for LAVI profiles across brain regions by shuffling region labels within each animal
% Regions are the repeated within-subject factor
% Obtains an observed F statistic, estimates permutation F distribution, computes permutation p-values, and then applies FDR correction.
% Post-hoc pairwise comparisons with Holm-Bonferroni correction for each significant frequency across animals
% Loops across all band-types for all states and compares regions (A1, BLA, HPC, PFC) within each state

clear; close all; clc;
rng(0);

% Data parameters
areas = {'A1','BLA','HPC','PFC'};
states = {'AW','QW','NREM','REM'};
animal_ids = {'r14','r16','r19','r20'};
condition = '1_habituation';

num_areas = numel(areas);
num_states = numel(states);
num_animals = numel(animal_ids);

% LAVI parameters
foi = logspace(log10(1), log10(40), 96);
num_frequencies = numel(foi);
alpha = 0.05;
n_perm = 10000;

DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';
save_path = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/LAVI profiles across regions';

if ~exist(save_path,'dir'), mkdir(save_path); end

% Region colours (A1 / BLA / HPC / PFC)
colors = {[0.1250 0.5625 0.9375], ...   % A1  – blue
          [1.0000 0.4980 0.0549], ...   % BLA – orange
          [0.1765 0.6275 0.1765], ...   % HPC – green
          [0.8392 0.1529 0.1608]};      % PFC – red

pair_names_cell = {'A1 vs BLA','A1 vs HPC','A1 vs PFC', ...
                   'BLA vs HPC','BLA vs PFC','HPC vs PFC'};
region_pairs = nchoosek(1:num_areas, 2);   % [6×2]
n_pairs      = size(region_pairs, 1);

% Initialise results CSV files
% omnibus_csv_file: permutation based RM-ANOVA results
% posthoc_csv_file: post-hoc pairwise comparison results
omnibus_csv_file = fullfile(save_path, 'omnibus_anova_lavi_across_regions.csv');
fid = fopen(omnibus_csv_file,'w');
fprintf(fid,'%s\n', strjoin({'State','Frequency_Hz','F_Observed', ...
    'P_permutation','P_FDR','Significant_FDR','Eta2p_Omnibus'}, ','));
fclose(fid);

posthoc_csv_file = fullfile(save_path, 'posthoc_results_lavi_across_regions.csv');
fid = fopen(posthoc_csv_file,'w');
fprintf(fid,'%s\n', strjoin({'Frequency_Hz','State','Region1','Region2', ...
    'N_Animals','Mean_Diff','SE_Diff','CI_Lower_95','CI_Upper_95', ...
    'T_Statistic','DF','Partial_Eta_Squared','P_raw','P_holm','Significant'}, ','));
fclose(fid);

% Pre-allocate storage: [states × areas × animals × frequencies]
lavi_data = cell(num_states, 1);

% Loads and averages LAVI profiles
for state_idx = 1:num_states
    target_state = states{state_idx};
    fprintf('\nLoading LAVI data for state: %s\n', target_state);

    lavi_state = cell(num_areas, num_animals);   % {area, animal}

    for area_idx = 1:num_areas
        for animal_idx = 1:num_animals
            animal_id = animal_ids{animal_idx};

            lavi_folder = fullfile([animal_id ' matrices'], ...
                [animal_id ' ' condition ' matrices'], areas{area_idx}, ...
                [areas{area_idx} ' ' target_state], ...
                ['LAVI processing (' areas{area_idx} ' ' target_state ')']);

            file_cut = fullfile(DATApath, lavi_folder, ...
                sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', ...
                animal_id, condition, areas{area_idx}, target_state));

            if ~isfile(file_cut)
                warning('Missing file: %s', file_cut);
                lavi_state{area_idx, animal_idx} = nan(1, num_frequencies);
                continue;
            end

            loaded = load(file_cut, 'LAVI_matrix');
            if ~isfield(loaded,'LAVI_matrix')
                warning('LAVI_matrix missing in: %s', file_cut);
                lavi_state{area_idx, animal_idx} = nan(1, num_frequencies);
                continue;
            end

            LAVI_matrix = loaded.LAVI_matrix; % [channels × freqs × epochs]
            [num_chan, ~, ~] = size(LAVI_matrix);

            % Averages across epochs then across channels (omits NaN) - produces one LAVI profile per animal per area per state
            chan_mean = nan(num_chan, num_frequencies);
            for ch = 1:num_chan
                ep_vals = squeeze(LAVI_matrix(ch, :, :)); % [freqs × epochs]
                chan_mean(ch, :) = mean(ep_vals, 2, 'omitnan')';
            end
            lavi_state{area_idx, animal_idx} = mean(chan_mean, 1, 'omitnan');
        end
    end

    lavi_data{state_idx} = lavi_state;
    fprintf('  Done (%d areas × %d animals loaded)\n', num_areas, num_animals);
end

% Stores per-state results for plotting
all_mean   = cell(num_states, 1); % {state}[areas × freqs]
all_sem    = cell(num_states, 1);
all_sigmat = cell(num_states, 1); % significance matrix [pairs × freqs]

for state_idx = 1:num_states
    target_state = states{state_idx};
    lavi_state   = lavi_data{state_idx};

    % Builds matrix [animals × areas × freqs]
    data_matrix = nan(num_animals, num_areas, num_frequencies);
    for a_idx = 1:num_areas
        for ani = 1:num_animals
            v = lavi_state{a_idx, ani};
            if ~isempty(v) && ~all(isnan(v))
                data_matrix(ani, a_idx, :) = v;
            end
        end
    end

    % Calculate grand mean and SEM across animals
    mean_by_area = nan(num_areas, num_frequencies);
    sem_by_area  = nan(num_areas, num_frequencies);
    for a_idx = 1:num_areas
        tmp = squeeze(data_matrix(:, a_idx, :)); % [animals × freqs]
        n_valid_animals = sum(~isnan(tmp), 1);
        mean_by_area(a_idx, :) = mean(tmp, 1, 'omitnan');
        sem_by_area(a_idx, :)  = std(tmp, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);
    end
    all_mean{state_idx} = mean_by_area;
    all_sem{state_idx}  = sem_by_area;

    % PERMUTATION RM-ANOVA (with FDR)
    % For each frequency computes a one-way repeated-measures ANOVA across regions (F statistic)
    % Then build a null distribution by repeatedly shuffling the region labels for each animal and recomputing the F statistic
    % The permutation p-value is how often the shuffled F statistic exceeds the observed F
    % FDR correction is then applied across frequencies
    fprintf('\n[%s] Running permutation RM-ANOVA (%d perms)...\n', target_state, n_perm);
    p_vals    = nan(num_frequencies, 1);
    F_obs     = nan(num_frequencies, 1);
    eta2p_obs = nan(num_frequencies, 1);

    for f = 1:num_frequencies
        Y = squeeze(data_matrix(:, :, f)); % [animals × areas]
        valid   = all(~isnan(Y), 2);
        Y       = Y(valid, :);
        n_valid = size(Y, 1);
        if n_valid < 3, continue; end

        grand_mean   = mean(Y, 'all');
        area_means   = mean(Y, 1);
        animal_means = mean(Y, 2);

        SS_total   = sum((Y - grand_mean).^2, 'all');
        SS_areas   = n_valid   * sum((area_means   - grand_mean).^2);
        SS_animals = num_areas * sum((animal_means - grand_mean).^2);
        SS_error   = SS_total - SS_areas - SS_animals;

        df_areas  = num_areas - 1;
        df_error  = (n_valid - 1) * (num_areas - 1);
        MS_areas  = SS_areas / df_areas;
        MS_error  = SS_error / df_error;

        if MS_error <= 0 || isnan(MS_error)
            F_obs(f) = Inf; p_vals(f) = 0; continue;
        end
        F_obs(f)     = MS_areas / MS_error;
        eta2p_obs(f) = SS_areas / (SS_areas + SS_error);

        % Permutation distribution
        perm_F = zeros(n_perm, 1);
        for perm = 1:n_perm
            Yp = Y;
            for a = 1:n_valid
                Yp(a,:) = Y(a, randperm(num_areas));
            end
            gm   = mean(Yp,'all');
            am_r = mean(Yp,1);
            am_a = mean(Yp,2);
            SSt  = sum((Yp-gm).^2,'all');
            SSas = n_valid   * sum((am_r-gm).^2);
            SSam = num_areas * sum((am_a-gm).^2);
            MSe  = (SSt - SSas - SSam) / df_error;
            if MSe > 0
                perm_F(perm) = (SSas / df_areas) / MSe;
            end
        end
        p_vals(f) = (sum(perm_F >= F_obs(f)) + 1) / (n_perm + 1);

        if mod(f,20)==0
            fprintf('  f %d/%d  p=%.4f\n', f, num_frequencies, p_vals(f));
        end
    end

    % FDR correction (Benjamini-Hochberg) - applied across all tested frequencies
    p_fdr     = mafdr(p_vals, 'BHFDR', true);
    sig_freqs = find(p_fdr < alpha);
    fprintf('[%s] Significant frequencies after FDR: %d / %d\n', ...
        target_state, numel(sig_freqs), sum(~isnan(p_vals)));

    % Writes omnibus ANOVA results to CSV
    fid = fopen(omnibus_csv_file,'a');
    for f = 1:num_frequencies
        if isnan(p_vals(f)), continue; end
        fprintf(fid,'%s,%.4f,%.4f,%.6f,%.6f,%d,%.4f\n', ...
            target_state, foi(f), F_obs(f), p_vals(f), p_fdr(f), p_fdr(f)<alpha, eta2p_obs(f));
    end
    fclose(fid);

    % POST-HOC PAIRWISE COMPARISONS (Holm-Bonferroni)
    % Only runs at frequencies where the omnibus ANOVA was significant following FDR
    % For each frequency every region pair is compared with a paired t-test on the animal-level difference
    % Then Holm-Bonferroni corrects the six pairwise p-values for multiple comparisons
    p_posthoc_holm = nan(n_pairs, num_frequencies);
    t_crit = tinv(1 - alpha/2, num_animals - 1);

    if ~isempty(sig_freqs)
        fprintf('[%s] Running post-hoc at %d significant frequencies...\n', ...
            target_state, numel(sig_freqs));

        for f_idx = sig_freqs'
            freq_hz = foi(f_idx);
            Y_ph    = squeeze(data_matrix(:, :, f_idx));   % [animals × areas]
            valid   = all(~isnan(Y_ph), 2);
            Y_ph    = Y_ph(valid, :);
            n_avail = size(Y_ph, 1);
            if n_avail < 3, continue; end

            p_raw      = nan(n_pairs,1);
            mean_diffs = nan(n_pairs,1);
            se_diffs   = nan(n_pairs,1);
            t_stats    = nan(n_pairs,1);
            eta2p      = nan(n_pairs,1);

            % Paired t-test for every region pair
            for i = 1:n_pairs
                r1 = region_pairs(i,1);  r2 = region_pairs(i,2);
                d  = Y_ph(:,r1) - Y_ph(:,r2);
                mean_diffs(i) = mean(d);
                se_diffs(i)   = std(d) / sqrt(numel(d));
                [~, p_raw(i), ~, stats] = ttest(d);
                t_stats(i) = stats.tstat;
                eta2p(i)   = t_stats(i)^2 / (t_stats(i)^2 + (numel(d)-1));
            end

            % Holm-Bonferroni correction
            [p_sorted, sort_idx] = sort(p_raw);
            holm_fac       = (n_pairs + 1 - (1:n_pairs)');
            p_holm_sorted  = min(cummax(p_sorted .* holm_fac), 1);
            p_holm         = nan(n_pairs,1);
            p_holm(sort_idx) = p_holm_sorted;
            p_posthoc_holm(:, f_idx) = p_holm;

            % Writes post-hoc results to CSV
            fid = fopen(posthoc_csv_file,'a');
            for i = 1:n_pairs
                ci_lo  = mean_diffs(i) - t_crit * se_diffs(i);
                ci_hi  = mean_diffs(i) + t_crit * se_diffs(i);
                is_sig = p_holm(i) < alpha;
                fprintf(fid,'%.4f,%s,%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.6f,%.6f,%d\n', ...
                    freq_hz, target_state, ...
                    areas{region_pairs(i,1)}, areas{region_pairs(i,2)}, ...
                    n_avail, mean_diffs(i), se_diffs(i), ci_lo, ci_hi, ...
                    t_stats(i), n_avail-1, eta2p(i), p_raw(i), p_holm(i), is_sig);
            end
            fclose(fid);
        end
    end

    % Accumulates three-value significance matrix for this state heatmap:
    % NaN = frequency not tested post-hoc (omnibus not significant)
    % 0 = post-hoc tested but not significant
    % 1 = significant after Holm correction
    sig_matrix = nan(n_pairs, num_frequencies);
    tested     = ~isnan(p_posthoc_holm);
    sig_matrix(tested) = double(p_posthoc_holm(tested) < alpha);
    all_sigmat{state_idx} = sig_matrix;

end

% Figure 1: 4-panel mean ± SEM line plots (one panel per state)
% Shows each region's mean occurrence line plot with shaded SEM band
fig_line = figure('Visible','off','Position',[50 50 1600 1000]);

for state_idx = 1:num_states
    ax = subplot(2, 2, state_idx);
    hold(ax,'on');
    set(ax,'XScale','log','Box','on');

    mean_by_area = all_mean{state_idx};
    sem_by_area  = all_sem{state_idx};
    leg_h = gobjects(num_areas,1);

    for a_idx = 1:num_areas
        x = foi(:);
        y = mean_by_area(a_idx,:)';
        e = sem_by_area(a_idx,:)';
        ok = ~isnan(y) & ~isnan(e);

        % SEM shaded bands
        fill(ax, [x(ok); flipud(x(ok))], ...
                 [y(ok)-e(ok); flipud(y(ok)+e(ok))], ...
             colors{a_idx}, 'FaceAlpha',0.18, 'EdgeColor','none', ...
             'HandleVisibility','off');

        % Mean line
        leg_h(a_idx) = semilogx(ax, foi, mean_by_area(a_idx,:), ...
            'LineWidth',2.2, 'Color',colors{a_idx}, ...
            'DisplayName',areas{a_idx});
    end

    set(ax, 'XTick',[1 2 4 8 10 20 40], ...
            'XTickLabel',{'1','2','4','8','10','20','40'}, ...
            'XLim',[0.9 42], 'FontSize',15);
    xlabel(ax,'Frequency (Hz)','FontSize',19);
    ylabel(ax,'LAVI','FontSize',19);
    title(ax, states{state_idx},'FontSize',20,'FontWeight','bold');
    % legend(leg_h,'Location','northeastoutside','FontSize',10); % COMMENTED OUT
    grid(ax,'on');
    hold(ax,'off');
end

lineplot_file = fullfile(save_path,'LAVI_grandavg_allStates.png');
print(fig_line,'-dpng','-r300', lineplot_file);
fprintf('\nLine-plot figure saved to: %s\n', lineplot_file);
close(fig_line);

% Figure 2: 4-panel significance heatmaps (one panel per state)
% Rows: region pairs
% Columns: frequencies
fig_heat = figure('Visible','off','Position',[50 50 1600 1000]);

for state_idx = 1:num_states
    ax = subplot(2, 2, state_idx);

    sig_mat  = all_sigmat{state_idx};
    plot_mat = sig_mat;
    plot_mat(isnan(sig_mat)) = -1;

    ytick_positions = 0.5:n_pairs-0.5;

    imagesc(ax, foi, 1:n_pairs, plot_mat);
    set(ax,'YDir','normal','XScale','log','TickLength',[0 0]);
    axis(ax,'tight');

    colormap(ax, [0.75 0.75 0.75;   % grey   = untested (omnibus not significant)
                  1.00 1.00 1.00;   % white  = tested but post-hoc not significant
                  0.00 0.00 0.00]); % black  = post-hoc significant
    caxis(ax,[-1 1]);

    set(ax, 'XTick',[1 2 4 8 10 20 40], ...
            'XTickLabel',{'1','2','4','8','10','20','40'}, ...
            'YTick',ytick_positions, ...
            'YTickLabel',pair_names_cell, ...
            'XLim',[0.9 42], ...
            'YLim',[0.5, n_pairs+0.5], ...
            'FontSize',15);
    xlabel(ax,'Frequency (Hz)','FontSize',18);
    title(ax, states{state_idx},'FontSize',20,'FontWeight','bold');
    grid(ax,'on'); box(ax,'on');
end

cb_ax = axes(fig_heat,'Position',[0.93 0.15 0.015 0.70],'Visible','off');
colormap(cb_ax,[0.75 0.75 0.75; 1 1 1; 0 0 0]);
caxis(cb_ax,[-1 1]);
cbar = colorbar(cb_ax,'Location','eastoutside');
set(cbar,'Ticks',[-1 0 1], ...
    'TickLabels',{'Untested (RM-ANOVA NS)','Tested (NS)','Significant (p<0.05)'}, ...
    'FontSize',11);

heatmap_file = fullfile(save_path,'LAVI_significance_heatmap_allStates.png');
print(fig_heat,'-dpng','-r300', heatmap_file);
fprintf('Heatmap figure saved to: %s\n', heatmap_file);
close(fig_heat);

fprintf('\n=== DONE ===\n');
fprintf('Omnibus ANOVA CSV : %s\n', omnibus_csv_file);
fprintf('Post-hoc CSV      : %s\n', posthoc_csv_file);
