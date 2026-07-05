% Finds per-frequency significant differences in rhythmicity across states
% Computes grand-averaged LAVI profiles across states per region
% Computes permutation repeated-measures ANOVA with FDR Correction for LAVI profiles across brain states by shuffling state labels within each animal
% States are the repeated within-subject factor
% Obtains an observed F statistic, estimates permutation F distribution, computes permutation p-values, and then applies FDR correction.
% Post-hoc pairwise comparisons with Holm-Bonferroni correction for each significant frequency across animals

clear; close all; clc;
rng(0);

% Data parameters
areas       = {'A1','BLA','HPC','PFC'};
states      = {'AW','QW','NREM','REM'};
animal_ids  = {'r14','r16','r19','r20'};
condition   = '1_habituation';

num_areas    = numel(areas);
num_states   = numel(states);
num_animals  = numel(animal_ids);

% LAVI parameters
foi             = logspace(log10(1), log10(40), 96);
num_frequencies = numel(foi);
alpha           = 0.05;
n_perm          = 10000;

DATApath  = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';
save_path = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/LAVI profiles across states';

if ~exist(save_path,'dir'), mkdir(save_path); end

% State colours (AW / QW / NREM / REM)
colors = {[0.80 0.60 0.00], ...   % AW  – amber
          [0.00 0.60 0.00], ...   % QW  – green
          [0.00 0.00 1.00], ...   % NREM – blue
          [1.00 0.00 0.00]};      % REM  – red

pair_names_cell = {'AW vs QW','AW vs NREM','AW vs REM', ...
                   'QW vs NREM','QW vs REM','NREM vs REM'};
state_pairs = nchoosek(1:num_states, 2);   % [6×2]
n_pairs     = size(state_pairs, 1);

% Initialise results CSV files
% omnibus_csv_file: permutation based RM-ANOVA results
% posthoc_csv_file: post-hoc pairwise comparison results
omnibus_csv_file = fullfile(save_path, 'omnibus_anova_lavi_across_states.csv');
fid = fopen(omnibus_csv_file,'w');
fprintf(fid,'%s\n', strjoin({'Area','Frequency_Hz','F_Observed', ...
    'P_permutation','P_FDR','Significant_FDR','Eta2p_Omnibus'}, ','));
fclose(fid);

posthoc_csv_file = fullfile(save_path, 'posthoc_results_lavi_across_states.csv');
fid = fopen(posthoc_csv_file,'w');
fprintf(fid,'%s\n', strjoin({'Frequency_Hz','Area','State1','State2', ...
    'N_Animals','Mean_Diff','SE_Diff','CI_Lower_95','CI_Upper_95', ...
    'T_Statistic','DF','Partial_Eta_Squared','P_raw','P_holm','Significant'}, ','));
fclose(fid);

% Pre-allocate storage: [areas × states × animals × frequencies]
lavi_data = cell(num_areas, 1);

% Loads and averages LAVI profiles
for area_idx = 1:num_areas
    area = areas{area_idx};
    fprintf('\nLoading LAVI data for area: %s\n', area);

    lavi_area = cell(num_states, num_animals);   % {state, animal}

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
                warning('Missing file: %s', file_cut);
                lavi_area{state_idx, animal_idx} = nan(1, num_frequencies);
                continue;
            end

            loaded = load(file_cut, 'LAVI_matrix');
            if ~isfield(loaded,'LAVI_matrix')
                warning('LAVI_matrix missing in: %s', file_cut);
                lavi_area{state_idx, animal_idx} = nan(1, num_frequencies);
                continue;
            end

            LAVI_matrix = loaded.LAVI_matrix; % [channels × freqs × epochs]
            [num_chan, num_freq_file, num_ep] = size(LAVI_matrix);

            % Averages across epochs then across channels (omits NaN) - produces one LAVI profile per animal per state per area
            chan_mean = nan(num_chan, num_frequencies);
            for ch = 1:num_chan
                ep_vals = squeeze(LAVI_matrix(ch, :, :)); % [freqs × epochs]
                chan_mean(ch, :) = mean(ep_vals, 2, 'omitnan')';
            end
            lavi_area{state_idx, animal_idx} = mean(chan_mean, 1, 'omitnan');
        end
    end

    lavi_data{area_idx} = lavi_area;
    fprintf('  Done (%d states × %d animals loaded)\n', num_states, num_animals);
end

% Stores per-area results for plotting
all_mean   = cell(num_areas, 1); % {area}[states × freqs]
all_sem    = cell(num_areas, 1);
all_sigmat = cell(num_areas, 1); % significance matrix [pairs × freqs]

for area_idx = 1:num_areas
    area      = areas{area_idx};
    lavi_area = lavi_data{area_idx};

    % Builds matrix [animals × states × freqs]
    data_matrix = nan(num_animals, num_states, num_frequencies);
    for s = 1:num_states
        for a = 1:num_animals
            v = lavi_area{s, a};
            if ~isempty(v) && ~all(isnan(v))
                data_matrix(a, s, :) = v;
            end
        end
    end

    % Calculate grand mean and SEM across animals
    mean_by_state = nan(num_states, num_frequencies);
    sem_by_state  = nan(num_states, num_frequencies);
    for s = 1:num_states
        tmp = squeeze(data_matrix(:, s, :)); % [animals × freqs]
        n_valid_animals = sum(~isnan(tmp), 1);
        mean_by_state(s, :) = mean(tmp, 1, 'omitnan');
        sem_by_state(s, :)  = std(tmp, 0, 1, 'omitnan') ./ sqrt(n_valid_animals);
    end
    all_mean{area_idx} = mean_by_state;
    all_sem{area_idx}  = sem_by_state;

    % PERMUTATION RM-ANOVA (with FDR)
    % For each frequency computes a one-way repeated-measures ANOVA across states (F statistic)
    % Then build a null distribution by repeatedly shuffling the state labels for each animal and recomputing the F statistic
    % The permutation p-value is how often the shuffled F statistic exceeds the observed F
    % FDR correction is then applied across frequencies
    fprintf('\n[%s] Running permutation RM-ANOVA (%d perms)...\n', area, n_perm);
    p_vals    = nan(num_frequencies, 1);
    F_obs     = nan(num_frequencies, 1);
    eta2p_obs = nan(num_frequencies, 1);

    for f = 1:num_frequencies
        Y = squeeze(data_matrix(:, :, f)); % [animals × states]
        valid   = all(~isnan(Y), 2);
        Y       = Y(valid, :);
        n_valid = size(Y, 1);
        if n_valid < 3, continue; end

        grand_mean   = mean(Y, 'all');
        state_means  = mean(Y, 1);
        animal_means = mean(Y, 2);

        SS_total   = sum((Y - grand_mean).^2, 'all');
        SS_states  = n_valid  * sum((state_means  - grand_mean).^2);
        SS_animals = num_states * sum((animal_means - grand_mean).^2);
        SS_error   = SS_total - SS_states - SS_animals;

        df_states = num_states - 1;
        df_error  = (n_valid - 1) * (num_states - 1);
        MS_states = SS_states / df_states;
        MS_error  = SS_error  / df_error;

        if MS_error <= 0 || isnan(MS_error)
            F_obs(f) = Inf; p_vals(f) = 0; continue;
        end
        F_obs(f)     = MS_states / MS_error;
        eta2p_obs(f) = SS_states / (SS_states + SS_error);

        % Permutation distribution
        perm_F = zeros(n_perm, 1);
        for perm = 1:n_perm
            Yp = Y;
            for a = 1:n_valid
                Yp(a,:) = Y(a, randperm(num_states));
            end
            gm   = mean(Yp,'all');
            sm   = mean(Yp,1);
            am   = mean(Yp,2);
            SSt  = sum((Yp-gm).^2,'all');
            SSts = n_valid * sum((sm-gm).^2);
            SSam = num_states * sum((am-gm).^2);
            MSe  = (SSt - SSts - SSam) / df_error;
            if MSe > 0
                perm_F(perm) = (SSts / df_states) / MSe;
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
        area, numel(sig_freqs), sum(~isnan(p_vals)));

    % Writes omnibus ANOVA results to CSV
    fid = fopen(omnibus_csv_file,'a');
    for f = 1:num_frequencies
        if isnan(p_vals(f)), continue; end
        fprintf(fid,'%s,%.4f,%.4f,%.6f,%.6f,%d,%.4f\n', ...
            area, foi(f), F_obs(f), p_vals(f), p_fdr(f), p_fdr(f)<alpha, eta2p_obs(f));
    end
    fclose(fid);

    % POST-HOC PAIRWISE COMPARISONS (Holm-Bonferroni)
    % Only runs at frequencies where the omnibus ANOVA was significant following FDR
    % For each frequency every state pair is compared with a paired t-test on the animal-level difference
    % Then Holm-Bonferroni corrects the six pairwise p-values for multiple comparisons
    p_posthoc_holm = nan(n_pairs, num_frequencies);
    t_crit = tinv(1 - alpha/2, num_animals - 1);

    if ~isempty(sig_freqs)
        fprintf('[%s] Running post-hoc at %d significant frequencies...\n', ...
            area, numel(sig_freqs));

        for f_idx = sig_freqs'
            freq_hz   = foi(f_idx);
            Y_ph      = squeeze(data_matrix(:, :, f_idx));   % [animals × states]
            valid     = all(~isnan(Y_ph), 2);
            Y_ph      = Y_ph(valid, :);
            n_avail   = size(Y_ph, 1);
            if n_avail < 3, continue; end

            p_raw      = nan(n_pairs,1);
            mean_diffs = nan(n_pairs,1);
            se_diffs   = nan(n_pairs,1);
            t_stats    = nan(n_pairs,1);
            eta2p      = nan(n_pairs,1);

            % Paired t-test for every state pair
            for i = 1:n_pairs
                s1 = state_pairs(i,1);  s2 = state_pairs(i,2);
                d  = Y_ph(:,s1) - Y_ph(:,s2);
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
                ci_lo = mean_diffs(i) - t_crit * se_diffs(i);
                ci_hi = mean_diffs(i) + t_crit * se_diffs(i);
                is_sig = p_holm(i) < alpha;
                fprintf(fid,'%.4f,%s,%s,%s,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.6f,%.6f,%d\n', ...
                    freq_hz, area, states{state_pairs(i,1)}, states{state_pairs(i,2)}, ...
                    n_avail, mean_diffs(i), se_diffs(i), ci_lo, ci_hi, ...
                    t_stats(i), n_avail-1, eta2p(i), p_raw(i), p_holm(i), is_sig);
            end
            fclose(fid);
        end
    end

    % Accumulates three-value significance matrix for this area heatmap:
    % NaN = frequency not tested post-hoc (omnibus not significant)
    % 0 = post-hoc tested but not significant
    % 1 = significant after Holm correction
    sig_matrix = nan(n_pairs, num_frequencies);
    tested     = ~isnan(p_posthoc_holm);
    sig_matrix(tested) = double(p_posthoc_holm(tested) < alpha);
    all_sigmat{area_idx} = sig_matrix;

end

% Figure 1: 4-panel mean ± SEM line plots (one panel per area)
% Shows each state's mean LAVI profile with shaded SEM band
fig_line = figure('Visible','off','Position',[50 50 1600 1000]);

for area_idx = 1:num_areas
    ax = subplot(2, 2, area_idx);
    hold(ax,'on');
    set(ax,'XScale','log','Box','on');

    mean_by_state = all_mean{area_idx};
    sem_by_state  = all_sem{area_idx};
    leg_h = gobjects(num_states,1);

    for s = 1:num_states
        x = foi(:);
        y = mean_by_state(s,:)';
        e = sem_by_state(s,:)';
        ok = ~isnan(y) & ~isnan(e);

        % SEM shaded bands
        fill(ax, [x(ok); flipud(x(ok))], ...
                 [y(ok)-e(ok); flipud(y(ok)+e(ok))], ...
             colors{s}, 'FaceAlpha',0.18, 'EdgeColor','none', ...
             'HandleVisibility','off');

        % Mean line
        leg_h(s) = semilogx(ax, foi, mean_by_state(s,:), ...
            'LineWidth',2.2, 'Color',colors{s}, ...
            'DisplayName',states{s});
    end

    set(ax, 'XTick',[1 2 4 8 10 20 40], ...
            'XTickLabel',{'1','2','4','8','10','20','40'}, ...
            'XLim',[0.9 42], 'FontSize',15);
    xlabel(ax,'Frequency (Hz)','FontSize',19);
    ylabel(ax,'LAVI','FontSize',19);
    title(ax, areas{area_idx},'FontSize',20,'FontWeight','bold');
    % legend(leg_h,'Location','northeastoutside','FontSize',10); % COMMENTED OUT
    grid(ax,'on');
    hold(ax,'off');
end

lineplot_file = fullfile(save_path,'LAVI_grandavg_allAreas.png');
print(fig_line,'-dpng','-r300', lineplot_file);
fprintf('\nLine-plot figure saved to: %s\n', lineplot_file);
close(fig_line);

% Figure 2: 4-panel significance heatmaps (one panel per area)
% Rows: state pairs
% Columns: frequencies
fig_heat = figure('Visible','off','Position',[50 50 1600 1000]);

for area_idx = 1:num_areas
    ax = subplot(2, 2, area_idx);

    sig_mat   = all_sigmat{area_idx};
    plot_mat  = sig_mat;
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
    grid(ax,'on'); box(ax,'on');
end

cb_ax = axes(fig_heat,'Position',[0.93 0.15 0.015 0.70],'Visible','off');
colormap(cb_ax,[0.75 0.75 0.75; 1 1 1; 0 0 0]);
caxis(cb_ax,[-1 1]);
cbar = colorbar(cb_ax,'Location','eastoutside');
set(cbar,'Ticks',[-1 0 1], ...
    'TickLabels',{'Untested (RM-ANOVA NS)','Tested (NS)','Significant (p<0.05)'}, ...
    'FontSize',11);

heatmap_file = fullfile(save_path,'LAVI_significance_heatmap_allAreas.png');
print(fig_heat,'-dpng','-r300', heatmap_file);
fprintf('Heatmap figure saved to: %s\n', heatmap_file);
close(fig_heat);

fprintf('\n=== DONE ===\n');
fprintf('Omnibus ANOVA CSV : %s\n', omnibus_csv_file);
fprintf('Post-hoc CSV      : %s\n', posthoc_csv_file);
