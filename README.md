# LAVI-ABBA-analysis
Pipeline for LAVI-ABBA analysis of rhythmicity and spectral organisation of electrophysiological data, based on the approach devised by Karvat et al. (2024).

This pipeline computes the rhythmicity profiles of LFP epoch data using LAVI and statistically detects rhythmically sustained and rhythmically transient bands using ABBA.
The pipeline also includes statistical analyses that identify significant differences in rhythmicity and spectral band distributions between states and regions.

Run the scripts in the following order:

1. Filter 50 Hz noise from epoch dataset: apply_50hz_filter
2. Compute LAVI profile and pink noise simulations: calculate_filtered_lavi_and_pink

3. Review LAVI and ABBA plots to identify unviable channels: plot_LAVI_and_ABBA_across_channels
4. Remove unviable channels: remove_channels
5. Remove unviable epochs per-channel using pink noise thresholds: remove_epochs_using_pink_thresh

6. Plot raw power (single plot per region with all states): raw_power_across_states
6. Plot raw power (single plot per state with all regions): raw_power_across_regions
7. Plot normalised power (single plot per region with all states): norm_power_across_states
7. Plot normalised power (single plot per state with all regions): norm_power_across_regions

8. Compute and plot significant differences in rhythmicity across states: LAVI_sig_perm_RM_ANOVA_across_states
8. Compute and plot significant differences in rhythmicity across regions: LAVI_sig_perm_RM_ANOVA_across_regions

9. Compute and plot significant differences in band presence across states: band_sig_perm_RM_ANOVA_across_states
9. Compute and plot significant differences in band presence across regions: band_sig_perm_RM_ANOVA_across_regions
10. Plot grand-averaged presence of all band-types (single subplot with all band-types): band_distributions

EXTRA: visualise non-averaged/averaged bands across epochs and averaged bands across channels
Review ABBA plots across epochs (non-averaged) - LAVI profile and ABBA bands: plot_ABBA_across_epochs
Review ABBA bands across channels (averaged across epochs): plot_ABBA_bands_across_channels
Review ABBA bands across epochs (non-averaged): plot_ABBA_bands_across_epochs
