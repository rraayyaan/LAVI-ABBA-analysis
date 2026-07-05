# LAVI-ABBA-analysis

Pipeline for LAVI-ABBA analysis of rhythmicity and spectral organisation of electrophysiological data, based on the approach devised by Karvat et al. (2024).

Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. C., & Landau, A. N. (2024). Universal rhythmic architecture uncovers distinct modes of neural dynamics. https://doi.org/10.1101/2024.12.05.627113

# The original LAVI toolbox

The LAVI toolbox is used to generate rhythmicity profile and automatic band detection of electrophysiological neural data. It is first developed and introduced in Karvat et al. (2024), Universal rhythmic architecture uncovers distinct modes of neural dynamics.  

LAVI quantifies the rhythmicity of the neural signal by comparing phase coherence at a fixed time lag.
ABBA statistically identifies rhythmically sustained and transient bands of electrophysiological activity.

# This implementation of the LAVI toolbox

This pipeline computes the rhythmicity profiles of LFP epoch data using LAVI and statistically detects rhythmically sustained and rhythmically transient bands using ABBA.

The pipeline also includes statistical analyses that identify significant differences in rhythmicity and spectral band distributions between states and regions.

THE PIPELINE USES THE FOLLOWING ORDER:

First complete the preprocessing and epoch creation with python using: LAVI-ABBA_data_preprocessing.ipynb

THEN using matlab:

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
- Review ABBA plots across epochs (non-averaged) - LAVI profile and ABBA bands: plot_ABBA_across_epochs
- Review ABBA bands across channels (averaged across epochs): plot_ABBA_bands_across_channels
- Review ABBA bands across epochs (non-averaged): plot_ABBA_bands_across_epochs

# References

When using this toolbox please cite the following publications:

**Core rhythmicity method (LAVI-ABBA)**

Karvat, G., Crespo-García, M., Vishne, G., Anderson, M., & Landau, A. N. (2024). Universal rhythmic architecture uncovers distinct modes of neural dynamics. bioRxiv. https://doi.org/10.1101/2024.12.05.627113

**Signal processing & preprocessing**

Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J.-M. (2011). FieldTrip: Open source software for advanced analysis of MEG, EEG, and invasive electrophysiological data. Computational Intelligence and Neuroscience, 2011, 156869. https://doi.org/10.1155/2011/156869

Siegle, J. H., López, A. C., Patel, Y. A., Abramov, K., Ohayon, S., & Voigts, J. (2017). Open Ephys: An open-source, plugin-based platform for multichannel electrophysiology. Journal of Neural Engineering, 14(4), 045003. https://doi.org/10.1088/1741-2552/aa5eea

Welch, P. D. (1967). The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms. IEEE Transactions on Audio and Electroacoustics, 15(2), 70–73. https://doi.org/10.1109/TAU.1967.1161901

**Surrogate data generation**

Venema, V. (2023). Surrogate time series and fields [Computer software]. MATLAB Central File Exchange. Retrieved January 17, 2023, from https://www.mathworks.com/matlabcentral/fileexchange/4783-surrogate-time-series-and-fields

Venema, V., Ament, F., & Simmer, C. (2006). A stochastic iterative amplitude adjusted Fourier transform algorithm with improved accuracy. Nonlinear Processes in Geophysics, 13(3), 321–328. https://doi.org/10.5194/npg-13-321-2006

**Statistical methods**

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B, 57(1), 289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6(2), 65–70.

Good, P. I. (2005). Permutation, Parametric, and Bootstrap Tests of Hypotheses (3rd ed.). Springer. https://doi.org/10.1007/b138696

**Software**

The MathWorks, Inc. (n.d.). MATLAB Curve Fitting Toolbox and Statistics and Machine Learning Toolbox. Natick, MA: The MathWorks, Inc.
Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., ... & van Mulbregt, P. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2
