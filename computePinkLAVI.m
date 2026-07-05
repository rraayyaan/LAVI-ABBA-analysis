function PINK = computePinkLAVI(cfg,data)
% Generates pink noise to estimate the significance level of ABBA detected spectral bands.

% Input:
% Data: N_channel x N_time array.

% cfg.Pink_reps : number of simulations created per channel. Default = 20.
% cfg.durs: the duration of each simulation. Default = duration of the data.

% Sets the defaults
if ~isfield(cfg,'Pink_reps'); cfg.Pink_reps = 100; end % number of repetitions of pink noise simulations
if ~isfield(cfg,'foi'); cfg.foi = 10.^(log10(0.5):0.025:log10(120)); end
if ~isfield(cfg,'fs'); cfg.fs = 1000; end
if ~isfield(cfg,'lag'); cfg.lag = 1.5; end
if ~isfield(cfg,'width'); cfg.width = 5; end
T = size(data,2)/cfg.fs;
if ~isfield(cfg,'durs'); cfg.durs = T; end % duration (in sec) of each simulation
if ~isfield(cfg,'thresholds'); cfg.thresholds = inf(size(data,1),1); end % default: no threshold per channel

if cfg.Pink_reps==0 || cfg.durs==0
    PINK = [];
else
    Pink_reps   = cfg.Pink_reps;
    durs        = cfg.durs;
    foi         = cfg.foi;
    fs          = cfg.fs;
    lag         = cfg.lag;
    width       = cfg.width;
    pmtr        = cfg; % to keep a copy of the original cfg
    w           = hanning(fs*2); % 2-sec window

    N.time      = size(data,2);
    N.chan      = size(data,1);
    N.freq      = length(foi);

    PINK = nan(Pink_reps, N.freq, N.chan); % dimord = rep_freq_chan
    tmax = min([durs, N.time]);
    data = data(:,1:floor(tmax*fs)); % take shorter duration then the original if requested
    N.time = size(data,2); % document the new number of samples

    % Defines a clean amplitude threshold (per channel) and remove peaks
    % Robust range of the background (excluding artifacts)
    thres = cfg.thresholds; % per channel threshold from cfg (default inf = no threshold)

    % For each channel makes a clean template by zeroing out big peaks
    data_clean = data; % copy
    for ch = 1:N.chan
        EEG = data(ch,:);

        % Finds samples where LFP exceeds threshold (in absolute amplitude)
        big = abs(EEG) > thres(ch);

        % Interpolates over those time points to keep the same length
        EEG_clean = EEG;
        if any(big)
            EEG_clean(big) = nan;
            EEG_clean = fillmissing(EEG_clean,'linear', 'SamplePoints',1:length(EEG_clean));
        end

        data_clean(ch,:) = EEG_clean;
    end

    % Overall timing
    t_.overall = tic;

    for ch = 1:size(data,1)
        t_.chan = tic;

        % Uses the clean version for surrogate spectrum generation
        EEG_clean = data_clean(ch,:);
        EEG_clean = EEG_clean - mean(EEG_clean); % demean

        % Computes surrogate spectrum coefficients with the clean data
        [coefsIntoSurr, fito] = get_pink_iafft_coefs_pow (EEG_clean,w,foi,fs);

        surrLAVI = nan(pmtr.Pink_reps,length(foi));
        rng = prctile(EEG_clean,99.9)-prctile(EEG_clean,0.01); % uses range if using random amplitude
        offset = prctile(EEG_clean,0.01); % values that are offset to the original LFP

        fprintf('Running PINK ANALYSIS Channel %d/%d', ch, N.chan);
        prev_len = 0;

        for ri = 1:pmtr.Pink_reps
            str = sprintf(' repeat %d/%d. So far analysis took %.1f seconds', ri, pmtr.Pink_reps, toc(t_.chan));

            fprintf(repmat('\b', 1, prev_len));
            fprintf(str);
            prev_len = length(str);

            % Runs the analysis
            pinkNoise = iaaft_loop_1d(coefsIntoSurr, sort(rand(size(EEG_clean))));

            % Optionally scale the surrogate to the clean channel's range (uncomment if needed)
            % pinkNoise = pinkNoise * rng + offset;

            for fi = 1:N.freq
                f = foi(fi);
                spectrum = waveletLight(pinkNoise, fs, f, width);
                surrLAVI(ri,fi) = compute_lavi(spectrum,fs,f,lag);
            end
        end

        fprintf('\n');

        PINK(:,:,ch) = surrLAVI;
    end

    fprintf('Total analysis took %.1f seconds\n', toc(t_.overall));
    fprintf('\n\n');
end
end
