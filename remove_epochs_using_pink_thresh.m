% Removes epochs from copy of LAVI/PINK .mat file based on mean pink noise on per-channel basis

clear; close all; clc;

% Set path
DATApath = '/Volumes/TOSHIBA_UO/Joao_LFP/raw dataset/Saved matrices/';

% Data parameters
animals = {'r14','r16','r19','r20'};
areas = {'A1','BLA','HPC','PFC'};
brain_states = {'AW','QW','REM','NREM'};
condition = '1_habituation';

% Set whether to overwrite original files (true) or generate a copy (false)
overwrite_original = false;

% LAVI parameters
foi = logspace(log10(1), log10(40), 96);
idx_0_6  = foi >= 0 & foi <= 6;
idx_6_40 = foi > 6 & foi <= 40;

range_0_6  = [0.35 0.45];
range_6_40 = [0.36 0.42];

% Loop over animals, brain areas, brain states:
for a = 1:numel(animals)
    animal_id = animals{a};

    for ar = 1:numel(areas)
        area = areas{ar};

        for bs = 1:numel(brain_states)
            brain_state = brain_states{bs};

            fprintf('\nProcessing: %s | %s | %s\n', animal_id, area, brain_state);

            % Builds input file path for this animal/area/state combination
            lavi_folder = fullfile([animal_id ' matrices'], ...
                [animal_id ' ' condition ' matrices'], area, [area ' ' brain_state], ...
                ['LAVI processing (' area ' ' brain_state ')']);

            inputFileName = fullfile(DATApath, lavi_folder, ...
                sprintf('%s_%s_%s_%s_lavi_pink.mat', animal_id, condition, area, brain_state));

            if overwrite_original
                outputFileName = inputFileName;
            else
                outputFileName = fullfile(DATApath, lavi_folder, ...
                    sprintf('%s_%s_%s_%s_lavi_pink_CUT.mat', animal_id, condition, area, brain_state));
            end

            if ~isfile(inputFileName)
                warning('Missing file: %s', inputFileName);
                continue;
            end

            try
                % Loads LAVI and PINK data matrices from file
                fprintf('Loading: %s\n', inputFileName);
                data = load(inputFileName);
                LAVI_matrix = data.LAVI_matrix;
                PINK_matrix = data.PINK_matrix;

                [nRep, nFreq, nChan, nEpoch] = size(PINK_matrix);
                fprintf('Data dimensions: Channels = %d | Frequencies = %d | Epochs = %d\n', ...
                    nChan, nFreq, nEpoch);

                % Averages pink noise across repetitions
                PINK_mean = mean(PINK_matrix, 1);
                PINK_mean = reshape(PINK_mean, [nFreq, nChan, nEpoch]);

                % Evaluates each channel/epoch against pink-noise:
                % Flags epoch-channel pairs whose mean pink-noise value falls outside either the acceptable range
                keep_mask = true(nChan, nEpoch);

                for ch = 1:nChan
                    for ep = 1:nEpoch
                        pink_spec = PINK_mean(:, ch, ep);

                        mean_0_6  = mean(pink_spec(idx_0_6), 'omitnan');
                        mean_6_40 = mean(pink_spec(idx_6_40), 'omitnan');

                        cond_0_6  = mean_0_6  >= range_0_6(1)  && mean_0_6  <= range_0_6(2);
                        cond_6_40 = mean_6_40 >= range_6_40(1) && mean_6_40 <= range_6_40(2);

                        if ~(cond_0_6 && cond_6_40)
                            keep_mask(ch, ep) = false;
                        end
                    end
                end

                % Applies exclusion mask
                % Sets NaN for rejected epoch-channel pairs
                for ch = 1:nChan
                    for ep = 1:nEpoch
                        if ~keep_mask(ch, ep)
                            LAVI_matrix(ch, :, ep) = NaN;
                            PINK_matrix(:, :, ch, ep) = NaN;
                        end
                    end
                end

                % Saves cleaned data
                save(outputFileName, 'LAVI_matrix', 'PINK_matrix', 'keep_mask');
                fprintf('Saved cleaned file: %s\n', outputFileName);
                fprintf('Excluded epoch-channel pairs set to NaN\n');

            catch ME
                warning('Failed for %s | %s | %s: %s', animal_id, area, brain_state, ME.message);
            end
        end
    end
end
