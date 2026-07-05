% Generates and saves significance limits based on different durations, sampling frequency, and aperiodic slope of pink noise power.
% Rationale: not to have to calculate every time for every session.
% With current parameteres, takes ~8 hours
% The parameters used to generate the table are saved in pmtrSIG

clear; close all; clc;

output_path = '/Users/rayan_1/Documents/MATLAB/LAVI (Rayan)'; % output folder
output_file = 'SIGLIM6'; % output file name. The original file is called SIGLIM

FS = 1000; % use 1000
DUR = 240; % 4-minute epochs
B = -3:0.1:-0.3; % aperiodic slope grid
a = 1;
reps = 100;

f = logspace(log10(1), log10(40), 96); % frequencies of interest
width = 5;
lag = 1.5;

N.freq = length(f);
N.dur = length(DUR);
N.fs = length(FS);
N.b = length(B);
SIGLIM = zeros(N.dur, N.fs, N.b,N.freq,2);
prev = fprintf(' ');
tt = tic;
for di = 1:N.dur
    T = DUR(di);
    for fi = 1:N.fs
        fs = FS(fi);
        n = round(fs*T);
        for bi = 1:N.b
            b = B(bi);
            coefsIntoSurr = get_pink_iafft_coefs_random_ap (n,fs,f,a,b);
            LAVI = zeros(reps, N.freq);            
            for ri = 1:reps
                if ri<10,r=['0' num2str(ri)]; else r=num2str(ri); end
                ttt = toc(tt);
                dak = floor(ttt/60); 
                sec = round(ttt-dak*60); if sec<10, sec = ['0' num2str(sec)]; else sec = num2str(sec); end
                dak = num2str(dak);
                prev = dispRMVprev (['Running T = ' num2str(T) '(' num2str(di) '/' num2str(N.dur) ')'...
                    ', fs = ' num2str(fs) '(' num2str(fi) '/' num2str(N.fs) ')'...
                    ', b = ' num2str(b) '(' num2str(bi) '/' num2str(N.b) ')'...
                    ', rep = ' r '/' num2str(reps)...
                    '. So far it took ' dak ':' sec ' minutes'],prev);
                in = rand(1,n);
                pink = iaaft_loop_1d(coefsIntoSurr, sort(in));
                cfg = [];
                cfg.foi = f;
                cfg.fs = fs;
                cfg.lag = lag;
                cfg.width = width;
                cfg.verbose = 0;
                LAVI(ri,:) = Prepare_LAVI(cfg,pink);
            end
            sig = [min(LAVI); max(LAVI)]; % the siginificance levels for EACH FREQUENCY
            SIGLIM(di,fi,bi,:,:) = sig';
        end
    end
end
disp('.');

% Prepare the output and save
[~,madeBy,~]=fileparts(matlab.desktop.editor.getActiveFilename);
pmtrSIG         = [];
pmtrSIG.B       = B;
pmtrSIG.DUR     = DUR;
pmtrSIG.FS      = FS;
pmtrSIG.f       = f;
pmtrSIG.dimord  = 'dur_fs_b_freq_min/max';
pmtrSIG.lag     = lag;
pmtrSIG.width   = width;
pmtrSIG.script  = madeBy;
save(fullfile(output_path, output_file),'SIGLIM', 'pmtrSIG');