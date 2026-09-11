%% breathing_freq_across_sessions_260629.m
% Per-session breathing frequency across a long multi-run recording, from
% the DLC tracking that is ALREADY ON DISK -- no video read, no SVD.
%
% For each run folder:
%   1. read the DLC csv, take the breathing dot (default dot4-x)
%   2. derive true fps from timestamps.csv (camera_timestamp_ns), NOT the
%      bogus nominal header fps (90) -- real rate is ~50.68 Hz here
%   3. high-pass + z-score the trace, detect breath peaks (findpeaks)
%   4. instantaneous freq = 1./ITI  (ITI = inter-peak interval, seconds)
%   5. summarise per run: mean +- std of 1/ITI
% Then plot the per-session mean +- std across recording time (wall clock).
%
% NEW standalone analysis -- does not touch any existing script.
% Outputs are written into folderPath alongside the data.
%
% RZ 2026-06-29

clear; clc;

%% ===================== USER PARAMS =====================
folderPath = 'Z:\Runqi\behavior\260505_breathing_wt\cam1';

BODYPART = 'dot4';   % DLC marker used as the breathing signal
AXIS     = '-x';     % 'x' | 'y' | '-x' | '-y' | 'magnitude'  (inverted: peaks = breath events)
LK_THR   = 0.50;     % likelihood threshold; below this -> treated as missing
N_DROP   = 30;       % drop first N frames of each run (stabilization)

HP_WIN_SEC = 1.0;    % moving-median high-pass window (s) to remove slow drift; 0 = off

% --- peak detection (per run, on z-scored trace) ---
MINPROM_FACTOR = 1.0;   % prominence = factor * std(breath_z)
MINDIST_SEC    = 0.10;  % min seconds between breaths (=> 10 Hz ceiling)
MINHEIGHT      = -Inf;  % min absolute peak height (z). -Inf = off
MINWIDTH_SEC   = 0;     % min peak width (s). 0 = off

BREATH_BAND_HZ = [1 14];  % plausible breath rate; 1./ITI outside this band is dropped from stats

TIME_AXIS = 'clock';   % 'clock' = elapsed hours from first run | 'run' = run number

EXAMPLE_RUN = 1;       % run number to draw a trace+peaks QC figure (0 = skip)
SAVE_OUT    = true;
%% =======================================================

assert(isfolder(folderPath), 'folderPath not found: %s', folderPath);

% Column layout in the DLC csv: each dot is [x y likelihood], dot1 starts at col 2.
dotNum = sscanf(BODYPART, 'dot%d');
assert(~isempty(dotNum), 'BODYPART must look like ''dot4''');
cx = 2 + 3*(dotNum-1);   cy = cx+1;   cl = cx+2;

% ---- enumerate run folders, sorted by run number ----
d = dir(fullfile(folderPath, 'cam1_*run*'));
d = d([d.isdir]);
runNum = nan(numel(d),1);
for i = 1:numel(d)
    tok = regexp(d(i).name, 'run(\d+)$', 'tokens', 'once');
    if ~isempty(tok), runNum(i) = str2double(tok{1}); end
end
keep = ~isnan(runNum);
d = d(keep);  runNum = runNum(keep);
[runNum, ord] = sort(runNum);  d = d(ord);
nRun = numel(d);
fprintf('Found %d run folders in %s\n', nRun, folderPath);
assert(nRun > 0, 'No run folders found.');

% ---- per-run accumulators ----
R = struct('run',{},'name',{},'fps',{},'recTime',{}, ...
           'nBreath',{},'meanFreq',{},'stdFreq',{},'medFreq',{}, ...
           'iqrFreq',{},'durSec',{},'nLowLk',{});

for i = 1:nRun
    rdir = fullfile(folderPath, d(i).name);

    % --- locate the DLC csv (the csv that is not timestamps.csv) ---
    csvs = dir(fullfile(rdir, '*.csv'));
    csvs = csvs(~strcmpi({csvs.name}, 'timestamps.csv'));
    dlcHit = csvs(contains({csvs.name}, 'DLC', 'IgnoreCase', true));
    if isempty(dlcHit)
        warning('run%03d: no DLC csv, skipped', runNum(i)); continue;
    end
    dlc_csv = fullfile(rdir, dlcHit(1).name);

    % --- true fps + recording start from timestamps.csv ---
    fps = NaN; recTime = NaT;
    tsf = fullfile(rdir, 'timestamps.csv');
    if isfile(tsf)
        Tts = readtable(tsf);
        if ismember('camera_timestamp_ns', Tts.Properties.VariableNames) && height(Tts) > 1
            fps = 1e9 / median(diff(Tts.camera_timestamp_ns));
        end
        if ismember('wall_time_s', Tts.Properties.VariableNames) && ~isempty(Tts.wall_time_s)
            recTime = datetime(Tts.wall_time_s(1), 'ConvertFrom', 'posixtime', 'TimeZone','UTC');
        end
    end
    if ~isfinite(fps)
        warning('run%03d: no usable timestamps, skipped', runNum(i)); continue;
    end
    if isnat(recTime)   % fall back to folder-name timestamp
        t2 = regexp(d(i).name, '_(\d{8})_(\d{6})_', 'tokens', 'once');
        if ~isempty(t2)
            recTime = datetime([t2{1} t2{2}], 'InputFormat','yyyyMMddHHmmss','TimeZone','UTC');
        end
    end

    % --- read DLC, build breathing signal ---
    raw = [];
    for attempt = 1:4   % network drive: retry transient read failures
        try
            raw = readmatrix(dlc_csv, 'NumHeaderLines', 3);
            break;
        catch ME
            if attempt == 4
                warning('run%03d: DLC read failed after retries (%s), skipped', ...
                        runNum(i), ME.message);
            else
                pause(0.5);
            end
        end
    end
    if isempty(raw), continue; end
    if size(raw,1) <= N_DROP || size(raw,2) < cl
        warning('run%03d: DLC csv too short/narrow, skipped', runNum(i)); continue;
    end
    raw(1:N_DROP, :) = [];

    xc = raw(:,cx);  yc = raw(:,cy);  pc = raw(:,cl);
    switch AXIS
        case 'x',         sig =  xc;
        case 'y',         sig =  yc;
        case '-x',        sig = -xc;
        case '-y',        sig = -yc;
        case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
        otherwise,        sig =  xc;
    end
    bad = pc < LK_THR;  nLow = sum(bad);  sig(bad) = NaN;
    if sum(~isnan(sig)) < 2
        warning('run%03d: signal all NaN, skipped', runNum(i)); continue;
    end
    sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');

    % --- detrend / high-pass / z-score ---
    sig = detrend(sig);
    if HP_WIN_SEC > 0
        w = max(3, round(HP_WIN_SEC * fps));
        sig = sig - movmedian(sig, w);
    end
    bz = (sig - mean(sig)) / std(sig);

    % --- peak detection -> ITI -> instantaneous frequency ---
    fp_args = {'MinPeakProminence', MINPROM_FACTOR*std(bz), ...
               'MinPeakDistance',   max(1, round(MINDIST_SEC*fps))};
    if isfinite(MINHEIGHT),  fp_args = [fp_args, {'MinPeakHeight', MINHEIGHT}]; end
    if MINWIDTH_SEC > 0,     fp_args = [fp_args, {'MinPeakWidth', round(MINWIDTH_SEC*fps)}]; end
    [~, pk] = findpeaks(bz, fp_args{:});

    ITI  = diff(pk) / fps;            % seconds
    freq = 1 ./ ITI;                  % Hz, instantaneous
    inb  = freq >= BREATH_BAND_HZ(1) & freq <= BREATH_BAND_HZ(2);
    freq = freq(inb);

    R(end+1) = struct( ...
        'run', runNum(i), 'name', d(i).name, 'fps', fps, 'recTime', recTime, ...
        'nBreath', numel(pk), 'meanFreq', mean(freq), 'stdFreq', std(freq), ...
        'medFreq', median(freq), 'iqrFreq', iqr(freq), ...
        'durSec', size(raw,1)/fps, 'nLowLk', nLow); %#ok<SAGROW>

    if mod(i, 20) == 0 || i == nRun
        fprintf('  processed %d/%d runs\n', i, nRun);
    end

    % --- optional QC trace for one run ---
    if runNum(i) == EXAMPLE_RUN
        tt = (0:numel(bz)-1)'/fps;
        figQC = figure('Name', sprintf('QC run%03d', runNum(i)), 'Color','w', ...
                       'Position',[120 120 1100 380]);
        plot(tt, bz, 'k'); hold on;
        plot(tt(pk), bz(pk), 'ro', 'MarkerFaceColor','r');
        xlabel('Time (s)'); ylabel('Breath (z)');
        title(sprintf('run%03d  fps=%.2f  %d peaks  mean 1/ITI = %.2f Hz', ...
              runNum(i), fps, numel(pk), mean(freq)));
        box off;
        if SAVE_OUT
            exportgraphics(figQC, fullfile(folderPath, ...
                sprintf('breathfreq_QC_run%03d.png', runNum(i))), 'Resolution',150);
        end
    end
end

nGood = numel(R);
fprintf('Summarised %d/%d runs with valid breathing stats.\n', nGood, nRun);
assert(nGood > 0, 'No runs produced valid stats -- check params.');

% ---- assemble vectors ----
rn   = [R.run]';
mF   = [R.meanFreq]';
sF   = [R.stdFreq]';
mdF  = [R.medFreq]';
recT = [R.recTime]';
nBr  = [R.nBreath]';

switch TIME_AXIS
    case 'clock'
        x  = hours(recT - min(recT));
        xl = 'Recording time (h from first run)';
    otherwise
        x  = rn;
        xl = 'Run number';
end

% ====================== MAIN FIGURE ======================
fig = figure('Name','Breathing frequency across sessions','Color','w', ...
             'Position',[100 100 1150 640]);

% (1) mean +- std per session vs time
ax1 = subplot(2,1,1);
errorbar(x, mF, sF, 'o-', 'Color',[0.10 0.45 0.80], ...
         'MarkerFaceColor',[0.10 0.45 0.80], 'CapSize',3, 'LineWidth',1.0);
hold on;
plot(x, mdF, '.', 'Color',[0.85 0.33 0.10], 'MarkerSize',8);
ylabel('Breath rate  1/ITI  (Hz)');
title(sprintf('%s  |  %s-%s  |  %d runs  |  mean\\pmstd per session (median = orange)', ...
      'breathing freq across sessions', BODYPART, AXIS, nGood), 'Interpreter','tex');
legend({'mean \pm std','median'}, 'Location','best'); legend boxoff;
grid on; box off;

% (2) breaths detected per session (coverage / QC)
ax2 = subplot(2,1,2);
bar(x, nBr, 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none');
ylabel('# breaths / run'); xlabel(xl);
grid on; box off;
linkaxes([ax1 ax2], 'x');

if SAVE_OUT
    stem = fullfile(folderPath, 'breathing_freq_across_sessions');
    exportgraphics(fig, [stem '.png'], 'Resolution',200);
    exportgraphics(fig, [stem '.pdf'], 'ContentType','vector', 'BackgroundColor','none');

    % tidy results table + mat
    Tsum = table(rn, recT, [R.fps]', nBr, mF, sF, mdF, [R.iqrFreq]', [R.durSec]', [R.nLowLk]', ...
        'VariableNames', {'run','recTime','fps','nBreath','meanFreq_Hz','stdFreq_Hz', ...
                          'medFreq_Hz','iqrFreq_Hz','durSec','nLowLk'});
    writetable(Tsum, [stem '.csv']);
    params = struct('BODYPART',BODYPART,'AXIS',AXIS,'LK_THR',LK_THR,'N_DROP',N_DROP, ...
        'HP_WIN_SEC',HP_WIN_SEC,'MINPROM_FACTOR',MINPROM_FACTOR,'MINDIST_SEC',MINDIST_SEC, ...
        'MINHEIGHT',MINHEIGHT,'MINWIDTH_SEC',MINWIDTH_SEC,'BREATH_BAND_HZ',BREATH_BAND_HZ);
    save([stem '.mat'], 'R', 'Tsum', 'params', 'folderPath');
    fprintf('Saved: %s .png/.pdf/.csv/.mat\n', stem);
end

fprintf('\nGrand mean breath rate across runs: %.2f Hz (range %.2f-%.2f)\n', ...
        mean(mF), min(mF), max(mF));
