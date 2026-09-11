% thermal_nostril_breath_lpsub.m
% Baseline removal by LOW-PASS SUBTRACTION (alternative to bandpass).
%   baseline = lowpass(raw, LP_CUT)        % the slow drift
%   detrended = raw - baseline             % breathing with baseline shift removed
% Reads the drawn-ROI traces (<stem>_nostrilROI.mat) resolved from the DLC csv.
%
%   FIG 1  per nostril: raw + LP baseline (top), detrended (bottom)
%   FIG 2  L+R averaged detrended (+ optional finer BP) + PSD
% Writes THE breathing signal (L+R averaged, LP-subtracted) to <stem>_breath.mat
% -- this is the canonical breathing file for later joint analysis with whisking.
%
% Zero-phase throughout (butter + filtfilt).
close all; clc; clear;

% Give ONLY the DLC csv; the _nostrilROI.mat is resolved via thermal_resolve_paths.
dlcCsv   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n4DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv";
dataRoot = "D:\260615_thermalNbasler";
LP_CUT   = 1;        % Hz, low-pass cutoff = the slow baseline to subtract off
INVERT   = true;       % inhale cools the nostril -> invert so inhale reads as a rise
FINE_BP  = [1 30];     % optional finer band of the detrended signal ([] = skip)

addpath(fileparts(mfilename('fullpath')));
P = thermal_resolve_paths(dlcCsv, dataRoot);
D = load(P.nostrilROI);
sides = {'L','R'};
fps = double(D.L.fps);
[lb,la] = butter(2, LP_CUT/(fps/2), 'low');     % low-pass for the baseline
invstr = ''; if INVERT, invstr = ' (inverted: inhale up)'; end

% per-side: baseline = LP(raw); detrended = raw - baseline
for i = 1:2
    sd  = sides{i};
    raw = double(D.(sd).trace);
    base = filtfilt(lb, la, raw);               % zero-phase low-pass
    detr = raw - base;
    if INVERT, detr = -detr; end
    D.(sd).baseline = base;
    D.(sd).detr = detr;
end

% ---- FIG 1: raw + baseline (top), detrended (bottom), per nostril ----
figure('Color','w','Position',[80 90 1150 620]);
for i = 1:2
    sd = sides{i}; raw = double(D.(sd).trace); base = D.(sd).baseline; detr = D.(sd).detr;
    t = (0:numel(raw)-1)/fps;
    subplot(2,2,i);
    plot(t, raw, '-'); hold on; plot(t, base, '-', 'LineWidth', 1.6); grid on;
    legend({'raw','LP baseline'}); ylabel('\circC');
    title(sprintf('%s raw + LP baseline (%g Hz)', D.(sd).side, LP_CUT));
    subplot(2,2,i+2);
    plot(t, detr, '-'); grid on; xlabel('s'); ylabel('\circC');
    title(sprintf('%s raw - baseline (detrended)', D.(sd).side));
end
sgtitle(sprintf('Baseline removal by LP-subtraction (cutoff %g Hz)%s', LP_CUT, invstr));

% ---- FIG 2: L+R averaged detrended (+ optional finer BP) + PSD ----
detrA = (D.L.detr + D.R.detr) / 2;
ta = (0:numel(detrA)-1)/fps;
figure('Color','w','Position',[140 130 1150 600]);
subplot(2,2,1); plot(ta, detrA, '-'); grid on; xlabel('s'); ylabel('\circC');
title('L+R avg detrended (raw - LP)');
subplot(2,2,2); psd_panel(detrA, fps, [LP_CUT 15], 'avg detrended PSD');

if ~isempty(FINE_BP)
    [fbb,faa] = butter(2, FINE_BP/(fps/2), 'bandpass');
    detrA_bp = filtfilt(fbb, faa, detrA);
    subplot(2,2,3); plot(ta, detrA_bp, '-'); grid on; xlabel('s'); ylabel('\circC');
    title(sprintf('L+R avg detrended, finer BP %g-%g Hz', FINE_BP(1), FINE_BP(2)));
    subplot(2,2,4); psd_panel(detrA_bp, fps, FINE_BP, 'avg finer PSD');
else
    detrA_bp = [];
end
sgtitle('Averaged nostril detrended (LP-subtraction)');

% ---- write THE breathing signal: L+R averaged, baseline-removed ----
B = struct();
B.breath    = detrA(:);     % PRIMARY breathing signal (L+R avg, raw - LP, inhale-up)
B.breath_bp = detrA_bp(:);  % finer-band version ([] if FINE_BP empty)
B.fps       = fps;
B.t         = ta(:);
B.method    = 'lpsub';
B.lp_cut    = LP_CUT;
B.fine_bp   = FINE_BP;
B.inverted  = INVERT;
B.animal    = P.animal; B.run = P.k;
B.src_csv   = char(dlcCsv); B.src_ats = P.ats;
save(P.breath, '-struct', 'B');
fprintf('saved %s (breathing signal: L+R avg, LP-subtracted)\n', P.breath);

% ---------- local function ----------
function psd_panel(x, fps, BND, ttl)
    nfft = 2^nextpow2(numel(x));
    P = abs(fft(x, nfft)).^2; fr = (0:nfft-1)*(fps/nfft);
    keep = fr <= min(40, fps/2);
    plot(fr(keep), P(keep), '-'); grid on; hold on;
    xline(BND(1),'r--'); xline(BND(2),'r--');
    inb = fr>=BND(1) & fr<=BND(2); [~,ip] = max(P(inb)); frb = fr(inb);
    xlabel('Hz'); ylabel('power'); title(sprintf('%s (peak %.2f Hz)', ttl, frb(ip)));
end
