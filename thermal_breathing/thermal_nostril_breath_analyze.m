% thermal_nostril_breath_analyze.mc
close all; clc; clear;

%%
% Analysis of the drawn-ROI nostril traces from thermal_nostril_roi_draw.m
% (<stem>_nostrilROI.mat). Bandpass + invert (inhale up), Hilbert instantaneous
% phase (L vs R), and the L+R averaged signal.
%
%   FIG 1  RAW deg-C ROI traces (L, R) + PSD
%   FIG 2  BANDPASSED + inverted traces (L, R) + PSD
%   FIG 3  Hilbert instantaneous phase, L vs R overlaid + phase difference + PLV
%   FIG 4  averaged (L+R)/2 ROI, raw + bandpassed + PSD
% Saves ONLY the final L+R averaged breathing signals (wide + finer band) to
% <stem>_breath.mat -- for later joint analysis with whisking. Nothing else saved.

% Give ONLY the DLC csv; the _nostrilROI.mat is resolved via thermal_resolve_paths.
dlcCsv    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n4DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv";
dataRoot  = "D:\260615_thermalNbasler";
BP        = [2 10];   % bandpass (Hz)
INVERT    = true;     % inhale cools the nostril -> invert so inhale reads as a rise
PHASE_WIN = 10;       % seconds shown initially on the phase figure (zoom for more)

addpath(fileparts(mfilename('fullpath')));
P = thermal_resolve_paths(dlcCsv, dataRoot);
roiPath = P.nostrilROI;
D = load(roiPath);
sides = {'L','R'};
fps = double(D.L.fps);
[bb,aa] = butter(2, BP/(fps/2), 'bandpass');

% per-side bandpass + invert
for i = 1:2
    sd = sides{i};
    tr = double(D.(sd).trace);
    trbp = filtfilt(bb, aa, tr - mean(tr,'omitnan'));
    if INVERT, trbp = -trbp; end
    D.(sd).trace_bp = trbp;
    D.(sd).phase = angle(hilbert(trbp));
end
D.bp = BP; D.inverted = INVERT;
invstr = ''; if INVERT, invstr = ' (inverted: inhale up)'; end

% ---- FIG 1: RAW ----
figure('Color','w','Position',[60 80 1150 600]);
for i = 1:2
    sd = sides{i}; tr = double(D.(sd).trace); t = (0:numel(tr)-1)/fps;
    subplot(2,2,i); plot(t, tr,'-'); grid on; xlabel('s'); ylabel('\circC');
    title(sprintf('%s RAW ROI temp (mean %.2f \\circC)', D.(sd).side, mean(tr,'omitnan')));
    subplot(2,2,i+2); psd_panel(tr - mean(tr,'omitnan'), fps, BP, sprintf('%s raw PSD', D.(sd).side));
end
sgtitle('Nostril ROI breathing — RAW (\circC)');

% ---- FIG 1b: WIDE bandpass 0.1-30 Hz (drift removed, near-raw) ----
WBP = [1 50];
[wb,wa] = butter(2, WBP/(fps/2), 'bandpass');
figure('Color','w','Position',[100 95 1150 600]);
for i = 1:2
    sd = sides{i}; tr = double(D.(sd).trace);
    trw = filtfilt(wb, wa, tr - mean(tr,'omitnan'));
    if INVERT, trw = -trw; end
    D.(sd).trace_wide = trw;
    t = (0:numel(trw)-1)/fps;
    subplot(2,2,i); plot(t, trw,'-'); grid on; xlabel('s'); ylabel('\circC (BP)');
    title(sprintf('%s %g-%g Hz', D.(sd).side, WBP(1), WBP(2)));
    subplot(2,2,i+2); psd_panel(trw, fps, WBP, sprintf('%s PSD', D.(sd).side));
end
sgtitle(sprintf('Nostril ROI — WIDE BANDPASS %g-%g Hz%s', WBP(1), WBP(2), invstr));

% ---- FIG 2: BANDPASSED ----
figure('Color','w','Position',[140 110 1150 600]);
for i = 1:2
    sd = sides{i}; trbp = D.(sd).trace_bp; t = (0:numel(trbp)-1)/fps;
    subplot(2,2,i); plot(t, trbp,'-'); grid on; xlabel('s'); ylabel('\circC (BP, inhale up)');
    title(sprintf('%s breathing %g-%g Hz', D.(sd).side, BP(1), BP(2)));
    subplot(2,2,i+2); psd_panel(trbp, fps, BP, sprintf('%s PSD', D.(sd).side));
end
sgtitle(sprintf('Nostril ROI breathing — BANDPASS %g-%g Hz%s', BP(1), BP(2), invstr));

% ---- FIG 3: Hilbert instantaneous phase, L vs R ----
phiL = D.L.phase; phiR = D.R.phase;
dphi = angle(exp(1i*(phiL - phiR))); plv = abs(mean(exp(1i*dphi)));
tp = (0:numel(phiL)-1)/fps;
figure('Color','w','Position',[200 140 1150 640]);
ax1 = subplot(3,1,1);
plot(ax1, tp, D.L.trace_bp,'-'); hold(ax1,'on'); plot(ax1, tp, D.R.trace_bp,'-'); grid(ax1,'on');
legend(ax1,{'L','R'}); ylabel(ax1,'\circC (BP)'); title(ax1,'bandpassed breathing (L vs R)');
ax2 = subplot(3,1,2);
plot(ax2, tp, phiL,'-.'); hold(ax2,'on'); plot(ax2, tp, phiR,'-.'); grid(ax2,'on');
ylim(ax2,[-pi pi]); yticks(ax2,[-pi 0 pi]); yticklabels(ax2,{'-\pi','0','\pi'});
ylabel(ax2,'inst. phase'); legend(ax2,{'L','R'});
title(ax2, sprintf('Hilbert instantaneous phase (L vs R)   PLV=%.2f', plv));
ax3 = subplot(3,1,3);
plot(ax3, tp, dphi,'-'); grid(ax3,'on');
ylim(ax3,[-pi pi]); yticks(ax3,[-pi 0 pi]); yticklabels(ax3,{'-\pi','0','\pi'});
xlabel(ax3,'s'); ylabel(ax3,'L - R phase');
title(ax3, sprintf('phase difference (mean %+.0f\\circ)', rad2deg(angle(mean(exp(1i*dphi))))));
linkaxes([ax1 ax2 ax3],'x'); xlim(ax3,[0 min(PHASE_WIN, tp(end))]);

% ---- FIG 4: averaged (L+R)/2 of the WIDE (0.5-30 Hz) signal + its finer BP ----
trA_wide = (D.L.trace_wide + D.R.trace_wide) / 2;   % average the wide (0.5-30) traces
trA_fine = filtfilt(bb, aa, trA_wide);              % then finer BP [BP] of that average
D.avg.trace_wide = trA_wide; D.avg.trace_bp = trA_fine; D.avg.fps = fps;
ta = (0:numel(trA_wide)-1)/fps;
figure('Color','w','Position',[260 160 1150 600]);
subplot(2,2,1); plot(ta, trA_wide,'-'); grid on; xlabel('s');
ylabel(sprintf('\\circC (%g-%g Hz)', WBP(1), WBP(2)));
title(sprintf('L+R avg WIDE %g-%g Hz', WBP(1), WBP(2)));
subplot(2,2,2); psd_panel(trA_wide, fps, WBP, 'avg wide PSD');
subplot(2,2,3); plot(ta, trA_fine,'-'); grid on; xlabel('s');
ylabel(sprintf('\\circC (%g-%g Hz)', BP(1), BP(2)));
title(sprintf('L+R avg finer BP %g-%g Hz', BP(1), BP(2)));
subplot(2,2,4); psd_panel(trA_fine, fps, BP, 'avg finer PSD');
sgtitle('Averaged nostril ROI (L+R)/2 — wide + finer');

% ---- FIG 5: L vs R linear regression (amplitude and phase), as in whisking ----
figure('Color','w','Position',[320 180 1050 480]);
subplot(1,2,1);
mdlA = fitlm(D.L.trace_bp, D.R.trace_bp);
plot(mdlA); legend(gca,'off');
xlabel('L (\circC, BP)'); ylabel('R (\circC, BP)');
title(sprintf('L vs R amplitude:  R^2=%.3f, p=%.3g', ...
    mdlA.Rsquared.Ordinary, mdlA.Coefficients.pValue(2)));
subplot(1,2,2);
mdlP = fitlm(phiL, phiR);                 % NOTE: phase is circular; linear fit as in whisking
plot(mdlP); legend(gca,'off');
xlabel('\phi_L (rad)'); ylabel('\phi_R (rad)');
title(sprintf('\\phi_L vs \\phi_R:  R^2=%.3f, p=%.3g', ...
    mdlP.Rsquared.Ordinary, mdlP.Coefficients.pValue(2)));
sgtitle('L vs R linear regression');
% ---- SAVE ONLY the final averaged breathing signals (for later joint analysis
%      with whisking). We do NOT re-save the full per-side/mask/phase struct. ----
B = struct();
B.wbp      = trA_wide(:);     % L+R averaged, WIDE band, inhale-up (deg C)
B.bp       = trA_fine(:);     % L+R averaged, finer band, inhale-up (deg C)
B.fps      = fps;
B.t        = ta(:);          % time vector (s)
B.wbp_band = WBP;
B.bp_band  = BP;
B.inverted = INVERT;
B.animal   = P.animal;
B.run      = P.k;
B.src_csv  = char(dlcCsv);
B.src_ats  = P.ats;
% NOTE: lpsub script owns the canonical _breath.mat; this bandpass version saves separately.
outPath = strrep(P.breath, '_breath.mat', '_breath_bpavg.mat');
save(outPath, '-struct', 'B');
fprintf('saved %s (bandpass-avg; canonical breathing = _breath.mat from lpsub)\n', outPath);

% ---------- local function ----------
function psd_panel(x, fps, BP, ttl)
    % averaged-FFT periodogram: |fft|^2 of detrended Hann-windowed overlapping
    % segments, averaged; band = +/-1 STD computed in LOG space (log-y).
    x = fillmissing(x(:),'linear');
    seg = max(min(round(4*fps), floor(numel(x)/8)), 16); nov = round(seg/2);
    nfft = 2^nextpow2(seg); w = hann(seg); step = max(seg-nov,1);
    st = 1:step:(numel(x)-seg+1); P = zeros(numel(st), nfft);
    for k = 1:numel(st)
        s = detrend(x(st(k):st(k)+seg-1)) .* w;
        P(k,:) = abs(fft(s,nfft)).^2;
    end
    fr = (0:nfft-1)*(fps/nfft); keep = fr <= min(40, fps/2); fr = fr(keep);
    Plog = log10(P(:,keep) + eps);
    ml = mean(Plog,1); sl = (size(Plog,1)>=2)*std(Plog,0,1);
    patch([fr fliplr(fr)], [10.^(ml+sl) fliplr(10.^(ml-sl))], [0.3 0.5 0.9], ...
        'FaceAlpha',0.2,'EdgeColor','none'); hold on;
    plot(fr, 10.^ml, '-','LineWidth',1.3); set(gca,'YScale','log'); grid on;
    xline(BP(1),'r--'); xline(BP(2),'r--');
    inb = fr>=BP(1) & fr<=BP(2); [~,ip] = max(ml(inb)); frb = fr(inb);
    xlabel('Hz'); ylabel('power'); title(sprintf('%s (peak %.2f Hz, %d segs)', ttl, frb(ip), numel(st)));
end
