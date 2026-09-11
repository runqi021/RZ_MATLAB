% C:\Users\Admin\Desktop\RZ_MATLAB\whisker_DLC_hilbert_RZ.m

close all; clc; clear;

csvPath = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos\5916296_whisk_n4DLC_Resnet50_260615_whiskJun17shuffle1_snapshot_best-140.csv";

%
% "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos\5916296_whisk_n2DLC_Resnet50_260615_whiskJun17shuffle1_snapshot_best-140.csv"
%
% "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos\5916297_whisk_n3DLC_Resnet50_260615_whiskJun17shuffle1_snapshot_best-140.csv"
% "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos\5916297_whisk_n4DLC_Resnet50_260615_whiskJun17shuffle1_snapshot_best-140.csv"
%% full
[folder, stem] = fileparts(csvPath);

M = dlc_gate_interp(csvPath, 0.6);     % cols: frame, vL1(x,y,lik), vL2, vR1, vR2  (1=base, 2=tip)    % lik<0.6 -> linear interp
fps = 400;
ts = fullfile(folder, 'timestamps.csv');
if isfile(ts), A = readmatrix(ts); fps = (size(A,1)-1)/((A(end,2)-A(1,2))/1e9); end
t = (0:size(M,1)-1)'/fps;

L0 = [M(:,2) -M(:,3)];  pL0 = M(:,4);     % vL1 = LEFT base
L1 = [M(:,5) -M(:,6)];  pL1 = M(:,7);     % vL2 = LEFT tip
R0 = [M(:,8) -M(:,9)]; pR0 = M(:,10);     % vR1 = RIGHT base
R1 = [M(:,11) -M(:,12)]; pR1 = M(:,13);   % vR2 = RIGHT tip
% --- sweep signal: base->tip ANGLE (atan2), now that tracking is consistently
% --- good. 'proj' (tip-only projection) kept as a fallback.
swMode = 'angle';                 % 'angle' (atan2 base->tip, deg) | 'proj' (tip-only, px)
switch swMode
    case 'angle'
        swL = whiskangle(L0, L1, true);  swR = whiskangle(R0, R1, false);  swUnit = 'deg';  % LEFT mirrored -> protraction + on both
    case 'proj'
        swL = projsweep(L1);       swR = projsweep(R1);        swUnit = 'px';
end

% --- baseline removal: subtract a 1 Hz low-pass (remove slow postural drift),
% --- same trick as breathing. swL/swR become the detrended whisking signal;
% --- swL_raw/swR_raw keep the raw angle if you want it.
% LP_CUT  = 1;                                  % Hz, baseline cutoff
% [lb,la] = butter(2, LP_CUT/(fps/2), 'low');
% swL_raw = swL;  swR_raw = swR;
% swL = swL - filtfilt(lb, la, fillmissing(swL, 'linear'));
% swR = swR - filtfilt(lb, la, fillmissing(swR, 'linear'));

%%
f = figure('Color','w','Position',[80 80 1100 600]);
subplot(2,1,1);
yyaxis left;  plot(t, swL, '-'); ylabel(sprintf('left sweep (%s)',swUnit));
yyaxis right; plot(t, pL1, 'k-'); ylabel('confidence'); ylim([0 1]);
title('LEFT tip (vL1)'); xlabel('s'); grid on;

subplot(2,1,2);
yyaxis left;  plot(t, swR, '-'); ylabel(sprintf('right sweep (%s)',swUnit));
yyaxis right; plot(t, pR1, 'k-'); ylabel('confidence'); ylim([0 1]);
title('RIGHT tip (vR1)'); xlabel('s'); grid on;

%%
figure;
plot(L0(:,1),L0(:,2), "."); hold on;
plot(L1(:,1),L1(:,2),"."); hold on;
plot(R0(:,1),R0(:,2),"."); hold on;
plot(R1(:,1),R1(:,2), "."); hold on;

figure;
plot(t, swL); hold on; plot(t, swR);

figure;
plot(swR, swL, '.');

%% xlim
t_epoc = [25 50];

f = figure('Color','w','Position',[80 80 1100 600]);

yyaxis left;  plot(t, swL, '-'); ylabel(sprintf('left sweep (%s)',swUnit));
yyaxis right;  plot(t, swR, '-'); ylabel(sprintf('right sweep (%s)',swUnit));
xlim(t_epoc);
title('L/R overlay'); xlabel('s'); grid on;

%% ===================== PSD : pick the whisking band (RUN THIS FIRST) =====================
% Read off where the whisking peak sits, then set lo/hi in the Hilbert
% section below. The band MUST bracket this peak or the Hilbert phase will
% stall/reverse (non-monotonic).
psd_view = [0 60];                 % PSD x-range to display (Hz)
psd_fmin = 1;                      % ignore < this for peak-find (drift)
psd_epoc = t_epoc;                 % PSD only on this window (s); = current t_epoc

ep  = t >= psd_epoc(1) & t <= psd_epoc(2);   % epoch sample mask
swL0 = fillmissing(swL(ep) - mean(swL(ep),'omitnan'), 'linear');
swR0 = fillmissing(swR(ep) - mean(swR(ep),'omitnan'), 'linear');

win  = min(round(2*fps), numel(swL0));   % ~2 s Welch segments (epoch is short)
nov  = round(win/2);
nfft = 2^nextpow2(win);
[PxL,Fp] = pwelch(swL0, hann(win), nov, nfft, fps);
[PxR,~ ] = pwelch(swR0, hann(win), nov, nfft, fps);

bm = Fp >= psd_fmin & Fp <= psd_view(2);
[~,iL] = max(PxL.*bm);  fpkL = Fp(iL);
[~,iR] = max(PxR.*bm);  fpkR = Fp(iR);

f = figure('Color','w','Position',[80 80 900 500]);
plot(Fp, PxL, '-', 'LineWidth',1.2); hold on;
plot(Fp, PxR, '-', 'LineWidth',1.2);
xline(fpkL,'--'); xline(fpkR,'--');
xlim(psd_view); xlabel('Hz'); ylabel(sprintf('PSD (%s^2/Hz)',swUnit)); grid on;
legend({'left','right', sprintf('L peak %.1f Hz',fpkL), ...
                         sprintf('R peak %.1f Hz',fpkR)});
title(sprintf('whisking PSD  %.0f-%.0f s — set bandpass lo/hi from the peak', psd_epoc(1), psd_epoc(2)));

%% ===================== Hilbert : instantaneous phase + envelope =====================
% EDIT lo/hi from the PSD peak above. Keep the band wide enough to preserve
% the whisk waveform shape, narrow enough to stay single-component.
lo = 3;   hi = 40;          % whisking band (Hz) -- EDIT from PSD
amp_thresh_frac = 0.20;     % mask phase where envelope < frac * 90th-pct envelope

[bb,aa] = butter(2, [lo hi]/(fps/2), 'bandpass');

% pin sign so protraction = positive on BOTH sides (aligns L/R despite the
% mirror-image pad geometry; in 'proj' mode also fixes the arbitrary eig sign).
% Needed for the in-phase vs anti-phase offset to be interpretable.
swLs = signpin(swL);  swRs = signpin(swR);

xL = filtfilt(bb,aa, fillmissing(swLs - mean(swLs,'omitnan'), 'linear'));  % zero-phase
xR = filtfilt(bb,aa, fillmissing(swRs - mean(swRs,'omitnan'), 'linear'));

zL = hilbert(xL);  phiL = angle(zL);  ampL = abs(zL);   % phase (-pi..pi), envelope
zR = hilbert(xR);  phiR = angle(zR);  ampR = abs(zR);

% mask phase to whisking bouts (phase is meaningless when not whisking)
thL = amp_thresh_frac * prctile(ampL,90);
thR = amp_thresh_frac * prctile(ampR,90);
phiLm = phiL;  phiLm(ampL < thL) = NaN;
phiRm = phiR;  phiRm(ampR < thR) = NaN;

% SANITY: instantaneous frequency should hover near the whisking peak and
% stay positive. Lots of negatives => band is wrong (usually lo too low).
ifL = diff(unwrap(phiL))*fps/(2*pi);
ifR = diff(unwrap(phiR))*fps/(2*pi);
fprintf('inst. freq median  L=%.1f Hz  R=%.1f Hz  (neg frac L=%.2f R=%.2f)\n', ...
    median(ifL), median(ifR), mean(ifL<0), mean(ifR<0));

f = figure('Color','w','Position',[80 80 1100 720]);
ax1 = subplot(3,1,1);
plot(t, xL); hold on; plot(t, ampL,'-'); plot(t,-ampL,'-');
ylabel(sprintf('left (%s)',swUnit)); grid on;
title(sprintf('LEFT  bandpass %.0f-%.0f Hz + envelope', lo, hi));
ax2 = subplot(3,1,2);
plot(t, xR); hold on; plot(t, ampR,'-'); plot(t,-ampR,'-');
ylabel(sprintf('right (%s)',swUnit)); grid on; title('RIGHT  bandpass + envelope');
ax3 = subplot(3,1,3);
plot(t, phiLm,'-.','MarkerSize',4); hold on; plot(t, phiRm,'-.','MarkerSize',4);
ylabel('phase (rad)'); xlabel('s'); ylim([-pi pi]); grid on;
yticks([-pi 0 pi]); yticklabels({'-\pi','0','\pi'});
legend(ax3,{'left','right'}); title(ax3,'instantaneous phase (masked to whisk bouts)');
linkaxes([ax1 ax2 ax3],'x'); xlim(ax3, t_epoc);

%% ===================== SYNCHRONIZATION : L vs R whisking =====================
% Phase-based sync is the right tool for two oscillators (correlation conflates
% amplitude with timing and collapses in low-amplitude periods).
% ALL phase metrics use whisking bouts only (both envelopes above threshold).
% NOTE on sign: projsweep's axis sign is arbitrary, so absolute in-phase vs
% anti-phase is only meaningful because signpin() pins protraction = positive
% on BOTH sides. PLV / coherence MAGNITUDE are robust regardless; the OFFSET
% is interpretable only under that protraction-positive convention.
sync_epoc = t_epoc;      % restrict to this window ([] = whole recording)
cohTW     = 4;           % Chronux time-bandwidth product
cohFpass  = [0 40];      % coherence frequency range (Hz)

inEp = true(size(t));
if ~isempty(sync_epoc), inEp = t>=sync_epoc(1) & t<=sync_epoc(2); end
bout = (ampL >= thL) & (ampR >= thR) & inEp;      % whisking-bout mask

% ---- (1) Phase Locking Value + mean phase offset (bout-masked) ----
dphi    = angle(exp(1i*(phiL - phiR)));           % wrapped L-R phase diff (-pi..pi)
PLV     = abs(mean(exp(1i*dphi(bout))));          % 0 = none, 1 = perfect locking
mu_dphi = angle(mean(exp(1i*dphi(bout))));        % mean L-R offset (rad)
fprintf('PLV = %.3f   mean L-R offset = %+.0f deg   (%d bout samples, %.0f%% of epoch)\n', ...
    PLV, rad2deg(mu_dphi), nnz(bout), 100*nnz(bout)/max(nnz(inEp),1));

% ---- (2) cross-correlation of bandpassed traces (timing/lag) ----
maxlag = round(0.5*fps);                          % +/- 0.5 s
a = xL(inEp) - mean(xL(inEp));
b = xR(inEp) - mean(xR(inEp));
[cc,lags]  = xcorr(a, b, maxlag, 'coeff');
[ccpk,ipk] = max(cc); lag_ms = 1000*lags(ipk)/fps;
fprintf('xcorr peak r = %.3f at lag = %+.1f ms (L relative to R)\n', ccpk, lag_ms);

% ---- (3) Chronux coherence (frequency-resolved, jackknife CI) ----
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'chronux_2_12')));
cp.Fs = fps; cp.tapers = [cohTW 2*cohTW-1]; cp.pad = 1; cp.fpass = cohFpass;
cp.err = [2 0.05]; cp.trialave = 0;
dA = swLs(inEp) - mean(swLs(inEp));
dB = swRs(inEp) - mean(swRs(inEp));
% NOTE: this Chronux copy prepends C12 (complex coherency) to the outputs
[~, Coh, cohPhi, ~, ~, ~, fC, confC, ~, Cerr] = coherencyc(dA, dB, cp);

% ---- figure ----
fig = figure('Color','w','Position',[80 80 1150 780]);

pax = subplot(2,2,1,polaraxes);
polarhistogram(pax, dphi(bout), 36, 'Normalization','probability', ...
    'FaceColor',[0.30 0.50 0.90], 'EdgeColor','none'); hold(pax,'on');
rmax = max(histcounts(dphi(bout), 36, 'Normalization','probability'));
polarplot(pax, [0 mu_dphi], [0 rmax], 'r-', 'LineWidth', 2);   % mean offset
title(pax, sprintf('L-R phase diff   PLV=%.2f, offset=%+.0f\\circ', PLV, rad2deg(mu_dphi)));

subplot(2,2,2);
plot(1000*lags/fps, cc, '-', 'LineWidth',1.1); hold on; xline(lag_ms,'r--');
xlabel('lag (ms), L rel. to R'); ylabel('xcorr (coeff)'); grid on; ylim([-1 1]);
title(sprintf('cross-corr   peak r=%.2f @ %+.0f ms', ccpk, lag_ms));

subplot(2,2,3);
plot(fC, Coh, '-', 'LineWidth',1.3); hold on;
if ~isempty(Cerr), plot(fC, Cerr(1,:), 'b--'); plot(fC, Cerr(2,:), 'b--'); end
yline(confC, 'k:', 'LineWidth',1);
xlim(cohFpass); ylim([0 1]); xlabel('Hz'); ylabel('coherence'); grid on;
title('L-R coherence (jackknife 95% CI; dotted = sig level)');

subplot(2,2,4);
sig = Coh > confC; phs = cohPhi; phs(~sig) = NaN;
plot(fC, rad2deg(phs), '.', 'MarkerSize',8);
xlim(cohFpass); ylim([-180 180]); yticks(-180:90:180);
xlabel('Hz'); ylabel('coh phase (deg)'); grid on;
title('coherence phase (only where C > sig)');

% %%
% plot(t, swL); hold on; plot(t, xL); xlim(t_epoc);
% %%
% plot(t, swR); hold on; plot(t, xR); xlim(t_epoc);
% %%
% plot(t, xL); hold on; plot(t, xR); xlim(t_epoc);
% %%
% plot(t, phiLm,'-','MarkerSize',4); hold on; plot(t, phiRm,'-','MarkerSize',4); xlim(t_epoc);

%%
t_epoc = [33 38];

figure;

yyaxis left;  plot(t, xL, '-'); ylabel(sprintf('left sweep (%s)',swUnit));
yyaxis right;  plot(t, xR, '-'); ylabel(sprintf('right sweep (%s)',swUnit));
xlim(t_epoc);
title('L/R overlay'); xlabel('s'); grid on;

%%
% on bandpassed traces, over the epoch (cleanest)
r_bp  = corr(xL(inEp), xR(inEp));
% on the raw demeaned sweeps, over the epoch
r_raw = corr(swL(inEp), swR(inEp));

%%
mdl = fitlm(swL, swR);

figure;
plot(mdl);
xlabel('xL');
ylabel('xR');

R2 = mdl.Rsquared.Ordinary;
p  = mdl.Coefficients.pValue(2);   % p-value for slope

title(sprintf('xR vs xL: R^2 = %.3f, p = %.3g', R2, p));

%%
mdl = fitlm(xL, xR);

figure;
plot(mdl);
xlabel('xL');
ylabel('xR');

R2 = mdl.Rsquared.Ordinary;
p  = mdl.Coefficients.pValue(2);   % p-value for slope

title(sprintf('xR vs xL: R^2 = %.3f, p = %.3g', R2, p));
%%
mdl = fitlm(phiL, phiR);

figure;
plot(mdl);
xlabel('xL');
ylabel('xR');

R2 = mdl.Rsquared.Ordinary;
p  = mdl.Coefficients.pValue(2);   % p-value for slope

title(sprintf('xR vs xL: R^2 = %.3f, p = %.3g', R2, p));
%%
function ang = whiskangle(base, tip, mirrorx)
% Whisker pointing angle (deg) = direction of the base->tip vector, UNWRAPPED
% (no +-180 seam jumps), not demeaned. With mirrorx=true the x-axis is reflected
% -- use it for the LEFT pad so that PROTRACTION IS POSITIVE ON BOTH SIDES (L and
% R read the same way instead of mirror-flipped) and the left stays off the
% +-180 seam (it points ~-x). Right pad: mirrorx=false.
if nargin < 3, mirrorx = false; end
d  = tip - base;                                    % [T x 2] base->tip vector
dx = d(:,1);  if mirrorx, dx = -dx; end             % reflect x for the left pad
ang = rad2deg(unwrap(atan2(d(:,2), dx)));           % angle (deg)
end

function s = projsweep(xy)
C = xy - mean(xy,1); [V,~] = eig(cov(C)); s = C*V(:,end);   % project on sweep axis
end

function s = signpin(s)
% pin sign so the sharper (protraction) excursions are positive -> reproducible phase
v = s(~isnan(s));
if skewness(v) < 0, s = -s; end
end
