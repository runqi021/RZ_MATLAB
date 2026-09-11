%% orofacial_crop_analyze.m
% STEP 2 of 2 -- analyze a cropped Basler ROI (extracted by
% orofacial_crop_extract.py) in space and time via PER-PIXEL dF/F.
%
% Per-pixel signal (SIGNAL param), all LINEAR (unlike rectified motion energy,
% which frequency-doubles):
%   'raw' = intensity minus its temporal mean  -> default; robust for BEHAVIOR
%           video (no division, so dark/background pixels stay well-behaved).
%   'dF'  = intensity minus sliding-median baseline (removes slow drift).
%   'dff' = (I-F0)/F0 -> ONLY for calcium ROIs. DIVIDES by baseline, so on a raw
%           video the near-black pixels (F0~0) explode to ~1e15 and swamp every
%           map/PSD/SVD. Do NOT use 'dff' on behavior video.
%
% Outputs five figure groups:
%   1. Spatial maps      : mean frame, per-pixel std, per-pixel peak
%   2. Frequency maps    : per-pixel dominant temporal freq + band-power maps
%   3. Time maps         : global trace, temporal PSD (Chronux), kymographs
%   4. SVD               : singular spectrum, spatial eigenimages, temporal modes
%   5. SVD traces        : top-K PC time series + movement-energy (DLC-free readout)
%
% Standalone. Run after the Python extractor has written crop_for_analysis.mat.

clear; close all; clc;

%% ============================= USER PARAMETERS =============================
% --- INPUTS (the Basler run folder) ---
RUN_DIR    = 'C:\Users\Admin\Desktop\260602_orofacial_ChAT_ChroME\cam1\cam1_20260602_180811_run001';
AVI_PATH   = fullfile(RUN_DIR, 'cam1_20260602_180811_run001.avi');
TS_PATH    = fullfile(RUN_DIR, 'timestamps.csv');          % absent -> header/FPS_OVERRIDE
OUT_MAT    = fullfile(RUN_DIR, 'crop_for_analysis.mat');   % cropped cube cache
PYTHON_EXE = 'C:\Program Files\Python314\python.exe';      % has cv2 + scipy

% --- EXTRACTION (Python reads the crop; runs only if OUT_MAT missing) ---
RE_EXTRACT      = true;  % true = re-run extraction/crop even if cached (new ROI/window)
EXTRACT_WIN     = [];     % time window (s): [] = full recording, else [start end]
TEMPORAL_STRIDE = 1;      % keep every Nth frame (1 = all)
SPATIAL_BIN     = 1;      % block-average bin (auto-raised to fit size budget)
FULL_FRAME      = false;   % true = whole frame (no crop / no ROI popup)
ROI_XYWH        = [];     % used only if FULL_FRAME=false. [] = interactive drag-box
FPS_OVERRIDE    = [];     % force fps (no timestamps); [] or 0 = auto (timestamps/header)

% --- ANALYSIS ---
T_CLIP         = [];      % optional FURTHER time clip [t0 t1] (s), [] = use all
% Signal per pixel:
%   'raw' = intensity minus its temporal mean  (robust for BEHAVIOR video)
%   'dF'  = intensity minus sliding-median baseline (linear, removes drift)
%   'dff' = (I-F0)/F0  -- ONLY for calcium ROIs; DIVIDES by baseline, so it
%           explodes on dark/background pixels of a raw video (do not use here)
SIGNAL         = 'raw';
BaselineWinSec = 15;      % sliding-median baseline window for 'dF'/'dff' (s)
DropFirstSec   = 0;       % drop first N s before dF/F (already clipped in Python)

TW             = 3;       % Chronux time-bandwidth product for temporal PSD
fBreath        = [0.5 5]; % breathing band (Hz) -- power map + PSD shading
fWhisk         = [8 15];  % whisking band (Hz)
fMap_search    = [0.3 15];% per-pixel dominant-frequency map: search range (Hz)

nSVD           = 20;      % SVD components to compute
nShowSVD       = 6;       % spatial eigenimages / temporal modes to display
nCombinePC     = 5;       % # top PCs combined into the movement-energy trace
SVD_MAX_BYTES  = 2e9;     % auto spatial-bin the SVD input to stay under this (double)
%% ===========================================================================

%% Path setup (repo root has the +helper package; chronux for mtspectrumc)
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
chronuxDir = fullfile(repoRoot, 'chronux_2_12');
if isfolder(chronuxDir)
    addpath(genpath(chronuxDir));
end

%% Extract cropped ROI via Python (skip if cached; RE_EXTRACT or delete OUT_MAT to redo)
if RE_EXTRACT && isfile(OUT_MAT)
    delete(OUT_MAT);
    fprintf('RE_EXTRACT: deleted cached crop -> will re-extract (pick a new ROI).\n');
end
if ~isfile(OUT_MAT)
    extractPy = fullfile(repoRoot, 'orofacial_crop_extract.py');
    assert(isfile(extractPy), 'Missing extractor: %s', extractPy);
    if isempty(EXTRACT_WIN)
        tStart = 0;
    else
        tStart = EXTRACT_WIN(1);
    end
    cmd = sprintf(['"%s" "%s" --avi "%s" --ts "%s" --out "%s" ' ...
                   '--t_start %g --temporal_stride %d --spatial_bin %d'], ...
        PYTHON_EXE, extractPy, AVI_PATH, TS_PATH, OUT_MAT, ...
        tStart, TEMPORAL_STRIDE, SPATIAL_BIN);
    if ~isempty(EXTRACT_WIN)
        cmd = [cmd sprintf(' --t_end %g', EXTRACT_WIN(2))];
    end
    if FULL_FRAME
        cmd = [cmd ' --fullframe'];
    elseif ~isempty(ROI_XYWH)
        cmd = [cmd sprintf(' --roi %d,%d,%d,%d', ROI_XYWH)];
    end
    if ~isempty(FPS_OVERRIDE) && FPS_OVERRIDE > 0
        cmd = [cmd sprintf(' --fps %g', FPS_OVERRIDE)];
    end
    fprintf('Extracting crop via Python (drag the ROI box when the window appears)...\n');
    status = system(cmd, '-echo');
    assert(status == 0 && isfile(OUT_MAT), ...
        'Extraction failed (status %d) -- see messages above.', status);
else
    fprintf('Reusing cached %s\n  (delete it to re-extract / pick a new ROI).\n', OUT_MAT);
end

%% Load extracted crop
S = load(OUT_MAT);
mov = S.mov;                         % uint8 [H, W, T]
t_s = double(S.t_s(:));              % [T] seconds (camera timestamps)
fps = double(S.fps_eff);
roi = S.roi_xywh(:)';
fprintf('Loaded %dx%d ROI, %d frames, fps=%.2f, bin=%g\n', ...
    size(mov,2), size(mov,1), size(mov,3), fps, double(S.spatial_bin));

%% Optional further time clip
if ~isempty(T_CLIP)
    keep = t_s >= T_CLIP(1) & t_s <= T_CLIP(2);
    mov  = mov(:,:,keep);
    t_s  = t_s(keep);
    fprintf('Clipped to [%.1f %.1f] s -> %d frames\n', T_CLIP(1), T_CLIP(2), numel(t_s));
end

[H, W, T] = size(mov);
Npix = H * W;

%% Per-pixel signal X [T2 x Npix] (one "ROI" per pixel)
% reshape [H,W,T] -> F [T x Npix] raw intensity, then form the chosen signal.
F = single(reshape(mov, Npix, T))';                   % [T x Npix] raw intensity
switch lower(SIGNAL)
    case 'raw'      % mean-subtracted intensity — no division, robust on dark pixels
        X = F - mean(F, 1);
        T2 = T;  t = t_s;
        sigLabel = 'intensity (mean-sub)';
    case 'df'       % sliding-median baseline SUBTRACTED (linear, no divide)
        win = round(BaselineWinSec * fps);
        if mod(win, 2) == 0, win = win + 1; end
        X = F - movmedian(F, win, 1);
        T2 = T;  t = t_s;
        sigLabel = '\DeltaF (baseline-sub)';
    case 'dff'      % (I-F0)/F0 — ONLY valid for calcium ROIs (divides by baseline)
        out = helper.dFF_RZ(F, 'FPS', fps, 'BaselineWinSec', BaselineWinSec, ...
                            'DropFirstSec', DropFirstSec, 'UseOASIS', false);
        X = out.dFF;  T2 = size(X, 1);  t = t_s(end-T2+1:end);
        sigLabel = 'dF/F';
    otherwise
        error('SIGNAL must be ''raw'', ''dF'', or ''dff''.');
end
fprintf('Signal: %s  (T=%d, N=%d)\n', sigLabel, T2, Npix);
clear F

%% Per-pixel summary maps
sdMap = reshape(std(X, 0, 1),      H, W);              % temporal std of signal
pkMap = reshape(max(abs(X), [], 1), H, W);             % peak |signal|
meanFrame = mean(double(mov), 3);                      % static structure

%% ===== FIGURE 1: spatial maps =====
figure('Name','1: spatial maps','Color','w','Position',[60 80 1200 420]);
subplot(1,3,1); imagesc(meanFrame); axis image off; colormap(gca,gray); colorbar;
title('mean frame (static structure)');
subplot(1,3,2); imagesc(sdMap); axis image off; colormap(gca,parula); colorbar;
title(sprintf('std(%s) per pixel  (where it moves)', sigLabel));
subplot(1,3,3); imagesc(pkMap); axis image off; colormap(gca,parula); colorbar;
title(sprintf('peak |%s| per pixel', sigLabel));
sgtitle(sprintf('Spatial heatmaps of per-pixel %s activity', sigLabel));

%% ===== FIGURE 2: spatial map of TEMPORAL frequency =====
% For each pixel, FFT its time course, then paint onto the FOV: (a) power in
% the breathing band, (b) power in the whisking band, (c) the DOMINANT temporal
% frequency per pixel. Shows WHERE each rhythm lives and at WHAT rate -- this is
% per-pixel temporal frequency, NOT the 2D spatial-frequency spectrum.
nf  = floor(T2/2) + 1;
fAx = (0:nf-1)' * (fps / T2);                          % temporal frequency axis (Hz)
win = hann(T2);
searchMask = fAx >= fMap_search(1) & fAx <= fMap_search(2);
breathMask = fAx >= fBreath(1)     & fAx <= fBreath(2);
whiskMask  = fAx >= fWhisk(1)      & fAx <= fWhisk(2);

domFreq   = zeros(Npix, 1, 'single');
breathPow = zeros(Npix, 1, 'single');
whiskPow  = zeros(Npix, 1, 'single');
chunk = 20000;                                          % pixels per FFT batch (memory)
for c0 = 1:chunk:Npix
    c1 = min(c0 + chunk - 1, Npix);
    Pf = abs(fft(X(:, c0:c1) .* win)).^2;              % [T2 x m] power
    Pf = Pf(1:nf, :);                                  % one-sided
    Ps = Pf; Ps(~searchMask, :) = 0;                   % restrict dominant-freq search
    [~, im] = max(Ps, [], 1);
    domFreq(c0:c1)   = fAx(im);
    breathPow(c0:c1) = sum(Pf(breathMask, :), 1);
    whiskPow(c0:c1)  = sum(Pf(whiskMask, :), 1);
end
breathPowMap = reshape(breathPow, H, W);
whiskPowMap  = reshape(whiskPow,  H, W);
domFreqMap   = reshape(domFreq,   H, W);
actMask = sdMap > prctile(sdMap(:), 75);               % only paint freq where it moves

figure('Name','2: frequency-distribution maps','Color','w','Position',[80 100 1200 760]);
subplot(2,2,1);
imagesc(breathPowMap); axis image off; colormap(gca,parula); colorbar;
title(sprintf('breathing-band power (%.1f-%.1f Hz)', fBreath(1), fBreath(2)));

subplot(2,2,2);
imagesc(whiskPowMap); axis image off; colormap(gca,parula); colorbar;
title(sprintf('%.0f-%.0f Hz band power', fWhisk(1), fWhisk(2)));

subplot(2,2,3);
hImg = imagesc(domFreqMap); axis image off; colormap(gca, jet);
set(hImg, 'AlphaData', actMask);                       % hide inactive background
set(gca, 'Color', [0.15 0.15 0.15]);
caxis(fMap_search); cb = colorbar; cb.Label.String = 'Hz';
title('dominant temporal frequency / pixel');

subplot(2,2,4);
histogram(domFreq(actMask(:)), 40); grid on;
xlabel('dominant frequency (Hz)'); ylabel('# active pixels'); xlim(fMap_search);
title('per-pixel dominant-frequency distribution');
sgtitle('Spatial distribution of temporal frequency across the FOV');

%% ===== FIGURE 3: time maps =====
gTrace = mean(X, 2);                                    % global mean signal over ROI
[fP, P] = tempPSD(gTrace, fps, TW);

% kymographs: collapse one spatial axis
cube  = reshape(X', H, W, T2);                          % [H, W, T2]
kymoX = squeeze(mean(cube, 1))';                        % [T2 x W]  (rows collapsed)
kymoY = squeeze(mean(cube, 2))';                        % [T2 x H]  (cols collapsed)

figure('Name','3: time maps','Color','w','Position',[100 60 1200 820]);
subplot(3,1,1);
plot(t, gTrace, 'LineWidth', 0.8); grid on; xlim([t(1) t(end)]);
xlabel('time (s)'); ylabel(sprintf('mean %s', sigLabel)); title('global movement trace');

subplot(3,2,3);
fhi  = min(25, fps/2);
pmax = max([P(fP <= fhi); eps]) * 1.05;
hb = shadeBand(fBreath, [0.85 0.92 1.0], pmax); hold on;
hw = shadeBand(fWhisk,  [1.0 0.9 0.85], pmax);
plot(fP, P, 'k', 'LineWidth', 1.4); grid on;
xlim([0 fhi]); ylim([0 pmax]); xlabel('frequency (Hz)'); ylabel('power');
title(sprintf('temporal PSD of %s (linear -> true frequencies)', sigLabel));
legend([hb hw], {'breathing','whisking'}, 'Location','northeast');

subplot(3,2,4); axis off;
text(0.0,0.5,sprintf(['Signal: %s (linear), so peaks sit at the\n' ...
    'REAL movement frequencies (no rectification\n' ...
    'frequency-doubling as motion energy would give).\n' ...
    'Breathing band: %g-%g Hz   Whisking band: %g-%g Hz'], ...
    sigLabel, fBreath(1),fBreath(2),fWhisk(1),fWhisk(2)),'FontSize',11);

subplot(3,2,5);
imagesc(t, 1:W, kymoX'); axis xy; colormap(gca,parula); colorbar;
xlabel('time (s)'); ylabel('column (x)'); title('kymograph: rows collapsed');
subplot(3,2,6);
imagesc(t, 1:H, kymoY'); axis xy; colormap(gca,parula); colorbar;
xlabel('time (s)'); ylabel('row (y)'); title('kymograph: cols collapsed');
sgtitle('Time maps: trace, temporal spectrum, space-time kymographs');

%% ===== FIGURE 4: SVD =====
% Decompose the signal into spatial eigenimages x temporal modes. Auto spatial-bin
% so the (double) matrix stays under SVD_MAX_BYTES.
sbin = 1;
while (Npix/sbin^2) * T2 * 8 > SVD_MAX_BYTES
    sbin = sbin + 1;
end
if sbin > 1
    cubeB = binCube(cube, sbin);
    fprintf('SVD: spatial-binned by %d -> %dx%d\n', sbin, size(cubeB,1), size(cubeB,2));
else
    cubeB = cube;
end
[Hb, Wb, ~] = size(cubeB);
Yd = reshape(cubeB, Hb*Wb, T2);                        % [pixels x time]
Yd = double(Yd) - mean(double(Yd), 2);                 % demean per pixel
[U, Sv, V] = svds(Yd, nSVD);
sv     = diag(Sv);
varExp = sv.^2 / sum(sv.^2);

figure('Name','4: SVD','Color','w','Position',[120 60 1280 820]);
subplot(3, nShowSVD, 1:2);
plot(100*varExp, 'o-', 'LineWidth',1.3); grid on;
xlabel('component'); ylabel('% variance'); title('singular spectrum');
subplot(3, nShowSVD, 3:nShowSVD);
plot(100*cumsum(varExp), 's-', 'LineWidth',1.3); grid on; ylim([0 100]);
xlabel('component'); ylabel('cumulative %'); title('cumulative variance explained');

for k = 1:nShowSVD
    comp = reshape(U(:,k), Hb, Wb);
    s = sign(sum(comp(:).^3) + eps);                   % stabilize sign for display
    subplot(3, nShowSVD, nShowSVD + k);
    imagesc(s*comp); axis image off; colormap(gca,parula);
    title(sprintf('mode %d (%.1f%%)', k, 100*varExp(k)));
    subplot(3, nShowSVD, 2*nShowSVD + k);
    plot(t, s*sv(k)*V(:,k), 'LineWidth', 0.7); grid on; xlim([t(1) t(end)]);
    xlabel('t (s)'); title(sprintf('mode %d time course', k));
end
sgtitle('SVD: spatial eigenimages (top) and their temporal modes (bottom)');

%% ===== FIGURE 5: SVD movement traces (DLC-free readout) =====
% Per-mode time course = sv_k * V(:,k). PC1 = dominant motion (= breathing here).
% You can't just ADD PC1..K (sign is arbitrary, modes are orthogonal = different
% motions). The meaningful combination is movement ENERGY = RMS across modes,
% which is sign-free but RECTIFIED (a pure oscillation appears at 2x its freq).
K   = min(nCombinePC, nSVD);
sgn = sign(sum(U(:,1:K).^3, 1) + eps);                 % stabilize each mode's sign
pcTrace    = (V(:,1:K) .* sgn) .* sv(1:K)';            % [T2 x K] signed, scaled PCs
moveEnergy = sqrt(sum((V(:,1:K) .* sv(1:K)').^2, 2));  % [T2 x 1] total movement
pcZ = (pcTrace - mean(pcTrace,1)) ./ std(pcTrace,0,1); % unit-variance per PC (display)

figure('Name','5: SVD movement traces','Color','w','Position',[140 60 1200 760]);
subplot(3,1,1); hold on;                                % (a) stacked top-K PCs
offStep = 6;
for k = K:-1:1
    plot(t, pcZ(:,k) + (k-1)*offStep, 'LineWidth', 0.7);
    text(t(1), (k-1)*offStep, sprintf(' PC%d (%.0f%%)', k, 100*varExp(k)), ...
        'FontWeight','bold', 'VerticalAlignment','bottom');
end
set(gca,'YTick',[]); grid on; xlim([t(1) t(end)]); xlabel('time (s)');
title(sprintf('top %d SVD temporal components (PC1 = dominant motion = breathing)', K));

subplot(3,1,2);                                         % (b) combined movement energy
plot(t, moveEnergy, 'k', 'LineWidth', 0.9); grid on; xlim([t(1) t(end)]);
xlabel('time (s)'); ylabel('movement energy');
title(sprintf('movement ENERGY = RMS of PC1..%d  (rectified -> doubles oscillation freq)', K));

subplot(3,1,3);                                         % (c) PSD: PC1 vs energy
[f1,P1] = tempPSD(pcTrace(:,1), fps, TW);
[fe,Pe] = tempPSD(moveEnergy,  fps, TW);
fhi2 = min(25, fps/2);
plot(f1, P1/max(P1), 'b', 'LineWidth',1.3); hold on;
plot(fe, Pe/max(Pe), 'r', 'LineWidth',1.3); grid on;
xlim([0 fhi2]); xlabel('frequency (Hz)'); ylabel('norm power');
legend('PC1 (signed: true freq)','energy (rectified: 2x)','Location','northeast');
title('PC1 gives the REAL breathing frequency; energy doubles it');
sgtitle('DLC-free breathing/movement readout from SVD');

%% Save results
[outDir, stem] = fileparts(OUT_MAT);
outFile = fullfile(outDir, [stem '_dFFanalysis.mat']);
save(outFile, 'sdMap', 'pkMap', 'meanFrame', 'gTrace', 'fP', 'P', ...
     'domFreqMap', 'breathPowMap', 'whiskPowMap', 'actMask', ...
     'U', 'sv', 'V', 'varExp', 'sbin', 'pcTrace', 'moveEnergy', 'nCombinePC', ...
     't', 'fps', 'roi', 'SIGNAL', 'BaselineWinSec', 'TW', '-v7.3');
fprintf('Saved analysis -> %s\n', outFile);

%% =============================== LOCAL FUNCTIONS ===============================
function [f, P] = tempPSD(x, Fs, TW)
% Temporal PSD: Chronux multitaper if available, else a simple Welch fallback.
    x = x(:) - mean(x(:));
    if exist('mtspectrumc', 'file') == 2
        params.Fs = Fs; params.pad = 0;
        params.tapers = [TW, 2*TW-1];
        params.fpass  = [0.1, min(25, Fs/2)];
        [P, f] = mtspectrumc(x, params);
        f = f(:); P = P(:);
    else
        seg = 2^floor(log2(numel(x)/8));                % ~8 segments
        seg = max(seg, 256);
        win = hann(seg); nov = floor(seg/2);
        step = seg - nov; nf = floor(seg/2)+1;
        acc = zeros(nf,1); cnt = 0;
        for s0 = 1:step:(numel(x)-seg+1)
            xw = x(s0:s0+seg-1) .* win;
            Xf = fft(xw); Pw = abs(Xf(1:nf)).^2;
            acc = acc + Pw; cnt = cnt + 1;
        end
        P = acc / max(cnt,1);
        f = (0:nf-1)' * (Fs/seg);
    end
end

function h = shadeBand(band, rgb, yhi)
% Shade a vertical frequency band [0 yhi] on the current axes; return handle.
    h = patch([band(1) band(2) band(2) band(1)], [0 0 yhi yhi], ...
              rgb, 'EdgeColor','none');
end

function cb = binCube(c, b)
% Block-mean spatial binning of an [H,W,T] cube by integer factor b.
    [H, W, T] = size(c);
    H2 = floor(H/b)*b; W2 = floor(W/b)*b;
    c = c(1:H2, 1:W2, :);
    c = reshape(c, b, H2/b, b, W2/b, T);
    cb = squeeze(mean(mean(c, 1), 3));                  % [H2/b, W2/b, T]
end
