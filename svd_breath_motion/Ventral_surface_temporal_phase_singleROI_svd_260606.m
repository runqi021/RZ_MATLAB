% Ventral_surface_temporal_phase_singleROI_svd_260606.m
% -----------------------------------------------------------------------
%  SINGLE-FILE temporal companion to chat_roi_quickview_260527.m.
%
%  Uses the EXACT piecewise-linear breath phase you already have
%  (piecewise_phase_local: inspiration ONSET = 0, breath PEAK = pi, linear
%  ramp in time between events) -- NOT a cosine fit, NOT Hilbert. The phase
%  analysis fed cos(phi) into Chronux coherence; here we use phi directly to
%  bin Ca-spikes and read out a preferred phase, then OVERLAY the coherence
%  result so the two methods can be compared on one polar axis.
%
%  Keeps the ENTIRE quickview analysis (same 7 panels, recomputed
%  identically) and ADDS three panels:
%
%    (A) AVG-PROJECTION crop of the ROI  [tile 10]
%        - average projection (AVG_*_MC_MC.tif, else mean of the MC_MC stack)
%        - ROI mask as a YELLOW boundary dilated +1 px (outline one pixel
%          outside the mask)
%        - radius = centroid -> FARTHEST boundary pixel
%        - square crop centered on centroid, half-side = radius + pad_um
%
%    (B) POLAR : Ca-spike probability vs piecewise breath phase, OVERLAID
%        with the coherence-method vector + jackknife error bar  [tile 11]
%          - gray rose       : P(spike | phase), occupancy-corrected,
%                               normalized to its max
%          - black arrow      : tuning resultant vector (mu, R) of that rose
%          - colored dot+bars : coherence method (th, r) with radial CI (Cerr)
%                               + phase-CI arc (1.96*phistd)
%          - dashed circle    : confC significance level
%        Convention matches the phase analysis: ONSET = 0, PEAK = pi.
%        0 marked RED (onset), pi marked SKY BLUE (peak) -- HSV 0/pi key.
%
%    (C) LINEAR twin of (B)  [tile 12]
%        Same occupancy-corrected spike-probability tuning curve as bars
%        across 0..2pi, with onset (0, red) + peak (pi, sky) reference lines
%        and the coherence preferred phase drawn as a vertical line, so the
%        direct-phase and cosine-coherence estimates line up side by side.
%
%  COHERENCE-METHOD overlay reproduces Ventral_surface_coherence_polar_260528
%  (cos(piecewise phase) x Ca-spike train, Chronux coherencyc, jackknife CI)
%  using that script's params (TW_coh, alpha_coh).
%
%  Inputs (in folderPath):
%     *_ch1_dFF.mat                      (dFF [T x N])
%     *DLC*breath_peak_data.mat          (breath waveform + PEAK idx)
%     *breath_insp_start_data.mat        (insp-start / foot idx)
%     ca_spike_data.mat                  (roi_spikes(roi).spike_train)
%     *_cpSAM_output.mat                 (maskL label image)
%     AVG_*_MC_MC.tif | *_preproc_MC_MC.tif  (avg projection / stack)
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m,
%                Image Processing TB (bwboundaries/imdilate).
% -----------------------------------------------------------------------

clear; close all; clc;

green = [0.2 0.7 0.2];
%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot=fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
folderPath = 'D:\Ventral_surface_summary\Vglut2_test\roi5_1400-1300-0_18lp_930_x4_512x512_6000f_00002';
roi        = 12;            % ROI index (== dFF column == maskL label)

nDrop      = 30;           % breath frames to toss (match calcium)
fallback_fps = 30;

% ---- quickview spectral params (PSD + breath x dFF coherence panels) ----
TW_spec    = 6;            % multitaper TW for PSD + dFF coherence
alpha_sig  = 0.01;

% ---- coherence-METHOD params (match Ventral_surface_coherence_polar) -----
TW_coh     = 4;            % multitaper TW for the spike-coherence overlay dot
alpha_coh  = 0.001;        % significance level for confC + jackknife err
minSpikes  = 2;            % require >= this many spikes to draw the dot
ca_lag_sec = 0.1;          % lead-comp: shift spikes EARLIER 3 frames @30Hz (0.1 s)

f_breath_search = [0.2 4]; % Hz, search band for breath PSD peak
fwhm_factor = 0.6;         % coherence/detection band = fwhm_factor x FWHM
min_bw      = 0.05;        % Hz, min band width
fmin        = 0.05;        % Hz, PSD/coherence lower bound
fmax        = 15;          % Hz, upper bound

% ---- phase-tuning params ----
nPhaseBins = 24;           % bins for the spike-probability rose / linear curve

% ---- avg-projection crop params ----
pad_um      = 20;          % square half-side = ROI radius + pad_um (um)
clip_pct    = [0.5 99.9];  % intensity clip percentiles for display
gamma_val   = 0.6;         % display gamma (<1 brightens midtones)
PixelSizeBase = 1.7778;    % um/px at zoom=1 (fallback if no pixelSize_um)
outlineLW   = 1.3;         % ROI boundary line width

% ---- trace panel ----
sortMode    = 'dt';        % 'postmean' | 'dt' | 'none'  (heatmap row order)
trace_xlim_sp = [];   % phase-sawtooth window (s); [] = full trace

% ---- colors ----
dffColor    = green;   % dF/F accent (quickview panels)
onsetCol    = [0.90 0.10 0.10];   % RED  : inspiration onset  (phase 0)
peakCol     = [0.35 0.75 1.00];   % SKY  : breath peak        (phase pi)

doSave      = true;
% =========================================================

lightDff = 0.30*dffColor + 0.70;
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
[~, recName] = fileparts(folderPath);

%% ===================== LOAD + ALIGN =====================
df  = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp  = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip  = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat')); end
sam = dir(fullfile(folderPath,'*_cpSAM_output.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No *DLC*breath_peak_data.mat in %s', folderPath);

[fps, sm] = detect_session_fps(folderPath, fallback_fps);
D   = load(fullfile(df(1).folder, df(1).name),'dFF');
BP  = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(D.dFF);
assert(roi>=1 && roi<=size(dff_all,2),'ROI %d out of range (1..%d)',roi,size(dff_all,2));

% pixel size (um/px) from RAW tiff metadata, with zoom / base fallback
px_um = NaN;
if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um>0
    px_um = sm.pixelSize_um;
elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor>0
    px_um = PixelSizeBase / sm.zoomFactor;
end

% breath waveform (toss nDrop, detrend, demean)
bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
nB = numel(BP.breath);

% PEAK event train (legacy "insp_onset_idx" == breath PEAK)
if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
    ev = double(BP.insp_onsets_train(:) ~= 0);
else
    ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
end
ev(1:min(nDrop,numel(ev))) = [];

% FOOT / inspiration-ONSET event train
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    ev_foot = zeros(nB,1); fi = round(IP.insp_start_idx(:));
    ev_foot(fi(fi>=1 & fi<=nB)) = 1;
    ev_foot(1:min(nDrop,numel(ev_foot))) = [];
else
    warning('No *breath_insp_start_data.mat -- onset overlays will be empty.');
    ev_foot = [];
end

% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay breath.
if contains(folderPath, fullfile('Vglut2','1124'))
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end

% spike train for this ROI
spk_train = [];
sp_file = fullfile(folderPath,'ca_spike_data.mat');
if isfile(sp_file)
    CA = load(sp_file,'roi_spikes');
    if isfield(CA,'roi_spikes') && roi <= numel(CA.roi_spikes)
        spk_train = double(CA.roi_spikes(roi).spike_train(:));
    end
end

% common length
T = min([size(dff_all,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
dff = dff_all(1:T, roi); bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); end
if ~isempty(spk_train)
    if numel(spk_train) < T, spk_train(end+1:T) = 0; end
    spk_train = spk_train(1:T);
else
    spk_train = zeros(T,1);
end
t = (0:T-1)'/fps;
fprintf('%s ROI%d | T=%d @%.3g Hz | %d peaks | %d onsets | %d spikes | %.4g um/px\n', ...
        recName, roi, T, fps, sum(ev), sum(ev_foot), sum(spk_train>0), px_um);

%% ===================== SPECTRA + BAND =====================
pB.Fs=fps; pB.tapers=[TW_spec,2*TW_spec-1]; pB.pad=0; pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=[2,alpha_sig];
[Sbw,fbw,SbwErr] = mtspectrumc(bw, pB);            Sbw=Sbw(:); fbw=fbw(:);
[Sdd,fdd,SddErr] = mtspectrumc(diff(dff)*fps, pB); Sdd=Sdd(:); fdd=fdd(:);

mm=fbw>=f_breath_search(1) & fbw<=f_breath_search(2);
[~,rl]=max(Sbw(mm)); ipk=find(mm,1)+rl-1; f_pk=fbw(ipk);
hh=Sbw(ipk)/2; lo=ipk; while lo>1&&Sbw(lo)>hh, lo=lo-1; end
hi=ipk;        while hi<numel(fbw)&&Sbw(hi)>hh, hi=hi+1; end
f_fwhm=[max(fbw(lo),f_breath_search(1)), min(fbw(hi),f_breath_search(2))];
bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
band=[max(f_pk-bwd/2,fmin), min(f_pk+bwd/2,fmax)];

%% ===================== COHERENCE (breath waveform x dF/F) =====================
pc.Fs=fps; pc.tapers=[TW_spec,2*TW_spec-1]; pc.pad=0; pc.fpass=[fmin,min(fmax,fps/2)]; pc.err=[2,alpha_sig];
[~,Cw,~,~,~,~,fcw,confCw,~,Cerrw] = coherencyc(bw, dff-mean(dff), pc);
fcw=fcw(:); Cw=Cw(:);
mb = fcw>=band(1) & fcw<=band(2); Cband = mean(Cw(mb));

%% ===================== PEAK-TRIGGERED  +  ONSET-TRIGGERED dF/F =====================
win = round(fps / f_pk);                 % +/- 1 breath period
tau = (-win:win)/fps;
on      = find(ev>0);        on      = on(on-win>=1 & on+win<=numel(dff));
foot_on = find(ev_foot>0);   if ~isempty(foot_on), foot_on = foot_on(foot_on-win>=1 & foot_on+win<=numel(dff)); end

E = zeros(numel(on), 2*win+1);
for k = 1:numel(on), E(k,:) = dff(on(k)-win : on(k)+win); end
dt_to_nearest_peak = nan(numel(on),1);
for k = 1:numel(on)
    others = on; others(k) = [];
    if isempty(others), continue; end
    [~, mi] = min(abs(others - on(k)));
    dt_to_nearest_peak(k) = (others(mi) - on(k)) / fps;
end
switch lower(sortMode)
    case 'postmean'
        key = mean(E(:, tau>=0), 2); [~,si] = sort(key,'descend');
        sortLbl = 'sorted by post-peak mean dF/F';
    case 'dt'
        [~,si] = sort(dt_to_nearest_peak,'ascend','MissingPlacement','last');
        sortLbl = 'sorted by dt to nearest peak';
    otherwise
        si = (1:numel(on))'; sortLbl = 'chronological';
end
Es = E(si,:);
cl = prctile(Es(:), [5 99.5]);
mu = mean(E,1); sd = std(E,0,1);

% --- spike-triggered dF/F (same +/- window as peak-triggered) ---
sp_trig = find(spk_train>0); sp_trig = sp_trig(sp_trig-win>=1 & sp_trig+win<=numel(dff));
if ~isempty(sp_trig)
    Esp = zeros(numel(sp_trig), 2*win+1);
    for k = 1:numel(sp_trig), Esp(k,:) = dff(sp_trig(k)-win : sp_trig(k)+win); end
    mu_sp = mean(Esp,1); sd_sp = std(Esp,0,1);
else
    Esp = zeros(0,2*win+1); mu_sp = []; sd_sp = [];
end

%% ===================== SPIKE FRAME INDICES (raster + onset-aligned time hist) =====================
sp_on_h = find(spk_train > 0);   % (onset-aligned time hist is built in tile_peri below)

%% ===================== EXACT BREATH PHASE + SPIKE TUNING =====================
% piecewise-linear phase (onset=0, peak=pi); NaN outside the event range.
peakAll = find(ev>0);
footAll = find(ev_foot>0);
phi_pw  = piecewise_phase_local(peakAll, footAll, T);

edgesP = linspace(0, 2*pi, nPhaseBins+1);
ctrsP  = (edgesP(1:end-1)+edgesP(2:end))/2;

% occupancy = time (frames) spent in each phase bin (insp/exp asymmetry ->
% non-uniform), used to turn spike counts into P(spike | phase).
valid = ~isnan(phi_pw);
occ   = histcounts(phi_pw(valid), edgesP);

% spikes, lead-compensated by ca_lag_sec (same as coherence method)
lagFrames = round(ca_lag_sec*fps);
spk_src = sp_on_h - lagFrames; spk_src = spk_src(spk_src>=1 & spk_src<=T);
spk_phi = phi_pw(spk_src); spk_phi = spk_phi(~isnan(spk_phi));
cntS    = histcounts(spk_phi, edgesP);

rate    = cntS ./ max(occ,1);              % P(spike | phase) per frame
rate_n  = rate / max(max(rate), eps);      % normalized for display [0,1]

% resultant vector of the (occupancy-corrected) tuning curve
wv  = rate(:); vec = sum(wv .* exp(1i*ctrsP(:)));
if sum(wv) > 0, muT = mod(angle(vec),2*pi); Rt = abs(vec)/sum(wv); else, muT = 0; Rt = 0; end

% --- per-cycle spike-phase distribution (matches Ventral_surface_spike_phase_polar) ---
phiW    = mod(phi_pw, 2*pi);                 % wrapped phase 0..2pi (display)
spkW    = mod(spk_phi, 2*pi);                % wrapped spike phase
cntW    = histcounts(spkW, edgesP);
phi_vv  = phi_pw(~isnan(phi_pw));
nCycles = (max(phi_vv) - min(phi_vv)) / (2*pi);   % breath cycles spanned
occW = histcounts(mod(phi_vv,2*pi), edgesP); Ob = mean(occW(occW>0));   % phase occupancy (dwell)
wOcc = ones(1,numel(occW)); wOcc(occW>0) = Ob ./ occW(occW>0);
pctC    = (100 * cntW / max(nCycles, eps)) .* wOcc;    % occupancy-compensated spk/cyc %
rmaxC   = max([pctC, 1]);
if ~isempty(spkW), muW = mod(angle(mean(exp(1i*spkW))),2*pi); else, muW = NaN; end

%% ===================== COHERENCE METHOD (cos(piecewise phase) x spikes) =====================
% reproduces Ventral_surface_coherence_polar_260528 for THIS ROI.
ref = cos(phi_pw); ref(isnan(ref)) = 0; ref = ref - mean(ref);

have_coh = sum(spk_train) >= minSpikes;
th_coh = NaN; r_coh = NaN; rlo = NaN; rhi = NaN; dphi = NaN; confC = NaN;
if have_coh
    pcoh.Fs=fps; pcoh.tapers=[TW_coh,2*TW_coh-1]; pcoh.pad=0;
    pcoh.fpass=band; pcoh.err=[2,alpha_coh];
    lag = round(ca_lag_sec*fps);
    stL = [spk_train(1+lag:end); zeros(lag,1)];   % lead-shift spikes earlier by ca_lag
    st  = stL - mean(stL);
    [~, Cmag, cphi, ~,~,~, fC, confC, phistd, Cerr] = coherencyc(ref, st, pcoh);
    fC = fC(:);
    mbc = fC>=band(1) & fC<=band(2); if ~any(mbc), mbc = true(size(fC)); end
    r_coh  = mean(Cmag(mbc));
    th_coh = angle(mean(exp(1i*(-cphi(mbc)))));
    rlo    = max(0, mean(Cerr(1,mbc)));
    rhi    = min(1, mean(Cerr(2,mbc)));
    dphi   = 1.96*mean(phistd(mbc));
end
th_coh2 = mod(th_coh, 2*pi);   % for the linear-axis guide

%% ===================== AVG PROJECTION + ROI MASK CROP =====================
have_proj = false; crop_img = []; bnd_crop = {}; barLen_pr = NaN; radius_um = NaN;
if ~isempty(sam) && isfinite(px_um)
    SAM = load(fullfile(sam(1).folder, sam(1).name), 'maskL');
    if isfield(SAM,'maskL')
        maskL = SAM.maskL;
        avgimg = read_avgproj_local(folderPath);
        if ~isempty(avgimg) && ~isequal(size(avgimg), size(maskL))
            warning('avg proj %s != maskL %s; using AVG_for_CP fallback.', ...
                    mat2str(size(avgimg)), mat2str(size(maskL)));
            fb = dir(fullfile(folderPath,'*_MC_MC_AVG_for_CP.tif'));
            if ~isempty(fb)
                V = tiffreadVolume(fullfile(fb(1).folder, fb(1).name));
                avgimg = mean(double(V),3);
            end
        end
        mask = (maskL == roi);
        if any(mask(:)) && ~isempty(avgimg) && isequal(size(avgimg),size(maskL))
            [H,W] = size(avgimg);
            [yy,xx] = find(mask);
            cx = mean(xx); cy = mean(yy);
            b0  = bwboundaries(mask,'noholes');
            ball = cat(1, b0{:});                      % [row col] = [y x]
            radius_px = max(hypot(ball(:,2)-cx, ball(:,1)-cy));
            half_px   = radius_px + pad_um/px_um;

            x0 = max(1, round(cx-half_px)); x1 = min(W, round(cx+half_px));
            y0 = max(1, round(cy-half_px)); y1 = min(H, round(cy+half_px));
            crop_raw = avgimg(y0:y1, x0:x1);

            loi = prctile(crop_raw(:), clip_pct(1));
            hii = prctile(crop_raw(:), clip_pct(2));
            crop_img = min(max((crop_raw-loi)/max(hii-loi,eps),0),1) .^ gamma_val;

            md = imdilate(mask, strel('square',3));    % +1 px outline
            bd = bwboundaries(md,'noholes');
            for k = 1:numel(bd)
                bd{k} = [bd{k}(:,1)-y0+1, bd{k}(:,2)-x0+1];
            end
            bnd_crop = bd;
            barLen_pr = max(1, round(20/px_um));        % 20 um scale bar (px)
            radius_um = radius_px*px_um;
            have_proj = true;
        end
    end
end

%% ===================== FIGURE =====================
% ---- TILE MAP (edit a tile number here to move that panel) -------------
nRows = 7;  nCols = 4;
tile_trace = 1;   span_trace = [2 4];   % rows 1-2 : dF/F + breath full trace
tile_saw   = 9;   span_saw   = [1 4];   % row 3    : phase sawtooth + events + spikes
tile_proj  = 13;                        % row 5    : avg-projection crop
tile_ptavg = 14;                        % row 4    : peak-triggered average
tile_psd   = 15;                        % row 4    : PSD
tile_cohSp = 16;                        % row 4    : breath x dF/F coherence
tile_pthm  = 17;                        % row 4    : peak-triggered heatmap
tile_peri  = 18;                        % row 5    : onset-aligned time hist
tile_phist = 19;                        % row 6    : linear 0..4pi

tile_epol  = 20;                        % row 6    : event-phase polar
tile_cpol  = 21;                        % row 6    : coherence polar
tile_stavg = 22;                        % row 6    : spike-triggered average
% ------------------------------------------------------------------------

fig = figure('Color','w','Name',sprintf('%s ROI%d temporal phase',recName,roi), ...
             'Units','normalized','Position',[0.03 0.03 0.9 0.8]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
title(tl, sprintf(['%s   ROI%d   |   fps=%.2f, breath peak %.2f Hz, band [%.2f %.2f] Hz, ', ...
      'coh-in-band=%.2f   |   RED=insp onset(0)  SKY=peak(\\pi)'], ...
      recName, roi, fps, f_pk, band(1), band(2), Cband), 'Interpreter','tex','FontWeight','bold');

% ===== tile_trace : dF/F + breath full trace =====
ax1 = nexttile(tl,tile_trace,span_trace);
yyaxis(ax1,'right');
plot(ax1, t, bw, '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
set(ax1,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax1,'breath');
yyaxis(ax1,'left');
plot(ax1, t, dff, '-','Color',dffColor,'LineWidth',0.8);
set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
xlim(ax1,[t(1) t(end)]); xlabel(ax1,'Time (s)'); box(ax1,'off');
%title(ax1,'dF/F (red) + breath (gray)  -- full trace');

% ===== tile_saw : phase sawtooth + events + spikes =====
if isempty(trace_xlim_sp), wsp=[t(1) t(end)]; else, wsp=[min(trace_xlim_sp) max(trace_xlim_sp)]; end
mwsp = t>=wsp(1) & t<=wsp(2);
bax = nexttile(tl,tile_saw,span_saw); hold(bax,'on');
plot(bax, t(mwsp), phiW(mwsp), 'k-','LineWidth',0.8);
set(bax,'YTick',[0 pi 2*pi],'YTickLabel',{'0','\pi','2\pi'}); ylabel(bax,'phase'); ylim(bax,[0 2*pi]);
yl = ylim(bax);
on_t = footAll/fps; pk_t = peakAll/fps;
on_t = on_t(on_t>=wsp(1)&on_t<=wsp(2)); pk_t = pk_t(pk_t>=wsp(1)&pk_t<=wsp(2));
for x=on_t(:)', plot(bax,[x x],yl,'-','Color',onsetCol,'LineWidth',0.8); end
for x=pk_t(:)', plot(bax,[x x],yl,'-','Color',peakCol,'LineWidth',0.8); end
spk_t = sp_on_h/fps; spk_t = spk_t(spk_t>=wsp(1)&spk_t<=wsp(2));
plot(bax, spk_t, (yl(1)+0.03*diff(yl))*ones(size(spk_t)),'|','Color','k','MarkerSize',6,'LineWidth',1);
xlim(bax,wsp); xlabel(bax,'Time (s)'); box(bax,'off');
%title(bax,'phase sawtooth (k) + onset (red) + peak (sky) + spikes (k ticks)');

% ===== tile_psd : PSD =====
ax2 = nexttile(tl,tile_psd); hold(ax2,'on');
%fill(ax2,[fbw;flipud(fbw)],10*log10([SbwErr(1,:)';flipud(SbwErr(2,:)')]),[0.85 0.85 0.85],'EdgeColor','none');
%fill(ax2,[fdd;flipud(fdd)],10*log10([SddErr(1,:)';flipud(SddErr(2,:)')]),lightDff,'EdgeColor','none');
plot(ax2, fbw,10*log10(Sbw),'Color',[0.5 0.5 0.5],'LineWidth',.8);
plot(ax2, fdd,10*log10(Sdd),'Color',dffColor,'LineWidth',0.8);
set(ax2,'XScale','log'); xlim(ax2,[fmin fmax]); xticks(ax2,[0.1 0.3 1 3 10]); set(ax2,'XMinorTick','off');
%add_band(ax2, band, f_pk); 
xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'power (dB)');
pbaspect(ax2,[1 1 1]);
%title(ax2,'PSD: breath (gray) + dF/F'' (color)');

% ===== tile_cohSp : breath x dF/F coherence =====
ax3 = nexttile(tl,tile_cohSp); hold(ax3,'on');
%fill(ax3,[fcw;flipud(fcw)],[Cerrw(1,:)';flipud(Cerrw(2,:)')],lightDff,'EdgeColor','none');
plot(ax3, fcw, Cw, 'Color',dffColor,'LineWidth',0.8);
yline(ax3, confCw,'k--','LineWidth',0.5);
set(ax3,'XScale','log'); xlim(ax3,[fmin fmax]); ylim(ax3,[0 1]); xticks(ax3,[0.1 0.3 1 3 10]); set(ax3,'XMinorTick','off');
add_band(ax3, band, f_pk); xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'coherence');
pbaspect(ax3,[1 1 1]);
%title(ax3, sprintf('breath x dF/F coherence (confC=%.2f)',confCw));

% ===== tile_ptavg : peak-triggered average =====
ax5 = nexttile(tl,tile_ptavg); hold(ax5,'on');
fill(ax5,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],lightDff,'EdgeColor','none');
plot(ax5, tau, mu, 'Color',dffColor,'LineWidth',1.5);
xline(ax5,0,'k--','LineWidth',0.8);
xlim(ax5,[tau(1) tau(end)]); xlabel(ax5,'time from breath peak (s)'); ylabel(ax5,'dF/F'); grid(ax5,'on');
pbaspect(ax5,[1 1 1]);
%title(ax5,'peak-triggered avg (\pm SD)');

% ===== tile_pthm : peak-triggered heatmap =====
ax4 = nexttile(tl,tile_pthm);
imagesc(ax4, tau, 1:size(Es,1), Es); axis(ax4,'tight');
colormap(ax4, flipud(gray(256))); caxis(ax4, cl);
hold(ax4,'on');
if ~isempty(foot_on)
    on_sorted = on(si); foot_x = []; foot_y = [];
    for k = 1:numel(on_sorted)
        p  = on_sorted(k);
        fw = foot_on(foot_on >= p-win & foot_on <= p+win);
        foot_x = [foot_x; (fw - p) / fps];        %#ok<AGROW>
        foot_y = [foot_y; repmat(k, numel(fw),1)]; %#ok<AGROW>
    end
    plot(ax4, foot_x, foot_y, '.', 'Color', onsetCol, 'MarkerSize', 1);
end
hold(ax4,'off');
set(ax4,'YDir','reverse'); xlabel(ax4,'time from breath peak (s)'); ylabel(ax4,'breath #');
pbaspect(ax4,[1 1 1]);
cb=colorbar(ax4); cb.Label.String='dF/F';
%title(ax4, sprintf('peak-triggered dF/F (n=%d, %s)', numel(on), sortLbl),'Interpreter','none');

% ===== tile_peri : ONSET-aligned time hist (spikes/cyc % + breath-peak distribution) =====
% 0 = inspiration ONSET; spikes on left axis (spk/cycle %, k), breath-PEAK times
% since the preceding onset on right axis (% of cycles, sky). Single cycle tiled to 1.5 period.
if ~isempty(footAll)
    timeWin = median(diff(footAll))/fps;                 % median breath period (s)
    if ~isfinite(timeWin) || timeWin<=0, timeWin = 1/max(f_pk,eps); end
    edgesTt = linspace(0, timeWin, nPhaseBins+1);
    ctrsTt  = (edgesTt(1:end-1)+edgesTt(2:end))/2;
    dtSpk   = time_from_prev_local(spk_src, footAll, T) / fps;   % spikes since onset
    dtPeak  = time_from_prev_local(peakAll, footAll, T) / fps;   % breath peak since onset
    spkH = 100*histcounts(dtSpk,  edgesTt)/max(nCycles,eps);     % spk/cycle %
    cP   = histcounts(dtPeak, edgesTt); pkH = 100*cP/max(sum(cP),1);  % peak distribution %
    nTileT = round(0.5*numel(ctrsTt));                           % tile out to 1.5 period
    ctrsX  = [ctrsTt, ctrsTt(1:nTileT)+timeWin];
    spkX   = [spkH, spkH(1:nTileT)];  pkX = [pkH, pkH(1:nTileT)];
    ax_h1 = nexttile(tl,tile_peri); hold(ax_h1,'on');
    yyaxis(ax_h1,'left');                                        % spikes
    bar(ax_h1, ctrsX, spkX, 1, 'FaceColor','k','FaceAlpha',0.85,'EdgeColor','none');
    set(ax_h1,'YColor','k'); ylim(ax_h1,[0 rmaxC*1.05]); ylabel(ax_h1,'spk/cyc %');
    yyaxis(ax_h1,'right');                                       % breath-peak distribution
    bar(ax_h1, ctrsX, pkX, 1, 'FaceColor',peakCol,'FaceAlpha',0.45,'EdgeColor','none');
    pkMax = max([pkX 1]); set(ax_h1,'YColor',peakCol*0.7); ylim(ax_h1,[0 pkMax*1.10]); ylabel(ax_h1,'peak %');
    xline(ax_h1, 0,       'Color',onsetCol,'LineWidth',1);       % insp ONSET at 0 (red)
    xline(ax_h1, timeWin, 'Color',onsetCol,'LineWidth',1);       % next onset at 1 period (red)
    xlim(ax_h1,[0 1.5*timeWin]);
    xlabel(ax_h1,'t from insp onset (s)');
    title(ax_h1, sprintf('%d spk, %.1f cyc, %.2f spk/cyc', ...
          sum(spk_train>0), nCycles, sum(spk_train>0)/max(nCycles,eps)),'Interpreter','none');
    pbaspect(ax_h1,[1 1 1]); box(ax_h1,'on');
end

% ===== tile_proj : avg-projection crop + yellow ROI boundary =====
ax_pr = nexttile(tl,tile_proj);
if have_proj
    imagesc(ax_pr, crop_img); colormap(ax_pr, gray(256)); caxis(ax_pr,[0 1]);
    axis(ax_pr,'image'); set(ax_pr,'XTick',[],'YTick',[]); hold(ax_pr,'on');
    for k = 1:numel(bnd_crop)
        plot(ax_pr, bnd_crop{k}(:,2), bnd_crop{k}(:,1), '-', 'Color',[1 1 0],'LineWidth',outlineLW);
    end
    [Hc,~] = size(crop_img);
    mgn = round(0.05*Hc); thk = max(2,round(0.02*Hc));
    rectangle(ax_pr,'Position',[mgn Hc-mgn-thk barLen_pr thk],'FaceColor','w','EdgeColor','none');
    hold(ax_pr,'off');
    %title(ax_pr, sprintf('avg proj + ROI (r=%.1f um, crop +/-%.0f um)', radius_um, pad_um),'Interpreter','none');
else
    axis(ax_pr,'off');
    text(ax_pr,0.5,0.5,'(no mask / avg proj)','Horizontal','center');
end

% ===== tile_epol : EVENT-PHASE polar (spikes/cycle %) =====
pax2 = polaraxes(tl); pax2.Layout.Tile = tile_epol; hold(pax2,'on');
polarplot(pax2, [ctrsP ctrsP(1)], [pctC pctC(1)], '-', ...
        'Color',[0.30 0.30 0.30],'LineWidth',1.4);
pax2.RLim=[0 rmaxC]; pax2.ThetaZeroLocation='right'; pax2.ThetaDir='counterclockwise';
pax2.RAxisLocation=180; pax2.FontSize=8; thetaticks(pax2,0:45:315);
title(pax2, sprintf('%.0f\\circ', rad2deg(muW)),'Interpreter','tex');

% ===== tile_cpol : COHERENCE polar (preferred phase + magnitude + jackknife CI) =====
cax = polaraxes(tl); cax.Layout.Tile = tile_cpol; hold(cax,'on');
if have_coh && isfinite(th_coh)
    thcc = linspace(0,2*pi,361);
    if isfinite(confC), polarplot(cax, thcc, confC*ones(size(thcc)),'k--','LineWidth',1); end
    polarplot(cax,[th_coh th_coh],[rlo rhi],'-','Color','k','LineWidth',1.8);        % magnitude CI
    if isfinite(dphi) && dphi>0
        arcc = linspace(th_coh-dphi, th_coh+dphi, 60);
        polarplot(cax, arcc, r_coh*ones(size(arcc)),'-','Color','k','LineWidth',1.8); % phase CI arc
    end
    polarplot(cax, th_coh, r_coh,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
end
cax.RLim=[0 1]; cax.ThetaZeroLocation='right'; cax.ThetaDir='counterclockwise';
cax.RAxisLocation=180; cax.FontSize=8; thetaticks(cax,0:45:315);
title(cax, sprintf('%.0f\\circ r=%.2f conf=%.2f', ...
      rad2deg(mod(th_coh,2*pi)), r_coh, confC),'Interpreter','tex');

% ===== tile_phist : LINEAR 0..2pi (single cycle), spikes/cycle % =====
lax2 = nexttile(tl,tile_phist); hold(lax2,'on');
bar(lax2, ctrsP, pctC, 1, 'FaceColor',[0.30 0.30 0.30],'FaceAlpha',0.75,'EdgeColor','none');
xline(lax2, 0,  'Color',onsetCol,'LineWidth',1);    % insp ONSET (red)
xline(lax2, pi, 'Color',peakCol, 'LineWidth',1);    % breath PEAK (sky)
xlim(lax2,[0 2*pi]); ylim(lax2,[0 rmaxC*1.05]);
set(lax2,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
xlabel(lax2,'breath phase (onset=0, peak=\pi)'); ylabel(lax2,'spikes/cycle (%)');
box(lax2,'on'); pbaspect(lax2,[1 1 1]);
title(lax2, sprintf('%.1f cycles, %.2f spk/cycle', ...
      nCycles, sum(cntW)/max(nCycles,eps)),'Interpreter','tex');

% ===== tile_stavg : spike-triggered average (same style as peak-triggered) =====
ax_sp = nexttile(tl,tile_stavg); hold(ax_sp,'on');
if ~isempty(mu_sp)
    fill(ax_sp,[tau fliplr(tau)],[mu_sp+sd_sp fliplr(mu_sp-sd_sp)],lightDff,'EdgeColor','none');
    plot(ax_sp, tau, mu_sp, 'Color',dffColor,'LineWidth',1.5);
    xline(ax_sp,0,'k--','LineWidth',0.8);
    xlim(ax_sp,[tau(1) tau(end)]);
end
xlabel(ax_sp,'time from spike (s)'); ylabel(ax_sp,'dF/F'); grid(ax_sp,'on');
pbaspect(ax_sp,[1 1 1]);
%title(ax_sp, sprintf('spike-triggered average (\\pm SD, n=%d)', size(Esp,1)));

%% ===================== SAVE =====================
if doSave
    base = fullfile(folderPath, sprintf('temporal_phase_svd_ROI%02d', roi));
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('Saved %s.png/.pdf\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function add_band(ax, band, f_pk)
    yl = ylim(ax);
    p = patch(ax, [band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.90 0.90 0.90], 'EdgeColor','none');
    xline(ax, f_pk, 'k--', 'LineWidth',1);
    uistack(p,'bottom'); ylim(ax, yl);
end

function avgimg = read_avgproj_local(folderPath)
% Prefer precomputed AVG_*_MC_MC.tif (single frame); else mean of MC_MC stack.
    avgimg = [];
    a = dir(fullfile(folderPath,'AVG_*_ch1_preproc_MC_MC.tif'));
    if isempty(a), a = dir(fullfile(folderPath,'AVG_*_MC_MC.tif')); end
    if ~isempty(a)
        V = tiffreadVolume(fullfile(a(1).folder, a(1).name));
        avgimg = mean(double(V),3); return;
    end
    s = dir(fullfile(folderPath,'*_ch1_preproc_MC_MC.tif'));
    s = s(~contains({s.name},'AVG','IgnoreCase',true));
    if ~isempty(s)
        V = tiffreadVolume(fullfile(s(1).folder, s(1).name));
        avgimg = mean(double(V),3);
    end
end

function dt = time_from_prev_local(frames, refIdx, T)
% time (frames) of each event since its PRECEDING reference event (0..period);
% NaN if before the first reference event.
dt = nan(numel(frames),1);
if isempty(refIdx) || isempty(frames), return; end
b = discretize(frames(:), [refIdx(:); T+1]);
ok = ~isnan(b); dt(ok) = frames(ok) - refIdx(b(ok));
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase reference: FEET at 0/2pi/..., PEAKS at pi/3pi/...
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
keep = true(size(events));
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases  = nan(size(events));
phi_cur = types(1) * pi;      % type=1 (peak) -> pi; type=0 (foot) -> 0
for i = 1:numel(events)
    phases(i) = phi_cur; phi_cur = phi_cur + pi;
end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
