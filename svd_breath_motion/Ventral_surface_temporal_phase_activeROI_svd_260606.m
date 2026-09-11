% Ventral_surface_temporal_phase_activeROI_svd_260606.m
% -----------------------------------------------------------------------
%  ACTIVE-ROI batch of Ventral_surface_temporal_phase.
%  Renders ONE temporal-phase figure per ACTIVE ROI (active = MORE THAN
%  minEvents spike events, i.e. nnz(spike_train>0) > minEvents), across all
%  ROIs that entered the coherence-polar analysis on the Ventral_surface_summary
%  dataset.
%
%  ROI list source:
%     <rootPath>\coherence_polar_svd_260606\coherence_polar_data.mat
%        labels{i} = 'group/date/recName/roi'  (every ROI used in coherence)
%  Each label is resolved to a recording folder by recursive search for
%  recName under <rootPath>\<group>.
%
%  Per ROI: full analysis + 12-panel figure (identical to the single-file
%  script) saved as:
%     <ROIfolder>\temporal_phase_svd_ROI##.png/.pdf
%     <rootPath>\temporal_phase_activeROI_svd_260606\active##_group_date_roi##.png/.pdf
%  Significant-ROI-only version: Ventral_surface_temporal_phase_sigROI_svd_260606.m
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m,
%                Image Processing TB (bwboundaries/imdilate).
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot=fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'temporal_phase_activeROI_svd_260606');

P = struct();
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.TW_coh       = 4;           P.alpha_coh    = 0.001;
P.minSpikes    = 2;           P.ca_lag_sec   = 0.1;    % lead-comp: spikes 3 frames earlier @30Hz (0.1 s)
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
P.nPhaseBins   = 24;
P.pad_um       = 20;          P.clip_pct     = [0.5 99.9];
P.gamma_val    = 0.6;         P.PixelSizeBase = 1.7778;  P.outlineLW = 1.3;
P.sortMode     = 'dt';        P.trace_xlim_sp = [];
P.dffColor     = [0.2 0.7 0.2];          % green accent
P.onsetCol     = [0.90 0.10 0.10];       % red  : onset (phase 0)
P.peakCol      = [0.35 0.75 1.00];       % sky  : peak  (phase pi)

minEvents  = 5;              % ACTIVE ROI = more than this many spike events (nnz>0)
doSave     = true;
closeAfter = true;           % close each figure after saving (avoid 11 open)
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ===================== ROI LIST (all coherence ROIs) =====================
assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData, 'labels');
fprintf('%d ROIs in coherence list | ACTIVE = more than %d events\n', numel(S.labels), minEvents);
caCache = containers.Map('KeyType','char','ValueType','any');   % per-folder roi_spikes

%% ===================== LOOP =====================
nActive = 0; nSkip = 0; nFail = 0;
for k = 1:numel(S.labels)
    lab   = S.labels{k};
    parts = regexp(lab, '/', 'split');
    if numel(parts) < 4
        warning('skip malformed label: %s', lab); nFail = nFail+1; continue;
    end
    group   = parts{1};
    recDate = parts{2};
    recName = strjoin(parts(3:end-1), '/');     % recName may (rarely) contain '/'
    roi     = str2double(parts{end});

    folderPath = resolve_folder(rootPath, group, recName);
    if isempty(folderPath)
        warning('skip (folder not found): %s', lab); nFail = nFail+1; continue;
    end

    % ---- activity gate: count calcium spike events for this ROI ----
    nEvents = count_events(folderPath, roi, caCache);
    if nEvents <= minEvents
        fprintf('[%d/%d] %s | ROI %d  skip (inactive, %d events)\n', ...
                k, numel(S.labels), lab, roi, nEvents);
        nSkip = nSkip + 1; continue;
    end

    nActive = nActive + 1;
    fprintf('\n[%d/%d] active#%d %s | ROI %d (%d events)\n   %s\n', ...
            k, numel(S.labels), nActive, lab, roi, nEvents, folderPath);
    try
        fig = make_temporal_phase_fig(folderPath, roi, recName, group, recDate, P);
        if doSave
            base = fullfile(folderPath, sprintf('temporal_phase_svd_ROI%02d', roi));
            exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
            exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            cbase = fullfile(outDir, sprintf('active%02d_%s_%s_roi%02d', nActive, group, recDate, roi));
            exportgraphics(fig, [cbase '.png'], 'Resolution',200, 'BackgroundColor','white');
            exportgraphics(fig, [cbase '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            fprintf('   saved ROI%02d (folder + %s)\n', roi, 'batch dir');
        end
        if closeAfter, close(fig); end
    catch ME
        warning('   ERROR %s ROI %d: %s', lab, roi, ME.message); nFail = nFail+1;
        if exist('fig','var') && ishandle(fig) && closeAfter, close(fig); end
    end
end
fprintf('\nDone. %d active figures rendered, %d skipped (inactive), %d failed. Batch dir: %s\n', ...
        nActive, nSkip, nFail, outDir);

%% ===================== ACTIVITY COUNT =====================
function nEvents = count_events(folderPath, roi, caCache)
% Number of calcium spike events for this ROI (nnz of its spike_train).
% roi_spikes is cached per folder to avoid reloading for each ROI.
    if ~isKey(caCache, folderPath)
        spf = fullfile(folderPath, 'ca_spike_data.mat');
        if isfile(spf)
            tmp = load(spf, 'roi_spikes'); caCache(folderPath) = tmp.roi_spikes;
        else
            caCache(folderPath) = [];
        end
    end
    rs = caCache(folderPath);
    if isempty(rs) || roi < 1 || roi > numel(rs)
        nEvents = 0;
    else
        nEvents = nnz(double(rs(roi).spike_train) > 0);
    end
end

%% ===================== FOLDER RESOLUTION =====================
function folderPath = resolve_folder(rootPath, group, recName)
% Recursively find a directory named recName under rootPath\group.
    folderPath = '';
    d = dir(fullfile(rootPath, group, '**', recName));
    d = d([d.isdir]);
    if isempty(d)
        d = dir(fullfile(rootPath, '**', recName)); d = d([d.isdir]);   % fallback: anywhere
    end
    if ~isempty(d)
        folderPath = fullfile(d(1).folder, d(1).name);
    end
end

%% ===================== PER-ROI FIGURE (port of single-file script) =====================
function fig = make_temporal_phase_fig(folderPath, roi, recName, group, recDate, P)
% Faithful port of Ventral_surface_temporal_phase_260601.m for ONE ROI.
lightDff = 0.30*P.dffColor + 0.70;
dffColor = P.dffColor; onsetCol = P.onsetCol; peakCol = P.peakCol;

%% ---- LOAD + ALIGN ----
df  = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp  = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
if isempty(bp), bp = dir(fullfile(folderPath,'breath_peak_pc1.mat')); end
ip  = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat')); end
sam = dir(fullfile(folderPath,'*_cpSAM_output.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No *DLC*breath_peak_data.mat in %s', folderPath);

[fps, sm] = detect_session_fps(folderPath, P.fallback_fps);
D   = load(fullfile(df(1).folder, df(1).name),'dFF');
BP  = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(D.dFF);
assert(roi>=1 && roi<=size(dff_all,2),'ROI %d out of range (1..%d)',roi,size(dff_all,2));

px_um = NaN;
if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um>0
    px_um = sm.pixelSize_um;
elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor>0
    px_um = P.PixelSizeBase / sm.zoomFactor;
end

bw = detrend(double(BP.breath(:))); bw(1:min(P.nDrop,numel(bw))) = []; bw = bw - mean(bw);
nB = numel(BP.breath);

if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
    ev = double(BP.insp_onsets_train(:) ~= 0);
else
    ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
end
ev(1:min(P.nDrop,numel(ev))) = [];

if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    ev_foot = zeros(nB,1); fi = round(IP.insp_start_idx(:));
    ev_foot(fi(fi>=1 & fi<=nB)) = 1; ev_foot(1:min(P.nDrop,numel(ev_foot))) = [];
else
    ev_foot = [];
end

% Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay breath.
if strcmpi(group,'Vglut2') && strcmp(recDate,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end

spk_train = [];
sp_file = fullfile(folderPath,'ca_spike_data.mat');
if isfile(sp_file)
    CA = load(sp_file,'roi_spikes');
    if isfield(CA,'roi_spikes') && roi <= numel(CA.roi_spikes)
        spk_train = double(CA.roi_spikes(roi).spike_train(:));
    end
end

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

%% ---- SPECTRA + BAND ----
pB.Fs=fps; pB.tapers=[P.TW_spec,2*P.TW_spec-1]; pB.pad=0; pB.fpass=[P.fmin,min(P.fmax,fps/2)]; pB.err=[2,P.alpha_sig];
[Sbw,fbw,SbwErr] = mtspectrumc(bw, pB);            Sbw=Sbw(:); fbw=fbw(:); %#ok<ASGLU>
[Sdd,fdd,SddErr] = mtspectrumc(diff(dff)*fps, pB); Sdd=Sdd(:); fdd=fdd(:); %#ok<ASGLU>

mm=fbw>=P.f_breath_search(1) & fbw<=P.f_breath_search(2);
[~,rl]=max(Sbw(mm)); ipk=find(mm,1)+rl-1; f_pk=fbw(ipk);
hh=Sbw(ipk)/2; lo=ipk; while lo>1&&Sbw(lo)>hh, lo=lo-1; end
hi=ipk;        while hi<numel(fbw)&&Sbw(hi)>hh, hi=hi+1; end
f_fwhm=[max(fbw(lo),P.f_breath_search(1)), min(fbw(hi),P.f_breath_search(2))];
bwd=max(diff(f_fwhm)*P.fwhm_factor, P.min_bw);
band=[max(f_pk-bwd/2,P.fmin), min(f_pk+bwd/2,P.fmax)];

%% ---- breath x dF/F coherence ----
pc.Fs=fps; pc.tapers=[P.TW_spec,2*P.TW_spec-1]; pc.pad=0; pc.fpass=[P.fmin,min(P.fmax,fps/2)]; pc.err=[2,P.alpha_sig];
[~,Cw,~,~,~,~,fcw,confCw,~,Cerrw] = coherencyc(bw, dff-mean(dff), pc);
fcw=fcw(:); Cw=Cw(:);
mb = fcw>=band(1) & fcw<=band(2); Cband = mean(Cw(mb));

%% ---- peak-triggered + spike-triggered dF/F ----
win = round(fps / f_pk); tau = (-win:win)/fps;
on      = find(ev>0);        on      = on(on-win>=1 & on+win<=numel(dff));
foot_on = find(ev_foot>0);   if ~isempty(foot_on), foot_on = foot_on(foot_on-win>=1 & foot_on+win<=numel(dff)); end

E = zeros(numel(on), 2*win+1);
for k = 1:numel(on), E(k,:) = dff(on(k)-win : on(k)+win); end
dt_to_nearest_peak = nan(numel(on),1);
for k = 1:numel(on)
    others = on; others(k) = [];
    if isempty(others), continue; end
    [~, mi] = min(abs(others - on(k))); dt_to_nearest_peak(k) = (others(mi) - on(k)) / fps;
end
switch lower(P.sortMode)
    case 'postmean', key = mean(E(:, tau>=0), 2); [~,si] = sort(key,'descend');
    case 'dt',       [~,si] = sort(dt_to_nearest_peak,'ascend','MissingPlacement','last');
    otherwise,       si = (1:numel(on))';
end
Es = E(si,:); cl = prctile(Es(:), [5 99.5]);
mu = mean(E,1); sd = std(E,0,1);

sp_trig = find(spk_train>0); sp_trig = sp_trig(sp_trig-win>=1 & sp_trig+win<=numel(dff));
if ~isempty(sp_trig)
    Esp = zeros(numel(sp_trig), 2*win+1);
    for k = 1:numel(sp_trig), Esp(k,:) = dff(sp_trig(k)-win : sp_trig(k)+win); end
    mu_sp = mean(Esp,1); sd_sp = std(Esp,0,1);
else
    Esp = zeros(0,2*win+1); mu_sp = []; sd_sp = [];
end

%% ---- spike frame indices (raster + onset-aligned time hist) ----
sp_on_h = find(spk_train > 0);   % (onset-aligned time hist is built in tile_peri below)

%% ---- exact breath phase + per-cycle spike distribution ----
peakAll = find(ev>0); footAll = find(ev_foot>0);
phi_pw  = piecewise_phase_local(peakAll, footAll, T);
edgesP = linspace(0, 2*pi, P.nPhaseBins+1); ctrsP = (edgesP(1:end-1)+edgesP(2:end))/2;

lagFrames = round(P.ca_lag_sec*fps);
spk_src = sp_on_h - lagFrames; spk_src = spk_src(spk_src>=1 & spk_src<=T);
spk_phi = phi_pw(spk_src); spk_phi = spk_phi(~isnan(spk_phi));

phiW    = mod(phi_pw, 2*pi);
spkW    = mod(spk_phi, 2*pi);
cntW    = histcounts(spkW, edgesP);
phi_vv  = phi_pw(~isnan(phi_pw));
nCycles = (max(phi_vv) - min(phi_vv)) / (2*pi);
occW = histcounts(mod(phi_vv,2*pi), edgesP); Ob = mean(occW(occW>0));   % phase occupancy (dwell)
wOcc = ones(1,numel(occW)); wOcc(occW>0) = Ob ./ occW(occW>0);
pctC    = (100 * cntW / max(nCycles, eps)) .* wOcc;    % occupancy-compensated spk/cyc %
rmaxC   = max([pctC, 1]);
if ~isempty(spkW), muW = mod(angle(mean(exp(1i*spkW))),2*pi); else, muW = NaN; end

%% ---- coherence method: BIT-FOR-BIT replica of Ventral_surface_coherence_polar_svd_260606 ----
% Independent of the TW_spec display band above: TW_coh PSD, band clamped to
% f_breath_search, Vglut2/1124 +1-frame fix, Tc=min(numel(bw),numel(spikes)).
have_coh = sum(spk_train) >= P.minSpikes;
th_coh = NaN; r_coh = NaN; rlo = NaN; rhi = NaN; dphi = NaN; confC = NaN;
if have_coh && ~isempty(ip)
    nD = P.nDrop;
    bwc = detrend(double(BP.breath(:))); bwc(1:min(nD,numel(bwc))) = []; bwc = bwc - mean(bwc);
    peakC = round(BP.insp_onset_idx(:)) - nD;     % breath PEAK (legacy name)
    footC = round(IP.insp_start_idx(:)) - nD;     % insp ONSET (foot)
    stkC  = double(CA.roi_spikes(roi).spike_train(:));
    if strcmpi(group,'Vglut2') && strcmp(recDate,'1124')   % rising-edge trigger fix
        peakC = peakC + 1; footC = footC + 1; bwc = [bwc(1); bwc(1:end-1)];
    end
    Tc = min(numel(bwc), numel(stkC));
    peakC = peakC(peakC>=1 & peakC<=Tc); footC = footC(footC>=1 & footC<=Tc);
    bwc = bwc(1:Tc); stc = stkC(1:Tc);
    if numel(peakC) >= 2 && numel(footC) >= 2
        % FWHM band from TW_coh PSD, clamped to f_breath_search
        pBc.Fs=fps; pBc.tapers=[P.TW_coh,2*P.TW_coh-1]; pBc.pad=0;
        pBc.fpass=[P.fmin,min(P.fmax,fps/2)]; pBc.err=0;
        [Sbc,fbc] = mtspectrumc(bwc, pBc); Sbc=Sbc(:); fbc=fbc(:);
        mc = fbc>=P.f_breath_search(1) & fbc<=P.f_breath_search(2);
        [~,rlc]=max(Sbc(mc)); ipc=find(mc,1)+rlc-1; f_pk_c=fbc(ipc);
        h2=Sbc(ipc)/2; loc=ipc; while loc>1&&Sbc(loc)>h2, loc=loc-1; end
        hic=ipc;        while hic<numel(fbc)&&Sbc(hic)>h2, hic=hic+1; end
        ff=[max(fbc(loc),P.f_breath_search(1)), min(fbc(hic),P.f_breath_search(2))];
        bwd2=max(diff(ff)*P.fwhm_factor, P.min_bw);
        band_c=[max(f_pk_c-bwd2/2,P.f_breath_search(1)), min(f_pk_c+bwd2/2,P.f_breath_search(2))];
        % coherency cos(phi) x spike train
        phiC = piecewise_phase_local(peakC, footC, Tc);
        refC = cos(phiC); refC(isnan(refC)) = 0; refC = refC - mean(refC);
        pcoh.Fs=fps; pcoh.tapers=[P.TW_coh,2*P.TW_coh-1]; pcoh.pad=0;
        pcoh.fpass=band_c; pcoh.err=[2,P.alpha_coh];
        lagC = round(P.ca_lag_sec*fps);
        stcL = [stc(1+lagC:end); zeros(lagC,1)];   % lead-shift spikes earlier
        [~, Cmag, cphi, ~,~,~, fC, confC, phistd, Cerr] = coherencyc(refC, stcL-mean(stcL), pcoh);
        fC = fC(:); mbc = fC>=band_c(1) & fC<=band_c(2); if ~any(mbc), mbc = true(size(fC)); end
        r_coh  = mean(Cmag(mbc));
        th_coh = angle(mean(exp(1i*(-cphi(mbc)))));
        rlo    = max(0, mean(Cerr(1,mbc))); rhi = min(1, mean(Cerr(2,mbc)));
        dphi   = 1.96*mean(phistd(mbc));
    end
end
fprintf('coh r=%.3f  th=%.0f deg  confC=%.3f\n', r_coh, rad2deg(mod(th_coh,2*pi)), confC);

%% ---- avg projection + ROI mask crop ----
have_proj = false; crop_img = []; bnd_crop = {}; barLen_pr = NaN;
if ~isempty(sam) && isfinite(px_um)
    SAM = load(fullfile(sam(1).folder, sam(1).name), 'maskL');
    if isfield(SAM,'maskL')
        maskL = SAM.maskL;
        avgimg = read_avgproj_local(folderPath);
        if ~isempty(avgimg) && ~isequal(size(avgimg), size(maskL))
            fb = dir(fullfile(folderPath,'*_MC_MC_AVG_for_CP.tif'));
            if ~isempty(fb)
                V = tiffreadVolume(fullfile(fb(1).folder, fb(1).name)); avgimg = mean(double(V),3);
            end
        end
        mask = (maskL == roi);
        if any(mask(:)) && ~isempty(avgimg) && isequal(size(avgimg),size(maskL))
            [Himg,Wimg] = size(avgimg);
            [yy,xx] = find(mask); cx = mean(xx); cy = mean(yy);
            b0  = bwboundaries(mask,'noholes'); ball = cat(1, b0{:});
            radius_px = max(hypot(ball(:,2)-cx, ball(:,1)-cy));
            half_px   = radius_px + P.pad_um/px_um;
            x0 = max(1, round(cx-half_px)); x1 = min(Wimg, round(cx+half_px));
            y0 = max(1, round(cy-half_px)); y1 = min(Himg, round(cy+half_px));
            crop_raw = avgimg(y0:y1, x0:x1);
            loi = prctile(crop_raw(:), P.clip_pct(1)); hii = prctile(crop_raw(:), P.clip_pct(2));
            crop_img = min(max((crop_raw-loi)/max(hii-loi,eps),0),1) .^ P.gamma_val;
            md = imdilate(mask, strel('square',3)); bd = bwboundaries(md,'noholes');
            for k = 1:numel(bd), bd{k} = [bd{k}(:,1)-y0+1, bd{k}(:,2)-x0+1]; end
            bnd_crop = bd; barLen_pr = max(1, round(20/px_um)); have_proj = true;
        end
    end
end

%% ===================== FIGURE =====================
% ---- TILE MAP (7x4) ----
nRows = 7;  nCols = 4;
tile_trace = 1;   span_trace = [2 4];
tile_saw   = 9;   span_saw   = [1 4];
tile_proj  = 13;  tile_ptavg = 14;  tile_psd  = 15;  tile_cohSp = 16;
tile_pthm  = 17;  tile_peri  = 18;  tile_phist = 19;
tile_epol  = 20;  tile_cpol  = 21;  tile_stavg = 22;

fig = figure('Color','w','Name',sprintf('%s ROI%d temporal phase',folderPath,roi), ...
             'Units','normalized','Position',[0.03 0.03 0.9 0.8]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
relPath = sprintf('%s/%s/%s', group, recDate, recName);   % short label, no C:\...
title(tl, sprintf(['%s   ROI%d   |   fps=%.2f, breath peak %.2f Hz'], ...
      relPath, roi, fps, f_pk), 'Interpreter','none','FontWeight','bold');

% trace
ax1 = nexttile(tl,tile_trace,span_trace);
yyaxis(ax1,'right'); plot(ax1, t, bw, '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
set(ax1,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax1,'breath');
yyaxis(ax1,'left'); plot(ax1, t, dff, '-','Color',dffColor,'LineWidth',0.8);
set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
xlim(ax1,[t(1) t(end)]); xlabel(ax1,'Time (s)'); box(ax1,'off');

% phase sawtooth
if isempty(P.trace_xlim_sp), wsp=[t(1) t(end)]; else, wsp=[min(P.trace_xlim_sp) max(P.trace_xlim_sp)]; end
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

% PSD
ax2 = nexttile(tl,tile_psd); hold(ax2,'on');
plot(ax2, fbw,10*log10(Sbw),'Color',[0.5 0.5 0.5],'LineWidth',.8);
plot(ax2, fdd,10*log10(Sdd),'Color',dffColor,'LineWidth',0.8);
set(ax2,'XScale','log'); xlim(ax2,[P.fmin P.fmax]); xticks(ax2,[0.1 0.3 1 3 10]); set(ax2,'XMinorTick','off');
xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'power (dB)'); pbaspect(ax2,[1 1 1]);

% breath x dF/F coherence
ax3 = nexttile(tl,tile_cohSp); hold(ax3,'on');
plot(ax3, fcw, Cw, 'Color',dffColor,'LineWidth',0.8);
yline(ax3, confCw,'k--','LineWidth',0.5);
set(ax3,'XScale','log'); xlim(ax3,[P.fmin P.fmax]); ylim(ax3,[0 1]); xticks(ax3,[0.1 0.3 1 3 10]); set(ax3,'XMinorTick','off');
add_band(ax3, band, f_pk); xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'coherence'); pbaspect(ax3,[1 1 1]);

% peak-triggered average
ax5 = nexttile(tl,tile_ptavg); hold(ax5,'on');
fill(ax5,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],lightDff,'EdgeColor','none');
plot(ax5, tau, mu, 'Color',dffColor,'LineWidth',1.5); xline(ax5,0,'k--','LineWidth',0.8);
xlim(ax5,[tau(1) tau(end)]); xlabel(ax5,'time from breath peak (s)'); ylabel(ax5,'dF/F'); grid(ax5,'on'); pbaspect(ax5,[1 1 1]);

% peak-triggered heatmap
ax4 = nexttile(tl,tile_pthm);
imagesc(ax4, tau, 1:size(Es,1), Es); axis(ax4,'tight');
colormap(ax4, flipud(gray(256))); caxis(ax4, cl); hold(ax4,'on');
if ~isempty(foot_on)
    on_sorted = on(si); foot_x = []; foot_y = [];
    for k = 1:numel(on_sorted)
        p  = on_sorted(k); fw = foot_on(foot_on >= p-win & foot_on <= p+win);
        foot_x = [foot_x; (fw - p) / fps]; foot_y = [foot_y; repmat(k, numel(fw),1)]; %#ok<AGROW>
    end
    plot(ax4, foot_x, foot_y, '.', 'Color', onsetCol, 'MarkerSize', 1);
end
hold(ax4,'off'); set(ax4,'YDir','reverse'); xlabel(ax4,'time from breath peak (s)'); ylabel(ax4,'breath #');
pbaspect(ax4,[1 1 1]); cb=colorbar(ax4); cb.Label.String='dF/F';

% ONSET-aligned time hist (spikes/cyc % + breath-peak distribution)
% 0 = inspiration ONSET; spikes left axis (spk/cycle %, k), breath-PEAK times
% since the preceding onset on right axis (% of cycles, sky). Single cycle tiled to 1.5 period.
if ~isempty(footAll)
    timeWin = median(diff(footAll))/fps;                 % median breath period (s)
    if ~isfinite(timeWin) || timeWin<=0, timeWin = 1/max(f_pk,eps); end
    edgesTt = linspace(0, timeWin, P.nPhaseBins+1);
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

% avg-projection crop
ax_pr = nexttile(tl,tile_proj);
if have_proj
    imagesc(ax_pr, crop_img); colormap(ax_pr, gray(256)); caxis(ax_pr,[0 1]);
    axis(ax_pr,'image'); set(ax_pr,'XTick',[],'YTick',[]); hold(ax_pr,'on');
    for k = 1:numel(bnd_crop)
        plot(ax_pr, bnd_crop{k}(:,2), bnd_crop{k}(:,1), '-', 'Color',[1 1 0],'LineWidth',P.outlineLW);
    end
    [Hc,~] = size(crop_img); mgn = round(0.05*Hc); thk = max(2,round(0.02*Hc));
    rectangle(ax_pr,'Position',[mgn Hc-mgn-thk barLen_pr thk],'FaceColor','w','EdgeColor','none');
    hold(ax_pr,'off');
else
    axis(ax_pr,'off'); text(ax_pr,0.5,0.5,'(no mask / avg proj)','Horizontal','center');
end

% event-phase polar (spikes/cycle %)
pax2 = polaraxes(tl); pax2.Layout.Tile = tile_epol; hold(pax2,'on');
polarplot(pax2, [ctrsP ctrsP(1)], [pctC pctC(1)], '-','Color',[0.30 0.30 0.30],'LineWidth',1.4);
pax2.RLim=[0 rmaxC]; pax2.ThetaZeroLocation='right'; pax2.ThetaDir='counterclockwise';
pax2.RAxisLocation=180; pax2.FontSize=8; thetaticks(pax2,0:45:315);
title(pax2, sprintf('%.0f\\circ', rad2deg(muW)),'Interpreter','tex');

% coherence polar
cax = polaraxes(tl); cax.Layout.Tile = tile_cpol; hold(cax,'on');
if have_coh && isfinite(th_coh)
    thcc = linspace(0,2*pi,361);
    if isfinite(confC), polarplot(cax, thcc, confC*ones(size(thcc)),'k--','LineWidth',1); end
    polarplot(cax,[th_coh th_coh],[rlo rhi],'-','Color','k','LineWidth',1.8);
    if isfinite(dphi) && dphi>0
        arcc = linspace(th_coh-dphi, th_coh+dphi, 60);
        polarplot(cax, arcc, r_coh*ones(size(arcc)),'-','Color','k','LineWidth',1.8);
    end
    polarplot(cax, th_coh, r_coh,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',6);
end
cax.RLim=[0 1]; cax.ThetaZeroLocation='right'; cax.ThetaDir='counterclockwise';
cax.RAxisLocation=180; cax.FontSize=8; thetaticks(cax,0:45:315);
title(cax, sprintf('%.0f\\circ r=%.2f conf=%.2f', rad2deg(mod(th_coh,2*pi)), r_coh, confC),'Interpreter','tex');

% linear 0..2pi (single cycle), spikes/cycle %
lax2 = nexttile(tl,tile_phist); hold(lax2,'on');
bar(lax2, ctrsP, pctC, 1, 'FaceColor',[0.30 0.30 0.30],'FaceAlpha',0.75,'EdgeColor','none');
xline(lax2, 0,  'Color',onsetCol,'LineWidth',1);    % insp ONSET (red)
xline(lax2, pi, 'Color',peakCol, 'LineWidth',1);    % breath PEAK (sky)
xlim(lax2,[0 2*pi]); ylim(lax2,[0 rmaxC*1.05]);
set(lax2,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
xlabel(lax2,'breath phase (onset=0, peak=\pi)'); ylabel(lax2,'spikes/cycle (%)');
box(lax2,'on'); pbaspect(lax2,[1 1 1]);
title(lax2, sprintf('%.1f cycles, %.2f spk/cycle', nCycles, sum(cntW)/max(nCycles,eps)),'Interpreter','tex');

% spike-triggered average
ax_sp = nexttile(tl,tile_stavg); hold(ax_sp,'on');
if ~isempty(mu_sp)
    fill(ax_sp,[tau fliplr(tau)],[mu_sp+sd_sp fliplr(mu_sp-sd_sp)],lightDff,'EdgeColor','none');
    plot(ax_sp, tau, mu_sp, 'Color',dffColor,'LineWidth',1.5); xline(ax_sp,0,'k--','LineWidth',0.8);
    xlim(ax_sp,[tau(1) tau(end)]);
end
xlabel(ax_sp,'time from spike (s)'); ylabel(ax_sp,'dF/F'); grid(ax_sp,'on'); pbaspect(ax_sp,[1 1 1]);
title(ax_sp, sprintf('spike-triggered average (\\pm SD, n=%d)', size(Esp,1)));
end

%% ===================== LOCAL HELPERS =====================
function add_band(ax, band, f_pk)
    yl = ylim(ax);
    p = patch(ax, [band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.90 0.90 0.90], 'EdgeColor','none');
    xline(ax, f_pk, 'k--', 'LineWidth',1); uistack(p,'bottom'); ylim(ax, yl);
end

function avgimg = read_avgproj_local(folderPath)
    avgimg = [];
    a = dir(fullfile(folderPath,'AVG_*_ch1_preproc_MC_MC.tif'));
    if isempty(a), a = dir(fullfile(folderPath,'AVG_*_MC_MC.tif')); end
    if ~isempty(a)
        V = tiffreadVolume(fullfile(a(1).folder, a(1).name)); avgimg = mean(double(V),3); return;
    end
    s = dir(fullfile(folderPath,'*_ch1_preproc_MC_MC.tif'));
    s = s(~contains({s.name},'AVG','IgnoreCase',true));
    if ~isempty(s)
        V = tiffreadVolume(fullfile(s(1).folder, s(1).name)); avgimg = mean(double(V),3);
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
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
[events, ord] = sort(events); types = types(ord);
keep = true(size(events));
for i = 2:numel(events), if types(i) == types(i-1), keep(i) = false; end, end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases = nan(size(events)); phi_cur = types(1) * pi;
for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
