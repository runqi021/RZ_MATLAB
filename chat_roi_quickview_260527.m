% chat_roi_quickview_260527.m
% -----------------------------------------------------------------------
%  QUICK single-ROI breath/calcium view. Set folderPath + roi, run.
%  One figure, five panels:
%    (1) dF/F (color) + breath (gray) full-trace overlay
%    (2) PSD : breath waveform + dF/F' (derivative), breath peak + band marked
%    (3) Coherence : breath waveform x dF/F, confC line + jackknife CI
%    (4) Breath-triggered dF/F heatmap, rows sorted (see sortMode)
%    (5) Breath-triggered average dF/F +/- SD
%
%  Inputs (in folderPath):  *_ch1_dFF.mat (dFF),  *DLC*breath_peak_data.mat
%  Alignment: toss first nDrop breath frames, truncate to common length
%             (breath cam is 2P-triggered, fps = imaging).
%
%  sortMode for the heatmap (panel 4):
%    'postmean' (default) : mean dF/F in the post-inspiration window (tau>=0),
%                           descending  -- strongest responders on top.
%    'dt'                 : inter-breath interval (time to NEXT onset),
%                           ascending  -- "dt nearest" ordering.
%    'none'               : chronological (breath order).
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
folderPath = 'D:\Ventral_surface_summary\ChAT\0521\cell1\roi5_7x_x-1200y200z-30_3000f_23lp_00001';
roi        = 1;            % ROI index (== dFF column)

%Vglut2
% D:\Ventral_surface_summary\Vglut2\1124\cell1\roi5_1400-1230-0_x4.4_15lp_6000f_00001  roi=5 

%ChAT
% D:\Ventral_surface_summary\ChAT\0521\cell1\roi5_7x_x-1200y200z-30_3000f_23lp_00001  roi=1
% D:\Ventral_surface_summary\ChAT\0523\cell1\roi1_4x_x-900y700z-15_6000f_13lp_00001   roi=1

nDrop      = 30;           % breath frames to toss (match calcium)
fallback_fps = 30;
TW_spec    = 6;            % multitaper TW for PSD + coherence
alpha_sig  = 0.01;

f_breath_search = [0.2 4]; % Hz, search band for breath PSD peak
fwhm_factor = 0.6;         % coherence/detection band = fwhm_factor x FWHM
min_bw      = 0.05;        % Hz, min band width
fmin        = 0.05;        % Hz, PSD/coherence lower bound
fmax        = 15;          % Hz, upper bound

trace_xlim  = [110 170];          % panel-1 window (s): [] = full trace; [t0 t1] or 1:10 (uses min..max)
sortMode    = 'dt';  % 'postmean' | 'dt' | 'none'   (see header)
pct_ylim_spk = [0 10];         % y-limit (%) for SPIKE histograms (panel 8 left + panel 9); [] = auto
pct_ylim_on  = [0 20];         % y-limit (%) for INSP-ONSET histogram (panel 8 right);    [] = auto

dffColor    = [0.85 0.10 0.10];   % single accent color used in ALL dF/F panels
doSave      = true;
% =========================================================

lightDff = 0.30*dffColor + 0.70;  % lightened version for CI/SD shaded fills
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
[~, recName] = fileparts(folderPath);

%% ===================== LOAD + ALIGN =====================
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp = dir(fullfile(folderPath,'*DLC*breath_peak_data.mat'));
ip = dir(fullfile(folderPath,'*DLC*breath_insp_start_data.mat'));
if isempty(ip), ip = dir(fullfile(folderPath,'*breath_insp_start_data.mat')); end
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No *DLC*breath_peak_data.mat in %s', folderPath);

fps = detect_session_fps(folderPath, fallback_fps);
D   = load(fullfile(df(1).folder, df(1).name),'dFF');
BP  = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(D.dFF);
assert(roi>=1 && roi<=size(dff_all,2),'ROI %d out of range (1..%d)',roi,size(dff_all,2));

bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
nB = numel(BP.breath);
if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
    ev = double(BP.insp_onsets_train(:) ~= 0);
else
    ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
end
ev(1:min(nDrop,numel(ev))) = [];

% --- insp-start (foot) train, same toss alignment ---
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    ev_foot = zeros(nB,1); fi = round(IP.insp_start_idx(:));
    ev_foot(fi(fi>=1 & fi<=nB)) = 1;
    ev_foot(1:min(nDrop,numel(ev_foot))) = [];
else
    warning('No *breath_insp_start_data.mat -- row 3 panels will be skipped.');
    ev_foot = [];
end

T   = min([size(dff_all,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
dff = dff_all(1:T, roi); bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); end
t   = (0:T-1)'/fps;
fprintf('%s ROI%d | T=%d @%.3g Hz | %d breaths | %d insp-starts\n', ...
        recName, roi, T, fps, sum(ev), sum(ev_foot));

%% ===================== SPECTRA + BAND =====================
pB.Fs=fps; pB.tapers=[TW_spec,2*TW_spec-1]; pB.pad=0; pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=[2,alpha_sig];
[Sbw,fbw,SbwErr] = mtspectrumc(bw, pB);            Sbw=Sbw(:); fbw=fbw(:);
[Sdd,fdd,SddErr] = mtspectrumc(diff(dff)*fps, pB); Sdd=Sdd(:); fdd=fdd(:);

mm=fbw>=f_breath_search(1) & fbw<=f_breath_search(2);
[~,rl]=max(Sbw(mm)); ip=find(mm,1)+rl-1; f_pk=fbw(ip);
hh=Sbw(ip)/2; lo=ip; while lo>1&&Sbw(lo)>hh, lo=lo-1; end
hi=ip;        while hi<numel(fbw)&&Sbw(hi)>hh, hi=hi+1; end
f_fwhm=[max(fbw(lo),f_breath_search(1)), min(fbw(hi),f_breath_search(2))];
bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
band=[max(f_pk-bwd/2,fmin), min(f_pk+bwd/2,fmax)];

%% ===================== COHERENCE (breath waveform x dF/F) =====================
pc.Fs=fps; pc.tapers=[TW_spec,2*TW_spec-1]; pc.pad=0; pc.fpass=[fmin,min(fmax,fps/2)]; pc.err=[2,alpha_sig];
[~,Cw,~,~,~,~,fcw,confCw,~,Cerrw] = coherencyc(bw, dff-mean(dff), pc);
fcw=fcw(:); Cw=Cw(:);
mb = fcw>=band(1) & fcw<=band(2); Cband = mean(Cw(mb));

%% ===================== PEAK-TRIGGERED  +  ONSET-TRIGGERED dF/F =====================
% Each heatmap is INDEPENDENT:
%   peak-triggered  -> triggered on every breath peak,  sorted by dt to next peak
%   onset-triggered -> triggered on every insp onset,   sorted by dt to next onset
win = round(fps / f_pk);                 % +/- 1 breath period
tau = (-win:win)/fps;
on      = find(ev>0);        on      = on(on-win>=1 & on+win<=numel(dff));
foot_on = find(ev_foot>0);   if ~isempty(foot_on), foot_on = foot_on(foot_on-win>=1 & foot_on+win<=numel(dff)); end

% --- peak-triggered ---
E = zeros(numel(on), 2*win+1);
for k = 1:numel(on), E(k,:) = dff(on(k)-win : on(k)+win); end
% signed dt from each peak to its NEAREST other peak (prev or next, whichever closer)
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
% dot overlay: nearest insp-start relative to each peak (signed)
dt_pk_to_foot = nan(numel(on),1);
if ~isempty(foot_on)
    for k = 1:numel(on)
        [~, mi] = min(abs(foot_on - on(k)));
        dt_pk_to_foot(k) = (foot_on(mi) - on(k)) / fps;
    end
end

% --- onset-triggered (independent set of rows, own sort) ---
Esf = []; clf_ = cl; mu_f = []; sd_f = []; si_f = []; dt_foot_to_peak = [];
if ~isempty(foot_on)
    Ef = zeros(numel(foot_on), 2*win+1);
    for k = 1:numel(foot_on), Ef(k,:) = dff(foot_on(k)-win : foot_on(k)+win); end
    dt_to_nearest_foot = nan(numel(foot_on),1);
    for k = 1:numel(foot_on)
        others = foot_on; others(k) = [];
        if isempty(others), continue; end
        [~, mi] = min(abs(others - foot_on(k)));
        dt_to_nearest_foot(k) = (others(mi) - foot_on(k)) / fps;
    end
    switch lower(sortMode)
        case 'postmean'
            key_f = mean(Ef(:, tau>=0), 2); [~,si_f] = sort(key_f,'descend');
        case 'dt'
            [~,si_f] = sort(dt_to_nearest_foot,'ascend','MissingPlacement','last');
        otherwise
            si_f = (1:numel(foot_on))';
    end
    Esf = Ef(si_f,:);
    clf_ = prctile(Esf(:), [5 99.5]);
    mu_f = mean(Ef,1); sd_f = std(Ef,0,1);
    % dot overlay: nearest peak relative to each onset (signed)
    dt_foot_to_peak = nan(numel(foot_on),1);
    for k = 1:numel(foot_on)
        [~, mi] = min(abs(on - foot_on(k)));
        dt_foot_to_peak(k) = (on(mi) - foot_on(k)) / fps;
    end
end
n_cyc = numel(foot_on);

%% ===================== SPIKE HISTOGRAMS AROUND BREATH PEAK (row 3) =====================
% Time histogram: for each peak, pool (spike-peak) AND (onset-peak) within +/- win.
% Phase histogram: piecewise phi(t) with FOOT=0, PEAK=pi; histogram spike phases.
have_hist = false;
sp_file = fullfile(folderPath,'ca_spike_data.mat');
if isfile(sp_file)
    CAh = load(sp_file,'roi_spikes');
    if isfield(CAh,'roi_spikes') && roi <= numel(CAh.roi_spikes)
        spk_train_h = double(CAh.roi_spikes(roi).spike_train(:));
        if numel(spk_train_h) < numel(dff), spk_train_h(end+1:numel(dff)) = 0; end
        spk_train_h = spk_train_h(1:numel(dff));
        sp_on_h = find(spk_train_h > 0);

        % per-peak relative spike times and onset times
        peri_t_spk    = [];
        peri_t_onset  = [];
        for k = 1:numel(on)
            p = on(k);
            spk_w = sp_on_h(sp_on_h >= p-win & sp_on_h <= p+win);
            peri_t_spk = [peri_t_spk; (spk_w - p) / fps];                 %#ok<AGROW>
            if ~isempty(foot_on)
                ft_w = foot_on(foot_on >= p-win & foot_on <= p+win);
                peri_t_onset = [peri_t_onset; (ft_w - p) / fps];          %#ok<AGROW>
            end
        end

        % piecewise phi(t): foot=0, peak=pi
        phi_h = piecewise_phase_local(on, foot_on, numel(dff));
        sp_phase_h = phi_h(sp_on_h);
        sp_phase_h = mod(sp_phase_h(~isnan(sp_phase_h)), 2*pi);

        have_hist = true;
    end
end

%% ===================== CA-SPIKE-TRIGGERED dF/F (legacy - kept for clarity, unused) =====================
sp_file = fullfile(folderPath,'ca_spike_data.mat');
have_sp = false;
if isfile(sp_file)
    CA = load(sp_file,'roi_spikes');
    if isfield(CA,'roi_spikes') && roi <= numel(CA.roi_spikes)
        spk_train = double(CA.roi_spikes(roi).spike_train(:));
        if numel(spk_train) < numel(dff), spk_train(end+1:numel(dff)) = 0; end
        spk_train = spk_train(1:numel(dff));
        sp_on = find(spk_train > 0);
        sp_on = sp_on(sp_on-win>=1 & sp_on+win<=numel(dff));
        if ~isempty(sp_on)
            Esp = zeros(numel(sp_on), 2*win+1);
            for k=1:numel(sp_on), Esp(k,:) = dff(sp_on(k)-win : sp_on(k)+win); end

            % dt to NEAREST insp-start and NEAREST breath peak per spike (s)
            dt_sp_to_foot = nan(numel(sp_on),1);
            dt_sp_to_peak = nan(numel(sp_on),1);
            if ~isempty(foot_on)
                for k = 1:numel(sp_on)
                    [~, mi] = min(abs(foot_on - sp_on(k)));
                    dt_sp_to_foot(k) = (foot_on(mi) - sp_on(k)) / fps;
                end
            end
            for k = 1:numel(sp_on)
                [~, mi] = min(abs(on - sp_on(k)));
                dt_sp_to_peak(k) = (on(mi) - sp_on(k)) / fps;
            end

            % sort by dt to nearest insp-start
            if ~isempty(foot_on)
                [~,si_sp] = sort(dt_sp_to_foot,'ascend','MissingPlacement','last');
            else
                [~,si_sp] = sort(dt_sp_to_peak,'ascend','MissingPlacement','last');
            end
            Esp_s = Esp(si_sp,:);
            cl_sp = prctile(Esp_s(:), [5 99.5]);
            have_sp = true;
        end
    end
end

%% (insp-onset-triggered E/Ef + sort were built jointly above)


%% ===================== FIGURE =====================
fig = figure('Color','w','Name',sprintf('%s ROI%d quickview',recName,roi), ...
             'Units','normalized','Position',[0.04 0.04 0.9 0.94]);
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s   ROI%d   |   fps=%.2f, breath peak %.2f Hz, band [%.2f %.2f] Hz, coh-in-band=%.2f', ...
      recName, roi, fps, f_pk, band(1), band(2), Cband), 'Interpreter','none','FontWeight','bold');

% (1) trace overlay (row 1, full width)
ax1 = nexttile(tl,1,[1 3]);
if isempty(trace_xlim), w=[t(1) t(end)]; else, w=[min(trace_xlim) max(trace_xlim)]; end
mw = t>=w(1) & t<=w(2);
yyaxis(ax1,'right');
plot(ax1, t(mw), bw(mw), '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
set(ax1,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax1,'breath');
yyaxis(ax1,'left');
plot(ax1, t(mw), dff(mw), '-','Color',dffColor,'LineWidth',0.8);
set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
xlim(ax1,w); xlabel(ax1,'Time (s)'); box(ax1,'off');
title(ax1,'dF/F (red) + breath (gray)');

% (2) PSD
ax2 = nexttile(tl,5); hold(ax2,'on');
fill(ax2,[fbw;flipud(fbw)],10*log10([SbwErr(1,:)';flipud(SbwErr(2,:)')]),[0.85 0.85 0.85],'EdgeColor','none');
fill(ax2,[fdd;flipud(fdd)],10*log10([SddErr(1,:)';flipud(SddErr(2,:)')]),lightDff,'EdgeColor','none');
plot(ax2, fbw,10*log10(Sbw),'Color',[0.5 0.5 0.5],'LineWidth',1.0);
plot(ax2, fdd,10*log10(Sdd),'Color',dffColor,'LineWidth',1.1);
set(ax2,'XScale','log'); xlim(ax2,[fmin fmax]); xticks(ax2,[0.1 0.3 1 3 10]); set(ax2,'XMinorTick','off');
add_band(ax2, band, f_pk); xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'power (dB)');
pbaspect(ax2,[1 1 1]);
title(ax2,'PSD: breath (gray) + dF/F'' (color)');

% (3) coherence
ax3 = nexttile(tl,6); hold(ax3,'on');
fill(ax3,[fcw;flipud(fcw)],[Cerrw(1,:)';flipud(Cerrw(2,:)')],lightDff,'EdgeColor','none');
plot(ax3, fcw, Cw, 'Color',dffColor,'LineWidth',1.1);
yline(ax3, confCw,'k--','LineWidth',1);
set(ax3,'XScale','log'); xlim(ax3,[fmin fmax]); ylim(ax3,[0 1]); xticks(ax3,[0.1 0.3 1 3 10]); set(ax3,'XMinorTick','off');
add_band(ax3, band, f_pk); xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'coherence');
pbaspect(ax3,[1 1 1]);
title(ax3, sprintf('breath x dF/F coherence (confC=%.2f)',confCw));

% (4) peak-triggered heatmap (sky-blue dots = ALL insp-starts in row's window)
skyCol  = [0.35 0.75 1.00];
pinkCol = [1.00 0.40 0.70];
ax4 = nexttile(tl,7);
imagesc(ax4, tau, 1:size(Es,1), Es); axis(ax4,'tight');
colormap(ax4, flipud(gray(256))); caxis(ax4, cl);
hold(ax4,'on');
% per row (sorted order si), pool every foot within +/- win of that peak
if ~isempty(foot_on)
    on_sorted = on(si);
    foot_x = []; foot_y = [];
    for k = 1:numel(on_sorted)
        p  = on_sorted(k);
        fw = foot_on(foot_on >= p-win & foot_on <= p+win);
        foot_x = [foot_x; (fw - p) / fps];        %#ok<AGROW>
        foot_y = [foot_y; repmat(k, numel(fw),1)]; %#ok<AGROW>
    end
    plot(ax4, foot_x, foot_y, '.', 'Color', skyCol, 'MarkerSize', 8);
end
hold(ax4,'off');
set(ax4,'YDir','reverse'); xlabel(ax4,'time from breath peak (s)'); ylabel(ax4,'breath #');
pbaspect(ax4,[1 1 1]);
cb=colorbar(ax4); cb.Label.String='dF/F';
title(ax4, sprintf('peak-triggered dF/F (n=%d, %s)', numel(on), sortLbl),'Interpreter','none');

% (5) peak-triggered average
ax5 = nexttile(tl,4); hold(ax5,'on');
fill(ax5,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],lightDff,'EdgeColor','none');
plot(ax5, tau, mu, 'Color',dffColor,'LineWidth',1.5);
xline(ax5,0,'k--','LineWidth',0.8);
xlim(ax5,[tau(1) tau(end)]); xlabel(ax5,'time from breath peak (s)'); ylabel(ax5,'dF/F'); grid(ax5,'on');
pbaspect(ax5,[1 1 1]);
title(ax5,'peak-triggered average (\pm SD)');


% ---- ROW 3 : SPIKE HISTOGRAMS (time around peak + phase) ----
if have_hist
    skyCol  = [0.35 0.75 1.00];
    pinkCol = [1.00 0.40 0.70];

    binT_sec = 0.033;   % 33 ms/bin (hardcoded)
    nHistP   = 24;
    edgesT   = -win/fps : binT_sec : win/fps;
    if edgesT(end) < win/fps, edgesT(end+1) = edgesT(end) + binT_sec; end
    ctrsT    = (edgesT(1:end-1)+edgesT(2:end))/2;
    edgesP = linspace(0, 2*pi, nHistP+1);
    ctrsP  = (edgesP(1:end-1)+edgesP(2:end))/2;
    ctrsP_tile = [ctrsP, ctrsP+2*pi];

    % time hist around breath peak: spikes (pink) + insp onsets (sky)
    ax_h1 = nexttile(tl,8); hold(ax_h1,'on');
    cnt_sp = histcounts(peri_t_spk,   edgesT);
    cnt_on = histcounts(peri_t_onset, edgesT);
    pct_sp = 100*cnt_sp / max(sum(cnt_sp),1);   % % of all peri-peak spikes
    pct_on = 100*cnt_on / max(sum(cnt_on),1);   % % of all peri-peak onsets
    yyaxis(ax_h1,'left');
    bar(ax_h1, ctrsT, pct_sp, 1, 'FaceColor', 'k', 'FaceAlpha',0.85, 'EdgeColor','none');
    ylabel(ax_h1,'spikes (%)'); set(ax_h1,'YColor',pinkCol);
    yyaxis(ax_h1,'right');
    bar(ax_h1, ctrsT, pct_on, 1, 'FaceColor', skyCol, 'FaceAlpha',0.55, 'EdgeColor','none');
    ylabel(ax_h1,'insp-onsets (%)'); set(ax_h1,'YColor',skyCol);
    if ~isempty(pct_ylim_on),  ylim(ax_h1, pct_ylim_on);  end   % right axis = onsets
    yyaxis(ax_h1,'left');
    if ~isempty(pct_ylim_spk), ylim(ax_h1, pct_ylim_spk); end   % left axis  = spikes
    xline(ax_h1, 0, 'Color', pinkCol, 'LineWidth', 1);
    xlim(ax_h1, [edgesT(1) edgesT(end)]);
    xlabel(ax_h1,'time from breath peak (s)');
    title(ax_h1, sprintf('peri-peak hist: spikes (pink, n=%d) + onsets (sky, n=%d)', numel(peri_t_spk), numel(peri_t_onset)),'Interpreter','none');
    pbaspect(ax_h1,[1 1 1]); box(ax_h1,'on');

    % phase hist: foot=0, peak=pi, tiled to [0, 4pi]
    ax_h2 = nexttile(tl,9); hold(ax_h2,'on');
    cnt_ph = histcounts(sp_phase_h, edgesP);
    pct_ph = 100*cnt_ph / max(sum(cnt_ph),1);   % % of all spikes per phase bin
    bar(ax_h2, ctrsP_tile, [pct_ph pct_ph], 1, 'FaceColor', 'k', 'FaceAlpha',0.85, 'EdgeColor','none');
    xline(ax_h2, 0,    'Color', skyCol, 'LineWidth', 1.0);   % onset @ 0
    xline(ax_h2, pi,   'Color', pinkCol,  'LineWidth', 1);            % peak  @ pi
    xline(ax_h2, 2*pi, 'Color', skyCol, 'LineWidth', 1.0);   % next onset
    xline(ax_h2, 3*pi, 'Color', pinkCol,  'LineWidth', 1);
    xlim(ax_h2, [0 4*pi]);
    if ~isempty(pct_ylim_spk), ylim(ax_h2, pct_ylim_spk); end
    set(ax_h2,'XTick',[0 pi 2*pi 3*pi 4*pi], 'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(ax_h2,'phase (onset = 0, peak = \pi)'); ylabel(ax_h2,'spikes (%)');
    title(ax_h2, sprintf('spike phase hist (n=%d)', numel(sp_phase_h)),'Interpreter','none');
    pbaspect(ax_h2,[1 1 1]); box(ax_h2,'on');
end

%% ===================== SAVE =====================
if doSave
    base = fullfile(folderPath, sprintf('quickview_ROI%02d', roi));
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

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase reference: FEET at 0/2pi/..., PEAKS at pi/3pi/...
% Linear ramps in time between consecutive events. NaN outside the range.
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
