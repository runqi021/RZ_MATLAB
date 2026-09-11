% Ventral_surface_temporal_phase_window_gui_svd_260607.m
% -----------------------------------------------------------------------
%  Interactive WINDOW-PICK version of Ventral_surface_temporal_phase_sigROI.
%
%  Produces the SAME 12-panel temporal-phase figure, but panels 1 (dF/F +
%  breath trace) and 2 (phase sawtooth) are clipped to a time window YOU pick
%  with a small GUI: a trace preview opens, you zoom/pan with the toolbar to
%  frame the window, then click "Use this window". Everything else in the
%  figure is unchanged (full-recording spectra / coherence / phase hist / etc.).
%
%  Pick the ROI by significant-ROI number (1..N) from the coherence data, or
%  set useSig=false and give folderPath/roi/group/recDate manually.
%
%  Conventions identical to the svd_260606 set (nDrop=30, ca_lag=0.1 s,
%  Vglut2/1124 +1-frame fix, detect_session_fps).
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
outDir   = fullfile(rootPath, 'representative_trace_nonIO');   % rendered figures saved here

selMode  = 'nonIO';         % 'nonIO' (default) | 'sig' | 'all' | 'manual'
startNo  = 1;               % which ROI in the list to open first (1..N)

% --- manual ROI (only used when selMode='manual') ---
folderPath_manual = 'D:\Ventral_surface_summary\Vglut2\0224\...';
roi_manual = 9; group_manual = 'Vglut2'; recDate_manual = '0224'; recName_manual = '';

% --- extra ROIs ALWAYS appended to the list (any selMode) ---
% each row: {folderPath, roi, group, recDate, recName}
% NOTE: group='Vglut2' & recDate='1124' triggers the +1-frame rising-edge fix.
extraROIs = {
  'D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing\roi5_1400-1300-0_18lp_930_x4_512x512_3000f_00001', 13, 'Vglut2', '1124', 'roi5_1400-1300-0_18lp_930_x4_512x512_3000f_00001'
};

doSave   = true;            % save the rendered figure (window tagged)

P = struct();
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.TW_coh       = 4;           P.alpha_coh    = 0.001;
P.minSpikes    = 2;           P.ca_lag_sec   = 0.1;
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
P.nPhaseBins   = 30;
P.trigWin_sec  = [];          % WAVEFORM-average HALF-window (s); [] = auto (1 breath period)
P.histWin_sec  = [];          % spike/peak HISTOGRAM HALF-window (s); [] = auto (1 breath period)
P.histBinFrames= 2;           % spike-hist bin width in FRAMES (1 = per-frame)
P.nShuffle     = 200;         % circular-shift shuffles for the spike-hist null (0 = off)
P.shiftMinCyc  = 3;           % min circular shift = this many breath periods (also >= window, >= 1 s)
P.pad_um       = 20;          P.clip_pct     = [0.5 99.9];
P.gamma_val    = 0.6;         P.PixelSizeBase = 1.7778;  P.outlineLW = 1.3;
P.sortMode     = 'dt';        P.trace_xlim_sp = [];     % set by the GUI
P.dffColor     = [0.2 0.7 0.2];
P.onsetCol     = [0.90 0.10 0.10];
P.peakCol      = [0.35 0.75 1.00];
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ===================== BUILD ROI LIST =====================
if strcmpi(selMode,'manual')
    rn = recName_manual; if isempty(rn), [~,rn] = fileparts(folderPath_manual); end
    items = struct('group',group_manual,'recDate',recDate_manual,'recName',rn, ...
        'roi',roi_manual,'r',NaN,'folderPath',folderPath_manual,'label','manual');
else
    assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
    S = load(cohData, 'PP','confC','labels','groups');
    ioGroup = find(strcmp(S.groups,'IO'));
    switch lower(selMode)
        case 'sig',   sel = find(S.PP.r >= S.confC);
        case 'nonio', sel = find(S.PP.colorIdx ~= ioGroup);
        case 'all',   sel = (1:numel(S.labels))';
        otherwise,    error('selMode must be nonIO | sig | all | manual');
    end
    assert(~isempty(sel), 'no ROIs selected for selMode=%s', selMode);
    items = struct('group',{},'recDate',{},'recName',{},'roi',{},'r',{},'folderPath',{},'label',{});
    for jj = 1:numel(sel)
        lab = S.labels{sel(jj)}; parts = regexp(lab,'/','split');
        if numel(parts)<4, continue; end
        items(end+1) = struct('group',parts{1},'recDate',parts{2}, ...
            'recName',strjoin(parts(3:end-1),'/'),'roi',str2double(parts{end}), ...
            'r',S.PP.r(sel(jj)),'folderPath','','label',lab); %#ok<SAGROW>
    end
    fprintf('selMode=%s -> %d ROIs\n', selMode, numel(items));
end

% append extra manual ROIs (any selMode)
for jj = 1:size(extraROIs,1)
    items(end+1) = struct('group',extraROIs{jj,3}, 'recDate',extraROIs{jj,4}, ...
        'recName',extraROIs{jj,5}, 'roi',extraROIs{jj,2}, 'r',NaN, ...
        'folderPath',extraROIs{jj,1}, ...
        'label',sprintf('extra:%s ROI%d', extraROIs{jj,3}, extraROIs{jj,2})); %#ok<SAGROW>
end
if ~isempty(extraROIs)
    fprintf('appended %d extra ROI(s); total %d\n', size(extraROIs,1), numel(items));
end

startIdx = min(max(startNo,1), numel(items));

%% ===================== LAUNCH GUI =====================
if doSave && ~isfolder(outDir), mkdir(outDir); end
launch_window_gui(items, startIdx, rootPath, P, doSave, outDir);
% Prev/Next switch ROI; type a window and click Render. Nothing renders until Render.

% =======================================================================
% =========================== LOCAL FUNCTIONS ===========================
% =======================================================================

function folderPath = resolve_folder(rootPath, group, recName)
    folderPath = '';
    d = dir(fullfile(rootPath, group, '**', recName)); d = d([d.isdir]);
    if isempty(d), d = dir(fullfile(rootPath, '**', recName)); d = d([d.isdir]); end
    if ~isempty(d), folderPath = fullfile(d(1).folder, d(1).name); end
end

% -----------------------------------------------------------------------
function [t, dff, bw, on_t] = load_trace(folderPath, roi, group, recDate, P)
% Minimal load/align for the preview trace (matches panel-1 alignment).
% on_t = inspiration-onset times (s) for counting breaths in a window.
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', folderPath);
fps = detect_session_fps(folderPath, P.fallback_fps);
D  = load(fullfile(df(1).folder, df(1).name),'dFF'); dff_all = double(D.dFF);
BP = load(fullfile(bp(1).folder, bp(1).name));
assert(roi>=1 && roi<=size(dff_all,2),'ROI %d out of range (1..%d)',roi,size(dff_all,2));
bw = detrend(double(BP.breath(:))); bw(1:min(P.nDrop,numel(bw))) = []; bw = bw - mean(bw);
isV1124 = strcmpi(group,'Vglut2') && strcmp(recDate,'1124');
if isV1124, bw = [bw(1); bw(1:end-1)]; end
T = min(size(dff_all,1), numel(bw));
dff = dff_all(1:T, roi); bw = bw(1:T); t = (0:T-1)'/fps;

on_t = [];
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    oi = round(IP.insp_start_idx(:)) - P.nDrop;     % to post-toss base
    if isV1124, oi = oi + 1; end                    % rising-edge: breath leads 1 frame
    oi = oi(oi>=1 & oi<=T);
    on_t = (oi-1)/fps;     % frame k -> t=(k-1)/fps, matches the preview trace
end
end

% -----------------------------------------------------------------------
function launch_window_gui(items, startIdx, rootPath, P, doSave, outDir)
% Persistent control panel over a LIST of ROIs. Prev/Next switch ROI; type a
% window (e.g. "0 30"), Apply to preview, Render to build the full figure.
% Nothing renders until Render is clicked.  Rendered figures saved into outDir.
curIdx = startIdx;
t = []; dff = []; bw = []; on_t = []; curFolder = '';   % per-ROI state (set by load_current)

f = figure('Color','w','Name','Temporal-phase window GUI', ...
           'Units','normalized','Position',[0.06 0.36 0.88 0.52]);
ax = axes('Parent',f,'Units','normalized','Position',[0.06 0.34 0.90 0.56]);

uicontrol(f,'Style','text','Units','normalized','Position',[0.06 0.135 0.40 0.06], ...
    'String','window (s):  start end   e.g.  0 30', ...
    'HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
uicontrol(f,'Style','text','Units','normalized','Position',[0.47 0.135 0.12 0.055], ...
    'String','wave win (s):', 'HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hTrig = uicontrol(f,'Style','edit','Units','normalized','Position',[0.595 0.135 0.07 0.06], ...
    'String','','FontSize',11);    % waveform-average half-window (blank = auto)
uicontrol(f,'Style','text','Units','normalized','Position',[0.69 0.135 0.12 0.055], ...
    'String','hist win (s):', 'HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hHist = uicontrol(f,'Style','edit','Units','normalized','Position',[0.815 0.135 0.07 0.06], ...
    'String','','FontSize',11);    % spike/peak histogram half-window (blank = auto)
hEdit = uicontrol(f,'Style','edit','Units','normalized','Position',[0.06 0.04 0.16 0.085], ...
    'String','0 30','FontSize',11,'Callback',@apply);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.23 0.04 0.075 0.085], ...
    'String','Apply','Callback',@apply);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.31 0.04 0.065 0.085], ...
    'String','Full','Callback',@(s,e) reset_full());
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.40 0.04 0.14 0.085], ...
    'String','Render figure','FontWeight','bold','Callback',@render);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.56 0.04 0.085 0.085], ...
    'String','<< Prev','Callback',@(s,e) step(-1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.65 0.04 0.085 0.085], ...
    'String','Next >>','Callback',@(s,e) step(1));
hStat = uicontrol(f,'Style','text','Units','normalized','Position',[0.75 0.04 0.22 0.085], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);

load_current();

    function load_current()
        it = items(curIdx);
        fp = it.folderPath;
        if isempty(fp), fp = resolve_folder(rootPath, it.group, it.recName); end
        if isempty(fp)
            set(hStat,'String', sprintf('folder not found: %s', it.label)); return;
        end
        curFolder = fp;
        try
            [t, dff, bw, on_t] = load_trace(fp, it.roi, it.group, it.recDate, P);
        catch ME
            set(hStat,'String', sprintf('load error: %s', ME.message)); return;
        end
        yyaxis(ax,'left');  cla(ax);
        yyaxis(ax,'right'); cla(ax);
        yyaxis(ax,'right'); plot(ax, t, bw, '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
        set(ax,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax,'breath');
        yyaxis(ax,'left');  plot(ax, t, dff, '-','Color',P.dffColor,'LineWidth',0.8);
        set(ax,'YColor','k'); ylabel(ax,'\DeltaF/F');
        xlim(ax,[t(1) t(end)]); xlabel(ax,'Time (s)'); box(ax,'off');
        if isnan(it.r), rstr = ''; else, rstr = sprintf('  r=%.3f', it.r); end
        title(ax, sprintf('[%d/%d]  %s/%s/%s  ROI%d%s', curIdx, numel(items), ...
              it.group, it.recDate, it.recName, it.roi, rstr),'Interpreter','none');
        apply([],[]);    % re-apply the current window box to the new trace
    end

    function step(d)
        curIdx = curIdx + d;
        if curIdx < 1,            curIdx = numel(items); end
        if curIdx > numel(items), curIdx = 1;            end
        load_current();
    end

    function apply(~,~)
        if isempty(t), return; end
        v = sscanf(strrep(get(hEdit,'String'), ',', ' '), '%f');
        if numel(v) >= 2
            a = max(min(v(1),v(2)), t(1));
            b = min(max(v(1),v(2)), t(end));
            if b > a
                xlim(ax, [a b]);
                mwsp = t>=a & t<=b;
                if nnz(mwsp) > 2
                    yyaxis(ax,'left');  ylim(ax, padlim(dff(mwsp)));
                    yyaxis(ax,'right'); ylim(ax, padlim(bw(mwsp)));
                    yyaxis(ax,'left');
                end
                nB = nnz(on_t>=a & on_t<=b);
                set(hStat,'String', sprintf('window %.2f-%.2f s  |  %d breaths', a, b, nB));
            else
                set(hStat,'String','end must be > start');
            end
        else
            set(hStat,'String','enter two numbers, e.g. 0 30');
        end
    end

    function reset_full()
        if isempty(t), return; end
        set(hEdit,'String', sprintf('%.0f %.0f', t(1), t(end)));
        xlim(ax,[t(1) t(end)]);
        yyaxis(ax,'left');  ylim(ax, padlim(dff));
        yyaxis(ax,'right'); ylim(ax, padlim(bw));
        yyaxis(ax,'left'); set(hStat,'String','full trace');
    end

    function render(~,~)
        if isempty(curFolder), set(hStat,'String','no ROI loaded'); return; end
        v = sscanf(strrep(get(hEdit,'String'), ',', ' '), '%f');
        if numel(v) >= 2
            wsp = [max(min(v(1),v(2)),t(1)), min(max(v(1),v(2)),t(end))];
        else
            wsp = [t(1) t(end)];
        end
        if wsp(2) <= wsp(1), wsp = [t(1) t(end)]; end
        it = items(curIdx);
        Pr = P; Pr.trace_xlim_sp = wsp;
        tv = sscanf(strrep(get(hTrig,'String'), ',', ' '), '%f');   % waveform-average half-window
        if ~isempty(tv) && isfinite(tv(1)) && tv(1) > 0
            Pr.trigWin_sec = tv(1);
        else
            Pr.trigWin_sec = [];
        end
        hv = sscanf(strrep(get(hHist,'String'), ',', ' '), '%f');   % spike/peak histogram half-window
        if ~isempty(hv) && isfinite(hv(1)) && hv(1) > 0
            Pr.histWin_sec = hv(1);
        else
            Pr.histWin_sec = [];
        end
        set(hStat,'String','rendering...'); drawnow;
        [figR, proj] = make_temporal_phase_fig(curFolder, it.roi, it.recName, it.group, it.recDate, Pr);
        if doSave
            tag  = sprintf('win%03d-%03d', round(wsp(1)), round(wsp(2)));
            stem = regexprep(sprintf('%s_%s_%s_ROI%02d', it.group, it.recDate, it.recName, it.roi), '[\\/:*?"<>|]', '_');
            base = fullfile(outDir, sprintf('%s_%s', stem, tag));
            exportgraphics(figR, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
            exportgraphics(figR, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            if proj.have_proj
                save_avgproj_png(proj, P, fullfile(outDir, sprintf('%s_avgproj.png', stem)));
            end
            set(hStat,'String', sprintf('saved -> %s.png (+ avgproj)', stem));
        else
            set(hStat,'String','rendered');
        end
    end

    function yl = padlim(y)
        lo = min(y); hi = max(y); pad = 0.05*max(hi-lo, eps); yl = [lo-pad hi+pad];
    end
end

% -----------------------------------------------------------------------
function save_avgproj_png(proj, P, outpng)
% Standalone, full-resolution PNG of the ROI avg-proj crop (outline + scalebar).
[H,W] = size(proj.crop_img);
fg = figure('Color','w','Visible','off','Units','pixels','Position',[100 100 max(W,60) max(H,60)]);
axp = axes('Parent',fg,'Position',[0 0 1 1]);
imagesc(axp, proj.crop_img); colormap(axp, gray(256)); caxis(axp,[0 1]);
axis(axp,'image','off'); hold(axp,'on');
for k = 1:numel(proj.bnd_crop)
    plot(axp, proj.bnd_crop{k}(:,2), proj.bnd_crop{k}(:,1), '-','Color',[1 1 0],'LineWidth',P.outlineLW);
end
mgn = round(0.05*H); thk = max(2, round(0.02*H));
rectangle(axp,'Position',[mgn H-mgn-thk proj.barLen_pr thk],'FaceColor','w','EdgeColor','none');
hold(axp,'off');
exportgraphics(axp, outpng, 'Resolution', 600, 'BackgroundColor','white');
close(fg);
end

% -----------------------------------------------------------------------
function [fig, proj] = make_temporal_phase_fig(folderPath, roi, recName, group, recDate, P)
% Port of Ventral_surface_temporal_phase_sigROI's figure; ONLY change vs that
% script: panels 1 & 2 are clipped to P.trace_xlim_sp (the GUI window).
% Second output `proj` carries the avg-proj crop so it can be saved separately
% at full resolution (the embedded panel gets rasterized in the combined figure).
lightDff = 0.30*P.dffColor + 0.70;
dffColor = P.dffColor; onsetCol = P.onsetCol; peakCol = P.peakCol;

%% ---- LOAD + ALIGN ----
df  = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp  = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip  = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
sam = dir(fullfile(folderPath,'*_cpSAM_output.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', folderPath);

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

% ---- WINDOW (panels 1 & 2) ----
if isempty(P.trace_xlim_sp)
    wsp = [t(1) t(end)];
else
    wsp = [max(min(P.trace_xlim_sp),t(1)), min(max(P.trace_xlim_sp),t(end))];
end
mwsp = t>=wsp(1) & t<=wsp(2);

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
mb = fcw>=band(1) & fcw<=band(2); Cband = mean(Cw(mb)); %#ok<NASGU>

%% ---- onset-triggered + spike-triggered dF/F ----
% Heatmap/average are triggered on inspiration ONSET (0 = onset); the breath
% PEAK is overlaid as the raster so it matches the onset-aligned panels.
if isfield(P,'trigWin_sec') && ~isempty(P.trigWin_sec) && P.trigWin_sec > 0
    win = max(1, round(P.trigWin_sec * fps));     % user-set half-window
else
    win = round(fps / f_pk);                      % auto: one breath period
end
tau = (-win:win)/fps;
on      = find(ev_foot>0);   on      = on(on-win>=1 & on+win<=numel(dff));     % trigger = insp ONSET
ovl_on  = find(ev>0);        if ~isempty(ovl_on), ovl_on = ovl_on(ovl_on-win>=1 & ovl_on+win<=numel(dff)); end  % overlay = breath PEAK

E = zeros(numel(on), 2*win+1);
for k = 1:numel(on), E(k,:) = dff(on(k)-win : on(k)+win); end
dt_to_nearest = nan(numel(on),1);
for k = 1:numel(on)
    others = on; others(k) = [];
    if isempty(others), continue; end
    [~, mi] = min(abs(others - on(k))); dt_to_nearest(k) = (others(mi) - on(k)) / fps;
end
switch lower(P.sortMode)
    case 'postmean', key = mean(E(:, tau>=0), 2); [~,si] = sort(key,'descend');
    case 'dt',       [~,si] = sort(dt_to_nearest,'ascend','MissingPlacement','last');
    otherwise,       si = (1:numel(on))';
end
Es = E(si,:); cl = prctile(Es(:), [0.5 99.5]);
mu = mean(E,1); sd = std(E,0,1);

% onset-triggered BREATH average (for the dF/F + breath overlay on ax5)
Eb = zeros(numel(on), 2*win+1);
for k = 1:numel(on), Eb(k,:) = bw(on(k)-win : on(k)+win); end
mu_b = mean(Eb,1);

% inspiration-PEAK-triggered dF/F average (trigger = breath peak)
if ~isempty(ovl_on)
    Ep = zeros(numel(ovl_on), 2*win+1);
    for k = 1:numel(ovl_on), Ep(k,:) = dff(ovl_on(k)-win : ovl_on(k)+win); end
    mu_p = mean(Ep,1); sd_p = std(Ep,0,1);
else
    Ep = zeros(0,2*win+1); mu_p = []; sd_p = [];
end

sp_trig = find(spk_train>0); sp_trig = sp_trig(sp_trig-win>=1 & sp_trig+win<=numel(dff));
if ~isempty(sp_trig)
    Esp = zeros(numel(sp_trig), 2*win+1);
    for k = 1:numel(sp_trig), Esp(k,:) = dff(sp_trig(k)-win : sp_trig(k)+win); end
    mu_sp = mean(Esp,1); sd_sp = std(Esp,0,1);
else
    Esp = zeros(0,2*win+1); mu_sp = []; sd_sp = [];
end

% breathing modulation score = 10^(peak onset-triggered dFF / peak spike-triggered dFF)
pk_on = max(mu);
if ~isempty(mu_sp), pk_sp = max(mu_sp); else, pk_sp = NaN; end
if isfinite(pk_sp) && pk_sp ~= 0
    BMS = 10^(pk_on / pk_sp);
else
    BMS = NaN;
end
fprintf('  breathing modulation score = %.3g   [10^(%.4g / %.4g)]\n', BMS, pk_on, pk_sp);

sp_on_h = find(spk_train > 0);

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
occW = histcounts(mod(phi_vv,2*pi), edgesP); Ob = mean(occW(occW>0));
wOcc = ones(1,numel(occW)); wOcc(occW>0) = Ob ./ occW(occW>0);
pctC    = (100 * cntW / max(nCycles, eps)) .* wOcc;
rmaxC   = max([pctC, 1]);
if ~isempty(spkW), muW = mod(angle(mean(exp(1i*spkW))),2*pi); else, muW = NaN; end

%% ---- coherence method: replica of Ventral_surface_coherence_polar_svd ----
have_coh = sum(spk_train) >= P.minSpikes;
th_coh = NaN; r_coh = NaN; rlo = NaN; rhi = NaN; dphi = NaN; confC = NaN;
if have_coh && ~isempty(ip)
    nD = P.nDrop;
    bwc = detrend(double(BP.breath(:))); bwc(1:min(nD,numel(bwc))) = []; bwc = bwc - mean(bwc);
    peakC = round(BP.insp_onset_idx(:)) - nD;
    footC = round(IP.insp_start_idx(:)) - nD;
    stkC  = double(CA.roi_spikes(roi).spike_train(:));
    if strcmpi(group,'Vglut2') && strcmp(recDate,'1124')
        peakC = peakC + 1; footC = footC + 1; bwc = [bwc(1); bwc(1:end-1)];
    end
    Tc = min(numel(bwc), numel(stkC));
    peakC = peakC(peakC>=1 & peakC<=Tc); footC = footC(footC>=1 & footC<=Tc);
    bwc = bwc(1:Tc); stc = stkC(1:Tc);
    if numel(peakC) >= 2 && numel(footC) >= 2
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
        phiC = piecewise_phase_local(peakC, footC, Tc);
        refC = cos(phiC); refC(isnan(refC)) = 0; refC = refC - mean(refC);
        pcoh.Fs=fps; pcoh.tapers=[P.TW_coh,2*P.TW_coh-1]; pcoh.pad=0;
        pcoh.fpass=band_c; pcoh.err=[2,P.alpha_coh];
        lagC = round(P.ca_lag_sec*fps);
        stcL = [stc(1+lagC:end); zeros(lagC,1)];
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
proj = struct('have_proj',have_proj,'crop_img',crop_img,'bnd_crop',{bnd_crop}, ...
              'barLen_pr',barLen_pr,'roi',roi);

%% ===================== FIGURE =====================
nRows = 7;  nCols = 4;
tile_trace = 1;   span_trace = [2 4];
tile_saw   = 9;   span_saw   = [1 4];
tile_proj  = 13;  tile_ptavg = 14;  tile_psd  = 15;  tile_cohSp = 16;
tile_pthm  = 17;  tile_peri  = 18;  tile_phist = 19;
tile_epol  = 20;  tile_cpol  = 21;  tile_stavg = 22;  tile_pktavg = 23;

fig = figure('Color','w','Name',sprintf('%s ROI%d temporal phase (window)',folderPath,roi), ...
             'Units','normalized','Position',[0.03 0.03 0.9 0.8]);
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
relPath = sprintf('%s/%s/%s', group, recDate, recName);
nBwin = nnz((footAll-1)/fps>=wsp(1) & (footAll-1)/fps<=wsp(2));   % breaths in the window
title(tl, sprintf(['%s   ROI%d   |   fps=%.2f, breath peak %.2f Hz   |   window %.1f-%.1f s (%d breaths)'], ...
      relPath, roi, fps, f_pk, wsp(1), wsp(2), nBwin), 'Interpreter','none','FontWeight','bold');

% trace (CLIPPED to window)
ax1 = nexttile(tl,tile_trace,span_trace);
yyaxis(ax1,'right'); plot(ax1, t, bw, '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
set(ax1,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax1,'breath');
if nnz(mwsp)>2, set_ylim_win(ax1, bw, mwsp); end
yyaxis(ax1,'left'); plot(ax1, t, dff, '-','Color',dffColor,'LineWidth',0.8);
set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
if nnz(mwsp)>2, set_ylim_win(ax1, dff, mwsp); end
xlim(ax1,wsp); xlabel(ax1,'Time (s)'); box(ax1,'off');

% phase sawtooth (CLIPPED to window)
bax = nexttile(tl,tile_saw,span_saw); hold(bax,'on');
% each cycle ramps 0 (onset) -> 2*pi (next onset); pin the top to 2*pi and break
% the line at the wrap so there is no down-stroke (red onset line marks the reset).
tw = t(mwsp); pw = phiW(mwsp); TT = []; PP = [];
for i = 1:numel(tw)
    if i>1 && (pw(i)-pw(i-1)) < -pi
        TT(end+1)=tw(i); PP(end+1)=2*pi;  %#ok<AGROW>  close previous ramp at 2*pi
        TT(end+1)=tw(i); PP(end+1)=NaN;   %#ok<AGROW>  break (no down-stroke)
    end
    TT(end+1)=tw(i); PP(end+1)=pw(i);     %#ok<AGROW>
end
plot(bax, TT, PP, 'k-','LineWidth',0.8);
set(bax,'YTick',[0 pi 2*pi],'YTickLabel',{'0','\pi','2\pi'}); ylabel(bax,'phase'); ylim(bax,[0 2*pi]);
yl = ylim(bax);
on_t = (footAll-1)/fps; pk_t = (peakAll-1)/fps;   % frame k -> t=(k-1)/fps, matches the sawtooth
on_t = on_t(on_t>=wsp(1)&on_t<=wsp(2)); pk_t = pk_t(pk_t>=wsp(1)&pk_t<=wsp(2));
for x=on_t(:)', plot(bax,[x x],yl,'-','Color',onsetCol,'LineWidth',0.8); end
for x=pk_t(:)', plot(bax,[x x],yl,'-','Color',peakCol,'LineWidth',0.8); end
spk_t = (spk_src-1)/fps; spk_t = spk_t(spk_t>=wsp(1)&spk_t<=wsp(2));   % lead-shifted (ca_lag), matches phase panels
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

% onset-triggered average: dF/F (left) + breath (right)
ax5 = nexttile(tl,tile_ptavg);
yyaxis(ax5,'left'); hold(ax5,'on');
fill(ax5,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],lightDff,'EdgeColor','none','HandleVisibility','off');
plot(ax5, tau, mu, '-', 'Color',dffColor,'LineWidth',1.5);
ylabel(ax5,'dF/F'); ax5.YAxis(1).Color = dffColor*0.8;
yyaxis(ax5,'right');
plot(ax5, tau, mu_b, '-', 'Color',[0.45 0.45 0.45],'LineWidth',1.2);
ylabel(ax5,'breath'); ax5.YAxis(2).Color = [0.45 0.45 0.45];
xline(ax5,0,'k--','LineWidth',0.8);
xlim(ax5,[tau(1) tau(end)]); xlabel(ax5,'time from insp onset (s)'); grid(ax5,'on'); pbaspect(ax5,[1 1 1]);
title(ax5, sprintf('breath mod score = %.3g', BMS));

% onset-triggered individual-cycle heatmap (0 = insp onset; peaks overlaid)
ax4 = nexttile(tl,tile_pthm);

imagesc(ax4, tau, 1:size(Es,1), Es); axis(ax4,'tight');
%colormap(ax4, flipud(gray(256))); caxis(ax4, cl); hold(ax4,'on');
colormap(ax4, gray(256)); caxis(ax4, cl); hold(ax4,'on');
if ~isempty(ovl_on)
    on_sorted = on(si); pk_x = []; pk_y = [];
    for k = 1:numel(on_sorted)
        p  = on_sorted(k); pw = ovl_on(ovl_on >= p-win & ovl_on <= p+win);
        pk_x = [pk_x; (pw - p) / fps]; pk_y = [pk_y; repmat(k, numel(pw),1)]; %#ok<AGROW>
    end
    plot(ax4, pk_x, pk_y, '.', 'Color', peakCol, 'MarkerSize', 1);
end
hold(ax4,'off'); set(ax4,'YDir','reverse'); xlabel(ax4,'time from insp onset (s)'); ylabel(ax4,'breath #');
pbaspect(ax4,[1 1 1]); cb=colorbar(ax4); cb.Label.String='dF/F';

% ONSET-triggered spike/peak histogram, CENTERED AT 0 (own window P.histWin_sec)
if isfield(P,'histWin_sec') && ~isempty(P.histWin_sec) && P.histWin_sec > 0
    winH = max(1, round(P.histWin_sec * fps));
else
    winH = round(fps / f_pk);                            % auto: one breath period
end
onH = find(ev_foot>0); onH = onH(onH-winH>=1 & onH+winH<=numel(dff));   % onsets valid for winH
if ~isempty(onH)
    Whist  = winH / fps;
    nbf    = 2; if isfield(P,'histBinFrames'), nbf = max(1,round(P.histBinFrames)); end
    Mb     = floor(winH/nbf);                            % # bins each side
    ctrsC  = (-Mb:Mb)*nbf/fps;                           % bin centers (0 = insp onset)
    edgesC = ((-Mb-0.5):(Mb+0.5))*nbf/fps;               % bin width = nbf frames
    spkS   = spk_src(:);
    dtS = []; dtP = [];
    for k = 1:numel(onH)
        o = onH(k);
        s = spkS(spkS >= o-winH & spkS <= o+winH);          dtS = [dtS; (s - o)/fps]; %#ok<AGROW>
        p = peakAll(peakAll >= o-winH & peakAll <= o+winH); dtP = [dtP; (p - o)/fps]; %#ok<AGROW>
    end
    spkH = 100*histcounts(dtS, edgesC)/max(numel(onH),1);   % spikes per cycle, %
    cP   = histcounts(dtP, edgesC); pkH = 100*cP/max(sum(cP),1);

    % ---- circular-shift shuffle null for the spike histogram ----
    nSh = 0; if isfield(P,'nShuffle'), nSh = P.nShuffle; end
    null_lo = []; null_hi = []; null_mu = [];
    if nSh > 0 && ~isempty(spkS)
        Tn = numel(dff);
        % min shift must exceed (a) the PSTH half-window and (b) several breath
        % periods, so shifted spikes are fully decoupled from inspiration phase.
        shiftMinCyc = 3; if isfield(P,'shiftMinCyc'), shiftMinCyc = P.shiftMinCyc; end
        minShift = max([winH+1, round(shiftMinCyc/max(f_pk,eps)*fps), round(1.0*fps)]);
        if Tn - minShift > minShift
            nb = numel(ctrsC); nullH = zeros(nSh, nb);
            for ii = 1:nSh
                sh = randi([minShift, Tn-minShift]);
                sS = mod(spkS - 1 + sh, Tn) + 1;            % circular shift in [1,Tn]
                cnt = zeros(1, nb);
                for k = 1:numel(onH)
                    o = onH(k);
                    cnt = cnt + histcounts((sS(sS>=o-winH & sS<=o+winH) - o)/fps, edgesC);
                end
                nullH(ii,:) = 100*cnt/max(numel(onH),1);
            end
            null_lo = prctile(nullH, 2.5, 1);
            null_hi = prctile(nullH, 97.5, 1);
            null_mu = mean(nullH, 1);
        end
    end

    ax_h1 = nexttile(tl,tile_peri); hold(ax_h1,'on');
    yyaxis(ax_h1,'left');
    if ~isempty(null_hi)
        fill(ax_h1, [ctrsC fliplr(ctrsC)], [null_hi fliplr(null_lo)], [0.6 0.6 0.6], ...
             'FaceAlpha',0.30, 'EdgeColor','none', 'HandleVisibility','off');   % 95% null band
        plot(ax_h1, ctrsC, null_mu, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',0.8);  % null mean
    end
    bar(ax_h1, ctrsC, spkH, 1, 'FaceColor','k','FaceAlpha',0.85,'EdgeColor','none');
    yTopL = max([spkH, null_hi, 1]) * 1.10;
    set(ax_h1,'YColor','k'); ylim(ax_h1,[0 yTopL]); ylabel(ax_h1,'spk/cyc %');
    yyaxis(ax_h1,'right');
    bar(ax_h1, ctrsC, pkH, 1, 'FaceColor',peakCol,'FaceAlpha',0.45,'EdgeColor','none');
    set(ax_h1,'YColor',peakCol*0.7); ylim(ax_h1,[0 max([pkH 1])*1.10]); ylabel(ax_h1,'peak %');
    xline(ax_h1, 0, 'Color',onsetCol,'LineWidth',1);
    xlim(ax_h1,[-Whist Whist]);
    xlabel(ax_h1,'time from insp onset (s)');
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

% event-phase polar
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

% linear 0..2pi
lax2 = nexttile(tl,tile_phist); hold(lax2,'on');
bar(lax2, ctrsP, pctC, 1, 'FaceColor',[0.30 0.30 0.30],'FaceAlpha',0.75,'EdgeColor','none');
xline(lax2, 0,  'Color',onsetCol,'LineWidth',1);
xline(lax2, pi, 'Color',peakCol, 'LineWidth',1);
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

% inspiration-peak-triggered average
ax_pk = nexttile(tl,tile_pktavg); hold(ax_pk,'on');
if ~isempty(mu_p)
    fill(ax_pk,[tau fliplr(tau)],[mu_p+sd_p fliplr(mu_p-sd_p)],lightDff,'EdgeColor','none');
    plot(ax_pk, tau, mu_p, 'Color',dffColor,'LineWidth',1.5);
    xlim(ax_pk,[tau(1) tau(end)]);
end
xline(ax_pk,0,'Color',peakCol,'LineWidth',0.8);
xlabel(ax_pk,'time from insp peak (s)'); ylabel(ax_pk,'dF/F'); grid(ax_pk,'on'); pbaspect(ax_pk,[1 1 1]);
title(ax_pk, sprintf('peak-triggered average (\\pm SD, n=%d)', size(Ep,1)));
end

%% ===================== LOCAL HELPERS =====================
function set_ylim_win(ax, y, mwsp)
    lo = min(y(mwsp)); hi = max(y(mwsp)); pad = 0.05*max(hi-lo, eps);
    ylim(ax, [lo-pad hi+pad]);
end

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
