function r = ca_recheck_load_obs_260901(o, P)
%CA_RECHECK_LOAD_OBS_260901  One (cell, recording) observation: dF/F + breath + events.
%
%   r = ca_recheck_load_obs_260901(o, P)
%   ca_recheck_load_obs_260901('clearcache')     % drop the folder cache
%
%   o : struct with fields
%         folder   recording folder (needs *_ch1_dFF.mat + breath_peak_pc1.mat)
%         roi      ROI index within THAT recording's dFF matrix
%         recName  label                      (optional)
%         group    ventral group label -- 'IO' is a SITE, not a genotype (optional)
%         recDate  MMDD                       (optional)
%   P : struct with .nDrop (default 30) and .fallback_fps (default 30)
%
%   THIS IS A DELIBERATE RE-IMPLEMENTATION of load_obs_local() inside
%   analysis_260727\coh_ca_breath\temporal_phase_cell_fig_260812.m, which is a
%   local function and so cannot be called from outside that file. Every
%   alignment step is kept identical -- nDrop on the breath side ONLY, the
%   Vglut2/1124 one-frame shift, and the T = min(dFF, breath, events) truncation
%   -- so a trace pulled by this loader is frame-for-frame the trace the per-cell
%   summary figure draws. The only thing dropped is the Chronux breath spectrum
%   (f_pk), which nothing here needs; IBI comes from the onset intervals.
%
%   ONE RECORDING SERVES MANY CELLS (2840 observations over 125 folders in the
%   archive), so the per-folder part -- the dF/F matrix, the breath trace, the
%   landmark trains, ca_spike_data -- is held in a persistent single-slot cache.
%   Calling this in folder order therefore costs one file read per folder rather
%   than one per observation. The slot is keyed on folder AND nDrop, so changing
%   nDrop mid-session cannot silently return a stale alignment.
%
%   Never throws for a missing ca_spike_data.mat -- that is reported as
%   r.has_spike = false, since the whole point of the re-check is that some
%   recordings were never curated.
%
%   Runqi Zhang / 2026-09-01

persistent CACHE

if nargin == 1 && (ischar(o) || isstring(o)) && strcmpi(o,'clearcache')
    CACHE = []; r = []; return
end

if nargin < 2 || isempty(P), P = struct(); end
if ~isfield(P,'nDrop'),        P.nDrop = 30;        end
if ~isfield(P,'fallback_fps'), P.fallback_fps = 30; end
if ~isfield(o,'group'),   o.group   = '';  end
if ~isfield(o,'recDate'), o.recDate = '';  end
if ~isfield(o,'recName'), o.recName = '';  end

key = sprintf('%s|%g', o.folder, P.nDrop);
if isempty(CACHE) || ~strcmp(CACHE.key, key)
    CACHE = local_load_folder(o, P);
    CACHE.key = key;
end
C = CACHE;

assert(o.roi>=1 && o.roi<=C.nROItot, 'ROI %d out of range (1..%d) in %s', ...
       o.roi, C.nROItot, o.folder);

dff = C.dff_all(1:C.T, o.roi);

Fraw = zeros(C.T,1);
if ~isempty(C.F_all) && o.roi <= size(C.F_all,2)
    Fraw = C.F_all(1:C.T, o.roi);
end

spk = zeros(C.T,1); has_spike = false;
if ~isempty(C.spk_all) && o.roi <= size(C.spk_all,2)
    s = C.spk_all(:, o.roi);
    if numel(s) < C.T, s(end+1:C.T) = 0; end
    spk = s(1:C.T); has_spike = C.has_spike;
end

r = struct( ...
    'folder',   o.folder, ...
    'recName',  o.recName, ...
    'roi',      o.roi, ...
    'group',    o.group, ...
    'recDate',  o.recDate, ...
    'fps',      C.fps, ...
    'T',        C.T, ...
    'dur_s',    C.T / C.fps, ...
    'px_um',    C.px_um, ...
    'dff',      dff, ...
    'F_roi',    Fraw, ...
    'breath',   C.breath, ...
    'spk',      spk, ...
    'spike_idx',find(spk>0), ...
    'peak',     C.peak, ...
    'foot',     C.foot, ...
    'IBI',      C.IBI, ...
    'nEv',      nnz(spk>0), ...
    'has_spike',has_spike, ...
    'spike_src',C.spike_src, ...
    'nROItot',  C.nROItot);
end

% =========================================================================
function C = local_load_folder(o, P)
%LOCAL_LOAD_FOLDER  Everything one recording folder contributes, ROI-independent.
df = dir(fullfile(o.folder,'*_ch1_dFF.mat'));
bp = dir(fullfile(o.folder,'breath_peak_pc1.mat'));
ip = dir(fullfile(o.folder,'breath_insp_start_pc1.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', o.folder);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', o.folder);

[fps, sm] = detect_session_fps(o.folder, P.fallback_fps);
D  = load(fullfile(df(1).folder, df(1).name));
BP = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(D.dFF);

% Raw F, the trace the dF/F was computed FROM, for the GUI's Raw F panel.
% F_roi in the dFF file is already frame-aligned with dFF (TossFrames removed),
% so it needs no shifting; dFFout.F_dff is the same matrix under the pipeline's
% other name, and the cpSAM F is the last resort.
F_all = [];
if     isfield(D,'F_roi'),  F_all = double(D.F_roi);
elseif isfield(D,'dFFout') && isfield(D.dFFout,'F_dff'), F_all = double(D.dFFout.F_dff);
else
    sam = dir(fullfile(o.folder,'*_cpSAM_output.mat'));
    if ~isempty(sam)
        SAMF = load(fullfile(sam(1).folder, sam(1).name),'F');
        if isfield(SAMF,'F'), F_all = double(SAMF.F); end
    end
end
if ~isempty(F_all) && size(F_all,1) < size(dff_all,1), F_all = []; end

px_um = NaN;
if isfield(sm,'pixelSize_um') && isfinite(sm.pixelSize_um) && sm.pixelSize_um>0
    px_um = sm.pixelSize_um;
elseif isfield(sm,'zoomFactor') && isfinite(sm.zoomFactor) && sm.zoomFactor>0
    px_um = 1.7778 / sm.zoomFactor;      % PixelSizeBase / zoom
end

% ---- breath waveform and its two landmark trains -------------------------
% nDrop is applied to the BREATH side only: the dF/F matrix already had
% TossFrames removed by the pipeline, so dropping again here would shift the
% calcium trace against breathing by nDrop frames.
bw = detrend(double(BP.breath(:)));
bw(1:min(P.nDrop,numel(bw))) = [];
bw = bw - mean(bw);

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
    ev_foot(fi(fi>=1 & fi<=nB)) = 1;
    ev_foot(1:min(P.nDrop,numel(ev_foot))) = [];
else
    ev_foot = [];
end

% Vglut2/1124 used a rising-edge trigger, so breath and events are one frame
% late. Gate on the GENOTYPE FOLDER, not the group label: that session's IO sites
% carry group 'IO' and would otherwise skip the fix.
gFix = o.group;
if isempty(gFix) || strcmpi(gFix,'IO'), gFix = local_detect_genotype(o.folder); end
if strcmpi(gFix,'Vglut2') && strcmp(o.recDate,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end

% ---- existing curated / detected events, ALL ROIs -------------------------
spk_all = []; has_spike = false; spike_src = '';
sp_file = fullfile(o.folder,'ca_spike_data.mat');
if isfile(sp_file)
    CA = load(sp_file,'roi_spikes');
    if isfield(CA,'roi_spikes') && ~isempty(CA.roi_spikes)
        rs = CA.roi_spikes;
        L  = max(arrayfun(@(x) numel(x.spike_train), rs));
        spk_all = zeros(L, numel(rs));
        for k = 1:numel(rs)
            v = double(rs(k).spike_train(:));
            spk_all(1:numel(v), k) = v;
        end
        has_spike = true; spike_src = sp_file;
    end
end

% ---- common length -------------------------------------------------------
T = min([size(dff_all,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); else, ev_foot = zeros(T,1); end

peak = find(ev>0);  foot = find(ev_foot>0);

IBI = NaN;
if numel(foot) >= 2, IBI = median(diff(foot))/fps; end
if (~isfinite(IBI) || IBI <= 0) && numel(peak) >= 2, IBI = median(diff(peak))/fps; end

C = struct('key','', 'fps',fps, 'T',T, 'px_um',px_um, 'dff_all',dff_all, ...
           'F_all',F_all, ...
           'breath',bw, 'peak',peak, 'foot',foot, 'IBI',IBI, ...
           'spk_all',spk_all, 'has_spike',has_spike, 'spike_src',spike_src, ...
           'nROItot',size(dff_all,2));
end

% -------------------------------------------------------------------------
function g = local_detect_genotype(fp)
%LOCAL_DETECT_GENOTYPE  Genotype folder name from an archive path.
g = '';
parts = regexp(fp, '[\\/]', 'split');
known = {'Vglut2_test','ChAT','Sst','Vgat','Vglut2','Sert'};
for k = numel(parts):-1:1
    hit = find(strcmpi(parts{k}, known), 1);
    if ~isempty(hit), g = known{hit}; return; end
end
end
