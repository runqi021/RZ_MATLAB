function cell_pool_260727()
%% cell_pool_260727  Pool every cell's observations into ONE lossless container.
% -----------------------------------------------------------------------
% STEP 2 of the cell-aware coherence pipeline. Still computes no coherence. It
% gathers, for every curated cell, all the raw material needed by ANY downstream
% analysis, once, in one file -- so later analyses never re-derive alignment and
% can never disagree about it.
%
% DESIGN RULE: nothing is lost and nothing is distorted.
%   * Breath data is stored ONCE PER RECORDING, not per observation, so multiple
%     ROIs from the same recording cannot end up with subtly different breath.
%   * Spike trains are stored SPARSE and EXACT (no thresholding, no rebinning,
%     no smoothing, no lead-shift). The 100 ms GCaMP lead-shift the foundation
%     applies is an ANALYSIS choice and is applied downstream, not baked in here.
%   * Recordings are NEVER resampled to a common rate. Each keeps its own fps and
%     its own frequency grid; pooling happens in the spectral domain later, which
%     is what makes that safe.
%   * The alignment truncation (T = min(breath, calcium)) is recorded frame by
%     frame -- n_frames_ca_orig, n_frames_breath_orig, T, and how many spikes fell
%     OUTSIDE the aligned window -- so any loss is visible instead of implied.
%   * Tossed observations are excluded from cells but still listed, with a reason.
%   * Recordings the link audit flagged as a suspect join are REFUSED, loudly,
%     rather than pooled on a join that may address the wrong neurons.
%
% Alignment (identical to coherence_polar_general_260727.m, on purpose):
%   fps      = detect_session_fps(recording folder)
%   nDrop    = 30 breath frames tossed up front to match the calcium's toss
%   T        = min(numel(breath), n calcium frames); breath truncated, NOT resampled
%   phase    = piecewise linear, foot = 0, peak = pi
%
% Output (into cfg.cellDir):  cell_pool.mat
%   pool.rec(k)   per recording: fps, T, breath waveform, phase, ref, band, f_pk
%   pool.obs(i)   per observation: cell_id, rec, roi_index, sparse spike train, dFF
%   pool.cells{c} observation indices belonging to cell c
%   pool.audit    everything that was dropped, and why
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
% These MUST match coherence_polar_general_260727.m or the pooled and per-ROI
% results stop being comparable. They are repeated here rather than imported so
% you can see them; if you change one, change it in both.
nDrop           = 30;        % breath frames tossed up front (align to Ca)
fallback_fps    = 30;
TW              = 4;         % multitaper TW (used for the breath PSD band here)
f_breath_search = [0.2 4];   % Hz, search band for the breath PSD peak
fwhm_factor     = 0.6;       % coherence band = fwhm_factor x FWHM
min_bw          = 0.05;      % Hz, minimum band width
fmin            = 0.05;      % Hz, PSD lower bound
fmax            = 15;        % Hz, PSD upper bound

INCLUDE_DFF     = true;      % also store each observation's dF/F trace (double).
                             %   true  = the container is self-sufficient
                             %   false = smaller file; traces reachable via dff_file
% =====================================================================

peak_file_name  = 'breath_peak_pc1.mat';
start_file_name = 'breath_insp_start_pc1.mat';

fprintf('\n=============== cell_pool_260727 ===============\n');
assert(isfile(cfg.linkFile), ['cell_link.mat not found:\n  %s\nRun cell_link_260727.m first.'], cfg.linkFile);
L = load(cfg.linkFile, 'link');  link = L.link;
obsT = link.obsT;  recT = link.recT;  cells = link.cells;
fprintf('link: %d observations, %d cells, %d recordings\n', ...
        height(obsT), numel(cells), height(recT));

% recordings the audit refuses to pool
badRec = [link.audit.rec_roi_count_mismatch(:); link.audit.rec_roi_index_overflow(:); ...
          link.audit.rec_labels_not_contiguous(:)];
badRec = unique(badRec);

%% ---- 1. per-recording breath + alignment ----
nRec = height(recT);
rec  = struct('name',{},'path',{},'fps',{},'T',{},'nCa',{},'nBreath',{}, ...
              'bw',{},'phi',{},'ref',{},'band',{},'f_pk',{}, ...
              'peak_idx',{},'foot_idx',{},'usable',{},'why',{});
for k = 1:nRec
    r = struct('name',recT.rec_name(k), 'path',recT.rec_path(k), 'fps',NaN, 'T',NaN, ...
               'nCa',NaN, 'nBreath',NaN, 'bw',[], 'phi',[], 'ref',[], 'band',[NaN NaN], ...
               'f_pk',NaN, 'peak_idx',[], 'foot_idx',[], 'usable',false, 'why',"");
    fp = char(recT.rec_path(k));
    if any(badRec == recT.rec_name(k))
        r.why = "REFUSED: link audit flagged this recording's ROI indexing";
        rec(end+1) = r; continue; %#ok<AGROW>
    end
    bpFile = fullfile(fp, peak_file_name);
    ipFile = fullfile(fp, start_file_name);
    if ~isfile(bpFile) || ~isfile(ipFile)
        r.why = "no breath_peak_pc1.mat / breath_insp_start_pc1.mat";
        rec(end+1) = r; continue; %#ok<AGROW>
    end
    try
        r.fps = detect_session_fps(fp, fallback_fps);
        CA = load(fullfile(fp,'ca_spike_data.mat'), 'roi_spikes');
        r.nCa = numel(CA.roi_spikes(1).spike_train);

        BP = load(bpFile);  IP = load(ipFile);
        bw = detrend(double(BP.breath(:)));
        r.nBreath = numel(bw);
        bw(1:min(nDrop,numel(bw))) = [];
        bw = bw - mean(bw);

        peak_idx = round(BP.insp_onset_idx(:)) - nDrop;
        foot_idx = round(IP.insp_start_idx(:)) - nDrop;

        T = min([numel(bw), r.nCa]);
        peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
        foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
        if numel(peak_idx) < 2 || numel(foot_idx) < 2
            r.why = "fewer than 2 breath peaks or feet in the aligned window";
            rec(end+1) = r; continue; %#ok<AGROW>
        end
        phi = piecewise_phase_local(peak_idx, foot_idx, T);
        ref = cos(phi);  ref(isnan(ref)) = 0;  ref = ref - mean(ref);
        bw  = bw(1:T);

        % breath PSD -> the coherence band for THIS recording (breathing rate
        % differs between recordings; each keeps its own band, nothing is forced
        % onto a shared band it does not occupy)
        pB.Fs = r.fps; pB.tapers = [TW, 2*TW-1]; pB.pad = 0;
        pB.fpass = [fmin, min(fmax, r.fps/2)]; pB.err = 0;
        [Sb, fb] = mtspectrumc(bw, pB);  Sb = Sb(:);  fb = fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        [~,rl] = max(Sb(m));  ipk = find(m,1)+rl-1;  f_pk = fb(ipk);
        h = Sb(ipk)/2;  lo = ipk;  while lo>1 && Sb(lo)>h, lo = lo-1; end
        hi = ipk;                  while hi<numel(fb) && Sb(hi)>h, hi = hi+1; end
        f_fwhm = [max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd  = max(diff(f_fwhm)*fwhm_factor, min_bw);
        band = [max(f_pk-bwd/2, f_breath_search(1)), min(f_pk+bwd/2, f_breath_search(2))];

        r.T = T;  r.bw = bw;  r.phi = phi;  r.ref = ref;
        r.band = band;  r.f_pk = f_pk;  r.peak_idx = peak_idx;  r.foot_idx = foot_idx;
        r.usable = true;  r.why = "";
        fprintf('  [%2d/%2d] %-52s fps %5.2f  T %5d  band [%.2f %.2f] Hz\n', ...
                k, nRec, r.name, r.fps, T, band(1), band(2));
    catch ME
        r.why = string(sprintf('ERROR: %s', ME.message));
        warning('  %s: %s', r.name, ME.message);
    end
    rec(end+1) = r; %#ok<AGROW>
end
recUsable = [rec.usable]';
fprintf('  usable recordings: %d of %d\n', nnz(recUsable), nRec);

%% ---- 2. per-observation spike trains (+ dF/F), aligned to that recording ----
nObs = height(obsT);
obs = struct('obs',{},'cell_id',{},'rec',{},'rec_name',{},'roi_index',{}, ...
             'status',{},'spikes',{},'n_spikes_full',{},'n_spikes_used',{}, ...
             'n_spikes_lost_to_alignment',{},'dff',{},'dff_file',{},'usable',{},'why',{});
dffCache = containers.Map('KeyType','char','ValueType','any');
for i = 1:nObs
    o = struct('obs',obsT.obs(i), 'cell_id',obsT.cell_id(i), 'rec',NaN, ...
               'rec_name',obsT.rec_name(i), 'roi_index',obsT.roi_index(i), ...
               'status',obsT.status(i), 'spikes',sparse(0,1), 'n_spikes_full',NaN, ...
               'n_spikes_used',NaN, 'n_spikes_lost_to_alignment',NaN, 'dff',[], ...
               'dff_file',"", 'usable',false, 'why',"");
    if o.status == "tossed"
        o.why = "tossed during curation";  obs(end+1) = o; continue; %#ok<AGROW>
    end
    k = find(recT.rec_path == obsT.rec_path(i), 1);
    if isempty(k)
        o.why = "no recording folder (matcher FOV without ca_spike_data.mat)";
        obs(end+1) = o; continue; %#ok<AGROW>
    end
    o.rec = k;
    if ~rec(k).usable
        o.why = "recording unusable: " + rec(k).why;  obs(end+1) = o; continue; %#ok<AGROW>
    end
    fp = char(recT.rec_path(k));
    try
        CA = load(fullfile(fp,'ca_spike_data.mat'), 'roi_spikes');
        ri = o.roi_index;
        if ri < 1 || ri > numel(CA.roi_spikes)
            o.why = sprintf("roi_index %d outside roi_spikes (n=%d)", ri, numel(CA.roi_spikes));
            obs(end+1) = o; continue; %#ok<AGROW>
        end
        st_full = double(CA.roi_spikes(ri).spike_train(:));
        T = rec(k).T;
        st = st_full(1:min(T,numel(st_full)));
        if numel(st) < T, st(end+1:T,1) = 0; end          % zero-pad, same as the foundation
        o.n_spikes_full = sum(st_full > 0);
        o.n_spikes_used = sum(st > 0);
        o.n_spikes_lost_to_alignment = o.n_spikes_full - o.n_spikes_used;
        o.spikes  = sparse(st);                            % EXACT values, no thresholding
        o.usable  = true;

        if INCLUDE_DFF
            if ~isKey(dffCache, fp)
                dd = dir(fullfile(fp, '*_dFF.mat'));
                if ~isempty(dd)
                    D = load(fullfile(dd(1).folder, dd(1).name), 'dFF');
                    dffCache(fp) = {D.dFF, string(fullfile(dd(1).folder, dd(1).name))};
                else
                    dffCache(fp) = {[], ""};
                end
            end
            cc = dffCache(fp);
            o.dff_file = cc{2};
            if ~isempty(cc{1}) && ri <= size(cc{1},2)
                d = cc{1}(:,ri);
                o.dff = d(1:min(T,numel(d)));
            end
        end
    catch ME
        o.why = string(sprintf('ERROR: %s', ME.message));
    end
    obs(end+1) = o; %#ok<AGROW>
end

%% ---- 3. rebuild the cells list over USABLE observations ----
usable   = [obs.usable]';
poolCells = cell(numel(cells),1);
for c = 1:numel(cells)
    ii = find(~isnan([obs.cell_id]') & [obs.cell_id]' == c & usable);
    poolCells{c} = ii(:);
end
poolSize = cellfun(@numel, poolCells);

%% ---- 4. AUDIT ----
a = struct();
a.n_obs               = nObs;
a.n_obs_usable        = nnz(usable);
a.n_obs_tossed        = nnz([obs.status]' == "tossed");
a.n_obs_dropped       = nObs - nnz(usable) - nnz([obs.status]' == "tossed");
a.n_rec               = nRec;
a.n_rec_usable        = nnz(recUsable);
a.n_cells             = numel(cells);
a.n_cells_with_data   = nnz(poolSize > 0);
a.n_cells_multi       = nnz(poolSize > 1);
a.n_spikes_total      = sum([obs(usable).n_spikes_used]);
a.n_spikes_lost_align = sum([obs(usable).n_spikes_lost_to_alignment]);
a.refused_recordings  = badRec;

fprintf('\n---- pooled ----\n');
fprintf('  observations : %d usable | %d tossed | %d dropped for other reasons\n', ...
        a.n_obs_usable, a.n_obs_tossed, a.n_obs_dropped);
fprintf('  recordings   : %d usable of %d\n', a.n_rec_usable, a.n_rec);
fprintf('  cells        : %d total | %d with >=1 usable observation | %d with >=2\n', ...
        a.n_cells, a.n_cells_with_data, a.n_cells_multi);
fprintf('  spikes       : %d pooled | %d fell outside the aligned window (%.2f%%)\n', ...
        a.n_spikes_total, a.n_spikes_lost_align, ...
        100*a.n_spikes_lost_align/max(a.n_spikes_total + a.n_spikes_lost_align, 1));
if a.n_obs_dropped > 0
    fprintf('\n  dropped observations, by reason:\n');
    why = [obs(~usable & [obs.status]' ~= "tossed").why]';
    u = unique(why);
    for q = 1:numel(u), fprintf('    %4d  %s\n', nnz(why==u(q)), u(q)); end
end
fprintf('\n  pooled cell size distribution: ');
for s = 1:max([poolSize;0]), if nnz(poolSize==s), fprintf('%dx:%d  ', s, nnz(poolSize==s)); end, end
fprintf('\n');

%% ---- 5. save ----
% rec and obs are struct ARRAYS -- they must be cell-wrapped or struct() would
% return a 1xN struct array instead of one scalar struct with those fields.
% Carry the link's observation table through verbatim. It holds identity
% (cell_id, cell_size, status) AND geometry (roi and cell coordinates), so every
% downstream analysis can attach a result to a place without re-joining anything.
pool = struct('obsT',obsT, 'recT',recT, ...
              'rec',{rec}, 'obs',{obs}, 'cells',{poolCells}, 'cellSize',poolSize, ...
              'params',struct('nDrop',nDrop,'TW',TW,'f_breath_search',f_breath_search, ...
                              'fwhm_factor',fwhm_factor,'min_bw',min_bw,'fmin',fmin, ...
                              'fmax',fmax,'fallback_fps',fallback_fps,'INCLUDE_DFF',INCLUDE_DFF), ...
              'audit',a, 'link_file',string(cfg.linkFile), 'cfg',cfg); %#ok<NASGU>
save(cfg.poolFile, 'pool', '-v7.3');
d = dir(cfg.poolFile);
fprintf('\nSaved %s (%.1f MB)\n', cfg.poolFile, d.bytes/1e6);
fprintf('Next: cell_coherence_pooled_260727.m\n');
end

%% ========================= helpers =========================
function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase: FEET at 0/2pi/..., PEAKS at pi/3pi/...; linear ramps
% between consecutive events; NaN outside the range.
% (Verbatim from coherence_polar_general_260727.m so the phase is identical.)
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
phi_cur = types(1) * pi;        % peak->pi, foot->0
for i = 1:numel(events)
    phases(i) = phi_cur; phi_cur = phi_cur + pi;
end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
