% polar_coh_rayleigh_260824.m
% -----------------------------------------------------------------------
%  FOUR polar figures, two scores x two scopes:
%
%     1  polar_coherence_new.png    coherence, 260824 vagotomised session only
%     2  polar_rayleigh_new.png     Rayleigh,  260824 vagotomised session only
%     3  polar_coherence_all.png    coherence, every dataset we have
%     4  polar_rayleigh_all.png     Rayleigh,  every dataset we have
%
%  WHAT EACH FIGURE PLOTS
%
%   COHERENCE     radius = band-averaged |C| at breath frequency
%                 angle  = coherence phase
%                 span   = JACKKNIFE error: the angular arc is th +/- 1.96*phistd
%                          (Chronux's jackknife phase s.d., averaged over the
%                          band), and for significant points the radial bar is
%                          the jackknife magnitude CI (Cerr).
%                 One dot per ROI-OBSERVATION (recording x ROI). Coherence is
%                 a per-recording spectral quantity; there is no honest way to
%                 pool it across recordings here, and taking each cell's BEST
%                 recording would bias r upward by max-selection. So this panel
%                 keeps its native unit and says so.
%
%   RAYLEIGH      radius = log Z   (Z = n_eff * Rbar^2, occupancy-corrected)
%                 angle  = CIRCULAR MEDIAN of the event phases
%                 span   = occupancy-weighted circular IQR about that median
%                 One dot per CELL, event phases pooled over all of the cell's
%                 recordings. Phase is the coordinate where pooling is legitimate:
%                 0 = insp onset and pi = insp peak mean the same landmark in
%                 every recording. Each recording keeps its own occupancy weights.
%
%   NOTE the Rayleigh dot is the MEDIAN here, not the mean direction. The 260808
%   figure drew the dot at the mean while the arc was a median-centred IQR, which
%   it flagged as a caveat (the two disagree by >30 deg for 61 of 445 cells).
%   Median + IQR about the median is internally consistent, and is what was asked.
%
%  BOTH scores read the SAME piecewise-linear breath phase, interpolated between
%  detected landmarks: insp ONSET (foot) -> 0, insp PEAK -> pi, next onset -> 2pi,
%  linear in time between them. Occupancy is ~7x denser in the slow expiratory
%  half, so the Rayleigh resultant is occupancy-normalised (weight = 1/frames in
%  that phase bin) and uses Kish's effective n. See the 260808 script's header for
%  the full argument; the maths below is copied from it unchanged.
%
%  SCOPES
%   "new"  = C:\260824_Vglut2-soma-g8s_vagotomized\phys, a FLAT session folder.
%            Cell identity comes from its own roi_match_out_260824\
%            roi_match_curated.csv -- the same curation the per-cell GUI-default
%            renders used, so cell 17 here is cell 17 there.
%   "all"  = that session PLUS the D:\Ventral_surface_summary archive
%            (ChAT/IO, Vglut2, Vgat, Sst, Sert), nested Group\Date\site\rec.
%            The vagotomised session is drawn as its OWN group in its own colour
%            rather than merged into Vglut2 -- it is a different preparation, and
%            burying it inside the Vglut2 cloud would hide exactly the comparison
%            it exists to make. Set includeNewInAll = false to leave it out.
%
%  OUTPUT  <outDir>\  four .png + four .pdf + two .csv + one .mat
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps,
%                coh_cfg_260727, cell_merge_overrides_260814
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));   % coh_cfg_260727
addpath(fullfile(repoRoot, 'analysis_260806'));                    % cell_merge_overrides

%% ===================== USER-EDITABLE PARAMETERS ======================
% ---- the new (vagotomised) session ----
newRoot       = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
newCuratedCsv = fullfile(newRoot, 'roi_match_out_260824', 'roi_match_curated.csv');
newGroup      = 'Vglut2-vagx';     % its own group label / colour
newDate       = '0824';

% ---- the existing archive ----
archiveRoot   = 'D:\Ventral_surface_summary';
scan_dirs     = {'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};

outDir        = fullfile(newRoot, 'polar_coh_rayleigh_260824');

includeNewInAll     = true;   % false -> figures 3/4 are the archive alone
applyMergeOverrides = true;   % fold in cell_merge_overrides_260814 (archive only)

% group label -> colour. IO is split out of ChAT: any archive folder whose
% cell-level dir name contains 'IO' counts as IO.
groups       = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert', newGroup};
group_colors = [0    0    0   ;     % IO           black (hollow marker)
                0.85 0.10 0.10;     % ChAT         red
                0.10 0.65 0.20;     % Vglut2       green
                0.10 0.30 0.85;     % Vgat         blue
                0.55 0.20 0.75;     % Sst          purple
                0.90 0.45 0.10;     % Sert         orange
                0.90 0.10 0.60];    % Vglut2-vagx  magenta

% cell registries for the archive sessions that have one. Sessions absent from
% this list fall back to one cell per ROI, which is a no-op for pooling.
cell_link_sources = { ...
    'Sert',   '0721', fullfile(archiveRoot, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(archiveRoot, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(archiveRoot, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0806', fullfile(archiveRoot, 'Sst',    '0806', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0807', fullfile(archiveRoot, 'Sst',    '0807', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0810', fullfile(archiveRoot, 'Vglut2', '0810', 'cell_pooled', 'cell_link.mat')
    };

% Recordings barred from anything pooled (different stage zero).
excludeRecordings = coh_cfg_260727().excludeRecordings;

% ---- analysis parameters: identical to Ventral_surface_polar_coh_vs_rayleigh_260808 ----
nDrop           = 30;        % breath frames tossed up front (align to Ca)
fallback_fps    = 30;
minSpikes       = 2;         % include ROI if it has >= this many events
TW              = 4;         % multitaper TW for coherence
alpha_sig       = 0.001;     % primary coherence level; also the jackknife CI level
alpha_sig2      = 0.05;      % secondary threshold (outer dashed circle)
ca_lag_sec      = 0.03;      % GCaMP8s lead: events shifted EARLIER before scoring

f_breath_search = [0.2 4];   % Hz, search band for the breath PSD peak
fwhm_factor     = 0.6;       % coherence band = fwhm_factor x FWHM
min_bw          = 0.05;      % Hz, minimum coherence band width
fmin            = 0.05;      % Hz, PSD lower bound
fmax            = 15;        % Hz, PSD upper bound

peak_file_name  = 'breath_peak_pc1.mat';
start_file_name = 'breath_insp_start_pc1.mat';

nPhaseBins      = 36;        % 10 deg bins for the occupancy histogram
alpha_ray       = 0.05;      % inner dashed circle -> log(-log(alpha))
alpha_ray2      = 0.001;     % outer dashed circle
nShuffleRay     = 500;       % circular-shift shuffles for the independent p
rayRLimPad      = 0.5;

% ---- how much error-bar ink ----
% Spans are drawn for EVERY point, graded by significance: faint for points
% inside the threshold circle, solid for points outside it. A non-significant
% point has a huge span almost by definition, so drawing them all at full
% strength would bury the tuned cells; dropping them entirely would hide that
% most of the cloud is unconstrained. Set arcAlphaNS = 0 to suppress them.
arcAlphaSig  = 0.75;   arcLWSig  = 1.3;
arcAlphaNS   = 0.16;   arcLWNS   = 0.6;
markerSize   = 6;

doSave = true;
% =====================================================================

logZ_crit  = log(-log(alpha_ray));      % alpha=0.05  -> Z=3.00 -> logZ=1.10
logZ_crit2 = log(-log(alpha_ray2));     % alpha=0.001 -> Z=6.91 -> logZ=1.93

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
% NVIDIA WGL throws inside exportgraphics ~1% of the time on this machine and
% silently writes a 9 KB blank. See reference_matlab_opengl_export_stubs.
try, opengl('software'); catch, end
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% =================== BUILD THE RECORDING JOB LIST ====================
% One job = one recording folder, tagged with its group, date and dataset scope.
% Both layouts (flat session, nested archive) collapse to this list, so the
% compute loop below is written once.
JOB = struct('folder',{},'group',{},'date',{},'recName',{},'scan',{},'isNew',{});

% ---- the new flat session ----
assert(isfile(newCuratedCsv), 'curation not found: %s', newCuratedCsv);
Tnew = readtable(newCuratedCsv);
reqCol = {'cell_id','fov_name','roi_index'};
for k = 1:numel(reqCol)
    assert(ismember(reqCol{k}, Tnew.Properties.VariableNames), ...
           '%s has no column "%s"', newCuratedCsv, reqCol{k});
end
newFovs = unique(Tnew.fov_name, 'stable');
nNoSpk = 0;
for k = 1:numel(newFovs)
    fp = fullfile(newRoot, newFovs{k});
    if ~isfile(fullfile(fp,'ca_spike_data.mat'))
        fprintf(2,'  [new] no ca_spike_data.mat, FOV dropped: %s\n', newFovs{k});
        nNoSpk = nNoSpk + 1; continue;
    end
    JOB(end+1) = struct('folder',fp,'group',newGroup,'date',newDate, ...
                        'recName',newFovs{k},'scan','NEW','isNew',true); %#ok<SAGROW>
end
fprintf('[new] %d curated FOVs, %d usable (%d without spike detection)\n', ...
        numel(newFovs), numel(JOB), nNoSpk);

% (fov,roi) -> cell key for the new session, straight from the curation
newCellMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:height(Tnew)
    fn = Tnew.fov_name{i};
    newCellMap(sprintf('%s|%d', fn, Tnew.roi_index(i))) = ...
        sprintf('%s/%s#c%d', newGroup, newDate, Tnew.cell_id(i));
end

% ---- the archive ----
nArchJob = 0;
if isfolder(archiveRoot)
    for sg = 1:numel(scan_dirs)
        gdir = fullfile(archiveRoot, scan_dirs{sg});
        if ~isfolder(gdir), warning('Folder missing: %s', gdir); continue; end
        allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
        for kk = 1:numel(allMat)
            fp = allMat(kk).folder;
            rn = folder_basename(fp);
            if any(strcmp(rn, excludeRecordings)), continue; end
            gname = scan_dirs{sg};
            if is_io_path(fp, gdir), gname = 'IO'; end
            JOB(end+1) = struct('folder',fp,'group',gname, ...
                                'date',date_from_path(fp,gdir),'recName',rn, ...
                                'scan',scan_dirs{sg},'isNew',false); %#ok<SAGROW>
            nArchJob = nArchJob + 1;
        end
    end
    fprintf('[archive] %d recordings across %s\n', nArchJob, strjoin(scan_dirs,', '));
else
    warning('Archive root missing: %s -- figures 3/4 will show the new session alone.', archiveRoot);
end

% ---- archive cell identity ----
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
for s = 1:size(cell_link_sources,1)
    pre = sprintf('%s/%s', cell_link_sources{s,1}, cell_link_sources{s,2});
    if ~isfile(cell_link_sources{s,3})
        warning('cell_link missing for %s -- that session stays per-ROI.', pre); continue;
    end
    Lk = load(cell_link_sources{s,3},'link');  Tk = Lk.link.obsT;
    for i = 1:height(Tk)
        kk = sprintf('%s/%s/%d', pre, Tk.rec_name(i), Tk.maskL_label(i));
        if isnan(Tk.cell_id(i)), tossed_set(kk) = true;
        else,                    cell_map(kk) = sprintf('%s#c%d', pre, Tk.cell_id(i));
        end
    end
    fprintf('  cell identity %-14s %4d masks -> %4d cells (%d tossed)\n', pre, height(Tk), ...
            numel(unique(Tk.cell_id(~isnan(Tk.cell_id)))), nnz(isnan(Tk.cell_id)));
end

%% ========================= SCAN AND COMPUTE ==========================
PP  = init_coll();
REC = struct('bidx',{},'occ',{},'T',{},'phiW',{});
OBS = struct('rec',{},'ev',{},'gi',{},'f_pk',{},'label',{},'cellKey',{}, ...
             'recName',{},'roi',{},'isNew',{});
labels = {};
confC  = NaN;
K_tap  = 2*TW - 1;
confC2 = sqrt(1 - alpha_sig2 ^ (1/(K_tap - 1)));      % Chronux closed form
edges  = linspace(0, 2*pi, nPhaseBins+1);

fprintf('\n=========== scanning %d recordings ===========\n', numel(JOB));
tScan = tic;
for j = 1:numel(JOB)
    J = JOB(j);
    try
        bpFile = fullfile(J.folder, peak_file_name);
        ipFile = fullfile(J.folder, start_file_name);
        if ~isfile(bpFile), fprintf('  skip (no %s): %s\n', peak_file_name, J.recName); continue; end
        if ~isfile(ipFile), fprintf('  skip (no %s): %s\n', start_file_name, J.recName); continue; end

        gi = find(strcmp(groups, J.group), 1);
        if isempty(gi), warning('No colour for group "%s" -- skipping %s', J.group, J.recName); continue; end

        fps  = detect_session_fps(J.folder, fallback_fps);
        CA   = load(fullfile(J.folder, 'ca_spike_data.mat'));
        nROI = numel(CA.roi_spikes);
        nCa  = numel(CA.roi_spikes(1).spike_train);

        BP = load(bpFile);     % SVD PC1 peaks (insp_onset_idx, breath)
        IP = load(ipFile);     % SVD PC1 feet  (insp_start_idx)
        bw = detrend(double(BP.breath(:)));
        bw(1:min(nDrop,numel(bw))) = [];
        bw = bw - mean(bw);

        peak_idx = round(BP.insp_onset_idx(:)) - nDrop;     % insp PEAK
        foot_idx = round(IP.insp_start_idx(:)) - nDrop;     % insp ONSET (foot)

        % Vglut2/1124 used a rising-edge trigger: breath sits 1 frame late.
        if strcmpi(J.scan,'Vglut2') && strcmp(J.date,'1124')
            peak_idx = peak_idx + 1;
            foot_idx = foot_idx + 1;
            bw = [bw(1); bw(1:end-1)];
        end

        T = min([numel(bw), nCa]);
        peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
        foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
        if numel(peak_idx) < 2 || numel(foot_idx) < 2
            fprintf('  skip (too few breath landmarks): %s\n', J.recName); continue;
        end
        phi = piecewise_phase_local(peak_idx, foot_idx, T);
        ref = cos(phi); ref(isnan(ref)) = 0;
        ref = ref - mean(ref(:));
        bw  = bw(1:T);

        % breath waveform PSD -> coherence band
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(bw, pB); Sb=Sb(:); fb=fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        [~,rl]=max(Sb(m)); ipk=find(m,1)+rl-1; f_pk=fb(ipk);
        h=Sb(ipk)/2; lo=ipk; while lo>1&&Sb(lo)>h, lo=lo-1; end
        hi=ipk;            while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

        pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
        pc.fpass=band; pc.err=[2,alpha_sig];        % err=[2 alpha] -> JACKKNIFE

        % occupancy: frames per phase bin, shared by every ROI in this recording
        phiW  = mod(phi, 2*pi);
        valid = ~isnan(phiW);
        bidx  = nan(T,1);
        bidx(valid) = discretize(phiW(valid), edges);
        occ   = accumarray(bidx(valid), 1, [nPhaseBins 1]);
        REC(end+1) = struct('bidx',bidx,'occ',occ,'T',T,'phiW',phiW); %#ok<SAGROW>
        recIdx = numel(REC);

        nInc = 0; nUncurated = 0;
        lag  = round(ca_lag_sec*fps);
        for rid = 1:nROI
            st = double(CA.roi_spikes(rid).spike_train(:));
            st = st(1:min(T,numel(st)));
            if numel(st)<T, st(end+1:T)=0; end
            if sum(st) < minSpikes, continue; end

            if J.isNew
                ckin = sprintf('%s|%d', J.recName, rid);
                if ~isKey(newCellMap, ckin), nUncurated = nUncurated + 1; continue; end
                ck = newCellMap(ckin);
            else
                ckin = sprintf('%s/%s/%s/%d', J.scan, J.date, J.recName, rid);
                if isKey(tossed_set, ckin), continue; end       % rejected in curation
                if isKey(cell_map, ckin), ck = cell_map(ckin);
                else,                     ck = sprintf('roi:%s', ckin); end
            end
            nInc = nInc + 1;
            stL = [st(1+lag:end); zeros(lag,1)];    % lead-shift earlier by lag frames
            [PP, confC] = add_coh(PP, ref, stL - mean(stL), pc, band, gi, confC, f_pk);

            lab = sprintf('%s/%s/%s/%d', J.group, J.date, J.recName, rid);
            OBS(end+1) = struct('rec',recIdx,'ev',find(stL > 0),'gi',gi, ...
                                'f_pk',f_pk,'label',lab,'cellKey',ck, ...
                                'recName',J.recName,'roi',rid,'isNew',J.isNew); %#ok<SAGROW>
            labels{end+1} = lab; %#ok<SAGROW>
        end
        if nUncurated > 0
            extra = sprintf('  (%d not in curation)', nUncurated);
        else
            extra = '';
        end
        fprintf('  [%3d/%3d] %-12s %-46s band [%.2f %.2f] Hz  %d/%d ROI%s\n', ...
                j, numel(JOB), J.group, J.recName, band(1), band(2), nInc, nROI, extra);
    catch ME
        warning('  ERROR %s: %s', J.recName, ME.message);
    end
end
if isempty(PP.r), error('No ROIs collected.'); end
fprintf('scan done in %.1f min: %d ROI-observations\n', toc(tScan)/60, numel(OBS));

%% ================ HAND-CURATED MERGES (archive only) =================
% cell_link's gate is deliberately tight, so genuine cross-FOV matches at
% different zoom or depth can fall outside it. Those are corrected by eye in
% cell_merge_overrides_260814 and applied here, keyed on (recName, roi) which
% survives any renumbering of the registry.
if applyMergeOverrides
    G = cell_merge_overrides_260814();
    keysAll = {OBS.cellKey};
    recroi  = arrayfun(@(o) sprintf('%s|%d', o.recName, o.roi), OBS, 'uni', 0);
    nDone = 0;
    for gi = 1:numel(G)
        grp = G{gi}; ks = {}; ok = true;
        for m = 1:numel(grp)
            kq = sprintf('%s|%d', grp{m}{1}, grp{m}{2});
            q  = find(strcmp(recroi, kq), 1);
            if isempty(q), ok = false; break; end
            ks{end+1} = keysAll{q}; %#ok<SAGROW>
        end
        if ~ok
            fprintf('  merge group %d: a member is not in this scan -- skipped\n', gi); continue;
        end
        ks = unique(ks);
        if numel(ks) < 2, continue; end          % already one cell
        for m = 2:numel(ks), keysAll(strcmp(keysAll, ks{m})) = ks(1); end
        nDone = nDone + 1;
    end
    [OBS.cellKey] = keysAll{:};
    fprintf('applied %d of %d hand-curated merges\n', nDone, numel(G));
end

%% ===================== POOL THE RAYLEIGH PER CELL ====================
RY = init_ray();
cellKeys = {OBS.cellKey};
[uCell, ~, obsOfCell] = unique(cellKeys, 'stable');
cellLabels = cell(numel(uCell),1);  cellNObs = zeros(numel(uCell),1);
cellIsNew  = false(numel(uCell),1);
fprintf('\npooling %d observations -> %d cells ...\n', numel(OBS), numel(uCell));
tPool = tic;
for c = 1:numel(uCell)
    m = find(obsOfCell == c);
    cellNObs(c)  = numel(m);
    cellIsNew(c) = OBS(m(1)).isNew;
    A = []; W = []; parts = struct('ev',{},'rec',{});
    for j = m(:)'
        r  = REC(OBS(j).rec);
        ev = OBS(j).ev;
        ev = ev(ev>=1 & ev<=r.T);
        ev = ev(~isnan(r.bidx(ev)));
        if isempty(ev), continue; end
        A = [A; r.phiW(ev(:))];                       %#ok<AGROW>
        W = [W; 1./max(r.occ(r.bidx(ev(:))),1)];      %#ok<AGROW>
        parts(end+1) = struct('ev',ev,'rec',OBS(j).rec); %#ok<AGROW>
    end
    RY = add_ray_pooled(RY, A, W, parts, REC, nShuffleRay, OBS(m(1)).gi, OBS(m(1)).f_pk);
    if cellNObs(c) == 1
        cellLabels{c} = OBS(m(1)).label;
    else
        cellLabels{c} = sprintf('%s  [+%d more rec]', OBS(m(1)).label, cellNObs(c)-1);
    end
end
fprintf('pooled in %.1f min (%d cells span >1 recording)\n', toc(tPool)/60, nnz(cellNObs>1));

obsIsNew = [OBS.isNew]';
obsGi    = [OBS.gi]';
gOfCell  = zeros(numel(uCell),1);
for c = 1:numel(uCell), gOfCell(c) = OBS(find(obsOfCell==c,1)).gi; end
for gi = 1:numel(groups)
    if ~any(obsGi == gi), continue; end
    fprintf('   %-12s %4d observations -> %4d cells\n', groups{gi}, ...
            nnz(obsGi==gi), nnz(gOfCell==gi));
end

%% ============================ FOUR FIGURES ===========================
scopes = struct( ...
    'tag',   {'new', 'all'}, ...
    'obsM',  {obsIsNew, true(size(obsIsNew))}, ...
    'celM',  {cellIsNew, true(size(cellIsNew))}, ...
    'name',  {sprintf('%s  %s (vagotomised)', newGroup, newDate), 'all datasets'});
if ~includeNewInAll
    scopes(2).obsM = ~obsIsNew;
    scopes(2).celM = ~cellIsNew;
    scopes(2).name = 'archive (vagotomised session excluded)';
end

figFiles = {};
for s = 1:numel(scopes)
    S  = scopes(s);
    om = S.obsM(:);  cm = S.celM(:);
    if ~any(om), fprintf(2,'scope "%s" is empty -- skipped\n', S.tag); continue; end

    % ---------- COHERENCE ----------
    f1 = figure('Color','w','Units','centimeters','Position',[2 2 19 16], ...
                'Name',sprintf('Coherence -- %s', S.name));
    set(f1,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    ax1 = polaraxes(f1,'Position',[0.06 0.04 0.62 0.70]);
    plot_coh_polar(ax1, PP, om, group_colors, confC, confC2, ...
                   arcAlphaSig, arcLWSig, arcAlphaNS, arcLWNS, markerSize);
    title(ax1, {'r = mean |C| in breath band,   \theta = coherence phase', ...
                'arc = jackknife 95% phase CI,   radial bar = jackknife |C| CI (significant only)'}, ...
          'Interpreter','tex');
    nObsG = arrayfun(@(g) nnz(om & obsGi==g), 1:numel(groups));
    nSigG = arrayfun(@(g) nnz(om & obsGi==g & PP.r(:)>=confC), 1:numel(groups));
    add_group_legend(f1, ax1, groups, group_colors, nObsG, nSigG, 'sig');
    sgtitle(f1, {sprintf('COHERENCE  |  %s', S.name), ...
        sprintf(['N = %d ROI-observations (recording \\times ROI) from %d cells   |   ' ...
                 'insp onset = 0,  insp peak = \\pi'], nnz(om), nnz(cm)), ...
        sprintf(['confC(\\alpha=%.3g) = %.2f \\rightarrow %d significant   ' ...
                 '(outer dashed: \\alpha=%.2g, confC = %.2f)   |   TW = %d'], ...
                 alpha_sig, confC, nnz(om & PP.r(:)>=confC), alpha_sig2, confC2, TW)}, ...
        'FontSize', 8.5);
    annotate_lag(f1, ca_lag_sec);
    figFiles{end+1} = save_fig(f1, outDir, sprintf('polar_coherence_%s', S.tag), doSave); %#ok<SAGROW>

    % ---------- RAYLEIGH ----------
    f2 = figure('Color','w','Units','centimeters','Position',[2 2 19 16], ...
                'Name',sprintf('Rayleigh -- %s', S.name));
    set(f2,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    ax2 = polaraxes(f2,'Position',[0.06 0.04 0.62 0.70]);
    plot_ray_polar(ax2, RY, cm, group_colors, logZ_crit, logZ_crit2, rayRLimPad, ...
                   arcAlphaSig, arcLWSig, arcAlphaNS, arcLWNS, markerSize);
    title(ax2, {'r = log Z,   \theta = circular MEDIAN event phase', ...
                'arc = occupancy-weighted circular IQR about that median'}, ...
          'Interpreter','tex');
    nCelG = arrayfun(@(g) nnz(cm & gOfCell==g), 1:numel(groups));
    nZ2G  = arrayfun(@(g) nnz(cm & gOfCell==g & RY.logZ(:)>=2), 1:numel(groups));
    add_group_legend(f2, ax2, groups, group_colors, nCelG, nZ2G, 'logZ>=2');
    sgtitle(f2, {sprintf('RAYLEIGH  |  %s', S.name), ...
        sprintf(['N = %d cells (%d ROI-observations pooled)   |   ' ...
                 'insp onset = 0,  insp peak = \\pi'], nnz(cm), nnz(om)), ...
        sprintf(['logZ crit(\\alpha=%.2g) = %.2f \\rightarrow %d significant ' ...
                 '(%d also pass the circular-shift shuffle)'], ...
                 alpha_ray, logZ_crit, nnz(cm & RY.logZ(:)>=logZ_crit), ...
                 nnz(cm & RY.logZ(:)>=logZ_crit & RY.p_sh(:)<alpha_ray))}, ...
        'FontSize', 8.5);
    annotate_lag(f2, ca_lag_sec);
    figFiles{end+1} = save_fig(f2, outDir, sprintf('polar_rayleigh_%s', S.tag), doSave); %#ok<SAGROW>

    fprintf('\n--- %s ---\n', S.name);
    fprintf('  coherence : %d / %d observations >= confC\n', nnz(om & PP.r(:)>=confC), nnz(om));
    fprintf('  Rayleigh  : %d / %d cells >= logZ crit,  %d also pass the shuffle\n', ...
            nnz(cm & RY.logZ(:)>=logZ_crit), nnz(cm), ...
            nnz(cm & RY.logZ(:)>=logZ_crit & RY.p_sh(:)<alpha_ray));
    fprintf('  median logZ = %.2f,  median |C| = %.3f\n', ...
            median(RY.logZ(cm)), median(PP.r(om)));
end

%% ============================== SAVE =================================
if doSave
    fid = fopen(fullfile(outDir,'polar_scores_percell.csv'),'w');
    fprintf(fid,['dataset,group,cell_key,n_recordings,first_label,n_events,n_eff,f_breath_Hz,' ...
                 'ray_Rbar,ray_Z,ray_logZ,ray_mean_deg,ray_med_deg,' ...
                 'ray_iqr_lo_deg,ray_iqr_hi_deg,ray_iqr_width_deg,' ...
                 'ray_p_analytic,ray_p_shuffle,ray_sig\n']);
    for c = 1:numel(uCell)
        if cellIsNew(c), dsTag = 'new'; else, dsTag = 'archive'; end
        fprintf(fid,'%s,%s,%s,%d,"%s",%d,%.2f,%.4f,%.6f,%.4f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.6g,%.4f,%d\n', ...
            dsTag, groups{RY.colorIdx(c)}, uCell{c}, ...
            cellNObs(c), cellLabels{c}, RY.nEv(c), RY.nEff(c), RY.f_pk(c), ...
            RY.Rbar(c), RY.Z(c), RY.logZ(c), rad2deg(RY.th(c)), rad2deg(RY.thMed(c)), ...
            rad2deg(RY.thMed(c)+RY.q25(c)), rad2deg(RY.thMed(c)+RY.q75(c)), ...
            rad2deg(RY.q75(c)-RY.q25(c)), RY.p_an(c), RY.p_sh(c), ...
            RY.logZ(c)>=logZ_crit);
    end
    fclose(fid);

    fid = fopen(fullfile(outDir,'polar_scores_perobs.csv'),'w');
    fprintf(fid,'dataset,group,label,cell_key,coh_r,coh_r_lo,coh_r_hi,coh_th_deg,coh_dphi_deg,coh_sig\n');
    for ii = 1:numel(labels)
        if OBS(ii).isNew, dsTag = 'new'; else, dsTag = 'archive'; end
        fprintf(fid,'%s,%s,%s,%s,%.6f,%.6f,%.6f,%.2f,%.2f,%d\n', ...
            dsTag, groups{PP.colorIdx(ii)}, labels{ii}, OBS(ii).cellKey, ...
            PP.r(ii), PP.rlo(ii), PP.rhi(ii), ...
            rad2deg(PP.th(ii)), rad2deg(PP.dphi(ii)), PP.r(ii)>=confC);
    end
    fclose(fid);

    save(fullfile(outDir,'polar_coh_rayleigh_260824_data.mat'), ...
         'PP','RY','labels','uCell','cellLabels','cellNObs','cellIsNew','obsIsNew', ...
         'obsOfCell','gOfCell','confC','confC2','groups','group_colors', ...
         'logZ_crit','logZ_crit2','TW','alpha_sig','alpha_sig2','alpha_ray','alpha_ray2', ...
         'nPhaseBins','nShuffleRay','minSpikes','nDrop','ca_lag_sec', ...
         'includeNewInAll','applyMergeOverrides','newRoot','archiveRoot');
    fprintf('\nSaved to %s\n', outDir);
    for k = 1:numel(figFiles), fprintf('   %s\n', figFiles{k}); end
end
fprintf('Done.\n');

%% ========================= LOCAL FUNCTIONS ==========================
function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[],'colorIdx',[],'f_pk',[]);
end

function C = init_ray()
    C = struct('th',[],'Rbar',[],'Z',[],'logZ',[],'p_an',[],'p_sh',[], ...
               'nEv',[],'nEff',[],'colorIdx',[],'f_pk',[],'thMed',[],'q25',[],'q75',[]);
end

function fn = save_fig(fh, outDir, stem, doSave)
    fn = fullfile(outDir, [stem '.png']);
    if ~doSave, return; end
    exportgraphics(fh, fn, 'Resolution',300, 'BackgroundColor','white');
    exportgraphics(fh, fullfile(outDir,[stem '.pdf']), 'ContentType','vector', ...
                   'BackgroundColor','white');
    d = dir(fn);
    if isempty(d) || d.bytes < 30e3     % OpenGL stub check
        fprintf(2,'  suspect PNG for %s -- retrying once\n', stem);
        exportgraphics(fh, fn, 'Resolution',300, 'BackgroundColor','white');
    end
end

function annotate_lag(fh, ca_lag_sec)
% The lead correction changes every reported angle, so it is stated ON the
% figure rather than left for a reader to guess.
    if ca_lag_sec == 0
        txt = 'NO GCaMP lead correction';
    else
        txt = sprintf('GCaMP lead correction: events shifted %.0f ms EARLIER', ca_lag_sec*1000);
    end
    annotation(fh, 'textbox', [0.02 0.005 0.96 0.045], 'String', txt, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'EdgeColor',[0.6 0.6 0.6],'BackgroundColor',[0.96 0.96 0.96], ...
        'FontSize',8.5,'FontWeight','bold','Margin',2,'Interpreter','tex');
end

function plot_coh_polar(ax, C, mask, group_colors, confC, confC2, aS, lwS, aN, lwN, ms)
% radius = |C|, angle = coherence phase, angular span = jackknife phase CI.
    hold(ax,'on');
    thc = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(confC2,1,360), '--','Color',[0.5 0.5 0.5],'LineWidth',0.8);
    polarplot(ax, thc, repmat(confC, 1,360), 'k--','LineWidth',1);
    text(ax, deg2rad(325), confC, sprintf('confC = %.2f', confC), ...
         'FontSize',7.5,'Color','k','HorizontalAlignment','center', ...
         'BackgroundColor','w','Margin',0.5);

    idx = find(mask(:)');
    % spans first, dots on top
    for k = idx
        ci = C.colorIdx(k);
        if isnan(C.th(k)) || isnan(C.r(k)) || ci<1 || ci>size(group_colors,1), continue; end
        col = group_colors(ci,:);
        sig = C.r(k) >= confC;
        if sig, aa = aS; lw = lwS; else, aa = aN; lw = lwN; end
        if aa <= 0, continue; end
        if ~isnan(C.dphi(k))
            % clamp: a jackknife CI wider than the circle carries no direction
            dp  = min(C.dphi(k), pi);
            arc = linspace(C.th(k)-dp, C.th(k)+dp, 40);
            polarplot(ax, arc, C.r(k)*ones(size(arc)), '-','Color',[col aa],'LineWidth',lw);
        end
        if sig
            polarplot(ax, [C.th(k) C.th(k)], [C.rlo(k) C.rhi(k)], ...
                      '-','Color',[col aa],'LineWidth',lw);
        end
    end
    for k = idx
        ci = C.colorIdx(k);
        if isnan(C.th(k)) || isnan(C.r(k)) || ci<1 || ci>size(group_colors,1), continue; end
        col = group_colors(ci,:);
        if ci == 1     % IO: hollow black
            polarplot(ax, C.th(k), C.r(k), 'o','MarkerFaceColor','none', ...
                      'MarkerEdgeColor','k','MarkerSize',ms,'LineWidth',0.8);
        else
            polarplot(ax, C.th(k), C.r(k), 'o','MarkerFaceColor',col, ...
                      'MarkerEdgeColor','k','MarkerSize',ms,'LineWidth',0.4);
        end
    end
    ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
end

function plot_ray_polar(ax, C, mask, group_colors, lzc, lzc2, pad, aS, lwS, aN, lwN, ms)
% radius = log Z, angle = circular MEDIAN, angular span = weighted circular IQR.
% log Z <= 0 (Z < 1) is clamped to the centre: a negative radius would be
% silently reflected by polarplot onto the OPPOSITE angle, i.e. the wrong phase.
    hold(ax,'on');
    idx  = find(mask(:)');
    rmax = max([C.logZ(idx); lzc2]) + pad;
    thc  = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(lzc2,1,360), '--','Color',[0.5 0.5 0.5],'LineWidth',0.8);
    polarplot(ax, thc, repmat(lzc, 1,360), 'k--','LineWidth',1);
    thLab = deg2rad(325);
    text(ax, thLab, lzc,  sprintf('\\alpha=%.2g', exp(-exp(lzc))), 'FontSize',7.5, ...
         'Color','k','HorizontalAlignment','center','BackgroundColor','w','Margin',0.5);
    text(ax, thLab, lzc2, sprintf('\\alpha=%.3g', exp(-exp(lzc2))), 'FontSize',7.5, ...
         'Color',[0.35 0.35 0.35],'HorizontalAlignment','center','BackgroundColor','w','Margin',0.5);

    for k = idx
        ci = C.colorIdx(k);
        if isnan(C.thMed(k)) || ci<1 || ci>size(group_colors,1), continue; end
        if isnan(C.q25(k)) || isnan(C.q75(k)), continue; end
        col = group_colors(ci,:);  if ci == 1, col = [0 0 0]; end
        if C.logZ(k) >= lzc, aa = aS; lw = lwS; else, aa = aN; lw = lwN; end
        if aa <= 0, continue; end
        r   = max(C.logZ(k), 0);
        arc = linspace(C.thMed(k)+C.q25(k), C.thMed(k)+C.q75(k), 40);
        polarplot(ax, arc, r*ones(size(arc)), '-','Color',[col aa],'LineWidth',lw);
    end
    for k = idx
        ci = C.colorIdx(k);
        if isnan(C.thMed(k)) || ci<1 || ci>size(group_colors,1), continue; end
        r   = max(C.logZ(k), 0);
        col = group_colors(ci,:);
        if ci == 1
            polarplot(ax, C.thMed(k), r, 'o','MarkerFaceColor','none', ...
                      'MarkerEdgeColor','k','MarkerSize',ms,'LineWidth',0.8);
        else
            polarplot(ax, C.thMed(k), r, 'o','MarkerFaceColor',col, ...
                      'MarkerEdgeColor','k','MarkerSize',ms,'LineWidth',0.4);
        end
    end
    ax.RLim=[0 rmax]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
end

function add_group_legend(fig, refAx, groups, group_colors, nTot, nSig, sigWord)
% Colour key: group -> hue, count, and how many cleared the threshold.
% Groups absent from this figure are omitted rather than listed with n = 0.
    pos = refAx.Position;
    legAx = axes(fig, 'Position', [pos(1)+pos(3)+0.01, pos(2)+0.30*pos(4), 0.22, 0.40*pos(4)]);
    hold(legAx,'on'); axis(legAx,'off');
    show = find(nTot > 0);
    xlim(legAx,[0 1]); ylim(legAx,[0 1]);
    if isempty(show), return; end
    y = linspace(0.90, 0.10, max(numel(show),2));
    for i = 1:numel(show)
        k = show(i);
        if k == 1
            plot(legAx, 0.08, y(i), 'o','MarkerFaceColor','none', ...
                 'MarkerEdgeColor','k','MarkerSize',7,'LineWidth',0.8);
        else
            plot(legAx, 0.08, y(i), 'o','MarkerFaceColor',group_colors(k,:), ...
                 'MarkerEdgeColor','k','MarkerSize',7,'LineWidth',0.4);
        end
        text(legAx, 0.22, y(i), sprintf('%s (n=%d, %s: %d)', groups{k}, nTot(k), sigWord, nSig(k)), ...
             'FontSize',8,'Interpreter','none');
    end
end

function [C, confC] = add_coh(C, x, y, pc, band, cidx, confC, f_pk)
% Chronux band-averaged magnitude + circular-mean phase + JACKKNIFE CIs.
% pc.err = [2 alpha] selects the jackknife estimator, so phistd and Cerr are
% both jackknife quantities -- that is what the figure's spans are.
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f  = f(:);
    mb = f>=band(1) & f<=band(2);
    if ~any(mb), mb = true(size(f)); end
    C.th(end+1,1)   = angle(mean(exp(1i*(-phi(mb)))));
    C.r(end+1,1)    = mean(Cxy(mb));
    C.rlo(end+1,1)  = max(0, mean(Cerr(1,mb)));
    C.rhi(end+1,1)  = min(1, mean(Cerr(2,mb)));
    C.dphi(end+1,1) = 1.96*mean(phistd(mb));
    C.colorIdx(end+1,1) = cidx;
    C.f_pk(end+1,1) = f_pk;
end

function C = add_ray_pooled(C, a, w, parts, REC, nSh, cidx, f_pk)
% Rayleigh for ONE CELL, pooled over every recording it appears in.
%   a  event phases (rad, 0 = insp onset, pi = insp peak)
%   w  event weights 1/(frames in that event's phase bin) IN ITS OWN recording
% The null shifts EACH recording's event train independently and re-pools, so it
% preserves both the per-recording occupancy and each train's own bursting.
    n = numel(a);
    if n < 1
        C.th(end+1,1)=NaN; C.Rbar(end+1,1)=0; C.Z(end+1,1)=0; C.logZ(end+1,1)=log(eps);
        C.p_an(end+1,1)=1; C.p_sh(end+1,1)=NaN; C.nEv(end+1,1)=0; C.nEff(end+1,1)=0;
        C.colorIdx(end+1,1)=cidx; C.f_pk(end+1,1)=f_pk;
        C.thMed(end+1,1)=NaN; C.q25(end+1,1)=NaN; C.q75(end+1,1)=NaN;
        return;
    end
    [th, Rbar, nEff] = wresultant(a, w);
    Z = nEff * Rbar^2;
    p_an = exp(-Z) * (1 + (2*Z - Z^2)/(4*nEff) ...
                        - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*nEff^2));
    p_an = min(max(p_an,0),1);

    if nSh > 0 && ~isempty(parts)
        Znull = zeros(nSh,1);
        for s = 1:nSh
            as = []; ws = [];
            for q = 1:numel(parts)
                r  = REC(parts(q).rec);
                sh = round(s * r.T / (nSh+1)) + q;      % a different shift per recording
                evs = mod(parts(q).ev - 1 + sh, r.T) + 1;
                evs = evs(~isnan(r.bidx(evs)));
                if isempty(evs), continue; end
                as = [as; r.phiW(evs(:))];                    %#ok<AGROW>
                ws = [ws; 1./max(r.occ(r.bidx(evs(:))),1)];   %#ok<AGROW>
            end
            if isempty(as), Znull(s) = 0; continue; end
            [~, Rs, nEs] = wresultant(as, ws);
            Znull(s) = nEs * Rs^2;
        end
        p_sh = (1 + nnz(Znull >= Z)) / (1 + nSh);
    else
        p_sh = NaN;
    end

    C.th(end+1,1)=th;  C.Rbar(end+1,1)=Rbar;  C.Z(end+1,1)=Z;
    C.logZ(end+1,1)=log(max(Z,eps));  C.p_an(end+1,1)=p_an;  C.p_sh(end+1,1)=p_sh;
    C.nEv(end+1,1)=n;  C.nEff(end+1,1)=nEff;
    C.colorIdx(end+1,1)=cidx;  C.f_pk(end+1,1)=f_pk;

    if n >= 2
        med = circ_median_w(a, w);
        u   = angle(exp(1i*(a - med)));         % offsets from the median, (-pi pi]
        C.thMed(end+1,1) = med;
        C.q25(end+1,1)   = wprctile(u, w, 25);
        C.q75(end+1,1)   = wprctile(u, w, 75);
    else
        C.thMed(end+1,1)=NaN; C.q25(end+1,1)=NaN; C.q75(end+1,1)=NaN;
    end
end

function [th, Rbar, nEff] = wresultant(a, w)
% Weighted circular resultant with Kish's effective sample size. The analytic
% Rayleigh Z = n*Rbar^2 assumes n equally weighted angles; feeding it a weighted
% Rbar with raw n gave p = 5e-79 where the shuffle said p = 0.05.
    S1 = sum(w);  S2 = sum(w.^2);
    if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
    nEff = S1^2 / max(S2, eps);
    v    = sum(w(:) .* exp(1i*a(:))) / S1;
    th   = angle(v);
    Rbar = min(abs(v), 1);
end

function m = circ_median_w(a, w)
% Weighted circular median: the direction minimising the total weighted angular
% distance to the data. Solved on a 0.5 deg grid -- there is no closed form on
% the circle, and the grid is far finer than the 10 deg occupancy histogram.
    g = linspace(-pi, pi, 721);  g(end) = [];
    d = abs(angle(exp(1i*(a(:) - g))));
    [~, k] = min(sum(w(:) .* d, 1));
    m = g(k);
end

function q = wprctile(v, w, p)
% Weighted percentile, midpoint rule on the cumulative weight.
    [v, o] = sort(v(:));  w = w(o);  w = w / sum(w);
    c = cumsum(w) - 0.5*w;
    if numel(v) < 2, q = v(1); return; end
    [c, iu] = unique(c);  v = v(iu);
    q = interp1(c, v, p/100, 'linear', 'extrap');
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase: FEET at 0/2pi/..., PEAKS at pi/3pi/..., linear in time
% between consecutive landmarks. NaN outside the first..last landmark.
    phi = nan(T,1);
    events = [peak_idx(:); foot_idx(:)];
    types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
    [events, ord] = sort(events);
    types = types(ord);
    keep = true(size(events));
    for i = 2:numel(events)
        if types(i) == types(i-1), keep(i) = false; end     % strict alternation
    end
    events = events(keep); types = types(keep);
    if numel(events) < 2, return; end
    phases  = nan(size(events));
    phi_cur = types(1) * pi;
    for i = 1:numel(events)
        phases(i) = phi_cur; phi_cur = phi_cur + pi;
    end
    for i = 1:numel(events)-1
        a = events(i); b = events(i+1);
        if a < 1 || b > T || b <= a, continue; end
        phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
    end
end

function d = date_from_path(folderPath, groupRoot)
% Path component immediately under groupRoot (e.g. '0721' under .../Sert/).
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    if isempty(parts), d=''; else, d = parts{1}; end
end

function tf = is_io_path(folderPath, groupRoot)
% True if the cell-level dir (right after the date dir) contains 'IO'.
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; rebuild the full last segment.
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end
