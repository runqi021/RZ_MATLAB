% Ventral_surface_polar_coh_vs_rayleigh_260729.m
% -----------------------------------------------------------------------
%  TWO polar panels of the SAME ROIs on the SAME phase, scored two ways:
%
%     PANEL A   coherence   radius = band-averaged |C|,  dashed confC
%     PANEL B   Rayleigh    radius = log Z,              dashed log(-ln alpha)
%
%  THE POINT IS THE PHASE. Both panels read the same piecewise-LINEAR phase,
%  interpolated between the detected landmarks:
%
%       insp ONSET (trough/foot) -> 0        insp PEAK -> pi
%       next ONSET -> 2pi, and so on; linear ramps in TIME between them.
%
%  Interpolating peak<->trough is what makes this usable when breathing is
%  irregular: every cycle is warped onto its own [0 2pi] regardless of how long
%  it lasted, so a long breath and a short breath contribute to the same phase
%  bin. That is the middle ground between assuming a fixed breathing frequency
%  (which irregular breathing breaks) and abandoning phase for absolute-time
%  PETH (breath_time_peth_260727.m). This script is the phase-domain half; the
%  PETH is the time-domain half.
%
%  WHY TWO SCORES. They fail differently, so agreement is informative:
%    - coherence asks "is there a linear relationship at breath frequency",
%      via cos(phi) and a multitaper cross-spectrum. It is a waveform measure
%      and can be carried by a few large events.
%    - Rayleigh asks "are the event phases non-uniform on the circle", using
%      only the phase AT each Ca event. No spectral step, no band choice.
%
%  OCCUPANCY CORRECTION (Rayleigh panel, essential). Phase advances faster
%  through the short half of the cycle, so frames are NOT spread evenly over
%  [0 2pi]. A cell firing uniformly in TIME would then look phase-modulated,
%  biased toward the slow (expiratory) half. So the Rayleigh resultant is built
%  from an OCCUPANCY-NORMALISED phase histogram: events per bin divided by
%  FRAMES per bin. The reported angle therefore stays in the raw linear-phase
%  convention -- 0 = onset and pi = peak for every recording, so angles remain
%  poolable and directly comparable with panel A. (Remapping phase by its own
%  ECDF also removes the bias, but then the peak lands at a different angle in
%  every recording and the pooled plot stops meaning one thing.)
%
%  Z = n * Rbar^2, p ~ exp(-Z) (Zar's correction applied), so the threshold is
%  n-INDEPENDENT: one dashed circle at log(-ln alpha) is valid for every ROI.
%  Because Z <= n, an ROI needs >= 4 events to reach alpha=0.05 at all -- the
%  statistic is self-limiting, so no extra event-count gate is imposed.
%  A circular-shift shuffle p is also computed (shifts preserve occupancy
%  exactly) and written to the CSV as an independent check on the analytic p.
%
%  Breath waveform + landmarks come from the DLC-FREE SVD PC1 pipeline:
%     breath waveform : BP.breath        from  breath_peak_pc1.mat
%     insp PEAK        : insp_onset_idx   from  breath_peak_pc1.mat
%     insp ONSET(foot) : insp_start_idx   from  breath_insp_start_pc1.mat
%
%  Layout assumed:
%     <rootPath>/<Group>/<Date>/<cell|IO>/<recording>/ca_spike_data.mat
%       + breath_peak_pc1.mat + breath_insp_start_pc1.mat
%
%  Outputs (under outDir):
%     polar_coh_vs_rayleigh.png / .pdf, polar_coh_vs_rayleigh_data.mat,
%     polar_scores.csv  (every ROI, both scores, both p-values)
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m
%  Runqi Zhang / 2026-07-29
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
% This script lives in <repoRoot>/analysis_260727/ -- add the repo root so
% detect_session_fps / Chronux / the +helper package resolve. (always addpath)
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(scriptDir);   % coh_cfg_260727, for the exclusion list

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath, 'polar_coh_vs_rayleigh_260729');

% group label -> color. IO is split out from ChAT: any folder under ChAT\
% whose cell-level dir name contains 'IO' counts as IO, the rest is true ChAT.
groups       = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};
group_colors = [0    0    0   ;     % IO      black
                0.85 0.10 0.10;     % ChAT    red
                0.10 0.65 0.20;     % Vglut2  green
                0.10 0.30 0.85;     % Vgat    blue
                0.55 0.20 0.75;     % Sst     purple
                0.90 0.45 0.10];    % Sert    orange, matching the cartoon map tint

% which top-level rootPath folders to scan, and how each maps to a group
scan_dirs = {'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};

% Recordings barred from anything spatial or pooled (different stage zero). Not
% present in the archive anyway, but honoured here so the two stay in step.
excludeRecordings = coh_cfg_260727().excludeRecordings;

% ---- POOL PER CELL, not per (recording, ROI) --------------------------------
% A cell re-imaged at a second Z or zoom appears once per recording. Scoring each
% appearance separately double-counts it AND makes significance recording-
% dependent: Sert cell 57 gave logZ [0.98 2.24 2.83], p_shuffle [0.166 0.006
% 0.020] -- non-significant in one recording, significant in the other two.
%
% Phase is the coordinate where pooling is legitimate. 0 = insp onset and
% pi = insp peak mean the same thing in every recording, so event phases can be
% concatenated across recordings; absolute time cannot. Each recording keeps its
% OWN occupancy weights (its cycle shape differs), and n_eff then accumulates,
% so a cell seen three times gets the statistical power of all three.
%
% Same registry as the cartoon map. Sessions with no matcher run fall back to one
% cell per ROI, so pooling is a no-op for them -- in practice this only merges
% Sert, where 68 ROI-observations collapse to 58 cells.
poolPerCell = true;
% Read from the ARCHIVE's own copies, not from the acquisition drives. The
% sources have already moved D: -> E: once and E: gets unplugged; a summary that
% silently loses its cell identities when a drive is absent is worse than one
% that cannot run at all. Copied in 2026-07-31 alongside each session's
% roi_match_out_*, so the archive is self-contained.
cell_link_sources = { ...
    'Sert',   '0721', fullfile(rootPath, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(rootPath, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(rootPath, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    };

nDrop           = 30;        % breath frames tossed up front (align to Ca)
fallback_fps    = 30;
minSpikes       = 2;         % include ROI if it has >= this many spikes
TW              = 4;         % multitaper TW for coherence
alpha_sig       = 0.001;     % primary significance level (jackknife err uses this)
alpha_sig2      = 0.05;      % secondary threshold (drawn as outer dashed circle)
% ---- GCaMP lead compensation: MEASURED, not inherited -----------------------
% Detection runs on raw dF/F (spike_params.detect_on = 'Raw dFF') and SNAPS TO A
% LOCAL MAX within snapWin = 0.2 s, so a detected "event" is the fluorescence
% PEAK, not the spike. Measured by aligning dF/F to the detected events:
%
%   cell                      rise starts   half-rise   detected event
%   Vglut2/1124 roi5 ROI5       -300 ms      -100 ms      0 (the peak)
%   Vglut2/0224 leftpFN ROI3    -133 ms       -67 ms      0 (the peak)
%   Sert/0721   ro5 ROI37       -167 ms       -67 ms      0 (the peak)
%
% The event sits exactly on the dF/F maximum every time, so it lags the underlying
% spike by the GCaMP rise. A lead correction IS justified, and the half-rise time
% (-67 to -100 ms) puts it right at the inherited 0.1 s.
%
% This was briefly set to 0 on the argument that a fixed time shift rotates fast and
% slow breathers by different amounts (120 deg at a 0.43 s cycle, 7 deg at 2.8 s).
% That reasoning was wrong: the lag is a fixed property of the INDICATOR, so a fixed
% time shift is the physically correct correction, and a larger phase consequence
% when breathing is fast is the true answer, not an artifact.
%
% Consequence worth knowing: for the fast breathers 100 ms is ~23% of a whole cycle,
% so a cell whose dF/F peaks just AFTER the breath peak can have fired just BEFORE
% it. Vglut2/0224 ROI3 is exactly that case -- dF/F event at 220 deg, spike estimate
% at 125 deg. The raw trace shows fluorescence; this shows the inferred spike.
%
% SET TO 30 ms (was 100 ms). The half-rise measurement above brackets the plausible
% range, but the half-rise is a LATE estimate of spike time -- the spike initiates
% the rise, so the true lag is nearer the foot than the midpoint, and 100 ms was
% over-correcting. 30 ms is roughly the GCaMP8s time-to-onset.
%
% NOTE this now differs from the temporal_phase* scripts, which still use 0.1, so
% their per-ROI figures will read ~20-70 deg later than this plot depending on
% breath rate. Not a disagreement about the data, just a different lag constant.
ca_lag_sec      = 0.03;

f_breath_search = [0.2 4];   % Hz, search band for breath PSD peak
fwhm_factor     = 0.6;       % coherence band = fwhm_factor x FWHM
min_bw          = 0.05;      % Hz, minimum coherence band width
fmin            = 0.05;      % Hz, PSD lower bound
fmax            = 15;        % Hz, PSD upper bound

% PC1 breath files written by the SVD peak/trough GUIs (fixed names)
peak_file_name  = 'breath_peak_pc1.mat';
start_file_name = 'breath_insp_start_pc1.mat';

% ---- panels ----
% The coherence panel is HIDDEN: log Z on the linear phase is the measure now.
% Coherence is still COMPUTED, so coh_r / coh_th stay in polar_scores.csv and the
% .mat -- it costs little and keeps the comparison available. Set true to show it.
showCoherencePanel = false;

% ---- Rayleigh panel ----
nPhaseBins      = 36;        % 10 deg bins for the occupancy-normalised histogram
alpha_ray       = 0.05;      % inner dashed circle  -> log(-log(alpha))
alpha_ray2      = 0.001;     % outer dashed circle
nShuffleRay     = 500;       % circular-shift shuffles for the independent p
rayRLimPad      = 0.5;       % headroom above max log Z on the radial axis

doSave          = true;
% =====================================================================
% Rayleigh thresholds. p ~ exp(-Z)  =>  Z_crit = -ln(alpha), independent of n,
% so these are single circles valid for every ROI.
%   alpha=0.05  -> Z=3.00 -> logZ=1.10        alpha=0.001 -> Z=6.91 -> logZ=1.93
logZ_crit  = log(-log(alpha_ray));
logZ_crit2 = log(-log(alpha_ray2));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ---- single collector across ALL groups / FOVs / ROIs ----
PP = init_coll();
RY = init_ray();
labels = {};
confC  = NaN;

% cell-identity lookup: '<scanDir>/<date>/<recName>/<maskL label>' -> cell key.
% Namespaced per session because cell ids restart at 1 in every matcher run.
% 'tossed' rows (cell_id NaN) are rejected curation and are dropped entirely.
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
if poolPerCell
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
end

% raw ingredients, kept so the pooling can be done after the scan
REC = struct('bidx',{},'occ',{},'T',{},'phiW',{});
OBS = struct('rec',{},'ev',{},'gi',{},'f_pk',{},'label',{},'cellKey',{});
% Chronux confC closed form: confC = sqrt(1 - alpha^(1/(K-1))), K = 2*TW-1
K_tap  = 2*TW - 1;
confC2 = sqrt(1 - alpha_sig2 ^ (1/(K_tap - 1)));
n_per_group = zeros(1, numel(groups));

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir)
        warning('Folder missing: %s', gdir); continue;
    end

    allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
    fprintf('\n=== [%s] %d recordings ===\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;
        recName    = folder_basename(folderPath);
        recDate    = date_from_path(folderPath, gdir);
        % assign group: any cell-level dir named "*IO*" -> IO (regardless of scan dir)
        gname = sname;
        if is_io_path(folderPath, gdir)
            gname = 'IO';
        end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), warning('No color for group "%s" -- skipping %s', gname, recName); continue; end
        if any(strcmp(recName, excludeRecordings))
            fprintf('  skip (excluded: different stage zero): %s\n', recName); continue;
        end
        try
            bpFile = fullfile(folderPath, peak_file_name);
            ipFile = fullfile(folderPath, start_file_name);
            if ~isfile(bpFile), fprintf('  skip (no %s): %s\n', peak_file_name, recName); continue; end
            if ~isfile(ipFile), fprintf('  skip (no %s): %s\n', start_file_name, recName); continue; end

            fps = detect_session_fps(folderPath, fallback_fps);
            CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));
            nROI = numel(CA.roi_spikes);
            nCa  = numel(CA.roi_spikes(1).spike_train);

            BP = load(bpFile);     % SVD PC1 peaks  (insp_onset_idx, breath)
            IP = load(ipFile);     % SVD PC1 feet   (insp_start_idx)
            bw = detrend(double(BP.breath(:)));
            bw(1:min(nDrop,numel(bw))) = [];
            bw = bw - mean(bw);

            % --- piecewise phase: foot=0, peak=pi  (post-toss frame index)
            peak_idx = round(BP.insp_onset_idx(:)) - nDrop;     % insp PEAK
            foot_idx = round(IP.insp_start_idx(:)) - nDrop;     % insp ONSET (foot/trough)

            % --- per-session timing fix: Vglut2/1124 used rising-edge trigger
            %     instead of falling-edge -> breath samples sit 1 frame late.
            %     Shift breath events +1 frame to realign with Ca.
            if strcmpi(sname,'Vglut2') && strcmp(recDate,'1124')
                peak_idx = peak_idx + 1;
                foot_idx = foot_idx + 1;
                bw = [bw(1); bw(1:end-1)];   % delay breath waveform by 1 frame
            end

            T = min([numel(bw), nCa]);
            peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
            foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
            if numel(peak_idx) < 2 || numel(foot_idx) < 2
                fprintf('  skip (too few events): %s\n', recName); continue;
            end
            phi = piecewise_phase_local(peak_idx, foot_idx, T);
            ref = cos(phi); ref(isnan(ref)) = 0;
            ref = ref - mean(ref(:));
            bw  = bw(1:T);

            %% breath waveform PSD -> coherence band
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
            pc.fpass=band; pc.err=[2,alpha_sig];

            % ---- occupancy: frames per phase bin, shared by every ROI here ----
            % phi is NaN outside the first..last landmark; those frames are not
            % part of any cycle and must be excluded from BOTH the event phases
            % and the occupancy, or the normalisation is against the wrong total.
            edges  = linspace(0, 2*pi, nPhaseBins+1);
            ctrs   = edges(1:end-1) + diff(edges)/2;
            phiW   = mod(phi, 2*pi);                  % phi accumulates past 2pi
            valid  = ~isnan(phiW);
            bidx   = nan(T,1);
            bidx(valid) = discretize(phiW(valid), edges);
            occ    = accumarray(bidx(valid), 1, [nPhaseBins 1]);   % FRAMES per bin
            REC(end+1) = struct('bidx',bidx,'occ',occ,'T',T,'phiW',phiW); %#ok<SAGROW>
            recIdx = numel(REC);

            nInc = 0;
            lag = round(ca_lag_sec*fps);          % GCaMP lead: shift spikes EARLIER before coherence
            for rid = 1:nROI
                st = double(CA.roi_spikes(rid).spike_train(:));
                st = st(1:min(T,numel(st)));
                if numel(st)<T, st(end+1:T)=0; end
                if sum(st) < minSpikes, continue; end
                % cell-identity key uses the SCAN DIR, not gname: the matcher ran
                % per session, before the IO folder override was applied.
                ckin = sprintf('%s/%s/%s/%d', sname, recDate, recName, rid);
                if isKey(tossed_set, ckin), continue; end     % rejected in curation
                nInc = nInc + 1;
                stL = [st(1+lag:end); zeros(lag,1)];   % lead-shift earlier by lag frames
                [PP, confC] = add_coh(PP, ref, stL - mean(stL), pc, band, gi, confC, f_pk);

                lab = sprintf('%s/%s/%s/%d', gname, recDate, recName, rid);
                if isKey(cell_map, ckin), ck = cell_map(ckin); else, ck = ['roi:' lab]; end
                OBS(end+1) = struct('rec',recIdx,'ev',find(stL > 0),'gi',gi, ...
                                    'f_pk',f_pk,'label',lab,'cellKey',ck); %#ok<SAGROW>
                labels{end+1} = lab; %#ok<SAGROW>
            end
            n_per_group(gi) = n_per_group(gi) + nInc;
            fprintf('  [%d] %-50s band [%.2f %.2f] Hz  %d/%d ROI\n', ...
                    kk, recName, band(1), band(2), nInc, nROI);
        catch ME
            warning('  ERROR %s: %s', recName, ME.message);
        end
    end
end

if isempty(PP.r), error('No ROIs collected.'); end

%% ============== POOL THE RAYLEIGH PER CELL ==============
% One entry per CELL. Event phases from all of that cell's recordings are
% concatenated, each carrying its own recording's occupancy weight; n_eff sums,
% so a cell seen three times is scored with the power of all three.
cellKeys = {OBS.cellKey};
[uCell, ~, obsOfCell] = unique(cellKeys, 'stable');
cellLabels = cell(numel(uCell),1);  cellNObs = zeros(numel(uCell),1);
for c = 1:numel(uCell)
    m  = find(obsOfCell == c);
    cellNObs(c) = numel(m);
    A = []; W = []; parts = struct('ev',{},'rec',{});
    for j = m(:)'
        r  = REC(OBS(j).rec);
        ev = OBS(j).ev;
        ev = ev(ev>=1 & ev<=r.T);
        ev = ev(~isnan(r.bidx(ev)));
        if isempty(ev), continue; end
        A = [A; r.phiW(ev(:))];                       %#ok<AGROW> event phases
        W = [W; 1./max(r.occ(r.bidx(ev(:))),1)];      %#ok<AGROW> per-recording weights
        parts(end+1) = struct('ev',ev,'rec',OBS(j).rec); %#ok<AGROW>
    end
    RY = add_ray_pooled(RY, A, W, parts, REC, ctrs, nPhaseBins, nShuffleRay, ...
                        OBS(m(1)).gi, OBS(m(1)).f_pk);
    if cellNObs(c) == 1
        cellLabels{c} = OBS(m(1)).label;
    else
        cellLabels{c} = sprintf('%s  [+%d more rec]', OBS(m(1)).label, cellNObs(c)-1);
    end
end
fprintf('\nPooled %d ROI-observations -> %d cells (%d cells span >1 recording)\n', ...
        numel(OBS), numel(uCell), nnz(cellNObs>1));
for gi = 1:numel(groups)
    mo = arrayfun(@(o) o.gi==gi, OBS);
    mc = arrayfun(@(c) OBS(find(obsOfCell==c,1)).gi==gi, 1:numel(uCell));
    if ~any(mo), continue; end
    fprintf('   %-8s %4d observations -> %4d cells\n', groups{gi}, nnz(mo), nnz(mc));
end
% n_per_group must now count CELLS, since that is what the panel plots
n_per_group(:) = 0;
for c = 1:numel(uCell)
    g = OBS(find(obsOfCell==c,1)).gi;  n_per_group(g) = n_per_group(g) + 1;
end

PP.th_raw = PP.th;
if ca_lag_sec == 0
    fprintf('No GCaMP lead correction: phases are those of the detected dF/F event (%d ROIs).\n', numel(PP.r));
else
    fprintf('GCaMP comp: events lead-shifted %.0f ms before scoring (%d ROIs).\n', ...
            ca_lag_sec*1000, numel(PP.r));
end

%% ============================ TWO POLAR PANELS =======================
% Arcs go on every cell OUTSIDE the alpha=0.05 circle, i.e. the same criterion the
% radius and the dashed circles already show. Previously they were restricted to
% cells that also passed the circular-shift shuffle, which was confusing and wrong:
%
%  - it was a knife-edge. Vglut2/1124 roi5 ROI5 sits at logZ 3.48, far outside both
%    circles, but p_shuffle = 0.056 vs 0.050 -- 28 versus 25 exceedances out of 500 --
%    so it lost its arc while staying plotted well outside the circle.
%  - worse, the shuffle is BIASED AGAINST rhythmic cells. A circular shift rotates
%    every event phase by nearly the same amount when breathing is regular, so it
%    moves the phase concentration rather than destroying it. That recording's breath
%    CV is 0.075 (Sert's is 0.110), which is why all six of its ROIs pile up at
%    p_shuffle 0.044-0.074 despite large logZ. The test is weakest exactly where the
%    modulation is cleanest.
%
% p_shuffle is still computed and written to polar_scores_percell.csv as a
% cross-check; it just no longer gates what the figure draws.
arcMask = RY.logZ >= logZ_crit;
if showCoherencePanel
    fig = figure('Color','w','Name','Ventral surface: coherence vs Rayleigh on linear phase', ...
                 'Units','centimeters','Position',[2 2 30 15]);
    set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    axA = polaraxes(fig,'Position',[0.045 0.06 0.36 0.76]);
    plot_panel(axA, PP, group_colors, confC, confC2, '');
    % plot_panel forces Interpreter 'none' (recording names elsewhere contain
    % underscores), which would print a literal \pi. Retitle with tex.
    title(axA, 'A   COHERENCE      r = mean|C| in band', 'Interpreter','tex');
    axB = polaraxes(fig,'Position',[0.505 0.06 0.36 0.76]);
    ttlB = ['B   RAYLEIGH      r = log Z      dot = circular median,' ...
            ' arc = IQR (shuffle-confirmed only)'];
else
    fig = figure('Color','w','Name','Ventral surface: Rayleigh on linear breath phase', ...
                 'Units','centimeters','Position',[2 2 19 16]);
    set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
    axB = polaraxes(fig,'Position',[0.06 0.05 0.66 0.78]);
    ttlB = ['r = log Z      dot = circular median,  arc = IQR' newline ...
            '(arcs on cells outside the \alpha = 0.05 circle)'];
end
plot_ray_panel(axB, RY, group_colors, logZ_crit, logZ_crit2, rayRLimPad, ttlB, arcMask);

add_group_legend(fig, axB, groups, group_colors, n_per_group);

nSigC = nnz(PP.r    >= confC);
nSigR = nnz(RY.logZ >= logZ_crit);
nSigS = nnz(RY.p_sh <  alpha_ray);
% The panel-B circle comes from the ANALYTIC threshold, which measures mildly
% liberal against the circular-shift null, so the shuffle-confirmed count is stated
% next to it rather than left implicit. Per-ROI shuffle p is in polar_scores.csv.
% Group counts live in the legend; repeating them here overflowed the title.
% Kept short: the narrow single-panel figure clips a long title. The full method
% is in the panel title, the axis labels and PROGRESS_260729.md.
if ca_lag_sec == 0
    lagTxt = 'NO GCaMP lead correction';
else
    lagTxt = sprintf('GCaMP lead correction: events shifted %.0f ms EARLIER', ca_lag_sec*1000);
end
l1 = sprintf('N = %d cells  (%d ROI-observations pooled)   |   insp onset = 0,  insp peak = \\pi', ...
             numel(uCell), numel(OBS));
l2 = sprintf('logZ crit(\\alpha=%.2g) = %.2f \\rightarrow %d sig   (%d also pass the shuffle)', ...
             alpha_ray, logZ_crit, nSigR, nnz(RY.logZ>=logZ_crit & RY.p_sh<alpha_ray));
if showCoherencePanel
    l2 = sprintf('A  confC(\\alpha=%.3f) = %.2f \\rightarrow %d sig        B  %s', ...
                 alpha_sig, confC, nSigC, l2);
end
sgtitle({l1, l2}, 'FontSize', 8.5);

% The lead correction changes every reported angle, so it is stated ON the figure
% rather than buried in the .mat -- a reader must not have to guess which
% convention a phase is in.
annotation(fig, 'textbox', [0.02 0.005 0.96 0.045], 'String', lagTxt, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'EdgeColor',[0.6 0.6 0.6], 'BackgroundColor',[0.96 0.96 0.96], ...
    'FontSize',8.5, 'FontWeight','bold', 'Margin',2, 'Interpreter','tex');
fprintf('  panel B: %d above the analytic circle, %d pass the shuffle, %d both\n', ...
        nSigR, nSigS, nnz(RY.logZ>=logZ_crit & RY.p_sh<alpha_ray));

%% ---- do the two scores agree? ----
% PP is per (recording, ROI); RY is now per CELL. Collapse coherence to the cell
% by taking its BEST recording -- "does any recording of this cell show
% coherence" -- so the two are comparable. Stated explicitly because it is a
% choice, not a pooled coherence (pooling that properly means summing
% cross-spectra, which is not what this script does).
sr = RY.logZ >= logZ_crit;
cohCellR  = zeros(numel(uCell),1);  cohCellTh = nan(numel(uCell),1);
for c = 1:numel(uCell)
    m = find(obsOfCell == c);
    [cohCellR(c), b] = max(PP.r(m));
    cohCellTh(c) = PP.th(m(b));
end
sc = cohCellR >= confC;
fprintf('\n=========== coherence vs Rayleigh agreement (per CELL) ===========\n');
fprintf('  significant in BOTH        : %d\n', nnz(sc & sr));
fprintf('  coherence only             : %d\n', nnz(sc & ~sr));
fprintf('  Rayleigh only              : %d\n', nnz(~sc & sr));
fprintf('  neither                    : %d\n', nnz(~sc & ~sr));
ok = ~isnan(cohCellTh) & ~isnan(RY.thMed);
if nnz(ok) > 2
    dth = angle(exp(1i*(cohCellTh(ok) - RY.thMed(ok))));
    fprintf('  phase offset (all cells)   : %+.2f rad over %d cells\n', ...
            angle(mean(exp(1i*dth))), nnz(ok));
    both = sc & sr & ok;
    if nnz(both) > 1
        d2 = angle(exp(1i*(cohCellTh(both) - RY.thMed(both))));
        fprintf('     restricted to the %d cells significant in both: %+.2f rad (sd %.2f)\n', ...
                nnz(both), angle(mean(exp(1i*d2))), std(d2));
    end
end
% analytic vs shuffle p. The dashed circle is drawn from the ANALYTIC threshold, so
% it is only honest if the analytic p tracks the shuffle. Report the comparison
% rather than assuming it.
okp = ~isnan(RY.p_sh);
fprintf('  Rayleigh p: analytic<%.2f = %d,  shuffle<%.2f = %d,  both = %d\n', ...
        alpha_ray, nnz(RY.p_an<alpha_ray), alpha_ray, nnz(RY.p_sh(okp)<alpha_ray), ...
        nnz(RY.p_an<alpha_ray & RY.p_sh<alpha_ray));
fprintf('  n_eff / n_events: median %.2f  (1.0 = uniform occupancy)\n', ...
        median(RY.nEff(RY.nEv>0) ./ RY.nEv(RY.nEv>0)));
fl = 1/(1+nShuffleRay);                       % shuffle p cannot go below this
comparable = okp & RY.p_sh > fl;              % exclude floored shuffle p
if nnz(comparable) > 3
    la = log10(max(RY.p_an(comparable), 1e-300));
    ls = log10(RY.p_sh(comparable));
    fprintf('  log10 p agreement over %d ROI (shuffle above its %.4f floor):\n', nnz(comparable), fl);
    fprintf('     median(analytic - shuffle) = %+.2f decades,  Spearman rho = %.3f\n', ...
            median(la - ls), corr(la, ls, 'type','Spearman'));
    if median(la - ls) < -1
        fprintf(2,'     WARNING analytic p is >10x optimistic -- trust ray_p_shuffle, not the circle.\n');
    end
end
% Benjamini-Hochberg across all ROIs, reported but NOT drawn (keeps the figure simple)
[~, bhSig] = bh_fdr(RY.p_an, 0.05);
fprintf('  Rayleigh after BH-FDR q=0.05 : %d / %d\n', nnz(bhSig), numel(bhSig));

%% ============ PRINT CELLS SIGNIFICANT BY EITHER SCORE =================
fprintf('\n============ cells significant by EITHER score ============\n');
fprintf('  %-9s %3s %5s %-6s %-6s %-7s %-8s %-8s  %s\n', ...
        'group','rec','nEv','r_coh','logZ','med_deg','p_an','p_shuf','cell');
either = find(sc | sr);
[~, ord] = sort(RY.logZ(either), 'descend');
for k = ord(:)'
    i = either(k);
    if sc(i) && sr(i), flag = 'BOTH'; elseif sc(i), flag = 'coh '; else, flag = 'ray '; end
    fprintf('  %-9s %3d %5d %6.3f %6.2f %+7.0f %8.1e %8.3f  [%s] %s\n', ...
            groups{RY.colorIdx(i)}, cellNObs(i), RY.nEv(i), cohCellR(i), RY.logZ(i), ...
            rad2deg(RY.thMed(i)), RY.p_an(i), RY.p_sh(i), flag, cellLabels{i});
end
fprintf('==========================================================\n\n');

if doSave
    % ---- per CELL: what the panel plots ----
    fid = fopen(fullfile(outDir,'polar_scores_percell.csv'),'w');
    fprintf(fid,['group,cell_key,n_recordings,first_label,n_events,n_eff,f_breath_Hz,' ...
                 'ray_Rbar,ray_Z,ray_logZ,ray_mean_deg,' ...
                 'ray_med_deg,ray_iqr_lo_deg,ray_iqr_hi_deg,ray_iqr_width_deg,' ...
                 'ray_p_analytic,ray_p_shuffle,ray_sig,ray_sig_BH,' ...
                 'coh_r_best_rec,coh_th_best_rec_deg,coh_sig\n']);
    for c = 1:numel(uCell)
        fprintf(fid,['%s,%s,%d,"%s",%d,%.2f,%.4f,%.6f,%.4f,%.4f,%.2f,' ...
                     '%.2f,%.2f,%.2f,%.2f,%.6g,%.4f,%d,%d,%.6f,%.2f,%d\n'], ...
            groups{RY.colorIdx(c)}, uCell{c}, cellNObs(c), cellLabels{c}, ...
            RY.nEv(c), RY.nEff(c), RY.f_pk(c), ...
            RY.Rbar(c), RY.Z(c), RY.logZ(c), rad2deg(RY.th(c)), ...
            rad2deg(RY.thMed(c)), rad2deg(RY.thMed(c)+RY.q25(c)), ...
            rad2deg(RY.thMed(c)+RY.q75(c)), rad2deg(RY.q75(c)-RY.q25(c)), ...
            RY.p_an(c), RY.p_sh(c), sr(c), bhSig(c), ...
            cohCellR(c), rad2deg(cohCellTh(c)), sc(c));
    end
    fclose(fid);
    fprintf('Saved polar_scores_percell.csv (%d cells) to %s\n', numel(uCell), outDir);

    % ---- per (recording, ROI): the unpooled coherence, kept for reference ----
    fid = fopen(fullfile(outDir,'polar_scores_perobs.csv'),'w');
    fprintf(fid,'group,date,fov,roi,label,cell_key,coh_r,coh_th_rad,coh_th_deg,coh_sig\n');
    for ii = 1:numel(labels)
        parts = regexp(labels{ii},'/','split');   % 'group/date/recName/rid'
        if numel(parts) >= 4
            date_part = parts{2};  fov_part = strjoin(parts(3:end-1),'/');  roi_part = parts{end};
        else
            date_part = '';  fov_part = labels{ii};  roi_part = '';
        end
        fprintf(fid,'%s,%s,"%s",%s,%s,%s,%.6f,%.6f,%.2f,%d\n', ...
            groups{PP.colorIdx(ii)}, date_part, fov_part, roi_part, labels{ii}, ...
            OBS(ii).cellKey, PP.r(ii), PP.th(ii), rad2deg(PP.th(ii)), PP.r(ii)>=confC);
    end
    fclose(fid);
    fprintf('Saved polar_scores_perobs.csv (%d observations) to %s\n', numel(labels), outDir);
end

if doSave
    exportgraphics(fig, fullfile(outDir,'polar_coh_vs_rayleigh.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'polar_coh_vs_rayleigh.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'polar_coh_vs_rayleigh_data.mat'), ...
         'PP','RY','labels','confC','confC2','groups','group_colors','n_per_group', ...
         'uCell','cellLabels','cellNObs','obsOfCell','cohCellR','cohCellTh', ...
         'sc','sr','bhSig','arcMask','poolPerCell','cell_link_sources', ...
         'TW','alpha_sig','alpha_sig2','alpha_ray','alpha_ray2', ...
         'logZ_crit','logZ_crit2','nPhaseBins','nShuffleRay', ...
         'f_breath_search','fwhm_factor','min_bw','minSpikes','nDrop','ca_lag_sec');
    fprintf('Saved polar_coh_vs_rayleigh.png/.pdf + .mat to %s\n', outDir);
end
fprintf('Done. %d ROIs across %d groups.\n', numel(PP.r), sum(n_per_group>0));

%% ========================= LOCAL FUNCTIONS ==========================
function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[],'colorIdx',[],'f_pk',[]);
end

function plot_ray_panel(ax, C, group_colors, lzc, lzc2, pad, ttl, arcMask)
% Same geometry and colour key as the coherence panel, so the two read together:
% angle = mean event phase (0 = insp onset, pi = insp peak), radius = log Z.
% log Z can be <= 0 when Z < 1; those ROIs are clamped to the centre rather than
% given a negative radius, which polarplot would silently reflect to the opposite
% angle and so place at the WRONG phase.
    hold(ax,'on');
    rmax = max([C.logZ(:); lzc2]) + pad;
    thc  = linspace(0,2*pi,360);
    % The two dashed circles are SIGNIFICANCE thresholds. Because p ~ exp(-Z), the
    % critical Z = -ln(alpha) does not depend on the event count, so one circle is
    % valid for every ROI -- which is the reason log Z works as a radius at all.
    % Labelled on the figure so the panel stands alone.
    polarplot(ax, thc, repmat(lzc2,1,360), '--','Color',[0.5 0.5 0.5],'LineWidth',0.8);
    polarplot(ax, thc, repmat(lzc, 1,360), 'k--','LineWidth',1);
    % White background box: the circles pass straight through the data cloud, so a
    % bare label is unreadable wherever it is placed.
    thLab = deg2rad(325);
    text(ax, thLab, lzc,  sprintf('\\alpha=%.2g', exp(-exp(lzc))), ...
         'FontSize',7.5,'Color','k','HorizontalAlignment','center', ...
         'BackgroundColor','w','Margin',0.5);
    text(ax, thLab, lzc2, sprintf('\\alpha=%.3g', exp(-exp(lzc2))), ...
         'FontSize',7.5,'Color',[0.35 0.35 0.35],'HorizontalAlignment','center', ...
         'BackgroundColor','w','Margin',0.5);

    % IQR arcs first, so the median dots sit on top of them.
    % Drawn ONLY for the shuffle-confirmed ROIs. An arc for all 368, or even for
    % the 78 above the liberal analytic circle, is mostly near-uniform cells whose
    % IQR approaches the full circle -- they pile into concentric rings of ink that
    % hide the tuned cells. Every ROI is still plotted as a dot.
    for k = 1:numel(C.logZ)
        ci = C.colorIdx(k);
        if isnan(C.thMed(k)) || ci<1 || ci>size(group_colors,1), continue; end
        if ~arcMask(k), continue; end
        r   = max(C.logZ(k), 0);
        col = group_colors(ci,:);
        if ci == 1, col = [0 0 0]; end
        arc = linspace(C.thMed(k)+C.q25(k), C.thMed(k)+C.q75(k), 40);
        polarplot(ax, arc, r*ones(size(arc)), '-', 'Color',[col 0.55], 'LineWidth',1.3);
    end
    for k = 1:numel(C.logZ)
        ci = C.colorIdx(k);
        if isnan(C.thMed(k)) || ci<1 || ci>size(group_colors,1), continue; end
        r = max(C.logZ(k), 0);
        col = group_colors(ci,:);
        if ci == 1     % IO: hollow black circle, matching panel A
            polarplot(ax, C.thMed(k), r, 'o', 'MarkerFaceColor','none', ...
                      'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.8);
        else
            polarplot(ax, C.thMed(k), r, 'o', 'MarkerFaceColor',col, ...
                      'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.4);
        end
    end
    ax.RLim=[0 rmax]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
    title(ax, ttl, 'Interpreter','tex');
end

function [pcrit, sig] = bh_fdr(p, q)
% Benjamini-Hochberg. Returns the critical p and the significance mask.
    p = p(:);  ok = ~isnan(p);
    sig = false(size(p));  pcrit = 0;
    ps = sort(p(ok));  m = numel(ps);
    if m == 0, return; end
    below = find(ps <= (1:m)'/m * q);
    if isempty(below), return; end
    pcrit = ps(max(below));
    sig(ok & p <= pcrit) = true;
end

function C = init_ray()
    C = struct('th',[],'Rbar',[],'Z',[],'logZ',[],'p_an',[],'p_sh',[], ...
               'nEv',[],'nEff',[],'colorIdx',[],'f_pk',[], ...
               'thMed',[],'q25',[],'q75',[]);
end

function m = circ_median_w(a, w)
% Weighted circular median: the direction minimising the total weighted angular
% distance to the data. Solved on a 0.5 deg grid -- a closed form does not exist
% on the circle, and the grid is exact enough next to a 10 deg histogram.
    g = linspace(-pi, pi, 721);  g(end) = [];
    d = abs(angle(exp(1i*(a(:) - g))));      % nEv x nGrid angular distances
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

function C = add_ray_pooled(C, a, w, parts, REC, ctrs, nB, nSh, cidx, f_pk)
% Rayleigh for ONE CELL, pooled over all recordings it appears in.
%
%   a  event phases   (rad, 0 = insp onset, pi = insp peak) from every recording
%   w  event weights  1/(frames in that event's phase bin) IN ITS OWN recording
%
% Pooling phases across recordings is legitimate precisely because phase is
% normalised -- 0 and pi mean the same landmark everywhere. Weights stay
% per-recording because each recording has its own cycle shape and occupancy.
%
% The null shifts EACH recording's event train independently and re-pools, so it
% preserves both the per-recording occupancy and the bursting of each train.
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
                sh = round(s * r.T / (nSh+1)) + q;         % different shift per recording
                evs = mod(parts(q).ev - 1 + sh, r.T) + 1;
                evs = evs(~isnan(r.bidx(evs)));
                if isempty(evs), continue; end
                as = [as; r.phiW(evs(:))];                       %#ok<AGROW>
                ws = [ws; 1./max(r.occ(r.bidx(evs(:))),1)];      %#ok<AGROW>
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
        u   = angle(exp(1i*(a - med)));
        C.thMed(end+1,1) = med;
        C.q25(end+1,1)   = wprctile(u, w, 25);
        C.q75(end+1,1)   = wprctile(u, w, 75);
    else
        C.thMed(end+1,1)=NaN; C.q25(end+1,1)=NaN; C.q75(end+1,1)=NaN;
    end
end

function [th, Rbar, nEff] = wresultant(a, w)
% Weighted circular resultant with Kish's effective sample size. Identical maths
% to ray_resultant but taking event lists rather than a binned histogram, so it
% can pool across recordings whose occupancy profiles differ.
    S1 = sum(w);  S2 = sum(w.^2);
    if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
    nEff = S1^2 / max(S2, eps);
    v    = sum(w(:) .* exp(1i*a(:))) / S1;
    th   = angle(v);
    Rbar = min(abs(v), 1);
end

function C = add_ray(C, ev, bidx, occ, ctrs, nB, T, nSh, cidx, f_pk, phiW) %#ok<DEFNU>
% Rayleigh on the OCCUPANCY-NORMALISED phase histogram of Ca event phases.
%
% Angle stays in raw linear-phase units (0 = insp onset, pi = insp peak) so it is
% comparable across recordings and with the coherence panel; the occupancy bias is
% removed in the STATISTIC instead, by dividing events-per-bin by frames-per-bin.
%
% Z = n*Rbar^2 with n = number of events, p ~ exp(-Z) with Zar's finite-n
% correction. The shuffle p uses CIRCULAR SHIFTS of the event train, which
% preserve the occupancy structure exactly, so it is a genuine independent check
% rather than a restatement of the same assumption.
    ev = ev(ev>=1 & ev<=T);
    ev = ev(~isnan(bidx(ev)));            % drop events outside any breath cycle
    n  = numel(ev);
    [th, Rbar, nEff] = ray_resultant(bidx(ev), occ, ctrs, nB);
    Z    = nEff * Rbar^2;
    % Zar (1999) eq 27.4: better than exp(-Z) alone at small n
    p_an = exp(-Z) * (1 + (2*Z - Z^2)/(4*nEff) ...
                        - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*nEff^2));
    p_an = min(max(p_an, 0), 1);

    % circular-shift null. Rigid shifts of the whole event train preserve BOTH the
    % occupancy profile and the train's own autocorrelation/bursting, neither of
    % which the analytic formula knows about -- so this is a genuinely independent
    % check, not a restatement of the same assumption.
    if n >= 1 && nSh > 0
        Znull = zeros(nSh,1);
        shifts = round(linspace(1, T, nSh+2));  shifts = shifts(2:end-1);
        for s = 1:nSh
            evs = mod(ev - 1 + shifts(s), T) + 1;
            evs = evs(~isnan(bidx(evs)));
            if isempty(evs), Znull(s) = 0; continue; end
            [~, Rs, nEs] = ray_resultant(bidx(evs), occ, ctrs, nB);
            Znull(s) = nEs * Rs^2;
        end
        p_sh = (1 + nnz(Znull >= Z)) / (1 + nSh);
    else
        p_sh = NaN;
    end
    C.nEff(end+1,1) = nEff;

    % ---- boxplot-style summary of the phase histogram ----
    % The dot is the circular MEDIAN and the arc the circular IQR, which is what
    % the per-ROI phase histogram actually shows. Both are OCCUPANCY-WEIGHTED
    % (w = 1/frames-in-bin) for the same reason the resultant is: without it the
    % median slides into the slow expiratory half where most frames live.
    % Median/IQR use the CONTINUOUS event phase, not the bin centre, so they are
    % not quantised to the 10 deg histogram grid.
    if n >= 2
        a = phiW(ev);  a = a(:);
        w = 1 ./ max(occ(bidx(ev)), 1);  w = w(:);
        med = circ_median_w(a, w);
        u   = angle(exp(1i*(a - med)));         % centred on the median, (-pi pi]
        C.thMed(end+1,1) = med;
        C.q25(end+1,1)   = wprctile(u, w, 25);  % offsets from the median, <= 0
        C.q75(end+1,1)   = wprctile(u, w, 75);  % offsets from the median, >= 0
    else
        C.thMed(end+1,1) = NaN;  C.q25(end+1,1) = NaN;  C.q75(end+1,1) = NaN;
    end

    C.th(end+1,1)   = th;
    C.Rbar(end+1,1) = Rbar;
    C.Z(end+1,1)    = Z;
    C.logZ(end+1,1) = log(max(Z, eps));
    C.p_an(end+1,1) = p_an;
    C.p_sh(end+1,1) = p_sh;
    C.nEv(end+1,1)  = n;
    C.colorIdx(end+1,1) = cidx;
    C.f_pk(end+1,1) = f_pk;
end

function [th, Rbar, nEff] = ray_resultant(evBins, occ, ctrs, nB)
% WEIGHTED Rayleigh resultant with an EFFECTIVE sample size.
%
% Occupancy correction = weight each event by 1/(frames in its phase bin), which
% is the same thing as taking the resultant of the events/frames rate profile.
% Measured occupancy here is ~7x denser in the slow expiratory half than in the
% fast inspiratory half (281-304 vs 27-47 frames per 10 deg bin), so those weights
% are very uneven -- a couple of events landing in a sparse inspiratory bin get
% ~7x leverage.
%
% That is exactly why Z must NOT use the raw event count. The analytic Rayleigh
% Z = n*Rbar^2 assumes n independent, EQUALLY weighted angles. Feeding it a
% weighted Rbar with raw n produced p = 5e-79 where the shuffle said p = 0.05.
% Kish's effective sample size is the standard correction:
%
%     nEff = (sum w)^2 / sum(w^2)
%
% which collapses to n when occupancy is uniform, and shrinks toward the number of
% well-sampled bins when it is not.
%
% Bins never visited (occ==0) carry no information and are dropped rather than
% divided by zero.
    th = NaN; Rbar = 0; nEff = 0;
    if isempty(evBins), return; end
    cnt  = accumarray(evBins(:), 1, [nB 1]);
    good = occ > 0;
    wsum = zeros(nB,1);                  % total weight per bin = cnt/occ
    wsum(good) = cnt(good) ./ occ(good);
    if ~any(wsum > 0), return; end
    S1 = sum(wsum);                                       % sum of event weights
    S2 = sum(cnt(good) ./ (occ(good).^2));                % sum of squared weights
    nEff = S1^2 / max(S2, eps);
    v    = sum(wsum .* exp(1i*ctrs(:))) / S1;
    th   = angle(v);
    Rbar = min(abs(v), 1);
end

function [C, confC] = add_coh(C, x, y, pc, band, cidx, confC, f_pk)
% Chronux band-averaged magnitude + circular-mean phase + jackknife CI.
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f = f(:);
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

function plot_panel(ax, C, group_colors, confC, confC2, ttl)
% theta-zero top, clockwise, dashed confC circles (primary + secondary),
% marker per ROI colored by group with magnitude-CI bar + phase-CI arc.
    hold(ax,'on');
    thc = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(confC, 1,360), 'k--','LineWidth',1);
    polarplot(ax, thc, repmat(confC2,1,360), '--','Color',[0.5 0.5 0.5],'LineWidth',0.8);

    for k = 1:numel(C.r)
        ci = C.colorIdx(k);
        if isnan(C.th(k)) || isnan(C.r(k)) || ci<1 || ci>size(group_colors,1), continue; end
        col = group_colors(ci,:);
        sig = C.r(k) >= confC;                  % above confC dashed circle
        if sig
            polarplot(ax, [C.th(k) C.th(k)], [C.rlo(k) C.rhi(k)], ...
                      '-', 'Color',col, 'LineWidth',1.0);
            if ~isnan(C.dphi(k))
                arc = linspace(C.th(k)-C.dphi(k), C.th(k)+C.dphi(k), 30);
                polarplot(ax, arc, C.r(k)*ones(size(arc)), ...
                          '-', 'Color',col, 'LineWidth',1.0);
            end
        end
        if ci == 1     % IO: hollow black circle
            polarplot(ax, C.th(k), C.r(k), 'o', ...
                      'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
                      'MarkerSize',6, 'LineWidth',0.8);
        else
            polarplot(ax, C.th(k), C.r(k), 'o', ...
                      'MarkerFaceColor',col, 'MarkerEdgeColor','k', ...
                      'MarkerSize',6, 'LineWidth',0.4);
        end
    end
    ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
    title(ax, ttl, 'Interpreter','none');
end

function add_group_legend(fig, refAx, groups, group_colors, n_per_group)
% Color-key legend showing group -> hue + ROI count.
    pos = refAx.Position;
    legAx = axes(fig, 'Position', [pos(1)+pos(3)+0.01, pos(2)+0.30*pos(4), 0.12, 0.40*pos(4)]);
    hold(legAx,'on'); axis(legAx,'off');
    N = numel(groups);
    y = linspace(0.90, 0.10, max(N,2));
    for k = 1:N
        if k == 1     % IO swatch: hollow black, matches plot
            plot(legAx, 0.08, y(k), 'o', 'MarkerFaceColor','none', ...
                 'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',0.8);
        else
            plot(legAx, 0.08, y(k), 'o', 'MarkerFaceColor',group_colors(k,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',0.4);
        end
        text(legAx, 0.25, y(k), sprintf('%s (n=%d)', groups{k}, n_per_group(k)), ...
             'FontSize',8, 'Interpreter','none');
    end
    xlim(legAx,[0 1]); ylim(legAx,[0 1]);
end

function d = date_from_path(folderPath, groupRoot)
% Path component immediately under groupRoot (e.g. '0521' under .../ChAT/).
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    if isempty(parts), d=''; else, d = parts{1}; end
end

function tf = is_io_path(folderPath, groupRoot)
% True if the cell-level dir (the one right after the date dir under
% groupRoot) contains 'IO' (case-insensitive). Layout:
%   <groupRoot>\<date>\<cell-or-IO>\<recording>\
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; this rebuilds the full
% last path segment for folders with dots in the name.
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase reference: FEET at 0/2pi/..., PEAKS at pi/3pi/...
% Linear ramps in time between consecutive events. NaN outside the range.
    phi = nan(T,1);
    events = [peak_idx(:); foot_idx(:)];
    types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
    [events, ord] = sort(events);
    types = types(ord);
    % strict alternation: drop adjacent same-type duplicates
    keep = true(size(events));
    for i = 2:numel(events)
        if types(i) == types(i-1), keep(i) = false; end
    end
    events = events(keep); types = types(keep);
    if numel(events) < 2, return; end
    phases  = nan(size(events));
    phi_cur = types(1) * pi;        % type=1 (peak)->pi; type=0 (foot)->0
    for i = 1:numel(events)
        phases(i) = phi_cur; phi_cur = phi_cur + pi;
    end
    for i = 1:numel(events)-1
        a = events(i); b = events(i+1);
        if a < 1 || b > T || b <= a, continue; end
        phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
    end
end
