function Ventral_surface_breath_time_summary_260808()
%% Ventral_surface_breath_time_summary_260731
% Population + per-cell Ca-event histograms in TIME and in LINEARIZED PHASE,
% for every ACTIVE cell in the ventral surface archive, one file set per
% genotype, with IO standing alone.
%
% Output: D:\Ventral_surface_summary\breath_time_summary\
%
% ------------------------------- LAYOUT ----------------------------------
% Three columns. Row 1 is the POOLED POPULATION:
%
%          ALL              SIG            NON-SIG
%   r1  [time|phase]    [time|phase]    [time|phase]     <- pooled over cells
%   ---------------------------------------------------
%   r2  [cell 1    ]    [cell 2    ]    [cell 3    ]     <- SIG cells first,
%   r3  [cell 4    ]    [cell 5    ]    [cell 6    ]        then NON-SIG
%   ...                                                     left->right,
%                                                           top->bottom
%
% Cells overflow onto further pages (`_p02`, `_p03`, ...); pages after the
% first carry cells only, no population row.
%
% Every slot is a PAIR of axes: the time histogram on the left, the phase
% histogram on the right, for the same cell (or the same pooled population).
%
% ---------------------------- COLOUR GRAMMAR -----------------------------
%   dark grey bars   pooled population
%   grey bars        one cell
%   RED              translucent histogram of the inspiration ONSETS,
%                    peak-triggered like everything else
%   SKY BLUE         inspiratory PEAK   -- vertical marker at t = 0
%   * + bold title   significant cell (and a heavier axes frame)
%                    -- deliberately NOT a red frame: red already means onset
%
% ----------------------------- TIME PANEL --------------------------------
% EVERYTHING IS PEAK-TRIGGERED. t = 0 is the inspiratory peak (sky blue line),
% one meaning, and the window runs +/- winIBI median breath cycles either side.
% Two histograms share that one trigger:
%
%   grey            Ca events of this cell (or of the pooled population)
%   translucent red BREATH ONSETS of the same recording(s)
%
% The red is the respiratory landmark itself, drawn as a distribution rather
% than a line, so you can see where inspiration begins relative to the peak AND
% how tightly it is locked. Read the grey against the red: that is the whole
% point of the panel. Nothing here is fitted or modelled.
%
% The red is NOT smoothed -- the width of that distribution IS the onset jitter,
% and smoothing would invent width. It is scaled to the grey panel's maximum,
% because onsets are locked to a few bins while Ca events spread over the cycle,
% so a shared percentage axis would let the red tower over and flatten the data.
% Only the red's POSITION and WIDTH mean anything, never its height.
% Set showOnsetOverlay = false for the Ca histogram alone.
%
% Window = +/- winIBI x THAT GENOTYPE'S median breath cycle. This is a
% deliberate choice and it has a cost. Vgat/0730 breathes at ~0.4 Hz, five
% times slower than everything else in the archive (median cycle 2.53 s vs
% 0.50 s for Vglut2), so a window fixed in seconds would show Vglut2 many
% cycles and Vgat less than one, and Vgat would read as flat. Scaling to each
% genotype's own IBI shows every genotype the same number of breaths -- but
% SECONDS ARE THEN NOT COMPARABLE BETWEEN GENOTYPES. The page title carries the
% median cycle so this is never read off by mistake. The phase panels are
% unaffected: phase is normalised per cycle, so those ARE comparable.
%
% ---------------------------- PHASE PANEL --------------------------------
% Linearized breath phase, piecewise-linear in time between landmarks:
% 0 = inspiration onset (red), pi = inspiratory peak (sky blue), 2pi = next
% onset. Each half-cycle gets exactly pi regardless of its duration.
%
% THE PHASE HISTOGRAM MUST BE OCCUPANCY-NORMALISED BEFORE IT IS TURNED INTO A
% PERCENTAGE. Inspiration (0..pi) is only ~15-21% of the cycle in TIME, so raw
% event counts per phase bin are ~7x denser in expiration for a completely
% unmodulated cell -- a raw percentage would show every cell as expiratory.
% So the bars are built in two steps:
%     rate(b) = events_in_bin(b) / (frames_in_bin(b) / fps)     events per second
%     pct(b)  = 100 * rate(b) / sum(rate)                       spk/cyc %
% Occupancy is accumulated per recording and applied before pooling, because
% each recording has its own cycle shape.
%
% -------------------------------- UNITS ----------------------------------
% Both panels are in SPK/CYC %: each bin as a percentage of that cell's events
% within one cycle, so the bars sum to 100 across the plotted cycle. A 7-event
% cell and a 315-event cell are then directly comparable. The absolute event
% count is in every panel title (n=...) and in the CSV.
%
% --------------------------- SIGNIFICANCE --------------------------------
% SIG = occupancy-weighted Rayleigh, analytic p < alpha_ray (0.001).
%
% The weighting is 1/(frames in that event's phase bin) in that event's OWN
% recording, and the sample size is Kish's effective n, nEff = (sum w)^2 /
% sum(w^2), NOT the raw event count. Using raw n with weighted Rbar returned
% p = 5e-79 where the shuffle said 0.05. Same maths as
% Ventral_surface_polar_coh_vs_rayleigh_260729.m, recomputed here rather than
% joined from its CSV so the label always matches the cell set actually drawn.
%
% NOTE ON THE ANALYTIC p. A circular-shift shuffle p has a floor of
% 1/(1+nShuffleRay); at 500 shuffles that is 0.002 and so it can never express
% alpha = 0.001. The analytic p (Zar 1999 eq 27.4 on nEff) can, and it is what
% decides the split. The shuffle p is still computed and written to the CSV as
% an independent cross-check. Raise nShuffleRay to >= 5000 to gate on it.
%
% -------------------------- ACTIVE CELL ----------------------------------
%     active cell <=> nnz(spike_train > 0) > activeMinEv   (i.e. >= 6 events)
% pooled over every recording the cell appears in. This is the archive's
% existing criterion (spike_trigger_dFF.m:34, temporal_phase_perROI.m:56).
% Deliberately looser than the >= 20 events the per-cell PETH test needs:
% pooling means a 6-event cell still contributes 6 real events, whereas a
% per-cell shuffle test on 6 events has no power at all. Cells below the gate
% are excluded entirely -- the user asked for the active cells.
%
% ---------------------------- CELL IDENTITY ------------------------------
% A MASK IS NOT A CELL. The same neuron re-imaged at a second Z or zoom gives
% one Cellpose mask per recording. Identity is resolved through each session's
% cell_link.mat, read from the ARCHIVE's own copies (the acquisition drives
% moved D: -> E: and E: gets unplugged; a summary that silently loses cell
% identities when a drive is absent is worse than one that refuses to run).
% Observations marked `tossed` in curation are dropped outright. Sessions with
% no matcher run fall back to one cell per ROI, so pooling is a no-op there.
%
% Runqi Zhang / 2026-07-31
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);           % detect_session_fps
addpath(scriptDir);
addpath(fullfile(scriptDir,'coh_ca_breath'));   % coh_cfg_260727, for the exclusion list

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath, 'breath_time_summary_260808');
% NEW folder on purpose: breath_time_summary is the reference set. This run
% differs only in that Sst\0806 and Sst\0807 were copied into the archive on
% 2026-08-08, so the Sst scan now also finds those 24 recordings. They have no
% cell_link entry, so each of their ROIs counts as its own cell (overcount).

% One file set per group. IO is split out of ChAT: any folder under a scanned
% genotype whose cell-level dir name contains 'IO' counts as IO, so IO stands
% alone and never pools with a genotype.
groups    = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};
scan_dirs = {'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};

% Cell identity, read from the archive's own copies. Sessions absent from this
% list fall back to one cell per ROI.
cell_link_sources = { ...
    'Sert',   '0721', fullfile(rootPath, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(rootPath, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(rootPath, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0806', fullfile(rootPath, 'Sst',    '0806', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0807', fullfile(rootPath, 'Sst',    '0807', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0810', fullfile(rootPath, 'Vglut2', '0810', 'cell_pooled', 'cell_link.mat')
    };

nDrop         = 30;      % breath frames tossed up front, to align with Ca
fallback_fps  = 30;
activeMinEv   = 5;       % ACTIVE = pooled nnz(spike_train>0) > this
ca_lag_sec    = 0.03;    % GCaMP lead: shift events EARLIER by this before use

% --- time panel ---
winIBI        = 2.0;     % window = +/- this many median breath cycles, around the PEAK
binWidth_s    = 0.050;   % fixed in SECONDS; must exceed the slowest frame period
smoothWidth_s = 0.150;   % fixed in SECONDS, so it means the same thing everywhere
ampFrac       = 0.20;    % trigger QC: amplitude only, long pauses never rejected
showOnsetOverlay = true; % overlay the ONSET-triggered histogram on the peak-centred axes
onsetAlpha       = 0.40; % its transparency, so the grey peak histogram stays readable

% --- phase panel ---
nPhaseBins    = 18;      % 20 deg bins.  WAS 36 (10 deg), changed 2026-08-11.
                         % This is a spike COUNT histogram, so its resolution is
                         % capped by the frame rate: on the fast sessions (0.53 s
                         % cycle at 30 fps) one frame spans 22.5 deg, so 16 bins is
                         % the hard limit and 36 was 2x oversampled -- the extra bins
                         % were interpolation-free empty resolution, not detail.
                         % The slow sessions (1.3-2.5 s cycles) could support 36, so
                         % this is a compromise set by the FASTEST session present.

% --- nulls ---
% nShuffle = 0 -> no grey band is drawn, bars only. Set it to e.g. 500 to get
% the circular-shift 95% band back behind the bars; everything else is unchanged.
nShuffle      = 0;
nShuffleRay   = 500;     % circular-shift shuffles for the Rayleigh cross-check p
alpha_ray     = 0.001;   % SIG <=> occupancy-weighted Rayleigh p_analytic < this
rngSeed       = 260731;

% --- figure ---
rowsPerPage   = 6;       % rows per page (page 1 spends one of them on the population)
colsPerRow    = 3;
doSave        = true;
% =====================================================================

excludeRecordings = coh_cfg_260727().excludeRecordings;   % different stage zero
rng(rngSeed);
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

COL_POP   = [0.25 0.25 0.25];    % pooled population bars, dark grey
COL_CELL  = [0.55 0.55 0.55];    % single-cell bars, grey
COL_ONSET = [0.85 0.10 0.10];    % RED  = inspiration onset (marker + overlay curve)
COL_PEAK  = [0.20 0.65 0.90];    % SKY  = inspiratory peak (marker at t = 0)
COL_BAND  = [0.82 0.82 0.82];    % shuffle band, only drawn when nShuffle > 0
COL_SHMU  = [0.45 0.45 0.45];
CL = struct('POP',COL_POP,'CELL',COL_CELL,'ONSET',COL_ONSET,'PEAK',COL_PEAK, ...
            'BAND',COL_BAND,'SHMU',COL_SHMU);

phEdges = linspace(0, 2*pi, nPhaseBins+1);
phCtrs  = phEdges(1:end-1) + diff(phEdges)/2;

fprintf('\n=========== Ventral_surface_breath_time_summary_260731 ===========\n');
fprintf('active = pooled >%d events | SIG = Rayleigh p_analytic < %.4g\n', ...
        activeMinEv, alpha_ray);
fprintf('time = PEAK-centred, +/-%.1f x genotype median IBI, %.0f ms bins, %.0f ms smoothing\n', ...
        winIBI, 1000*binWidth_s, 1000*smoothWidth_s);
fprintf('phase = %d bins, occupancy-normalised; both panels in spk/cyc %%\n\n', nPhaseBins);

%% ---- cell-identity lookup ----------------------------------------------
% Key '<scanDir>/<date>/<recName>/<maskL label>' -> cell key. Namespaced per
% session because cell ids restart at 1 in every matcher run, and cell_id is
% positional: splitting one dataset differently renumbers everything.
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

%% ---- SCAN: one pass over every recording in the archive -----------------
REC = struct('bidx',{},'occ',{},'T',{},'phiW',{},'fps',{}, ...
             'onTrig',{},'pkTrig',{},'opLag',{},'gi',{});
OBS = struct('rec',{},'ev',{},'gi',{},'label',{},'cellKey',{},'nEv',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), warning('Folder missing: %s', gdir); continue; end

    allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
    fprintf('\n=== [%s] %d recordings with spikes ===\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;
        recName    = folder_basename(folderPath);
        recDate    = date_from_path(folderPath, gdir);

        gname = sname;
        if is_io_path(folderPath, gdir), gname = 'IO'; end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), warning('No group for "%s" -- skipping %s', gname, recName); continue; end
        if any(strcmp(recName, excludeRecordings))
            fprintf('  skip (excluded: different stage zero): %s\n', recName); continue;
        end
        try
            bpFile = fullfile(folderPath, 'breath_peak_pc1.mat');
            ipFile = fullfile(folderPath, 'breath_insp_start_pc1.mat');
            if ~isfile(bpFile) || ~isfile(ipFile)
                fprintf('  skip (no breath peak/onset): %s\n', recName); continue;
            end

            fps  = detect_session_fps(folderPath, fallback_fps);
            CA   = load(fullfile(folderPath, 'ca_spike_data.mat'), 'roi_spikes');
            nROI = numel(CA.roi_spikes);
            nCa  = numel(CA.roi_spikes(1).spike_train);

            BP = load(bpFile);   % SVD PC1 peaks -> insp_onset_idx IS THE PEAK
            IP = load(ipFile);   % SVD PC1 feet  -> insp_start_idx IS THE ONSET
            bw = detrend(double(BP.breath(:)));
            bw(1:min(nDrop,numel(bw))) = [];
            bw = bw - mean(bw);

            peak_idx = round(BP.insp_onset_idx(:)) - nDrop;   % insp PEAK  (misleading field name)
            foot_idx = round(IP.insp_start_idx(:)) - nDrop;   % insp ONSET (foot/trough)

            % per-session timing fix: Vglut2/1124 used a rising-edge trigger, so
            % its breath samples sit one frame late relative to Ca.
            if strcmpi(sname,'Vglut2') && strcmp(recDate,'1124')
                peak_idx = peak_idx + 1;
                foot_idx = foot_idx + 1;
                bw = [bw(1); bw(1:end-1)];
            end

            T = min([numel(bw), nCa]);
            peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
            foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
            if numel(peak_idx) < 2 || numel(foot_idx) < 2
                fprintf('  skip (too few breath landmarks): %s\n', recName); continue;
            end
            bw = bw(1:T);

            % ---- linearized phase, and occupancy in FRAMES per phase bin ----
            % phi is NaN outside the first..last landmark; those frames are in no
            % cycle and must leave BOTH the event phases and the occupancy, or the
            % normalisation is against the wrong total.
            phi   = piecewise_phase_local(peak_idx, foot_idx, T);
            phiW  = mod(phi, 2*pi);
            valid = ~isnan(phiW);
            bidx  = nan(T,1);
            bidx(valid) = discretize(phiW(valid), phEdges);
            occ   = accumarray(bidx(valid), 1, [nPhaseBins 1]);

            % ---- trigger QC: amplitude only ----
            % Robust z, peak-minus-foot amplitude of each cycle, compared with THIS
            % recording's own median. A global median would conflate between-recording
            % brightness with within-recording pauses. Amplitude only -- a long
            % respiratory pause is a real breath and is never rejected.
            bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
            ft  = sort(foot_idx);
            amp = nan(numel(ft)-1,1);  pkf = nan(numel(ft)-1,1);
            for i = 1:numel(ft)-1
                q = peak_idx(peak_idx>ft(i) & peak_idx<ft(i+1));
                if ~isempty(q), amp(i) = bwz(q(1)) - bwz(ft(i)); pkf(i) = q(1); end
            end
            good   = amp > ampFrac*median(amp,'omitnan');
            onTrig = ft(good);
            pkPair = pkf(good);
            % onset-to-peak interval, from the SAME accepted pairs -- this is where
            % the red onset marker goes on the peak-centred axis.
            opLag  = median((pkPair - onTrig)/fps, 'omitnan');
            pkTrig = pkPair(~isnan(pkPair));

            REC(end+1) = struct('bidx',bidx,'occ',occ,'T',T,'phiW',phiW,'fps',fps, ...
                                'onTrig',onTrig(:),'pkTrig',pkTrig(:), ...
                                'opLag',opLag,'gi',gi); %#ok<AGROW>
            recIdx = numel(REC);

            % ---- ROIs ----
            % NO spike gate here. A cell's zero-event recording still contributes
            % triggers and exposure; dropping it inflates the cell's rate. The
            % ACTIVE gate is applied later, pooled over the whole cell.
            lag  = round(ca_lag_sec*fps);
            nInc = 0;
            for rid = 1:nROI
                st = double(CA.roi_spikes(rid).spike_train(:));
                st = st(1:min(T,numel(st)));
                if numel(st) < T, st(end+1:T,1) = 0; end %#ok<AGROW>

                % the matcher ran per session, BEFORE the IO folder override, so the
                % identity key uses the SCAN DIR while the display label uses gname.
                ckin = sprintf('%s/%s/%s/%d', sname, recDate, recName, rid);
                if isKey(tossed_set, ckin), continue; end   % rejected in curation: not a cell
                stL = [st(1+lag:end); zeros(lag,1)];        % lead-shift earlier

                lab = sprintf('%s/%s/%s/%d', gname, recDate, recName, rid);
                if isKey(cell_map, ckin), ck = cell_map(ckin); else, ck = ['roi:' lab]; end
                OBS(end+1) = struct('rec',recIdx,'ev',find(stL>0),'gi',gi, ...
                                    'label',lab,'cellKey',ck,'nEv',nnz(stL>0)); %#ok<AGROW>
                nInc = nInc + 1;
            end
            fprintf('  [%2d] %-52s %s  fps %.2f  %d ROI\n', kk, recName, gname, fps, nInc);
        catch ME
            warning('  ERROR %s: %s', recName, ME.message);
        end
    end
end
if isempty(OBS), error('No observations collected from %s', rootPath); end

%% ---- CELLS: pool observations, apply the ACTIVE gate --------------------
cellKeys = {OBS.cellKey};
[uCell, ~, obsOfCell] = unique(cellKeys, 'stable');
nCellAll = numel(uCell);

CELL = struct('key',{},'gi',{},'obs',{},'nEv',{},'nObs',{},'label',{});
for c = 1:nCellAll
    m  = find(obsOfCell == c);
    ne = sum([OBS(m).nEv]);
    if ne <= activeMinEv, continue; end          % ACTIVE gate, pooled
    if numel(m) == 1
        lab = OBS(m(1)).label;
    else
        lab = sprintf('%s [+%d rec]', OBS(m(1)).label, numel(m)-1);
    end
    CELL(end+1) = struct('key',uCell{c},'gi',OBS(m(1)).gi,'obs',m(:)', ...
                         'nEv',ne,'nObs',numel(m),'label',lab); %#ok<AGROW>
end
fprintf('\n%d ROI-observations -> %d cells -> %d ACTIVE (>%d events)\n', ...
        numel(OBS), nCellAll, numel(CELL), activeMinEv);
if isempty(CELL), error('No active cells.'); end

%% ---- per-genotype median breath cycle and onset-to-peak lag -------------
% Both drive the time axis: the cycle sets the window, the lag places the red
% onset marker. Computed over the ACTIVE cells' recordings only, so they match
% the data actually drawn.
medCyc = nan(1, numel(groups));
medOP  = nan(1, numel(groups));
for gi = 1:numel(groups)
    rk = unique([OBS([CELL([CELL.gi]==gi).obs]).rec]);
    cyc = [];  lags = [];
    for k = rk
        f = sort(REC(k).onTrig);
        if numel(f) > 1, cyc = [cyc; diff(f)/REC(k).fps]; end %#ok<AGROW>
        if isfinite(REC(k).opLag), lags = [lags; REC(k).opLag]; end %#ok<AGROW>
    end
    if ~isempty(cyc),  medCyc(gi) = median(cyc);  end
    if ~isempty(lags), medOP(gi)  = median(lags); end
end

%% ---- per-cell histograms and Rayleigh -----------------------------------
S = struct('cnt_pk',{},'cnt_bon',{},'sh_mu',{},'sh_var',{}, ...
           'ph_rate',{},'ph_shmu',{},'ph_shvar',{}, ...
           'Rbar',{},'Z',{},'logZ',{},'nEff',{},'th',{},'p_an',{},'p_sh',{},'sig',{}, ...
           'nTrigOn',{},'nTrigPk',{},'nEvUsed',{});

fprintf('\ncomputing per-cell histograms (%d cells)...\n', numel(CELL));
tStart = tic;
for c = 1:numel(CELL)
    gi = CELL(c).gi;
    [tEdges, ~, smB] = time_axis(medCyc(gi), winIBI, binWidth_s, smoothWidth_s);
    S(c) = cell_stats(CELL(c), OBS, REC, tEdges, smB, nPhaseBins, ...
                      nShuffle, nShuffleRay, alpha_ray); %#ok<AGROW>
    if mod(c, 50) == 0
        fprintf('  %d/%d  (%.0f s)\n', c, numel(CELL), toc(tStart));
    end
end
fprintf('  done in %.0f s\n', toc(tStart));

fprintf('\nSIG (Rayleigh p_analytic < %.4g): %d of %d active cells\n', ...
        alpha_ray, nnz([S.sig]), numel(CELL));

%% ---- figures, one file set per group ------------------------------------
if doSave && ~isfolder(outDir), mkdir(outDir); end
rowsCSV = {};

for gi = 1:numel(groups)
    idx = find([CELL.gi] == gi);
    if isempty(idx)
        fprintf('\n[%s] no active cells -- skipped\n', groups{gi}); continue;
    end
    sig = [S(idx).sig];
    % SIG block first, then NON-SIG; within each, strongest modulation first so
    % the leading panels are the informative ones.
    [~, oS] = sort([S(idx(sig)).logZ],  'descend');
    [~, oN] = sort([S(idx(~sig)).logZ], 'descend');
    iSig = idx(sig);   iSig = iSig(oS);
    iNon = idx(~sig);  iNon = iNon(oN);
    ord  = [iSig, iNon];

    [tEdges, tCtrs, smB] = time_axis(medCyc(gi), winIBI, binWidth_s, smoothWidth_s);

    POP = struct('name',{},'S',{},'n',{});
    POP(1) = struct('name','ALL',     'S',pool_cells(S(idx)),  'n',numel(idx));
    POP(2) = struct('name','SIG',     'S',pool_cells(S(iSig)), 'n',numel(iSig));
    POP(3) = struct('name','NON-SIG', 'S',pool_cells(S(iNon)), 'n',numel(iNon));

    slotsP1 = colsPerRow * (rowsPerPage - 1);   % page 1 spends row 1 on the population
    slotsPn = colsPerRow * rowsPerPage;
    nPage   = 1 + ceil(max(numel(ord) - slotsP1, 0) / slotsPn);

    fprintf('\n[%s] %d active cells (%d sig / %d non-sig), cycle %.3f s, onset->peak %.3f s -> %d page(s)\n', ...
            groups{gi}, numel(idx), numel(iSig), numel(iNon), medCyc(gi), medOP(gi), nPage);

    for pg = 1:nPage
        if pg == 1
            take = ord(1:min(slotsP1, numel(ord)));
            nRow = 1 + ceil(numel(take)/colsPerRow);
        else
            a = slotsP1 + (pg-2)*slotsPn + 1;
            take = ord(a : min(a+slotsPn-1, numel(ord)));
            nRow = ceil(numel(take)/colsPerRow);
        end
        nRow = max(nRow, 1);

        hf = figure('Color','w','Visible','off','Units','centimeters', ...
                    'Position',[1 1 34 3.35*nRow + 2.6]);
        set(hf,'DefaultAxesFontSize',7);

        yTop = 0.915; yBot = 0.050;
        rowH = (yTop - yBot)/nRow;
        % 6 sub-panels across, so the gaps carry the tick numbers: subGap
        % separates a pair's time and phase axes, colGap separates the columns.
        colW = 0.288; colGap = 0.028; X0 = 0.050;
        subGap = 0.050; subW = (colW - subGap)/2;

        rowsUsed = 0;
        % ---- row 1 of page 1: the pooled population ----
        if pg == 1
            for q = 1:3
                x0 = X0 + (q-1)*(colW + colGap);
                y0 = yTop - 1*rowH + rowH*0.30;
                axT = axes('Parent',hf,'Position',[x0,              y0, subW, rowH*0.60]); %#ok<LAXES>
                axP = axes('Parent',hf,'Position',[x0+subW+subGap,  y0, subW, rowH*0.60]); %#ok<LAXES>
                if POP(q).n == 0
                    axis(axT,'off'); axis(axP,'off');
                    text(axT, 0.5, 0.5, sprintf('%s: no cells', POP(q).name), ...
                         'Units','normalized','HorizontalAlignment','center','FontSize',8);
                    continue;
                end
                ttl = sprintf('%s  |  %s POP  |  %d cells, n=%d', ...
                              groups{gi}, POP(q).name, POP(q).n, POP(q).S.nEvUsed);
                draw_time( axT, tCtrs, POP(q).S, smB, COL_POP, COL_ONSET, COL_PEAK, ...
                           COL_BAND, COL_SHMU, medCyc(gi), tEdges, ttl, ...
                           true, q==1, showOnsetOverlay, onsetAlpha);
                draw_phase(axP, phCtrs, POP(q).S, COL_POP, COL_ONSET, COL_PEAK, ...
                           COL_BAND, COL_SHMU, true, q==1);
            end
            rowsUsed = 1;
        end

        % ---- cells: 3 per row, left->right, top->bottom ----
        for k = 1:numel(take)
            c   = take(k);
            r   = rowsUsed + ceil(k/colsPerRow);
            q   = mod(k-1, colsPerRow) + 1;
            x0  = X0 + (q-1)*(colW + colGap);
            y0  = yTop - r*rowH + rowH*0.30;
            axT = axes('Parent',hf,'Position',[x0,             y0, subW, rowH*0.60]); %#ok<LAXES>
            axP = axes('Parent',hf,'Position',[x0+subW+subGap, y0, subW, rowH*0.60]); %#ok<LAXES>

            isLast = (r == nRow);
            mark = '';
            if S(c).sig, mark = ' *'; end
            % the group is already in the page title, so the label drops it
            ttl = sprintf('%s%s   n=%d  logZ %.2f  p=%.1e', ...
                          disp_label(CELL(c).label), mark, S(c).nEvUsed, S(c).logZ, S(c).p_an);
            draw_time( axT, tCtrs, S(c), smB, COL_CELL, COL_ONSET, COL_PEAK, ...
                       COL_BAND, COL_SHMU, medCyc(gi), tEdges, ttl, ...
                       isLast, q==1, showOnsetOverlay, onsetAlpha);
            draw_phase(axP, phCtrs, S(c), COL_CELL, COL_ONSET, COL_PEAK, ...
                       COL_BAND, COL_SHMU, isLast, q==1);
            if S(c).sig
                % Significant cells get a HEAVY BLACK frame and a bold title, not a
                % coloured one: red is already spoken for as the inspiration-onset
                % marker, and a red frame beside a red landmark line is unreadable.
                % The sig block is also sorted to the front, so ordering carries most
                % of the findability and the frame only has to confirm it.
                set([axT axP], 'LineWidth',1.6);
                set(get(axT,'Title'), 'FontWeight','bold');
            end
        end

        sgtitle({sprintf(['%s   |   page %d/%d   |   %d active cells (%d sig, %d non-sig)   ' ...
                          '|   spikes/cycle %%   |   window +/-%.1f x median cycle (%.3f s)'], ...
                         groups{gi}, pg, nPage, numel(idx), numel(iSig), numel(iNon), ...
                         winIBI, medCyc(gi)), ...
                 sprintf(['LEFT = time, ALL PEAK-TRIGGERED (t=0 = insp peak, sky line):   ' ...
                          'grey = Ca events,   translucent red = BREATH ONSETS (scaled to panel)   |   ' ...
                          'RIGHT = linearized phase, occupancy-normalised']), ...
                 sprintf('* , bold title and heavy frame = significant, occupancy-weighted Rayleigh p_analytic < %.4g (Kish n_eff)', alpha_ray)}, ...
                'Interpreter','none','FontSize',8);

        if doSave
            base = fullfile(outDir, sprintf('%s_p%02d', groups{gi}, pg));
            exportgraphics(hf, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
            exportgraphics(hf, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            fprintf('   saved %s.png/.pdf  (%d rows, %d cells)\n', ...
                    folder_basename(base), nRow, numel(take));
        end
        close(hf);
    end

    %% ---- population-only figure: one row, ALL | SIG | NON-SIG ----------
    % The same three pooled panels as row 1 of page 1, on their own and drawn
    % large enough to read. Nothing is recomputed -- these are the identical
    % POP structs, so this figure can never disagree with the per-cell pages.
    hp = pop_figure(POP, tCtrs, tEdges, smB, phCtrs, medCyc(gi), groups{gi}, ...
                    numel(idx), numel(iSig), numel(iNon), winIBI, alpha_ray, CL, ...
                    showOnsetOverlay, onsetAlpha, '');
    if doSave
        base = fullfile(outDir, sprintf('%s_population', groups{gi}));
        exportgraphics(hp, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
        exportgraphics(hp, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        fprintf('   saved %s.png/.pdf\n', folder_basename(base));
    end
    close(hp);

    for k = 1:numel(ord)
        c = ord(k);
        rowsCSV(end+1,:) = {groups{gi}, CELL(c).key, CELL(c).label, CELL(c).nObs, ...
                            S(c).nEvUsed, S(c).nEff, S(c).Rbar, S(c).Z, S(c).logZ, ...
                            rad2deg(S(c).th), S(c).p_an, S(c).p_sh, double(S(c).sig), ...
                            medCyc(gi), medOP(gi), k}; %#ok<AGROW>
    end
end

%% ---- GRAND population: every active cell, every animal, every genotype ---
% PHASE POOLING HERE IS LEGITIMATE, ABSOLUTE-TIME POOLING IS NOT. Phase is
% normalised per cycle, so 0 and pi mean the same landmark in every animal and
% the phase panel is a genuine pooled result. Seconds are not: Vgat/0730
% breathes at 2.53 s per cycle and Vglut2 at 0.50 s, a 5x spread, so one second
% axis puts different fractions of the respiratory cycle into the same bin. The
% grand TIME panel is drawn on +/-winIBI x the median cycle over ALL recordings
% and is therefore dominated by whichever genotypes breathe near that median --
% read it as a summary of mixed cycle lengths, NOT as a latency measurement.
% IO is pooled in here too; it still has its own per-group figure.
allCyc = [];
for k = unique([OBS([CELL.obs]).rec])
    f = sort(REC(k).onTrig);
    if numel(f) > 1, allCyc = [allCyc; diff(f)/REC(k).fps]; end %#ok<AGROW>
end
grandCyc = median(allCyc);
[gEdges, gCtrs, gSmB] = time_axis(grandCyc, winIBI, binWidth_s, smoothWidth_s);

fprintf('\n[ALL GENOTYPES] pooling %d active cells onto a common %.3f s cycle\n', ...
        numel(CELL), grandCyc);

% Recompute each cell's TIME histogram on the common grid. The phase histogram
% and the Rayleigh are grid-independent and come out identical, so nShuffleRay
% is 0 here -- the sig labels are taken from the first pass, never re-derived.
SG = S(1);
for c = 1:numel(CELL)
    SG(c) = cell_stats(CELL(c), OBS, REC, gEdges, gSmB, nPhaseBins, ...
                       nShuffle, 0, alpha_ray); %#ok<AGROW>
end
sigAll = [S.sig];
POPG(1) = struct('name','ALL',     'S',pool_cells(SG),          'n',numel(SG));
POPG(2) = struct('name','SIG',     'S',pool_cells(SG(sigAll)),  'n',nnz(sigAll));
POPG(3) = struct('name','NON-SIG', 'S',pool_cells(SG(~sigAll)), 'n',nnz(~sigAll));

hp = pop_figure(POPG, gCtrs, gEdges, gSmB, phCtrs, grandCyc, 'ALL GENOTYPES', ...
                numel(SG), nnz(sigAll), nnz(~sigAll), winIBI, alpha_ray, CL, ...
                showOnsetOverlay, onsetAlpha, ...
                'phase pools validly across animals; the TIME panel mixes 0.50-2.53 s cycles');
if doSave
    base = fullfile(outDir, 'ALL_GENOTYPES_population');
    exportgraphics(hp, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(hp, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('   saved ALL_GENOTYPES_population.png/.pdf\n');
end
close(hp);

%% ---- CSV ---------------------------------------------------------------
if doSave && ~isempty(rowsCSV)
    Tc = cell2table(rowsCSV, 'VariableNames', ...
        {'group','cell_key','label','n_recordings','n_events','n_eff','ray_Rbar', ...
         'ray_Z','ray_logZ','ray_mean_deg','ray_p_analytic','ray_p_shuffle','sig', ...
         'genotype_median_cycle_s','genotype_onset_to_peak_s','panel_order'});
    writetable(Tc, fullfile(outDir,'breath_time_summary_cells.csv'));
    fprintf('\nWrote breath_time_summary_cells.csv (%d active cells)\n', height(Tc));
end

fprintf('\nOutput: %s\n', outDir);
fprintf('==================================================================\n\n');
end

% =========================================================================
%                              HELPERS
% =========================================================================

function hp = pop_figure(POP, tCtrs, tEdges, smB, phCtrs, medCyc, name, ...
                         nAll, nSig, nNon, winIBI, alphaRay, CL, showOverlay, ...
                         onsetAlpha, note)
% One row of pooled panels: ALL | SIG | NON-SIG, each as a [time | phase] pair.
% Shared by the per-genotype figures and the grand all-animal figure so the two
% can never be drawn by different code.
hp = figure('Color','w','Visible','off','Units','centimeters', ...
            'Position',[1 1 34 8.4]);
set(hp,'DefaultAxesFontSize',8);
colW = 0.288; colGap = 0.028; X0 = 0.050;
subGap = 0.050; subW = (colW - subGap)/2;
for q = 1:3
    x0  = X0 + (q-1)*(colW + colGap);
    axT = axes('Parent',hp,'Position',[x0,             0.19, subW, 0.55]); %#ok<LAXES>
    axP = axes('Parent',hp,'Position',[x0+subW+subGap, 0.19, subW, 0.55]); %#ok<LAXES>
    if POP(q).n == 0
        axis(axT,'off'); axis(axP,'off');
        text(axT, 0.5, 0.5, sprintf('%s: no cells', POP(q).name), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',9);
        continue;
    end
    ttl = sprintf('%s  |  %s  |  %d cells, n=%d events', ...
                  name, POP(q).name, POP(q).n, POP(q).S.nEvUsed);
    draw_time( axT, tCtrs, POP(q).S, smB, CL.POP, CL.ONSET, CL.PEAK, ...
               CL.BAND, CL.SHMU, medCyc, tEdges, ttl, true, q==1, ...
               showOverlay, onsetAlpha);
    draw_phase(axP, phCtrs, POP(q).S, CL.POP, CL.ONSET, CL.PEAK, ...
               CL.BAND, CL.SHMU, true, q==1);
end
% the panel fonts are sized for a 6-row page; on a single row they can breathe
axAll = findall(hp,'Type','axes');
set(axAll, 'FontSize',8);
for a = axAll(:)'
    set(get(a,'Title'),  'FontSize',8.5);
    set(get(a,'XLabel'), 'FontSize',8.5);
    set(get(a,'YLabel'), 'FontSize',8.5);
end
lines = {sprintf(['%s POPULATION   |   %d active cells (%d sig, %d non-sig)   |   ' ...
                  'spikes/cycle %%   |   window +/-%.1f x median cycle (%.3f s)'], ...
                 name, nAll, nSig, nNon, winIBI, medCyc), ...
         sprintf(['LEFT = time, ALL PEAK-TRIGGERED (t=0 = insp peak, sky line):   ' ...
                  'dark grey = Ca events,   translucent red = BREATH ONSETS (scaled to panel)   |   ' ...
                  'RIGHT = linearized phase, occupancy-normalised']), ...
         sprintf('SIG = occupancy-weighted Rayleigh p_analytic < %.4g (Kish n_eff)', alphaRay)};
if ~isempty(note), lines{end+1} = note; end
sgtitle(lines, 'Interpreter','none','FontSize',9);
end

function [edges, ctrs, smB] = time_axis(medCyc, winIBI, binWidth_s, smoothWidth_s)
% Symmetric window of +/- winIBI median cycles about the PEAK, snapped to whole
% bins so that t = 0 sits on a bin EDGE and the histogram is not half a bin out.
if ~isfinite(medCyc) || medCyc <= 0, medCyc = 0.5; end
nHalf = max(round(winIBI*medCyc/binWidth_s), 4);
edges = (-nHalf:nHalf) * binWidth_s;
ctrs  = edges(1:end-1) + binWidth_s/2;
smB   = max(1, round(smoothWidth_s/binWidth_s));
end

function S = cell_stats(C, OBS, REC, tEdges, smB, nPh, nSh, nShRay, alphaRay)
% Every histogram and statistic for ONE cell, pooled over its observations.
% Phase binning is NOT redone here: REC(k).bidx already maps every frame to its
% phase bin, computed once per recording at scan time with the same edges.
nB = numel(tEdges) - 1;
S = struct('cnt_pk',zeros(1,nB),'cnt_bon',zeros(1,nB), ...
           'sh_mu',zeros(1,nB),'sh_var',zeros(1,nB), ...
           'ph_rate',zeros(1,nPh),'ph_shmu',zeros(1,nPh),'ph_shvar',zeros(1,nPh), ...
           'Rbar',0,'Z',0,'logZ',log(eps),'nEff',0,'th',NaN,'p_an',1,'p_sh',NaN, ...
           'sig',false,'nTrigOn',0,'nTrigPk',0,'nEvUsed',0);

phCntObs = zeros(1,nPh);  occSec = zeros(1,nPh);
phSH     = zeros(nSh, nPh);
A = []; W = []; parts = struct('ev',{},'rec',{});

for j = C.obs
    o  = OBS(j);
    r  = REC(o.rec);
    T  = r.T;  fs = r.fps;
    ev = zeros(T,1);
    e  = o.ev(o.ev>=1 & o.ev<=T);
    ev(e) = 1;

    % ---- TIME: EVERYTHING IS PEAK-TRIGGERED ----
    % which = 1: this cell's Ca events, relative to the inspiratory peak.
    % which = 2: the BREATH ONSETS of the same recording(s), relative to the
    %            same peak. That second histogram is the respiratory landmark
    %            itself -- where inspiration begins around the peak and how
    %            tightly -- so the Ca events can be read against the breath on
    %            one time base with one meaning for t = 0.
    onTrain = zeros(T,1);
    oi = r.onTrig(r.onTrig>=1 & r.onTrig<=T);
    onTrain(oi) = 1;
    for which = 1:2
        if which == 1, srcTrain = ev; else, srcTrain = onTrain; end
        [cnt, ccf, lagIdx, binOf, nTrig] = one_peth(srcTrain, r.pkTrig, T, fs, tEdges, nB);
        if which == 1
            S.cnt_pk  = S.cnt_pk + cnt;   S.nTrigPk = S.nTrigPk + nTrig;
            if nTrig > 0 && nSh > 0
                % Circular shift of the EVENT train: preserves event count and burst
                % structure exactly. Only boundary-safe triggers are in the trigger
                % vector, so the circular correlogram IS the linear one over the
                % plotted lag range, and the shift stays legitimate.
                sh = zeros(nSh, nB);
                for s = 1:nSh
                    d = randi(T) - 1;
                    sh(s,:) = accumarray(binOf, ccf(mod(lagIdx-1+d,T)+1), [nB 1])';
                end
                sh = movmean(sh, smB, 2);   % smooth the null exactly as the bars are smoothed
                S.sh_mu  = S.sh_mu  + mean(sh,1);
                S.sh_var = S.sh_var + var(sh,0,1);
            end
        else
            S.cnt_bon = S.cnt_bon + cnt;  S.nTrigOn = S.nTrigOn + numel(oi);
        end
    end

    % ---- PHASE: counts, occupancy, weights ----
    e2 = e(~isnan(r.bidx(e)));
    if ~isempty(e2)
        b = r.bidx(e2);
        phCntObs = phCntObs + accumarray(b, 1, [nPh 1])';
        A = [A; r.phiW(e2(:))];                     %#ok<AGROW>
        W = [W; 1./max(r.occ(b(:)),1)];             %#ok<AGROW>
        parts(end+1) = struct('ev',e2,'rec',o.rec); %#ok<AGROW>
    end
    occSec = occSec + (r.occ(:)') / fs;             % occupancy in SECONDS
    S.nEvUsed = S.nEvUsed + numel(e2);

    if nSh > 0 && ~isempty(e2)
        for s = 1:nSh
            d  = randi(T) - 1;
            es = mod(e2 - 1 + d, T) + 1;
            es = es(~isnan(r.bidx(es)));
            if isempty(es), continue; end
            phSH(s,:) = phSH(s,:) + accumarray(r.bidx(es), 1, [nPh 1])';
        end
    end
end

% occupancy-normalised RATE: events per SECOND spent in that phase bin. Raw
% counts here would be ~7x denser in expiration for an unmodulated cell. The
% conversion to spk/cyc %% happens at draw time, AFTER this correction.
den        = max(occSec, eps);
S.ph_rate  = phCntObs ./ den;
if nSh > 0
    phSHr      = phSH ./ den;
    S.ph_shmu  = mean(phSHr, 1);
    S.ph_shvar = var(phSHr, 0, 1);
end

% ---- occupancy-weighted Rayleigh, Kish effective n ----
if ~isempty(A)
    [th, Rbar, nEff] = wresultant(A, W);
    Z = nEff * Rbar^2;
    % Zar (1999) eq 27.4, evaluated on nEff -- NOT on the raw event count.
    p_an = exp(-Z) * (1 + (2*Z - Z^2)/(4*nEff) ...
                        - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4)/(288*nEff^2));
    S.th = th; S.Rbar = Rbar; S.nEff = nEff; S.Z = Z;
    S.logZ = log(max(Z,eps));
    S.p_an = min(max(p_an,0),1);
    S.sig  = S.p_an < alphaRay;

    % independent cross-check only; its floor is 1/(1+nShRay), so at 500 it
    % cannot express alpha = 0.001 and does NOT decide the split.
    if nShRay > 0 && ~isempty(parts)
        Zn = zeros(nShRay,1);
        for s = 1:nShRay
            as = []; ws = [];
            for q = 1:numel(parts)
                r  = REC(parts(q).rec);
                d  = round(s * r.T / (nShRay+1)) + q;   % a different shift per recording
                es = mod(parts(q).ev - 1 + d, r.T) + 1;
                es = es(~isnan(r.bidx(es)));
                if isempty(es), continue; end
                as = [as; r.phiW(es(:))];               %#ok<AGROW>
                ws = [ws; 1./max(r.occ(r.bidx(es(:))),1)]; %#ok<AGROW>
            end
            if isempty(as), continue; end
            [~, Rs, nEs] = wresultant(as, ws);
            Zn(s) = nEs * Rs^2;
        end
        S.p_sh = (1 + nnz(Zn >= Z)) / (1 + nShRay);
    end
end
end

function [cnt, ccf, lagIdx, binOf, nTrig] = one_peth(ev, f, T, fs, edges, nB)
% Trigger-aligned event count by FFT cross-correlation.
%   ccf(m+1) = sum_t trig(t) * ev(t+m)
% Only BOUNDARY-SAFE triggers enter `trig`, which makes the circular
% correlation exactly the linear one over the plotted lag range -- while still
% letting the shuffle use a circular shift with no wraparound artefact.
lo = ceil(edges(1)*fs);  hi = floor(edges(end)*fs);
f  = f(f+lo>=1 & f+hi<=T);
nTrig = numel(f);
cnt = zeros(1,nB);
m   = (lo:hi)';  bo = discretize(m/fs, edges);
kp  = ~isnan(bo);  m = m(kp);  bo = bo(kp);
lagIdx = mod(m,T)+1;  binOf = bo;
if nTrig < 1, ccf = zeros(T,1); return; end
trig = zeros(T,1);  trig(f) = 1;
ccf  = real(ifft(conj(fft(trig)) .* fft(ev)));
cnt  = accumarray(bo, ccf(lagIdx), [nB 1])';
end

function P = pool_cells(Sset)
% Pool a set of per-cell results into one population. Counts and rates add;
% VARIANCES add because the cells' shuffles are independent -- percentiles do
% NOT, so a pooled band (when drawn) is a 1.96-sigma normal approximation
% rather than an empirical quantile. The percentage conversion happens at draw
% time on the pooled totals, so a cell contributes in proportion to its events.
if isempty(Sset)
    P = struct('cnt_pk',0,'cnt_bon',0,'sh_mu',0,'sh_var',0, ...
               'ph_rate',0,'ph_shmu',0,'ph_shvar',0,'nEvUsed',0,'sig',false, ...
               'logZ',NaN,'p_an',NaN,'nTrigOn',0,'nTrigPk',0);
    return;
end
P = struct('cnt_pk',  sum(cat(1,Sset.cnt_pk),1), ...
           'cnt_bon', sum(cat(1,Sset.cnt_bon),1), ...
           'sh_mu',   sum(cat(1,Sset.sh_mu),1), ...
           'sh_var',  sum(cat(1,Sset.sh_var),1), ...
           'ph_rate', sum(cat(1,Sset.ph_rate),1), ...
           'ph_shmu', sum(cat(1,Sset.ph_shmu),1), ...
           'ph_shvar',sum(cat(1,Sset.ph_shvar),1), ...
           'nEvUsed', sum([Sset.nEvUsed]), 'sig', false, ...
           'logZ', NaN, 'p_an', NaN, ...
           'nTrigOn', sum([Sset.nTrigOn]), 'nTrigPk', sum([Sset.nTrigPk]));
end

function draw_time(ax, ctrs, S, smB, colBar, colOn, colPk, colBand, colShMu, ...
                   medCyc, edges, ttl, isLast, showY, showOverlay, onsetAlpha)
% Peak-centred: t = 0 IS the inspiratory peak (sky blue). The onset sits one
% onset-to-peak interval earlier (red). Bars are spk/cyc %.
hold(ax,'on'); box(ax,'on');
pk  = movmean(S.cnt_pk, smB);
sc  = 100 / max(sum(S.cnt_pk), eps);
pk  = pk * sc;
top = max([pk, eps]);

% optional shuffle band, scaled by the same factor as the bars
if any(S.sh_var > 0)
    sd = sqrt(S.sh_var) * sc;
    mu = S.sh_mu * sc;
    fill(ax, [ctrs fliplr(ctrs)], [mu+1.96*sd fliplr(max(mu-1.96*sd,0))], ...
         colBand, 'EdgeColor','none');
    plot(ax, ctrs, mu, '-', 'Color',colShMu, 'LineWidth',0.8);
    top = max([top, mu+1.96*sd]);
end

bar(ax, ctrs, pk, 1, 'FaceColor',colBar, 'EdgeColor','none');

% BREATH ONSET histogram, on the SAME peak trigger, in translucent red. This is
% the respiratory landmark itself: where inspiration begins relative to the peak,
% and how much that varies. Two deliberate choices:
%   - NOT smoothed. The width of this distribution IS the onset jitter, and
%     150 ms of smoothing would invent width that is not there.
%   - Scaled to the Ca panel's own maximum, not left in its own percentage.
%     Onsets are locked to within a few bins while Ca events are spread over the
%     whole cycle, so a shared percentage axis would make the red spike tower
%     over the grey and flatten the data the panel is actually about. Only the
%     POSITION and WIDTH of the red carry meaning, never its height.
if showOverlay && sum(S.cnt_bon) > 0
    bo = S.cnt_bon / max(S.cnt_bon) * top;
    bar(ax, ctrs, bo, 1, 'FaceColor',colOn, 'EdgeColor','none', 'FaceAlpha',onsetAlpha);
end

% Only the PEAK gets a marker line, at t = 0, because that is the trigger. The
% onset needs no line: the red histogram shows the whole onset distribution,
% which a single median line could only misrepresent as exact. Neighbouring
% cycles are marked faintly so a side-peak can be located.
faint = @(c) 1 - 0.55*(1 - c);
for k = -3:3
    xk = k*medCyc;
    if abs(xk) > edges(end), continue; end
    if k == 0, xline(ax, 0,  '-', 'Color',colPk,        'LineWidth',1.4);
    else,      xline(ax, xk, ':', 'Color',faint(colPk), 'LineWidth',0.6);
    end
end

xlim(ax, [edges(1) edges(end)]);
ylim(ax, [0 top*1.12]);
if isLast, xlabel(ax, 't from insp peak (s)', 'FontSize',6.5); else, set(ax,'XTickLabel',[]); end
% ylabel only on the leftmost column: at six sub-panels across, a label on every
% axes runs into the neighbouring column's tick numbers. The tick VALUES stay on
% every panel, which is what carries the scale.
if showY, ylabel(ax, 'spk/cyc %', 'FontSize',6.5); end
if ~isempty(ttl)
    % Left-aligned and anchored to the axes' left edge so the title runs across
    % the whole time+phase PAIR instead of being centred on (and overflowing)
    % the narrow time axes alone.
    title(ax, ttl, 'Interpreter','none', 'FontSize',6, 'FontWeight','normal', ...
          'HorizontalAlignment','left', 'Units','normalized', 'Position',[0 1.03 0]);
end
set(ax,'TickDir','out','FontSize',6);
end

function draw_phase(ax, ctrs, S, colBar, colOn, colPk, colBand, colShMu, isLast, showY)
% Occupancy-normalised rate, then rescaled to spk/cyc %. The occupancy step
% must come first or an unmodulated cell reads as expiratory.
hold(ax,'on'); box(ax,'on');
sc  = 100 / max(sum(S.ph_rate), eps);
p   = S.ph_rate * sc;
top = max([p, eps]);

if any(S.ph_shvar > 0)
    sd = sqrt(S.ph_shvar) * sc;
    mu = S.ph_shmu * sc;
    fill(ax, [ctrs fliplr(ctrs)], [mu+1.96*sd fliplr(max(mu-1.96*sd,0))], ...
         colBand, 'EdgeColor','none');
    plot(ax, ctrs, mu, '-', 'Color',colShMu, 'LineWidth',0.8);
    top = max([top, mu+1.96*sd]);
end

bar(ax, ctrs, p, 1, 'FaceColor',colBar, 'EdgeColor','none');
xline(ax, 0,      '-', 'Color',colOn, 'LineWidth',1.4);   % onset
xline(ax, 2*pi,   '-', 'Color',colOn, 'LineWidth',1.4);   % next onset
xline(ax, pi,     '-', 'Color',colPk, 'LineWidth',1.4);   % peak
xlim(ax, [0 2*pi]);
ylim(ax, [0 top*1.12]);
set(ax,'XTick',[0 pi/2 pi 3*pi/2 2*pi]);
if isLast
    set(ax,'XTickLabel',{'0','','\pi','','2\pi'});
    xlabel(ax, 'breath phase', 'FontSize',6.5);
else
    set(ax,'XTickLabel',[]);
end
if showY, ylabel(ax, 'spk/cyc %', 'FontSize',6.5); end
set(ax,'TickDir','out','FontSize',6);
end

function [th, Rbar, nEff] = wresultant(a, w)
% Weighted circular resultant with Kish's effective sample size. Feeding a
% weighted Rbar the RAW event count returned p = 5e-79 where the shuffle said
% 0.05 -- nEff is the correction that makes the analytic p usable.
S1 = sum(w);  S2 = sum(w.^2);
if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
nEff = S1^2 / max(S2, eps);
v    = sum(w(:) .* exp(1i*a(:))) / S1;
th   = angle(v);
Rbar = min(abs(v), 1);
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase: FEET at 0/2pi/..., PEAKS at pi/3pi/..., linear ramps
% in time between consecutive landmarks. NaN outside the landmark range.
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
keep = true(size(events));                     % strict alternation
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases  = nan(size(events));
phi_cur = types(1) * pi;                       % peak -> pi, foot -> 0
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
rel = strrep(folderPath, groupRoot, '');
rel = regexprep(rel, '^[\\/]+', '');
parts = regexp(rel, '[\\/]', 'split');
if isempty(parts), d = ''; else, d = parts{1}; end
end

function tf = is_io_path(folderPath, groupRoot)
% True if the cell-level dir (right after the date dir) contains 'IO'.
%   <groupRoot>\<date>\<cell-or-IO>\<recording>\
rel = strrep(folderPath, groupRoot, '');
rel = regexprep(rel, '^[\\/]+', '');
parts = regexp(rel, '[\\/]', 'split');
tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; this rebuilds the full last
% path segment for folders with dots in the name (e.g. "1.7x", "15.5lp").
p = char(p);
while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
[~,n,e] = fileparts(p);
name = [n e];
end

function s = disp_label(x)
% Panel label for one cell. The group is already in the page title, so it is
% dropped, leaving '<date>/<recording>/<roi>'. The recording name is the long
% part and is elided in the MIDDLE -- its head (roiN_zoom) and tail (the
% _0000N run index) are what distinguish two recordings of the same site.
x = char(x);
p = regexp(x, '/', 'split');
if numel(p) >= 4
    rec = p{end-1};
    if numel(rec) > 26, rec = [rec(1:13) '..' rec(end-10:end)]; end
    s = sprintf('%s/%s/%s', p{end-2}, rec, p{end});
else
    s = x;
end
if numel(s) > 46, s = [s(1:22) '..' s(end-21:end)]; end
end
