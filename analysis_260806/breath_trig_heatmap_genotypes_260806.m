% breath_trig_heatmap_genotypes_260806.m
% -----------------------------------------------------------------------
%  The pooled breath-triggered figure (breath_trig_heatmap_pooled_260806.m),
%  extended across GENOTYPES: the new 260806 Sst session plus every recording
%  in D:\Ventral_surface_summary under Vglut2, Vgat, Sst, ChAT and Sert.
%
%  Identical panels and identical rules -- no sorting, no selection, every ROI
%  of every recording that has both a breath trigger and a dFF file.  The only
%  change is the colour code and the row order: ROIs are grouped by GENOTYPE
%  (then recording, then native ROI order), and colour now means genotype.
%
%      row 1   breath trace, mean per genotype        | breath power spectra
%      row 2   all ROI dF/F traces (5% bar) | dF/F heatmap | z heatmap | z power
%
%  Pooling across genotypes has the same rate problem as pooling across FOVs,
%  only worse: breath rate and frame rate both differ between sessions.  So, as
%  in the pooled script, every PETH is built on its own axis and resampled onto
%  a common axis in BREATH CYCLES, then labelled in seconds at the pooled median
%  IBI.  Only the trigger and the cycle structure are aligned; no amplitude is
%  rescaled.
%
%  Output: <outRoot>\breath_trig_heatmap_GENOTYPES_<TRIGGER>.png / .pdf / .mat
%
%  Runqi Zhang / 2026-08-06

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
% Each group = one row block of the figure, in this order.
GROUPS = { ...
    'Vglut2',     'D:\Ventral_surface_summary\Vglut2'     ; ...
    'Vgat',       'D:\Ventral_surface_summary\Vgat'       ; ...
    'Sst',        'D:\Ventral_surface_summary\Sst'        ; ...
    'ChAT',       'D:\Ventral_surface_summary\ChAT'       ; ...
    'Sert',       'D:\Ventral_surface_summary\Sert'       };
% The 260806/260807 Sst sessions were copied into Sst\0806 and Sst\0807 on
% 2026-08-08 (archive_sst_260806_260807), so they now arrive through the Sst root.
% They must NOT also be listed as C: paths -- that counted every one of their ROIs
% twice.  They are hand-curated for activity but are NOT in the polar analysis, so
% none of their ROIs can carry a Rayleigh logZ: they land in 'active'/'inactive'.

outRoot    = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806';

% Rayleigh phase locking, recycled from the polar analysis (per CELL -> its ROIs).
% The measure is RY.logZ -- the polar figure's own radial coordinate, log(Rayleigh Z).
% Heatmap bands are cut at logZ = 3, 2, 1.  For reference those are p_analytic of
% about 2e-9, 6e-4 and 0.07, and the polar figure's dashed circles sit at
% logZ = 1.097 (p .05) and 1.933 (p .001).
rayFile = 'D:\Ventral_surface_summary\polar_coh_vs_rayleigh_260808\polar_coh_vs_rayleigh_data.mat';
% REPOINTED 2026-08-11 from the 260729 run. That one predates Sst 0806/0807,
% Vglut2 0810 and the Vgat 120 um depth bar, so it could only label 389 of the
% ROIs here and none of the new sessions could ever land in a logZ band.
zCuts   = [3 2 1];         % white lines here, highest first
zLock   = 1;               % spectra: RED at or above this logZ

TRIGGER    = 'peak';       % 'peak' | 'onset'
% Windows are HARD-FIXED in seconds (not derived from the breath rate).  Sessions
% here breathe at 0.4-2 Hz, so the +/-1 cycle peaks WILL smear across sessions --
% that is the price of a fixed axis and it is deliberate.
tWin_top   = 10;           % +/- seconds, top breath trace
tWin_panel = 5;            % +/- seconds, row 2 heatmaps
nDrop      = 30;           % calcium TossFrames
prcLim     = [0.5 99.5];   % percentile clamp for the heatmap colour scales
nGrid      = 241;          % samples across a panel window on the common cycle axis
nPhaseBins = 120;          % INTERPOLATION grid for the dF/F phase panels (Pcyc/PZcyc).
                           % This is resampling a continuous trace, not counting, so it
                           % is not limited by the frame rate and can stay fine.
nSpkPhBins = 18;           % HISTOGRAM bins for the spike-phase distribution (SpkPh).
                           % 18 = 20 deg. Set by the sampling, not by taste: on the
                           % fast sessions (0.53 s cycle at 30 fps) one frame already
                           % spans 22.5 deg, so 16 bins is the hard resolution limit
                           % and anything finer is empty resolution.
                           % bins per breath cycle for the PHASE version.
                           % Phase is piecewise: inspiration onset = 0, inspiratory
                           % peak = pi, next onset = 2pi, linear within each half.
                           % Onset comes from breath_insp_start_pc1.mat; a recording
                           % without that file gets time panels but no phase.
nCycShow   = 20;           % heatmap rows per ROI: this many cycles, drawn at RANDOM
                           %   from that ROI's boundary-safe cycles (all of them if it
                           %   has fewer).  Inf = every cycle.
                           %   20 is safe: the worst recording of the 91 still has 24
                           %   boundary-safe cycles (median 76), so no ROI is short.
rngSeed    = 0;            % fixed so the same 10 cycles come back on a re-run

fpass      = [0.05 14];    % Hz.  NOTE the multitaper half-bandwidth below is 0.25 Hz,
                           % so everything under ~0.25 Hz is inside one resolution
                           % element -- it is drawn, but it is not resolved.
W_hz       = 0.25;         % multitaper half-bandwidth (smoothing), Hz
                           % nothing is normalised: breath PC1 and the calcium
                           % traces are both z-scored, so both spectra are z^2/Hz

stackGain    = 1.0;        % spacing of the dF/F traces
scaleBarDFF  = 0.05;       % 5% dF/F calibration bar
minSpikes    = 0;          % an ROI counts as ACTIVE if nnz(spike_train>0) > minSpikes.
                           %   0 = any detected event at all.  5 = archive ACTIVE criterion.
                           % A recording with no ca_spike_data.mat was curated by hand as
                           % having no activity, so those ROIs count as INACTIVE (grey) --
                           % that is a curation result, not missing data.
fallback_fps = 30;
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(rngSeed);              % reproducible cycle subsampling

%% ---- common axes ----
tt = linspace(-tWin_panel, tWin_panel, nGrid);        % SECONDS, common to all sessions
tb = linspace(-tWin_top,   tWin_top,   3*nGrid);
fg = logspace(log10(fpass(1)), log10(fpass(2)), 600).';

nGrp   = size(GROUPS,1);
M_all  = []; Mz_all = []; grpOf = []; Sc_all = [];
spkOf  = [];                      % per ROI: 1 = active (detected events), 0 = inactive
roiLbl = {};                      % per ROI: 'Genotype/MMDD/FOV/roiIndex' (polar-figure key)
C_cyc  = {};  Z_cyc = {};  cycGrp = {};  cycRoi = {};   % EVERY cycle of every ROI (row = one cycle)
P_cyc  = {};  PZ_cyc = {};                             % the same cycles, resampled in PHASE
phGrid   = linspace(0, 2*pi, nPhaseBins+1); phGrid(end) = [];   % dF/F interp grid
phEdges  = linspace(0, 2*pi, nSpkPhBins+1);                     % SPIKE histogram edges
phSpkCtr = phEdges(1:end-1) + diff(phEdges(1:2))/2;
SpkPh   = [];   % [nROI x nPhaseBins] spikes per cycle in each phase bin, in PERCENT
B_grp  = cell(nGrp,1); Sb_grp = cell(nGrp,1);
ibi_all = []; nRec_grp = zeros(nGrp,1); nTrig_grp = zeros(nGrp,1);

for g = 1:nGrp
    gname = GROUPS{g,1}; groot = GROUPS{g,2};
    if ~isfolder(groot), warning('%s: %s not found, skipped.', gname, groot); continue; end

    % every folder holding a trigger file AND a dFF file, at any depth
    hits = dir(fullfile(groot, '**', 'breath_peak_pc1.mat'));
    fps_list = {};
    for h = 1:numel(hits)
        fp = hits(h).folder;
        df = dir(fullfile(fp,'*_ch1_dFF.mat'));
        if isempty(df), continue; end

        D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
        dFF = double(D.dFF);
        switch lower(TRIGGER)
            case 'peak'
                B = load(fullfile(fp,'breath_peak_pc1.mat'));      trig = round(B.insp_onset_idx(:));
            case 'onset'
                ipf = fullfile(fp,'breath_insp_start_pc1.mat');
                if ~isfile(ipf), continue; end
                B = load(ipf);                                     trig = round(B.insp_start_idx(:));
        end

        fps    = detect_session_fps(fp, fallback_fps);
        breath = B.breath(:); breath = breath(nDrop+1:end);
        T      = min(size(dFF,1), numel(breath));
        if T < 100, continue; end
        dFF = dFF(1:T,:); breath = breath(1:T);
        trig = trig - nDrop;
        % Vglut2/1124: rising-edge 2P trigger, breath leads calcium by one frame
        if contains(fp, fullfile('Vglut2','1124')), trig = trig + 1; end

        nROI = size(dFF,2);
        Z    = (dFF - mean(dFF,1)) ./ std(dFF,0,1);

        % ACTIVE = has detected calcium events.  No ca_spike_data.mat means the
        % recording was curated as having no activity, so those ROIs are inactive (0).
        spk = zeros(nROI,1);  spkTrains = [];
        caf = fullfile(fp,'ca_spike_data.mat');
        if isfile(caf)
            try
                CA  = load(caf,'roi_spikes');
                nEv = arrayfun(@(r) nnz(r.spike_train>0), CA.roi_spikes); nEv = nEv(:);
                if numel(nEv) >= nROI
                    spk = double(nEv(1:nROI) > minSpikes);
                    spkTrains = CA.roi_spikes;
                end
            catch
            end
        end
        spkOf = [spkOf; spk]; %#ok<AGROW>

        % polar-figure key WITHOUT the group prefix: MMDD/FOV/roi.
        % The prefix has to go: the polar analysis calls IO-site recordings 'IO/...'
        % while they live under ChAT\ and Vglut2\ here, so a Genotype-prefixed key
        % silently missed every IO cell.
        pp   = strsplit(fp, filesep);
        iGen = find(strcmpi(pp, GROUPS{g,1}), 1, 'last');
        if ~isempty(iGen) && numel(pp) >= iGen+2
            stem = sprintf('%s/%s', pp{iGen+1}, pp{end});
        else
            stem = sprintf('?/%s', pp{end});
        end
        roiLbl = [roiLbl; arrayfun(@(r) sprintf('%s/%d', stem, r), (1:nROI).', 'uni',0)]; %#ok<AGROW>

        ibi   = median(diff(sort(trig)))/fps;
        if ~(ibi > 0), continue; end
        win   = round(tWin_panel*fps);          % fixed seconds -> frames, per fps
        win_b = round(tWin_top  *fps);
        tg    = trig(trig-win   >= 1 & trig+win   <= T);
        tg_b  = trig(trig-win_b >= 1 & trig+win_b <= T);
        if isempty(tg) || isempty(tg_b), continue; end
        idx   = tg(:).'   + (-win:win).';
        idx_b = tg_b(:).' + (-win_b:win_b).';
        tl_u  = (-win:win)    /fps;             % local axis in SECONDS
        tl_ub = (-win_b:win_b)/fps;

        % ---- breath PHASE for the same cycles (onset 0, peak pi, next onset 2pi) ----
        phTarget = []; phRow = [];
        ipf = fullfile(fp,'breath_insp_start_pc1.mat');
        if isfile(ipf)
            IP  = load(ipf, 'insp_start_idx');
            ons = round(IP.insp_start_idx(:)) - nDrop;
            if contains(fp, fullfile('Vglut2','1124')), ons = ons + 1; end
            ons = ons(ons >= 1 & ons <= T);
            phi = piecewise_phase(round(trig), ons, T);       % cumulative, NaN outside
            kOf = discretize(tg, ons);                        % which cycle each trigger is in
            good = ~isnan(kOf) & ~isnan(phi(min(max(tg,1),T)));
            phRow = find(good);
            phTarget = (kOf(good).'-1)*2*pi + phGrid.';       % [nPhaseBins x nGoodTrig]
            okPhi = ~isnan(phi);
        end

        % ---- spikes per cycle in each phase bin, in PERCENT, per ROI ----
        % Denominator is the number of COMPLETE breath cycles the recording has, so
        % the unit is "spikes per cycle in this bin x 100" and is comparable between
        % ROIs and recordings regardless of how long the recording ran.
        Sph = nan(nROI, nSpkPhBins);
        if ~isempty(phTarget) && ~isempty(spkTrains)
            nCycRec = numel(unique(floor(phi(okPhi)/(2*pi))));
            if nCycRec > 0
                for r = 1:nROI
                    st = double(spkTrains(r).spike_train(:));
                    if numel(st) < T, st(end+1:T) = 0; end
                    sIdx = find(st(1:T) > 0);
                    sPh  = mod(phi(sIdx), 2*pi);
                    sPh  = sPh(~isnan(sPh));
                    Sph(r,:) = 100 * histcounts(sPh, phEdges) / nCycRec;
                end
            end
        end
        SpkPh = [SpkPh; Sph]; %#ok<AGROW>

        Mk = nan(nROI, nGrid); Mzk = nan(nROI, nGrid);
        for r = 1:nROI
            tr = dFF(:,r); zr = Z(:,r);
            Ei = interp1(tl_u, tr(idx), tt).';        % [nTrig x nGrid] cycles, dF/F
            Zi = interp1(tl_u, zr(idx), tt).';        %                 cycles, z
            % the SAME cycles in phase (NaN row where the cycle has no onset pair)
            Pi = nan(size(Ei,1), nPhaseBins);  PZi = Pi;
            if ~isempty(phTarget)
                Pi(phRow,:)  = reshape(interp1(phi(okPhi), tr(okPhi), phTarget(:)), nPhaseBins, []).';
                PZi(phRow,:) = reshape(interp1(phi(okPhi), zr(okPhi), phTarget(:)), nPhaseBins, []).';
            end
            % keep nCycShow cycles at random (all of them if it has fewer)
            nC = size(Ei,1);
            if nC > nCycShow
                sel = sort(randperm(nC, nCycShow));
                Ei = Ei(sel,:); Zi = Zi(sel,:); Pi = Pi(sel,:); PZi = PZi(sel,:);
            end
            C_cyc{end+1}  = single(Ei);              %#ok<SAGROW>
            Z_cyc{end+1}  = single(Zi);              %#ok<SAGROW>
            P_cyc{end+1}  = single(Pi);              %#ok<SAGROW>
            PZ_cyc{end+1} = single(PZi);             %#ok<SAGROW>
            cycGrp{end+1} = repmat(g, size(Ei,1), 1);%#ok<SAGROW>
            cycRoi{end+1} = repmat(numel(spkOf)-nROI+r, size(Ei,1), 1); %#ok<SAGROW>
            Mk(r,:)  = interp1(tl_u, mean(tr(idx),2), tt, 'linear', NaN);
            Mzk(r,:) = interp1(tl_u, mean(zr(idx),2), tt, 'linear', NaN);
        end
        M_all  = [M_all;  Mk];   %#ok<AGROW>
        Mz_all = [Mz_all; Mzk];  %#ok<AGROW>
        grpOf  = [grpOf;  repmat(g, nROI, 1)]; %#ok<AGROW>

        B_grp{g}(end+1,:) = interp1(tl_ub, mean(breath(idx_b),2), tb, 'linear', NaN);

        TW = max(2, round(W_hz * T/fps)); Kt = 2*TW - 1;
        params = struct('Fs',fps, 'tapers',[TW Kt], 'fpass',fpass, 'pad',0, 'trialave',0, 'err',0);
        [Sb, fb] = mtspectrumc(detrend(breath), params);
        [Sc, fc] = mtspectrumc(detrend(Z),      params);
        Sb_grp{g}(:,end+1) = interp1(fb, Sb, fg, 'linear', NaN);
        Sc_all = [Sc_all, interp1(fc, Sc, fg, 'linear', NaN)]; %#ok<AGROW>

        ibi_all(end+1) = ibi;           %#ok<AGROW>
        nRec_grp(g)  = nRec_grp(g) + 1;
        nTrig_grp(g) = nTrig_grp(g) + numel(tg);
        fps_list{end+1} = sprintf('%.0f', fps); %#ok<AGROW>
    end
    fprintf('  %-12s %2d recordings | %4d ROIs | %5d triggers | fps {%s}\n', ...
            gname, nRec_grp(g), nnz(grpOf==g), nTrig_grp(g), strjoin(unique(fps_list), ','));
end

ibi_pool = median(ibi_all);
nROIall = size(M_all,1);
Ccyc   = vertcat(C_cyc{:});     clear C_cyc
Zcyc   = vertcat(Z_cyc{:});     clear Z_cyc
Pcyc   = vertcat(P_cyc{:});     clear P_cyc
PZcyc  = vertcat(PZ_cyc{:});    clear PZ_cyc
cycGrp = vertcat(cycGrp{:});
cycRoi = vertcat(cycRoi{:});
fprintf('  phase: %d of %d cycles resampled (%.0f%%)\n', ...
        nnz(~all(isnan(Pcyc),2)), size(Pcyc,1), 100*nnz(~all(isnan(Pcyc),2))/size(Pcyc,1));
fprintf('POOLED across genotypes: %d ROIs, %d recordings, %d triggers, median IBI %.3f s (%.2f Hz)\n', ...
        nROIall, sum(nRec_grp), sum(nTrig_grp), ibi_pool, 1/ibi_pool);
fprintf('  single-cycle matrix: %d rows (ROI x cycle) x %d columns\n', size(Ccyc,1), size(Ccyc,2));

%% ---- Rayleigh phase locking, per ROI ----
% RY is per CELL; obsOfCell maps each ROI observation to its cell, so an ROI inherits
% the phase-locking verdict of the cell it belongs to.
rayZ = nan(nROIall,1);                       % log(Rayleigh Z), per ROI
if isfile(rayFile)
    R = load(rayFile, 'RY','labels','obsOfCell');
    zObs = nan(numel(R.labels),1);
    ok   = R.obsOfCell >= 1 & R.obsOfCell <= numel(R.RY.logZ);
    zObs(ok) = R.RY.logZ(R.obsOfCell(ok));
    key = cellfun(@(s) strjoin(subsref(regexp(s,'/','split'), substruct('()',{2:4})),'/'), ...
                  R.labels(:), 'uni', 0);     % drop the group prefix -> MMDD/FOV/roi
    [uk, ia] = unique(key);
    if numel(uk) < numel(key)
        warning('%d duplicate polar keys after dropping the group prefix.', numel(key)-numel(uk));
    end
    M = containers.Map(uk, num2cell(zObs(ia)));
    for i = 1:nROIall
        if isKey(M, roiLbl{i}), rayZ(i) = M(roiLbl{i}); end
    end
    fprintf('  Rayleigh: %d of %d ROIs matched (of %d polar observations); logZ>=3 %d, >=2 %d, >=1 %d\n', ...
            nnz(~isnan(rayZ)), nROIall, numel(R.labels), ...
            nnz(rayZ>=3), nnz(rayZ>=2), nnz(rayZ>=1));
else
    warning('Rayleigh file not found: %s -- no ROI will be red.', rayFile);
end

%% ---- ROI class: 3 = phase-locked, 2 = active, 1 = inactive ----
cls = ones(nROIall,1);
cls(spkOf == 1)     = 2;
cls(rayZ >= zLock)  = 3;
isLock = cls == 3;  isAct = cls == 2;  isOff = cls == 1;

% ---- heatmap row order ----
% Top to bottom: ACTIVE ROIs sorted by Rayleigh p (smallest first), then active
% ROIs that were never tested, then INACTIVE ROIs in random order.  White lines
% are drawn where p crosses 0.001, 0.005 and 0.05.  Cycles of one ROI always stay
% together and in time order -- only whole ROIs are moved.
% YDir is 'normal' so the LAST row is the top: build the order bottom-up.
% White lines at logZ = 3, 2, 1 and at the active|inactive boundary.
nTier = numel(zCuts) + 2;
pTier = nTier*ones(nROIall,1);             % last tier = inactive
pTier(cls >= 2) = nTier - 1;               % active: logZ < min(zCuts) or never tested
for k = numel(zCuts):-1:1
    pTier(cls >= 2 & rayZ >= zCuts(k)) = k;
end
tierName = [arrayfun(@(z) sprintf('logZ>=%g',z), zCuts, 'uni',0), {'active','inactive'}];

zKey = rayZ;  zKey(isnan(zKey)) = -Inf;    % untested sort below every tested ROI
rp   = randperm(nROIall).';                % random order inside the inactive block
[~, roiOrder] = sortrows([-pTier, zKey, (pTier==nTier).*rp]);   % bottom-up
% -> roiOrder(1) is the BOTTOM row, roiOrder(end) the TOP (smallest p)

rank = zeros(nROIall,1); rank(roiOrder) = 1:nROIall;
[~, ordC] = sortrows([rank(cycRoi), (1:numel(cycRoi)).']);   % keep cycles together
Ccyc = Ccyc(ordC,:);  Zcyc = Zcyc(ordC,:);
Pcyc = Pcyc(ordC,:);  PZcyc = PZcyc(ordC,:);
cycTier = pTier(cycRoi(ordC));

cbnd = find(diff(cycTier) ~= 0) + 0.5;
cctr = arrayfun(@(t) mean(find(cycTier==t)), 1:nTier);
clbl = arrayfun(@(t) sprintf('%s (%d)', tierName{t}, nnz(cycTier==t)), 1:nTier, 'uni',0);
keep = ~isnan(cctr);  cctr = cctr(keep);  clbl = clbl(keep);
[cctr, iSort] = sort(cctr);  clbl = clbl(iSort);     % YTick must increase
fprintf('  heatmap tiers (ROIs):');
for t = 1:nTier, fprintf(' %s %d |', tierName{t}, nnz(pTier==t)); end
fprintf('\n');

%% ---- figure ----
grpCol = lines(nGrp);
bnd    = find(diff(grpOf) ~= 0) + 0.5;                       % genotype boundaries
ctr    = arrayfun(@(g) mean(find(grpOf==g)), 1:nGrp);        % row centre per group
lbl    = arrayfun(@(g) sprintf('%s (%d)', GROUPS{g,1}, nnz(grpOf==g)), 1:nGrp, 'uni',0);

fh = figure('Color','w','Position',[30 30 1500 990]);
tl = tiledlayout(fh, 5, 3, 'TileSpacing','compact','Padding','compact');

% ---- row 1 ----
axb = nexttile(tl, 1, [1 2]); hold(axb,'on');
for g = 1:nGrp
    if isempty(B_grp{g}), continue; end
    plot(axb, tb, mean(B_grp{g},1,'omitnan'), '-', 'Color', grpCol(g,:), 'LineWidth', 1.4);
end
xline(axb, 0, 'r-', 'LineWidth', 1);
for e = [-1 1]*tWin_panel, xline(axb, e, ':', 'Color',[0.5 0.5 0.5]); end
xlim(axb, [tb(1) tb(end)]); box(axb,'off');
set(axb, 'TickDir','out', 'XTick', -tWin_top:2:tWin_top);
xlabel(axb, sprintf('time from inspiratory %s (s)', TRIGGER));
ylabel(axb, 'breath (PC1, z)');
title(axb, sprintf('%s-triggered breath, +/-%g s, mean per genotype', TRIGGER, tWin_top));
legend(axb, GROUPS(:,1), 'Location','eastoutside', 'Box','off', 'Interpreter','none', 'FontSize',8);

axs1 = nexttile(tl, 3, [1 1]); hold(axs1,'on');
for g = 1:nGrp
    if isempty(Sb_grp{g}), continue; end
    loglog(axs1, fg, Sb_grp{g}, '-', 'Color',[grpCol(g,:) 0.5], 'LineWidth',0.7);
end
finish_spec(axs1, fpass, 1/ibi_pool, 'power (z^2/Hz)', 'breath (PC1) power, per recording');

% ---- row 2: SINGLE CYCLES only, no trial averaging, sig on top ----
if isfinite(nCycShow)
    cycTtl = sprintf('%d ROIs x %d random cycles = %d', nROIall, nCycShow, size(Ccyc,1));
else
    cycTtl = sprintf('%d ROIs x all cycles = %d', nROIall, size(Ccyc,1));
end
draw_heatmap(nexttile(tl, 4, [4 1]), tt, Ccyc, prcLim, cbnd, cctr, clbl, ...
             sprintf('single-cycle dF/F, %s', cycTtl), 'dF/F', TRIGGER, gray(256));
draw_heatmap(nexttile(tl, 5, [4 1]), tt, Zcyc, prcLim, cbnd, cctr, clbl, ...
             sprintf('single-cycle z-scored dF/F, %s', cycTtl), 'z (SD)', TRIGGER, parula(256));

% per-ROI spectra, coloured by whether that ROI has detected calcium events
% RED = Rayleigh logZ >= zLock | GREEN = active | GREY = inactive
% spectra: ACTIVE = green, INACTIVE = red.  Two classes only -- phase locking is
% not colour-coded here, so this panel says nothing about logZ.
axs2 = nexttile(tl, 6, [4 1]); hold(axs2,'on');
anyAct = isAct | isLock;               % locked ROIs are active by definition
specAlpha = 0.6;
for i = find(~anyAct).'
    loglog(axs2, fg, Sc_all(:,i), '-', 'Color',[0.85 0.1 0.1 specAlpha], 'LineWidth',0.4);
end
for i = find(anyAct).'
    loglog(axs2, fg, Sc_all(:,i), '-', 'Color',[0 0.65 0.2 specAlpha], 'LineWidth',0.5);
end
finish_spec(axs2, fpass, 1/ibi_pool, 'power (SD^2/Hz)', ...
            sprintf('z dF/F power: %d active (green), %d inactive (red)', ...
                    nnz(anyAct), nnz(~anyAct)));
fprintf('  spike labels: %d ROIs with events, %d without, %d with no ca_spike_data.mat\n', ...
        nnz(spkOf==1), nnz(spkOf==0), nnz(isnan(spkOf)));

title(tl, sprintf('breath-triggered dF/F across genotypes   |   %d recordings, %d ROIs, no sorting', ...
      sum(nRec_grp), nROIall), 'Interpreter','none','FontWeight','bold');

%% ---- save ----
if doSave
    if ~isfolder(outRoot), mkdir(outRoot); end
    base = fullfile(outRoot, sprintf('breath_trig_heatmap_GENOTYPES_%s', lower(TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], '-v7.3', 'Ccyc','Zcyc','Pcyc','PZcyc','phGrid','phEdges','phSpkCtr','nPhaseBins','nSpkPhBins','SpkPh', ...
                        'cycRoi','cycTier','ordC','tierName','cls', ...
                        'M_all','Mz_all','grpOf','spkOf','minSpikes','rayZ','zCuts','zLock','roiLbl','cycGrp','GROUPS','B_grp','Sb_grp','Sc_all','fg', ...
                        'tt','tb','ibi_all','ibi_pool','nRec_grp','nTrig_grp', ...
                        'TRIGGER','tWin_panel','tWin_top','nDrop','prcLim','fpass','W_hz');
    fprintf('saved %s.{png,pdf,mat}\n', base);
end

%% ---- local ----
function phi = piecewise_phase(peakIdx, onsetIdx, T)
%PIECEWISE_PHASE  Cumulative breath phase: onset = 0, peak = pi, next onset = 2pi.
% Linear inside each half separately, NOT a uniform ramp between onsets -- the
% inspiratory half is far shorter than the expiratory one, so a uniform ramp would
% put the peak at the wrong phase and smear every phase-locked cell.
% Cumulative (cycle k spans (k-1)*2pi .. k*2pi) so one interp1 can serve all cycles.
phi = nan(T,1);
peakIdx  = sort(peakIdx(:));
onsetIdx = sort(onsetIdx(:));
for k = 1:numel(onsetIdx)-1
    o1 = onsetIdx(k); o2 = onsetIdx(k+1);
    if o2 <= o1 + 1, continue; end
    p = peakIdx(peakIdx > o1 & peakIdx < o2);
    if isempty(p), continue; end            % no peak in this cycle -> leave NaN
    p = p(1);
    if p <= o1 || p >= o2, continue; end
    base = (k-1)*2*pi;
    phi(o1:p) = base + pi*(0:(p-o1)).'/(p-o1);
    phi(p:o2) = base + pi + pi*(0:(o2-p)).'/(o2-p);
end
end

function draw_heatmap(ax, tt, A, prcLim, bnd, ctr, lbl, ttl, cbl, TRIGGER, cmap)
c = prctile(A(:), prcLim);
if ~(c(2) > c(1)), c = [min(A(:)) max(A(:))+eps]; end
imagesc(ax, tt, 1:size(A,1), A);
set(ax, 'CLim', c, 'YDir','normal', 'TickDir','out', 'Layer','top');
colormap(ax, cmap); hold(ax,'on');
xline(ax, 0, 'r-', 'LineWidth', 1);
for b = bnd(:).', plot(ax, [tt(1) tt(end)], [b b], 'w-', 'LineWidth', 1.2); end
set(ax, 'YTick', ctr, 'YTickLabel', lbl, 'TickLabelInterpreter','none', 'FontSize',8);
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
title(ax, ttl);
cb = colorbar(ax); cb.Label.String = cbl;
end

function draw_stack(ax, tt, A, colPerRow, gain, barDFF, TRIGGER, bnd, lbl, ctr)
nR   = size(A,1);
step = gain * median(range(A, 2), 'omitnan');
if ~(step > 0), step = 1; end
hold(ax,'on');
for k = 1:nR
    plot(ax, tt, A(k,:) - mean(A(k,:),'omitnan') + (k-1)*step, '-', ...
         'Color', colPerRow(k,:), 'LineWidth', 0.4);
end
xline(ax, 0, 'k-', 'LineWidth', 1);
xlim(ax, [tt(1) tt(end)]); ylim(ax, [-step, (nR-1)*step + step]);
for b = bnd(:).', yline(ax, (b-1)*step, '-', 'Color',[0.75 0.75 0.75]); end
set(ax, 'YTick', (ctr-1)*step, 'YTickLabel', lbl, 'TickDir','out', 'Box','off', ...
        'FontSize',8, 'TickLabelInterpreter','none');
xb = tt(end) - 0.02*(tt(end)-tt(1)); yb = -0.6*step;
plot(ax, [xb xb], yb + [0 barDFF], 'k-', 'LineWidth', 2.5);
text(ax, xb - 0.02*(tt(end)-tt(1)), yb + barDFF/2, sprintf('%g%% dF/F', 100*barDFF), ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontSize', 9);
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
title(ax, sprintf('mean dF/F, all %d ROIs (common scale)', nR));
end

function finish_spec(ax, fpass, f_breath, ylb, ttl)
hold(ax,'on');
xline(ax, f_breath, 'r--', 'LineWidth', 1);
set(ax, 'XScale','log', 'YScale','log', 'TickDir','out', 'Box','off');
xlim(ax, fpass);
xt = [0.05 0.1 0.2 0.5 1 2 3 5 7 10 14]; xt = xt(xt>=fpass(1) & xt<=fpass(2));
set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt));
xlabel(ax, 'frequency (Hz)'); ylabel(ax, ylb); title(ax, ttl);
end
