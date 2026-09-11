% breath_amp_vs_dff_260825.m
% -----------------------------------------------------------------------
%  GRADED per-cycle relationship between BREATH AMPLITUDE and Ca RESPONSE,
%  one XY scatter per cell.
%
%      X_i = breath amplitude of cycle i   (robust-z peak-to-trough of PC1)
%      Y_i = Ca response in cycle i        (baseline-subtracted dF/F peak in a
%                                           FIXED window after inspiration onset)
%
%  THE QUESTION: does a bigger breath go with a bigger response IN THIS CELL?
%  The unit is the breath CYCLE and the comparison is WITHIN a recording, which
%  is what makes it clean: inside one recording the laser power, depth, zoom and
%  the breath video's arbitrary PC1 scale are all constant, and a rank
%  correlation is invariant to any monotonic per-recording rescaling anyway. So
%  none of those confound the result. Amplitudes are NEVER pooled across
%  recordings -- only rank-transformed values are.
%
%  ---------------------------------------------------------------------
%  THREE THINGS THAT WOULD FAKE THE EFFECT, AND WHAT IS DONE ABOUT EACH
%
%  1. CYCLE DURATION.  If the response window were "the cycle", a longer cycle
%     would give the dF/F peak more time to accumulate, so Y would rise with
%     duration -- and duration correlates with amplitude (more so after
%     vagotomy). That alone produces a positive X-Y slope with no neural content.
%     => the window is FIXED in absolute seconds (respWin_sec), identical for
%     every cycle, matching the absolute-time convention of analysis_260727.
%
%  2. GCaMP CARRY-OVER.  GCaMP8s decays with tau ~0.45 s against ~1 s cycles, so
%     a raw peak carries a large tail from the PREVIOUS cycle -- that is
%     carry-over, not response.
%     => Y subtracts a short pre-onset baseline (basePre_sec), and a second
%        response measure, the peak of d(dF/F)/dt, is computed alongside it.
%        A decay tail contributes NEGATIVE slope, so the derivative version is
%        nearly immune to carry-over. If the two disagree, believe neither
%        until you know why.
%
%  3. THE BREATH-LOCKED OPTICAL ARTIFACT.  This is the one that matters. A
%     bigger breath displaces the animal more, so it perturbs the optics more,
%     so apparent dF/F rises -- in EVERY ROI at once, within a single recording,
%     with no imaging-parameter confound anywhere. In deep Vgat the top SVD mode
%     tracks chest breathing at r = 0.96 and survives motion correction.
%     => a COMMON-MODE control: the FOV population mean response is regressed
%        out (rho_partial_pop). An optical modulation is shared by every ROI in
%        the field, so partialling the population mean removes it; a cell with
%        genuine private tuning survives.
%        THIS CONTROL IS CONSERVATIVE, NOT DECISIVE. Real co-modulated neural
%        activity is also common-mode and is removed too, so rho_partial_pop is
%        a LOWER bound on the neural effect. The decisive test is a background
%        annulus around each soma, which needs the registered movie
%        (*_ch1_preproc_MC_MC.tif) and so cannot run off the archive -- see
%        "NOT DONE HERE" below.
%
%  ALSO REPORTED: rho_partial_ibi. Breath amplitude and rate co-vary (big
%  breaths come in fast bouts), so an apparent amplitude effect can be a rate
%  effect wearing a different hat. Points on the scatter are coloured by IBI so
%  this is visible rather than merely tested.
%
%  ---------------------------------------------------------------------
%  ALL CYCLES ARE USED, INCLUDING FAILED BREATHS.  The registry pipeline gates
%  cycles at ampFrac = 0.20 x median amplitude and discards the rest. That gate
%  is NOT applied here: the small-amplitude cycles are precisely the low end of
%  the X axis, and dropping them would truncate the range this analysis exists
%  to measure. Breath landmarks are therefore recomputed from the pc1 files
%  rather than taken from REC.trig, which is already gated.
%
%  SIGNIFICANCE: circular-shift null. The dF/F trace of each observation is
%  circularly shifted against the breath by at least shiftMinCyc cycles and rho
%  recomputed, nShuffle times. This preserves the Ca autocorrelation and the
%  breath statistics, which a naive Spearman p-value does not -- with ~400
%  serially correlated cycles that p would be wildly overconfident.
%
%  POOLING ACROSS A CELL'S RECORDINGS: X and Y are converted to percentile ranks
%  WITHIN each recording and then concatenated, so every recording contributes on
%  a common [0 1] scale and no absolute amplitude ever crosses a recording
%  boundary. Per-observation rho is stored too, for checking.
%
%  CELL IDENTITY comes from the registry (event_latency_data.mat), so a cell id
%  here is the same pooled id as in per-cell-summary_active_260812 and in
%  cell_id_map_0824.csv. That also means this runs on the ACTIVE set only
%  (registry gate: pooled > 5 events).
%
%  NOT DONE HERE: the background-annulus artifact control. It requires the
%  registered movie, which is not in the archive (only the acquisition drives
%  hold *_ch1_preproc_MC_MC.tif). rho_partial_pop is the archive-only stand-in.
%
%  OUTPUT  <outDir>\
%     breath_amp_vs_dff_cells.csv     one row per cell, every statistic
%     breath_amp_vs_dff_cycles.mat    per-cycle X/Y/IBI, for re-plotting
%     population_summary.png/.pdf     rho distributions + the common-mode check
%     cells_page%02d.png              paged per-cell scatter gallery
%
%  Dependencies: detect_session_fps.m, Statistics TB (tiedrank, prctile)
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(fullfile(repoRoot, 'analysis_260806'));

%% ===================== USER-EDITABLE =====================
rootPath   = 'D:\Ventral_surface_summary';
regFile    = fullfile(rootPath, 'event_latency_260811', 'event_latency_data.mat');
outDir     = fullfile(rootPath, 'breath_amp_vs_dff_260825');

% Restrict to a subset of cells. [] = every cell in the registry.
% e.g. 282:297 for the 260824 vagotomised session alone.
onlyCells  = 282:297;
% Restrict by group name, {} = all. e.g. {'Vglut2'}
onlyGroups = {};

% ---- the two axes -----------------------------------------------------
respWin_sec = 0.50;    % FIXED response window after inspiration onset (s).
                       % Must stay inside one cycle: the breath band across this
                       % archive is ~0.65-1.2 Hz, so a cycle is ~0.8-1.5 s.
basePre_sec = 0.10;    % pre-onset baseline averaged over this (s)

nDrop        = 30;     % breath frames tossed up front, to align with Ca
fallback_fps = 30;

minCycles    = 20;     % need this many usable cycles in an observation to score it
minCyclesCell= 40;     % and this many pooled across the cell

% ---- null -------------------------------------------------------------
nShuffle     = 1000;
shiftMinCyc  = 3;      % reject shifts closer than this many cycles to no-shift
alphaSig     = 0.05;

% ---- plotting ---------------------------------------------------------
nBinX        = 8;      % bins for the running median overlaid on each scatter
rowsPerPage  = 5;
colsPerRow   = 6;
maxPanels    = 300;    % cap on gallery panels (0 = no gallery)
dotAlpha     = 0.30;
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
% NVIDIA WGL throws inside exportgraphics ~1% of the time here and writes a 9 KB
% blank. See reference_matlab_opengl_export_stubs.
try, opengl('software'); catch, end
if doSave && ~isfolder(outDir), mkdir(outDir); end
rng(260825);

%% ===================== LOAD THE REGISTRY =============================
assert(isfile(regFile), 'registry not found: %s', regFile);
D = load(regFile, 'CELL','OBS','REC','groups');
CELL = D.CELL; OBS = D.OBS; REC = D.REC; groups = D.groups;
fprintf('\n=========== breath_amp_vs_dff_260825 ===========\n');
fprintf('registry : %d cells, %d observations, %d recordings\n', ...
        numel(CELL), numel(OBS), numel(REC));

% hand-curated merges, same as every other consumer
obsOfCell = pooled_obs_260814(CELL, OBS);
cellIdx = find(~cellfun(@isempty, obsOfCell));
if ~isempty(onlyCells), cellIdx = cellIdx(ismember(cellIdx, onlyCells(:)')); end
if ~isempty(onlyGroups)
    keep = false(size(cellIdx));
    for k = 1:numel(cellIdx)
        keep(k) = any(strcmp(groups{CELL(cellIdx(k)).gi}, onlyGroups));
    end
    cellIdx = cellIdx(keep);
end
fprintf('scoring  : %d cells\n', numel(cellIdx));
fprintf('window   : peak of dF/F over [onset, onset+%.2f s] minus mean over %.0f ms before\n', ...
        respWin_sec, basePre_sec*1000);

%% ============ PER-RECORDING CYCLE TABLE + dF/F ACCESSORS =============
% Everything a recording contributes is computed ONCE and cached: the cycle
% landmarks, the per-cycle breath amplitude and IBI, and the circularly-wrapped
% running max / running mean of every ROI's dF/F. With those, evaluating a cycle
% response at ANY circular shift is a single indexing operation, which is what
% makes 1000 shuffles per cell affordable.
recCache = cell(numel(REC),1);
needRec  = unique(arrayfun(@(o) OBS(o).rec, [obsOfCell{cellIdx}]));
fprintf('\npreparing %d recordings ...\n', numel(needRec));
tPrep = tic; nSkipRec = 0;
for r = needRec(:)'
    C = prep_recording(REC(r).folder, nDrop, fallback_fps, respWin_sec, basePre_sec);
    if isempty(C)
        nSkipRec = nSkipRec + 1;
    end
    recCache{r} = C;
end
fprintf('  done in %.1f min (%d unusable)\n', toc(tPrep)/60, nSkipRec);

%% ========================= SCORE EACH CELL ===========================
R = struct('cell',{},'group',{},'label',{},'nObs',{},'nCyc',{}, ...
           'rho',{},'p',{},'rho_deriv',{},'rho_pIBI',{},'rho_pPop',{}, ...
           'rhoObs',{},'slopeBin',{},'X',{},'Y',{},'IBI',{});
fprintf('\nscoring %d cells (%d shuffles each) ...\n', numel(cellIdx), nShuffle);
tScore = tic;
for k = 1:numel(cellIdx)
    c   = cellIdx(k);
    obs = obsOfCell{c};

    Xp = []; Yp = []; Yd = []; Ip = []; Pp = [];    % pooled percentile ranks
    Xr = []; Yr = []; Ir = [];                      % pooled RAW, for the scatter
    rhoObs = []; parts = struct('rec',{},'roi',{},'cyc',{});
    for o = obs(:)'
        r = OBS(o).rec;  Cc = recCache{r};
        if isempty(Cc), continue; end
        roi = roi_from_label(OBS(o).label);
        if isnan(roi) || roi < 1 || roi > Cc.nROI, continue; end

        [y, yd] = cycle_response(Cc, roi, 0);
        ok = ~isnan(y) & ~isnan(Cc.amp(:)) & ~isnan(Cc.ibi(:));
        if nnz(ok) < minCycles, continue; end

        x  = Cc.amp(ok);  yy = y(ok);  ydd = yd(ok);  ib = Cc.ibi(ok);
        pop = Cc.popY(ok);                      % FOV common mode, this recording

        rhoObs(end+1,1) = spear(x, yy); %#ok<SAGROW>

        % within-recording percentile ranks, so recordings pool on a common
        % scale without any absolute amplitude crossing a recording boundary
        Xp = [Xp; pct(x)];    Yp = [Yp; pct(yy)];   Yd = [Yd; pct(ydd)]; %#ok<AGROW>
        Ip = [Ip; pct(ib)];   Pp = [Pp; pct(pop)]; %#ok<AGROW>
        Xr = [Xr; x];         Yr = [Yr; yy];        Ir = [Ir; ib]; %#ok<AGROW>
        parts(end+1) = struct('rec',r,'roi',roi,'cyc',find(ok)); %#ok<SAGROW>
    end
    if numel(Xp) < minCyclesCell || isempty(parts), continue; end

    rho   = spear(Xp, Yp);
    rhoD  = spear(Xp, Yd);
    rhoI  = partial_spear(Xp, Yp, Ip);
    rhoP  = partial_spear(Xp, Yp, Pp);

    % ---- circular-shift null -------------------------------------------
    % Each observation is shifted INDEPENDENTLY and the cell re-pooled, so the
    % null preserves each recording's own breath structure and Ca bursting.
    rnull = nan(nShuffle,1);
    for s = 1:nShuffle
        Xs = []; Ys = [];
        for q = 1:numel(parts)
            Cc = recCache{parts(q).rec};
            sh = draw_shift(Cc.T, shiftMinCyc*Cc.medIBIfr);
            ys = cycle_response(Cc, parts(q).roi, sh);
            ys = ys(parts(q).cyc);
            xs = Cc.amp(parts(q).cyc);
            g  = ~isnan(ys);
            if nnz(g) < 3, continue; end
            Xs = [Xs; pct(xs(g))]; Ys = [Ys; pct(ys(g))]; %#ok<AGROW>
        end
        if numel(Xs) > 3, rnull(s) = spear(Xs, Ys); end
    end
    nv = rnull(~isnan(rnull));
    if isempty(nv), pShift = NaN;
    else,           pShift = (1 + nnz(abs(nv) >= abs(rho))) / (1 + numel(nv));
    end

    R(end+1) = struct('cell',c, 'group',groups{CELL(c).gi}, 'label',CELL(c).label, ...
        'nObs',numel(parts), 'nCyc',numel(Xp), 'rho',rho, 'p',pShift, ...
        'rho_deriv',rhoD, 'rho_pIBI',rhoI, 'rho_pPop',rhoP, ...
        'rhoObs',rhoObs.', 'slopeBin',{bin_median(Xp, Yp, nBinX)}, ...
        'X',Xr.', 'Y',Yr.', 'IBI',Ir.'); %#ok<SAGROW>

    if mod(k,25) == 0
        fprintf('  %d/%d  (%.1f min)\n', k, numel(cellIdx), toc(tScore)/60);
    end
end
fprintf('  scored %d cells in %.1f min\n', numel(R), toc(tScore)/60);
assert(~isempty(R), 'no cell had enough usable cycles');

%% ========================= REPORT ====================================
rho  = [R.rho]';   pS = [R.p]';   rhoP = [R.rho_pPop]';  rhoI = [R.rho_pIBI]';
sig  = pS < alphaSig;
fprintf('\n--------------------------------------------------\n');
fprintf('%d cells scored, median rho = %+.3f\n', numel(R), median(rho));
fprintf('significant (circular-shift p<%.2g): %d  (%d positive, %d negative)\n', ...
        alphaSig, nnz(sig), nnz(sig & rho>0), nnz(sig & rho<0));
fprintf('sign test on all rho vs 0: p = %.3g\n', signrank_safe(rho));
fprintf('\nafter removing the FOV common mode (the artifact stand-in):\n');
fprintf('   median rho %+.3f -> %+.3f   (%d of %d sig cells keep |rho|>0.1)\n', ...
        median(rho), median(rhoP), nnz(sig & abs(rhoP)>0.1), nnz(sig));
fprintf('after removing IBI:  median rho %+.3f -> %+.3f\n', median(rho), median(rhoI));
fprintf('derivative response: median rho %+.3f\n', median([R.rho_deriv]'));
for gi = 1:numel(groups)
    m = strcmp({R.group}, groups{gi});
    if ~any(m), continue; end
    fprintf('   %-8s n=%3d  median rho %+.3f (partial-pop %+.3f)  %d sig\n', ...
            groups{gi}, nnz(m), median(rho(m)), median(rhoP(m)), nnz(sig(m)));
end
if median(abs(rhoP)) < 0.5*median(abs(rho))
    fprintf(2,['\nWARNING the effect halves when the FOV common mode is removed.\n' ...
               'That is the signature of a field-wide optical modulation. The\n' ...
               'background-annulus test on the registered movie is needed before\n' ...
               'any of this is called neural.\n']);
end

%% ========================= FIGURES ===================================
figure_population(R, groups, alphaSig, outDir, doSave);
if maxPanels > 0
    figure_gallery(R, rowsPerPage, colsPerRow, maxPanels, dotAlpha, alphaSig, outDir, doSave);
end

%% ========================= SAVE ======================================
if doSave
    fid = fopen(fullfile(outDir,'breath_amp_vs_dff_cells.csv'),'w');
    fprintf(fid,['cell,group,n_obs,n_cycles,rho,p_shift,sig,' ...
                 'rho_deriv,rho_partial_ibi,rho_partial_pop,label\n']);
    for i = 1:numel(R)
        fprintf(fid,'%d,%s,%d,%d,%.4f,%.5f,%d,%.4f,%.4f,%.4f,"%s"\n', ...
            R(i).cell, R(i).group, R(i).nObs, R(i).nCyc, R(i).rho, R(i).p, ...
            R(i).p < alphaSig, R(i).rho_deriv, R(i).rho_pIBI, R(i).rho_pPop, R(i).label);
    end
    fclose(fid);
    save(fullfile(outDir,'breath_amp_vs_dff_cycles.mat'), 'R', ...
         'respWin_sec','basePre_sec','nDrop','nShuffle','shiftMinCyc', ...
         'minCycles','minCyclesCell','alphaSig','regFile','-v7.3');
    fprintf('\nSaved to %s\n', outDir);
end
fprintf('Done.\n');

% =======================================================================
% =========================== LOCAL FUNCTIONS ===========================
% =======================================================================
function C = prep_recording(folderPath, nDrop, fallback_fps, respWin_sec, basePre_sec)
%PREP_RECORDING  Cycle table + circularly-wrapped response accessors for one FOV.
%  Returns [] if the recording cannot be used.
C = [];
try
    df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
    bp = fullfile(folderPath,'breath_peak_pc1.mat');
    ip = fullfile(folderPath,'breath_insp_start_pc1.mat');
    if isempty(df) || ~isfile(bp) || ~isfile(ip), return; end

    fps = detect_session_fps(folderPath, fallback_fps);
    Dd  = load(fullfile(df(1).folder, df(1).name),'dFF');
    dff = double(Dd.dFF);

    BP = load(bp);  IP = load(ip);
    bw = detrend(double(BP.breath(:)));
    bw(1:min(nDrop,numel(bw))) = [];  bw = bw - mean(bw);
    peak_idx = round(BP.insp_onset_idx(:)) - nDrop;
    foot_idx = round(IP.insp_start_idx(:)) - nDrop;

    T = min(size(dff,1), numel(bw));
    if T < 100, return; end
    dff = dff(1:T,:);  bw = bw(1:T);
    peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
    foot_idx = sort(foot_idx(foot_idx>=1 & foot_idx<=T));
    if numel(foot_idx) < 5, return; end

    % ---- one row per cycle: foot_i -> peak_i -> foot_{i+1} --------------
    % Amplitude is peak MINUS foot on a robust-z breath trace, never the raw
    % peak: the PC1 baseline drifts, so an absolute peak is not an amplitude.
    % The z-scaling is per recording and monotonic, so it cannot affect a rank
    % correlation -- it only makes the axis readable.
    bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
    nC  = numel(foot_idx) - 1;
    amp = nan(nC,1); ibi = nan(nC,1); on = nan(nC,1);
    for i = 1:nC
        q = peak_idx(peak_idx > foot_idx(i) & peak_idx < foot_idx(i+1));
        if isempty(q), continue; end
        amp(i) = bwz(q(1)) - bwz(foot_idx(i));
        ibi(i) = (foot_idx(i+1) - foot_idx(i)) / fps;
        on(i)  = foot_idx(i);
    end
    good = ~isnan(amp) & ~isnan(on);
    if nnz(good) < 5, return; end

    w = max(2, round(respWin_sec*fps));      % response window, in frames
    b = max(1, round(basePre_sec*fps));      % baseline window, in frames

    % ---- circularly-wrapped running max / running mean -------------------
    % Doubling the trace makes every window circular, so the observed data and
    % every circular shift are evaluated by exactly the same rule -- no edge
    % handling that differs between the two.
    x2   = [dff; dff];
    Mx   = movmax(x2, [0 w-1], 1);                     % Mx(t) = max over [t, t+w-1]
    xpre = [x2(end,:); x2(1:end-1,:)];                 % xpre(t) = x(t-1)
    Bx   = movmean(xpre, [b-1 0], 1);                  % Bx(t) = mean over [t-b, t-1]
    dx2  = [diff(x2,1,1); zeros(1,size(x2,2))] * fps;  % d(dF/F)/dt
    Mdx  = movmax(dx2, [0 w-1], 1);

    C = struct('folder',folderPath,'fps',fps,'T',T,'nROI',size(dff,2), ...
               'amp',amp(good),'ibi',ibi(good),'on',on(good), ...
               'Mx',Mx,'Bx',Bx,'Mdx',Mdx, ...
               'medIBIfr',median(ibi(good))*fps, 'popY',[]);

    % ---- FOV common mode -------------------------------------------------
    % Mean response across EVERY ROI in the field, each z-scored first so one
    % bright ROI cannot dominate. A breath-locked optical modulation is shared
    % by the whole field, so this carries it.
    Yall = nan(numel(C.amp), C.nROI);
    for j = 1:C.nROI
        Yall(:,j) = cycle_response(C, j, 0);
    end
    Z = (Yall - mean(Yall,1,'omitnan')) ./ max(std(Yall,0,1,'omitnan'), eps);
    C.popY = mean(Z, 2, 'omitnan');
catch
    C = [];
end
end

% -----------------------------------------------------------------------
function [y, yd] = cycle_response(C, roi, shiftFrames)
%CYCLE_RESPONSE  Baseline-subtracted dF/F peak per cycle, at a circular shift.
%  Shifting the TRACE by s is identical to reading every cycle window s frames
%  later, so a shuffle costs one indexing operation instead of a re-filter.
idx = mod(C.on - 1 + shiftFrames, C.T) + 1;
y   = C.Mx(idx, roi) - C.Bx(idx, roi);
if nargout > 1, yd = C.Mdx(idx, roi); end
end

% -----------------------------------------------------------------------
function s = draw_shift(T, minFr)
%DRAW_SHIFT  Uniform circular shift excluding +/- minFr around no-shift.
if 2*minFr >= T, s = randi(T) - 1; return; end
s = round(minFr + (T - 2*minFr)*rand);
end

% -----------------------------------------------------------------------
function r = spear(a, b)
%SPEAR  Spearman rho, tie-corrected, without the Statistics toolbox's corr().
a = a(:); b = b(:);
m = ~isnan(a) & ~isnan(b);
if nnz(m) < 4, r = NaN; return; end
ra = tiedrank(a(m)); rb = tiedrank(b(m));
ra = ra - mean(ra);  rb = rb - mean(rb);
d  = sqrt(sum(ra.^2) * sum(rb.^2));
if d <= 0, r = NaN; else, r = sum(ra.*rb) / d; end
end

% -----------------------------------------------------------------------
function r = partial_spear(a, b, z)
%PARTIAL_SPEAR  Spearman partial correlation of a,b controlling for z.
%  Ranks first, then residualise linearly -- the standard rank-partial.
a = a(:); b = b(:); z = z(:);
m = ~isnan(a) & ~isnan(b) & ~isnan(z);
if nnz(m) < 5, r = NaN; return; end
ra = tiedrank(a(m)); rb = tiedrank(b(m)); rz = tiedrank(z(m));
X  = [ones(nnz(m),1), rz];
ea = ra - X*(X\ra);
eb = rb - X*(X\rb);
d  = sqrt(sum(ea.^2) * sum(eb.^2));
if d <= 0, r = NaN; else, r = sum(ea.*eb) / d; end
end

% -----------------------------------------------------------------------
function p = pct(v)
%PCT  Percentile rank in [0 1] within this vector. Monotonic, so it cannot
%  change a rank correlation -- it only puts recordings on a common scale.
v = v(:);
p = nan(size(v));
m = ~isnan(v);
if nnz(m) < 2, p(m) = 0.5; return; end
p(m) = (tiedrank(v(m)) - 0.5) / nnz(m);
end

% -----------------------------------------------------------------------
function B = bin_median(x, y, nBin)
%BIN_MEDIAN  Running median of y in nBin equal-count bins of x.
%  Equal COUNT, not equal width: amplitude distributions are skewed, and equal
%  width leaves the top bins with two points each.
x = x(:); y = y(:);
m = ~isnan(x) & ~isnan(y);  x = x(m); y = y(m);
B = struct('xc',[],'ym',[],'ylo',[],'yhi',[]);
if numel(x) < nBin*3, return; end
e = prctile(x, linspace(0,100,nBin+1));
e(1) = -inf; e(end) = inf;
for i = 1:nBin
    s = x >= e(i) & x < e(i+1);
    if nnz(s) < 3, continue; end
    B.xc(end+1,1)  = median(x(s));
    B.ym(end+1,1)  = median(y(s));
    B.ylo(end+1,1) = prctile(y(s), 25);
    B.yhi(end+1,1) = prctile(y(s), 75);
end
end

% -----------------------------------------------------------------------
function roi = roi_from_label(lab)
%ROI_FROM_LABEL  Trailing ROI index of a 'Group/Date/rec/roi' label.
p = regexp(lab, '/', 'split');
if numel(p) < 2, roi = NaN; else, roi = str2double(p{end}); end
end

% -----------------------------------------------------------------------
function p = signrank_safe(v)
%SIGNRANK_SAFE  Wilcoxon signed rank vs 0, or NaN if unavailable.
v = v(~isnan(v));
if numel(v) < 6, p = NaN; return; end
try, p = signrank(v); catch, p = NaN; end
end

% -----------------------------------------------------------------------
function figure_population(R, groups, alphaSig, outDir, doSave)
%FIGURE_POPULATION  rho distribution, the common-mode check, and rho vs n.
rho = [R.rho]'; pS = [R.p]'; rhoP = [R.rho_pPop]'; rhoD = [R.rho_deriv]';
sig = pS < alphaSig;
f = figure('Color','w','Units','centimeters','Position',[2 2 30 9]);
set(f,'DefaultAxesFontSize',8);

ax = subplot(1,3,1); hold(ax,'on');
edges = -1:0.05:1;
histogram(ax, rho, edges, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','none');
histogram(ax, rho(sig), edges, 'FaceColor',[0.85 0.20 0.10], 'EdgeColor','none');
xline(ax, 0, 'k-'); xline(ax, median(rho), 'k--','LineWidth',1);
xlabel(ax,'Spearman \rho  (breath amplitude vs \DeltaF/F)'); ylabel(ax,'cells');
title(ax, sprintf('%d cells, %d sig (p<%.2g)\nmedian \\rho = %+.3f', ...
      numel(rho), nnz(sig), alphaSig, median(rho)), 'FontWeight','normal');
box(ax,'off');

% The decisive panel: how much of rho survives removing the field common mode.
% Points on the diagonal are private to the cell; points collapsing to the
% x-axis were field-wide, i.e. the optical-artifact signature.
ax = subplot(1,3,2); hold(ax,'on');
plot(ax, [-1 1], [-1 1], '-','Color',[0.8 0.8 0.8]);
plot(ax, [-1 1], [0 0], '-','Color',[0.8 0.8 0.8]);
scatter(ax, rho(~sig), rhoP(~sig), 14, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha',0.5);
scatter(ax, rho(sig),  rhoP(sig),  18, [0.85 0.20 0.10], 'filled', 'MarkerFaceAlpha',0.8);
axis(ax,[-1 1 -1 1]); axis(ax,'square');
xlabel(ax,'\rho'); ylabel(ax,'\rho after removing FOV common mode');
title(ax, sprintf('common-mode check\nmedian %+.3f \\rightarrow %+.3f', ...
      median(rho), median(rhoP)), 'FontWeight','normal');
box(ax,'off');

ax = subplot(1,3,3); hold(ax,'on');
gi = 0; xt = []; xl = {};
for g = 1:numel(groups)
    m = strcmp({R.group}, groups{g});
    if ~any(m), continue; end
    gi = gi + 1;
    v = rho(m);
    scatter(ax, gi + 0.12*randn(nnz(m),1), v, 12, [0.6 0.6 0.6], 'filled', ...
            'MarkerFaceAlpha',0.6);
    plot(ax, gi + [-0.25 0.25], [median(v) median(v)], 'k-','LineWidth',2);
    xt(end+1) = gi; xl{end+1} = sprintf('%s (%d)', groups{g}, nnz(m)); %#ok<AGROW>
end
yline(ax, 0, 'k--');
set(ax,'XTick',xt,'XTickLabel',xl,'XTickLabelRotation',30);
xlim(ax,[0.4 gi+0.6]); ylabel(ax,'\rho');
title(ax,'by group','FontWeight','normal'); box(ax,'off');

sgtitle(sprintf(['breath amplitude vs \\DeltaF/F, per cycle   |   ' ...
    'derivative response median \\rho = %+.3f'], median(rhoD)), 'FontSize',9);
if doSave
    exportgraphics(f, fullfile(outDir,'population_summary.png'), 'Resolution',300, ...
                   'BackgroundColor','white');
    exportgraphics(f, fullfile(outDir,'population_summary.pdf'), 'ContentType','vector', ...
                   'BackgroundColor','white');
end
close(f);
end

% -----------------------------------------------------------------------
function figure_gallery(R, rowsPerPage, colsPerRow, maxPanels, dotAlpha, alphaSig, outDir, doSave)
%FIGURE_GALLERY  One XY scatter per cell, strongest |rho| first, paged.
%  Points are coloured by IBI so the rate confound is VISIBLE: if the trend runs
%  along the colour gradient, what is being measured is breathing rate.
[~, ord] = sort(abs([R.rho]), 'descend');
ord = ord(1:min(maxPanels, numel(ord)));
per = rowsPerPage*colsPerRow;
nPage = ceil(numel(ord)/per);
for pg = 1:nPage
    take = ord((pg-1)*per + 1 : min(pg*per, numel(ord)));
    f = figure('Color','w','Visible','off','Units','centimeters', ...
               'Position',[1 1 4.6*colsPerRow, 4.0*rowsPerPage + 1.2]);
    set(f,'DefaultAxesFontSize',6.5);
    for k = 1:numel(take)
        i  = take(k);
        ax = subplot(rowsPerPage, colsPerRow, k); hold(ax,'on');
        x = R(i).X(:); y = R(i).Y(:); ib = R(i).IBI(:);
        scatter(ax, x, y, 6, ib, 'filled', 'MarkerFaceAlpha', dotAlpha);
        colormap(ax, parula);
        B = R(i).slopeBin;
        if ~isempty(B.xc)
            % the binned median is on percentile-rank axes; rescale it onto the
            % raw axes so it can be drawn over the scatter
            xr = prctile(x, 100*B.xc);  yr = prctile(y, 100*B.ym);
            plot(ax, xr, yr, '-o','Color',[0.85 0.1 0.1],'LineWidth',1.2, ...
                 'MarkerSize',2.5,'MarkerFaceColor',[0.85 0.1 0.1]);
        end
        mark = ''; if R(i).p < alphaSig, mark = ' *'; end
        title(ax, sprintf('cell %d %s%s\n\\rho %+.2f (pop %+.2f) n=%d', ...
              R(i).cell, R(i).group, mark, R(i).rho, R(i).rho_pPop, R(i).nCyc), ...
              'FontWeight','normal','FontSize',6.5);
        if mod(k-1,colsPerRow) == 0, ylabel(ax,'\DeltaF/F'); end
        if k > numel(take)-colsPerRow, xlabel(ax,'breath amp (z)'); end
        box(ax,'off');
    end
    if doSave
        fn = fullfile(outDir, sprintf('cells_page%02d.png', pg));
        exportgraphics(f, fn, 'Resolution',200, 'BackgroundColor','white');
        d = dir(fn);
        if isempty(d) || d.bytes < 30e3      % OpenGL stub guard
            exportgraphics(f, fn, 'Resolution',200, 'BackgroundColor','white');
        end
    end
    close(f);
end
fprintf('gallery  : %d page(s) of per-cell scatters\n', nPage);
end
