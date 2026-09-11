function phase_rayleigh_260727()
%% phase_rayleigh_260727  Respiratory phase modulation of calcium events, per cell.
%  Rayleigh statistic on the occupancy-corrected phase distribution, after
%  Sirota / Karalis, adapted to a preparation where the published thresholds
%  do not transfer.
% -----------------------------------------------------------------------
% This replaces coherence as the primary breathing-modulation score. Coherence
% asks whether respiration and calcium share power at a stable FREQUENCY;
% it weakens whenever breathing rate drifts, even for a cell that fires at
% exactly the same point of every breath. This asks the latter question directly.
%
% ================== WHY OCCUPANCY CORRECTION IS MANDATORY HERE ==============
% Phase is cycle-interpolated, not Hilbert: foot = 0, peak = pi, next foot = 2pi
% (already computed in cell_pool.mat, so it cannot drift from the coherence
% pipeline). That definition is robust to breathing-rate change -- which is the
% reason to prefer it here -- but inspiration and expiration have very unequal
% durations, so phase is NOT uniformly occupied in time.
%
% MEASURED ON THIS DATASET: only 15% of each cycle is spent in [0,pi), which is
% half the phase axis; the 24-bin occupancy max/min ratio averages 12x and
% reaches 22x. An entirely unmodulated cell, firing uniformly in TIME, therefore
% produces a large and highly significant raw Rayleigh vector pointing into
% expiration. Uncorrected phase statistics on this data measure breathing
% asymmetry, not neural tuning. The ECDF correction is load-bearing, not a
% refinement.
%
%   psi = 2*pi*F(phi),  F = empirical CDF of phase over ALL valid imaging frames
%
% F is built PER RECORDING (occupancy differs between recordings: the fraction
% of time in [0,pi) ranges 0.10 to 0.24 here) and applied before pooling.
%
% ========================= WHAT IS COMPUTED PER CELL ========================
%   n, R = |sum exp(i*psi)|, rbar = R/n
%   Z    = R^2/n                      (Rayleigh)
%   logZ = log(Z)                     the modulation score, variance-stabilised
%   preferred phase -- TWO of them, see below
%   modulation depth  (max-min)/(max+min) of the occupancy-normalised rate profile
%   cycle reliability, both denominators:
%       recruitment = breaths with an event near preferred phase / ALL breaths
%       precision   = breaths with an event near preferred phase / breaths with ANY event
%
% ------------------------- ON PREFERRED PHASE -------------------------------
% Sirota reports preferred phase from RAW phases, because the ECDF distorts the
% biological phase axis. That reasoning holds at mild asymmetry. At the 12x
% occupancy ratio measured here it does not: the raw resultant is dragged toward
% expiration for EVERY cell, modulated or not. So this reports both, and treats
% the occupancy-normalised one as primary:
%   th_rate  PRIMARY. Circular mean of the rate profile (events per bin divided
%            by time spent in that bin). A true rate maximum, on the undistorted
%            biological axis.
%   th_raw   Sirota's convention, kept for comparison. Expect it to sit later in
%            the cycle than th_rate.
%
% ---------------------------- ON SIGNIFICANCE -------------------------------
% Sirota excluded point processes with fewer than 200 events because Z is
% sample-size biased. On this dataset only 15 cells reach 200, so that threshold
% cannot be copied and the asymptotic Rayleigh p-value cannot be trusted either.
% Significance therefore comes ENTIRELY from a cycle-preserving shuffle, which is
% valid at any n because the null is built from the cell's own event train.
%
% TWO nulls are run, because they fail in opposite directions and disagreement
% between them is itself informative:
%   'breath'  circular shift by a whole number of complete breaths (as proposed).
%             Preserves the event train's autocorrelation AND its alignment to
%             cycle boundaries. CAVEAT: when breathing is very regular, a whole-
%             breath shift moves events by nearly a whole cycle, so it barely
%             changes phase and the null becomes conservative. The diagnostic
%             null_phase_shift_rad reports how far the null actually moved the
%             phases -- if it is small, this null has little power and the
%             'uniform' p-value is the honest one.
%   'uniform' circular shift by a uniform random offset. Fully decorrelates phase,
%             preserves count and autocorrelation exactly. Cannot be degenerate.
% Both p-values are reported per cell; primary = 'breath', as requested.
% Benjamini-Hochberg FDR across tested cells.
%
% INCLUSION: pooled event rate >= 1 event/min (79 cells here). Every cell with
% >=1 event is still written to the CSV with its counts and tested = false, so
% nothing is silently dropped.
%
% Input : cell_pool.mat (from cell_pool_260727.m) -- phase, events, breaths, fps
% Output: <phys>\analysis_260727\phase_rayleigh\
%           phase_rayleigh_data.mat, phase_rayleigh_cells.csv,
%           phase_occupancy_qc.png/.pdf
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir), 'coh_ca_breath'));   % so coh_cfg_260727 resolves

cfg = coh_cfg_260727();
if ~isfolder(cfg.phaseDir), mkdir(cfg.phaseDir); end

%% ===================== USER-EDITABLE PARAMETERS ======================
minRate_perMin = 1;        % INCLUSION: pooled events per minute required to be tested
nBins          = 24;       % phase bins for rate profile / modulation depth
relWin_rad     = pi/4;     % +/- window around preferred phase for cycle reliability
nShuffle       = 1000;     % shuffle iterations per cell (per null)
fdr_q          = 0.05;     % Benjamini-Hochberg level
rngSeed        = 260727;   % reproducible nulls
ecdfGridN      = 3600;     % resolution of the per-recording empirical phase CDF
doSave         = true;
% =====================================================================

rng(rngSeed);
fprintf('\n============== phase_rayleigh_260727 ==============\n');
assert(isfile(cfg.poolFile), ['cell_pool.mat not found:\n  %s\nRun cell_pool_260727.m first.'], cfg.poolFile);
P = load(cfg.poolFile,'pool');  pool = P.pool;
rec = pool.rec;  obs = pool.obs;  cells = pool.cells;
fprintf('pool: %d cells, %d usable observations, %d recordings\n', ...
        numel(cells), pool.audit.n_obs_usable, numel(rec));

%% ---- 1. per-recording phase ECDF + occupancy ----
nRec = numel(rec);
E = struct('grid',{},'cdf',{},'occ',{},'occ_frac_first_half',{},'occ_ratio',{}, ...
           'nValid',{},'nBreath',{},'usable',{});
binEdges = linspace(0, 2*pi, nBins+1);
binCtrs  = (binEdges(1:end-1) + binEdges(2:end))/2;
for k = 1:nRec
    e = struct('grid',[],'cdf',[],'occ',zeros(1,nBins),'occ_frac_first_half',NaN, ...
               'occ_ratio',NaN,'nValid',0,'nBreath',0,'usable',false);
    if rec(k).usable
        ph = mod(rec(k).phi, 2*pi);
        ph = ph(~isnan(ph));
        if numel(ph) > 10
            g = linspace(0, 2*pi, ecdfGridN+1);
            c = histcounts(ph, g);
            e.grid = g;
            e.cdf  = [0, cumsum(c)/sum(c)];        % F(g(j)) , F(0)=0, F(2pi)=1
            e.occ  = histcounts(ph, binEdges) / rec(k).fps;    % SECONDS per bin
            e.occ_frac_first_half = mean(ph < pi);
            e.occ_ratio = max(e.occ)/max(min(e.occ), eps);
            e.nValid = numel(ph);
            e.nBreath = max(numel(rec(k).foot_idx)-1, 0);
            e.usable = true;
        end
    end
    E(k) = e; %#ok<AGROW>
end
uk = [E.usable];
fprintf('ECDF built for %d recordings | time in [0,pi): %.3f (uniform=0.500) | occupancy max/min %.1fx\n', ...
        nnz(uk), mean([E(uk).occ_frac_first_half]), mean([E(uk).occ_ratio]));

%% ---- 2. gather each cell's events, and the per-observation material for nulls ----
nCell = numel(cells);
Cell = struct('cell_id',{},'n_obs',{},'n_events',{},'dur_s',{},'rate_perMin',{}, ...
              'n_breaths',{},'tested',{},'why',{}, ...
              'phi_raw',{},'psi',{},'rate_prof',{},'occ_prof',{}, ...
              'R',{},'rbar',{},'Z',{},'logZ',{},'p_rayleigh_asym',{}, ...
              'th_rate',{},'th_raw',{},'mod_depth',{},'mod_depth_cos',{},'th_cos',{},'r2_cos',{}, ...
              'recruitment',{},'precision',{},'n_breaths_with_event',{}, ...
              'p_breath',{},'p_uniform',{},'logZ_null_p95_breath',{},'logZ_null_p95_uniform',{}, ...
              'null_phase_shift_rad',{},'null_logZ_breath',{},'null_logZ_uniform',{},'rec_names',{});

for c = 1:nCell
    ii = cells{c};
    S = struct('cell_id',c,'n_obs',0,'n_events',0,'dur_s',0,'rate_perMin',0, ...
               'n_breaths',0,'tested',false,'why',"", ...
               'phi_raw',[],'psi',[],'rate_prof',nan(1,nBins),'occ_prof',zeros(1,nBins), ...
               'R',NaN,'rbar',NaN,'Z',NaN,'logZ',NaN,'p_rayleigh_asym',NaN, ...
               'th_rate',NaN,'th_raw',NaN,'mod_depth',NaN,'mod_depth_cos',NaN, ...
               'th_cos',NaN,'r2_cos',NaN, ...
               'recruitment',NaN,'precision',NaN,'n_breaths_with_event',0, ...
               'p_breath',NaN,'p_uniform',NaN,'logZ_null_p95_breath',NaN, ...
               'logZ_null_p95_uniform',NaN,'null_phase_shift_rad',NaN, ...
               'null_logZ_breath',single([]),'null_logZ_uniform',single([]), ...
               'rec_names',strings(0,1));
    if isempty(ii), Cell(c) = S; continue; end %#ok<AGROW>

    phiAll = []; psiAll = []; evCount = zeros(1,nBins); occSum = zeros(1,nBins);
    for q = 1:numel(ii)
        o = obs(ii(q));
        if ~o.usable, continue; end
        k = o.rec;  if ~E(k).usable, continue; end
        S.n_obs = S.n_obs + 1;
        S.rec_names(end+1,1) = rec(k).name;
        S.dur_s     = S.dur_s + rec(k).T/rec(k).fps;
        S.n_breaths = S.n_breaths + E(k).nBreath;
        occSum = occSum + E(k).occ;

        fr = find(o.spikes);                       % event frames (train is binary)
        ph = mod(rec(k).phi(fr), 2*pi);
        ph = ph(~isnan(ph));
        if isempty(ph), continue; end
        phiAll = [phiAll; ph(:)]; %#ok<AGROW>
        psiAll = [psiAll; ecdf_apply(E(k), ph(:))]; %#ok<AGROW>
        evCount = evCount + histcounts(ph, binEdges);
    end
    S.n_events    = numel(phiAll);
    S.rate_perMin = S.n_events / max(S.dur_s, eps) * 60;
    S.occ_prof    = occSum;
    S.rate_prof   = evCount ./ max(occSum, eps);      % events per SECOND in each phase bin
    S.phi_raw = phiAll;  S.psi = psiAll;

    if S.n_events == 0
        S.why = "no events";
    elseif S.rate_perMin < minRate_perMin
        S.why = sprintf("rate %.2f/min below the %.2f/min inclusion threshold", S.rate_perMin, minRate_perMin);
    else
        S.tested = true;
    end
    Cell(c) = S; %#ok<AGROW>
end

hasEv  = [Cell.n_events] > 0;
tested = [Cell.tested];
fprintf('cells: %d with >=1 event | %d tested (rate >= %g/min)\n', nnz(hasEv), nnz(tested), minRate_perMin);

%% ---- 3. statistics + shuffle nulls for the tested cells ----
testedIdx = find(tested);
fprintf('shuffling (%d iterations x 2 nulls x %d cells)...\n', nShuffle, numel(testedIdx));
for tt = 1:numel(testedIdx)
    c = testedIdx(tt);
    S = Cell(c);

    [S.R, S.rbar, S.Z, S.logZ] = rayleigh_stats(S.psi);
    S.p_rayleigh_asym = exp(-S.Z);                     % asymptotic; reported, NOT used
    S.th_raw  = angle(sum(exp(1i*S.phi_raw)));
    S.th_rate = angle(sum(S.rate_prof .* exp(1i*binCtrs)));
    mx = max(S.rate_prof); mn = min(S.rate_prof);
    S.mod_depth = (mx - mn)/max(mx + mn, eps);
    % NOTE on mod_depth: (max-min)/(max+min) saturates at exactly 1 as soon as ANY
    % bin is empty, which at a median of 70 events over 24 bins is nearly every
    % cell -- so it carries almost no information here. mod_depth_cos below is the
    % usable effect size; mod_depth is kept only for continuity.
    [S.mod_depth_cos, S.th_cos, S.r2_cos] = cosine_depth(S.rate_prof, binCtrs, S.occ_prof);

    % ---- cycle reliability, both denominators, around th_rate ----
    [S.recruitment, S.precision, S.n_breaths_with_event] = ...
        cycle_reliability(cells{c}, obs, rec, E, S.th_rate, relWin_rad);

    % ---- nulls ----
    [nb, du] = deal(nan(nShuffle,1));
    shiftRad = nan(nShuffle,1);
    for s = 1:nShuffle
        [psiB, dphi] = shuffled_psi(cells{c}, obs, rec, E, 'breath');
        [psiU, ~   ] = shuffled_psi(cells{c}, obs, rec, E, 'uniform');
        if ~isempty(psiB), [~,~,~,nb(s)] = rayleigh_stats(psiB); end
        if ~isempty(psiU), [~,~,~,du(s)] = rayleigh_stats(psiU); end
        shiftRad(s) = dphi;
    end
    nb = nb(isfinite(nb));  du = du(isfinite(du));
    S.p_breath  = (1 + nnz(nb >= S.logZ)) / (1 + numel(nb));   % +1: never report p = 0
    S.p_uniform = (1 + nnz(du >= S.logZ)) / (1 + numel(du));
    S.logZ_null_p95_breath  = prctile(nb, 95);
    S.logZ_null_p95_uniform = prctile(du, 95);
    S.null_phase_shift_rad  = mean(shiftRad, 'omitnan');
    S.null_logZ_breath      = single(nb);    % kept so the per-cell QC can show the
    S.null_logZ_uniform     = single(du);    %   actual null, not just its 95th pct

    Cell(c) = S;
    if mod(tt, 10) == 0 || tt == numel(testedIdx)
        fprintf('  %3d/%3d cells\n', tt, numel(testedIdx));
    end
end

%% ---- 4. FDR over tested cells ----
sigB = false(1,nCell);  sigU = false(1,nCell);
qB = nan(1,nCell);      qU = nan(1,nCell);
if ~isempty(testedIdx)
    [sigB(testedIdx), qB(testedIdx)] = bh_fdr([Cell(testedIdx).p_breath],  fdr_q);
    [sigU(testedIdx), qU(testedIdx)] = bh_fdr([Cell(testedIdx).p_uniform], fdr_q);
end

fprintf('\n---- results ----\n');
fprintf('  tested cells                    : %d\n', numel(testedIdx));
fprintf('  significant, breath-shift null  : %d  (BH q<%.2f)\n', nnz(sigB), fdr_q);
fprintf('  significant, uniform-shift null : %d  (BH q<%.2f)\n', nnz(sigU), fdr_q);
if ~isempty(testedIdx)
    msr = mean([Cell(testedIdx).null_phase_shift_rad], 'omitnan');
    fprintf('  breath-shift null moved phases by %.2f rad on average (pi/2=1.57 would be a full decorrelation)\n', msr);
    if msr < 0.5
        fprintf(['  NOTE: the breath-shift null barely moves phase on this data, so it is\n' ...
                 '        conservative -- treat the uniform-shift p-value as the honest one.\n']);
    end
    fprintf('  logZ over tested cells: median %.2f  max %.2f\n', ...
            median([Cell(testedIdx).logZ]), max([Cell(testedIdx).logZ]));
    dth = angdiff_local([Cell(testedIdx).th_rate], [Cell(testedIdx).th_raw]);
    fprintf('  preferred phase: rate-normalised vs raw differ by median %.0f deg (occupancy bias)\n', ...
            rad2deg(median(abs(dth))));
end

%% ---- 5. occupancy QC figure ----
fq = figure('Color','w','Name','phase occupancy QC','Units','centimeters','Position',[2 2 30 10]);
ax1 = subplot(1,3,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
for k = find(uk)
    plot(ax1, binCtrs, E(k).occ/sum(E(k).occ), '-', 'Color',[.72 .72 .72], 'LineWidth',0.6);
end
Ok = cell2mat(arrayfun(@(e) e.occ/sum(e.occ), E(uk), 'UniformOutput',false)');
plot(ax1, binCtrs, mean(Ok,1), 'r-', 'LineWidth',2);
yline(ax1, 1/nBins, 'k--','LineWidth',1);
xlim(ax1,[0 2*pi]); set(ax1,'XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
xlabel(ax1,'breath phase (foot=0, peak=\pi)'); ylabel(ax1,'fraction of TIME');
title(ax1,{'phase is not uniformly occupied','dashed = uniform; red = mean across recordings'},'FontSize',8);

ax2 = subplot(1,3,2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for k = find(uk), plot(ax2, E(k).grid, E(k).cdf, '-', 'Color',[.72 .72 .72], 'LineWidth',0.6); end
plot(ax2, [0 2*pi],[0 1],'k--','LineWidth',1);
xlim(ax2,[0 2*pi]); ylim(ax2,[0 1]);
set(ax2,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
xlabel(ax2,'raw phase \phi'); ylabel(ax2,'F(\phi)');
title(ax2,{'the ECDF transform \psi = 2\pi F(\phi)','dashed = the identity it would be if occupancy were uniform'},'FontSize',8);

ax3 = subplot(1,3,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
if any(tested)
    allRaw = vertcat(Cell(tested).phi_raw);
    allPsi = vertcat(Cell(tested).psi);
    histogram(ax3, allRaw, binEdges, 'Normalization','probability', ...
              'FaceColor',[.85 .33 .10],'FaceAlpha',0.55, 'DisplayName','raw \phi');
    histogram(ax3, allPsi, binEdges, 'Normalization','probability', ...
              'FaceColor',[0 .45 .74],'FaceAlpha',0.55, 'DisplayName','corrected \psi');
end
yline(ax3, 1/nBins, 'k--','LineWidth',1,'HandleVisibility','off');
xlim(ax3,[0 2*pi]); set(ax3,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
xlabel(ax3,'phase'); ylabel(ax3,'fraction of ALL events'); legend(ax3,'Location','best');
title(ax3,{'pooled event phases before and after correction','raw piles into expiration because that is where the time is'},'FontSize',8);
sgtitle(sprintf('phase occupancy QC  |  %d recordings, %.1f%% of each cycle in [0,\\pi), occupancy ratio %.1fx', ...
        nnz(uk), 100*mean([E(uk).occ_frac_first_half]), mean([E(uk).occ_ratio])));

%% ---- 6. save ----
if doSave
    exportgraphics(fq, fullfile(cfg.phaseDir,'phase_occupancy_qc.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fq, fullfile(cfg.phaseDir,'phase_occupancy_qc.pdf'),'ContentType','vector','BackgroundColor','white');

    rep = find(hasEv);
    T = table([Cell(rep).cell_id]', [Cell(rep).n_obs]', [Cell(rep).n_events]', ...
              [Cell(rep).dur_s]', [Cell(rep).rate_perMin]', [Cell(rep).n_breaths]', ...
              [Cell(rep).tested]', ...
              [Cell(rep).logZ]', [Cell(rep).Z]', [Cell(rep).rbar]', [Cell(rep).R]', ...
              rad2deg([Cell(rep).th_rate]'), rad2deg([Cell(rep).th_raw]'), ...
              [Cell(rep).mod_depth]', [Cell(rep).mod_depth_cos]', rad2deg([Cell(rep).th_cos]'), ...
              [Cell(rep).r2_cos]', [Cell(rep).recruitment]', [Cell(rep).precision]', ...
              [Cell(rep).n_breaths_with_event]', ...
              [Cell(rep).p_breath]', qB(rep)', sigB(rep)', ...
              [Cell(rep).p_uniform]', qU(rep)', sigU(rep)', ...
              [Cell(rep).p_rayleigh_asym]', [Cell(rep).null_phase_shift_rad]', ...
              arrayfun(@(s) strjoin(cellstr(s.rec_names),'|'), Cell(rep), 'UniformOutput',false)', ...
              [Cell(rep).why]', ...
        'VariableNames', {'cell_id','n_obs','n_events','duration_s','rate_per_min','n_breaths', ...
                          'tested','logZ','Z','rbar','R','pref_phase_rate_deg','pref_phase_raw_deg', ...
                          'mod_depth','mod_depth_cos','pref_phase_cos_deg','r2_cos', ...
                          'recruitment','precision','n_breaths_with_event', ...
                          'p_breathshift','q_breathshift','sig_breathshift', ...
                          'p_uniformshift','q_uniformshift','sig_uniformshift', ...
                          'p_rayleigh_asymptotic','null_phase_shift_rad','recordings','excluded_reason'});
    writetable(T, fullfile(cfg.phaseDir,'phase_rayleigh_cells.csv'));

    params = struct('minRate_perMin',minRate_perMin,'nBins',nBins,'relWin_rad',relWin_rad, ...
                    'nShuffle',nShuffle,'fdr_q',fdr_q,'rngSeed',rngSeed,'ecdfGridN',ecdfGridN); %#ok<NASGU>
    save(fullfile(cfg.phaseDir,'phase_rayleigh_data.mat'), ...
         'Cell','E','binEdges','binCtrs','sigB','sigU','qB','qU','params','cfg','-v7.3');
    fprintf('\nSaved phase_rayleigh_data.mat / _cells.csv / phase_occupancy_qc.png+pdf to\n  %s\n', cfg.phaseDir);
end
fprintf('Next: phase_rayleigh_polar_260727.m\n');
end

%% ========================= helpers =========================
function psi = ecdf_apply(e, phi)
% psi = 2*pi*F(phi) using this recording's empirical phase-occupancy CDF.
% The grid is uniform on [0,2pi], so this is direct index arithmetic rather than
% interp1 -- it is called ~2 million times inside the shuffle loop.
p = mod(phi(:), 2*pi);
n = numel(e.cdf) - 1;                       % = ecdfGridN
x = p/(2*pi)*n;
i0 = min(floor(x), n-1);                    % 0-based cell index
f  = x - i0;
c  = e.cdf(:);
psi = 2*pi * ((1-f).*c(i0+1) + f.*c(i0+2));
psi = mod(psi, 2*pi);
end

function [depth, mu, r2] = cosine_depth(rate, ctrs, occ)
% Effect size from a cosine fit to the occupancy-normalised rate profile:
%       rate(theta) ~ a + b*cos(theta - mu),   depth = b/a
% Robust to empty bins, which is exactly where (max-min)/(max+min) fails at these
% event counts. Matches the a + b*cos(phase - mu) convention already used
% elsewhere in this project.
%
% WEIGHTED by occupancy: for a Poisson rate estimated as N/T, var ~ lambda/T, so
% precision scales with the time spent in the bin. Inspiration bins hold ~1/7 of
% the time here and are correspondingly noisier, so they are weighted down. That
% is statistically correct, not a preference for expiration -- the fitted mu is
% still free to point anywhere.
depth = NaN; mu = NaN; r2 = NaN;
rate = rate(:); ctrs = ctrs(:); w = occ(:);
ok = isfinite(rate) & isfinite(w) & w > 0;
if nnz(ok) < 4, return; end
rate = rate(ok); ctrs = ctrs(ok); w = w(ok);
X  = [ones(numel(ctrs),1), cos(ctrs), sin(ctrs)];
W  = sqrt(w/sum(w));
beta = (X .* W) \ (rate .* W);
a = beta(1);  b = hypot(beta(2), beta(3));
mu = atan2(beta(3), beta(2));
if a > 0, depth = b/a; end
fit  = X*beta;
mw   = sum(w.*rate)/sum(w);
ssr  = sum(w.*(rate - fit).^2);
sst  = sum(w.*(rate - mw).^2);
if sst > 0, r2 = 1 - ssr/sst; end
end

function [R, rbar, Z, logZ] = rayleigh_stats(psi)
n = numel(psi);
if n == 0, [R,rbar,Z,logZ] = deal(NaN); return; end
C = sum(cos(psi));  S = sum(sin(psi));
R = hypot(C, S);  rbar = R/n;  Z = R^2/n;
logZ = log(max(Z, realmin));
end

function [psi, meanShiftRad] = shuffled_psi(ii, obs, rec, E, mode)
% One null realisation: circularly shift each observation's event train, then
% re-read the phases at the shifted frames. Shifting the INDEX rather than the
% train keeps this O(n events) instead of O(T).
psi = [];  d_ph = [];
for q = 1:numel(ii)
    o = obs(ii(q));
    if ~o.usable, continue; end
    k = o.rec;  if ~E(k).usable, continue; end
    T = rec(k).T;
    fr = find(o.spikes);
    if isempty(fr), continue; end
    switch mode
        case 'breath'
            % shift by a whole number of complete breaths
            f = sort(rec(k).foot_idx(:));
            if numel(f) < 3, d = randi(T) - 1; else
                j = randi(numel(f)-1) + 1;  d = f(j) - f(1);
            end
        otherwise                              % 'uniform'
            d = randi(T) - 1;
    end
    fr2 = mod(fr - d - 1, T) + 1;
    p1 = mod(rec(k).phi(fr),  2*pi);
    p2 = mod(rec(k).phi(fr2), 2*pi);
    ok = ~isnan(p2);
    if ~any(ok), continue; end
    psi = [psi; ecdf_apply(E(k), p2(ok))]; %#ok<AGROW>
    both = ok & ~isnan(p1);
    if any(both), d_ph = [d_ph; abs(angdiff_local(p2(both), p1(both)))]; end %#ok<AGROW>
end
meanShiftRad = mean(d_ph, 'omitnan');
end

function [recruit, precis, nBrEv] = cycle_reliability(ii, obs, rec, E, thPref, win)
% recruitment = breaths with an event within +/-win of thPref / ALL breaths
% precision   = breaths with an event within +/-win of thPref / breaths with ANY event
nBrTot = 0; nBrEv = 0; nBrNear = 0;
for q = 1:numel(ii)
    o = obs(ii(q));
    if ~o.usable, continue; end
    k = o.rec;  if ~E(k).usable, continue; end
    f = sort(rec(k).foot_idx(:));
    if numel(f) < 2, continue; end
    nBrTot = nBrTot + numel(f) - 1;
    fr = find(o.spikes);
    if isempty(fr), continue; end
    ph = mod(rec(k).phi(fr), 2*pi);
    b  = discretize(fr, f);                       % which breath each event fell in
    ok = ~isnan(b) & ~isnan(ph);
    if ~any(ok), continue; end
    b = b(ok); ph = ph(ok);
    near = abs(angdiff_local(ph, thPref)) <= win;
    nBrEv   = nBrEv   + numel(unique(b));
    nBrNear = nBrNear + numel(unique(b(near)));
end
recruit = nBrNear / max(nBrTot, 1);
precis  = nBrNear / max(nBrEv, 1);
if nBrEv == 0, precis = NaN; end
end

function d = angdiff_local(a, b)
% signed smallest angular difference a - b, wrapped to [-pi, pi]
d = angle(exp(1i*(a(:) - b(:))));
d = reshape(d, size(a));
end

function [sig, q] = bh_fdr(p, level)
% Benjamini-Hochberg. Returns the rejection mask and the adjusted q-values.
p = p(:); n = numel(p);
[ps, ord] = sort(p);
qs = min(1, ps .* n ./ (1:n)');
for i = n-1:-1:1, qs(i) = min(qs(i), qs(i+1)); end
q = nan(n,1);  q(ord) = qs;
sig = q <= level;
sig = reshape(sig, 1, []);  q = reshape(q, 1, []);
end
