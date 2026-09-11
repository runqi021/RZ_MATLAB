function breath_time_peth_260727(triggerIn)
%% breath_time_peth_260727  PRIMARY breathing-modulation analysis.
%  Absolute-time inspiration-triggered PETH, scored against a circular-shift null.
% -----------------------------------------------------------------------
% The question this answers:
%
%   Does this cell's event activity have a reproducible temporal relationship to
%   inspiration onset, compared with shuffled versions of the SAME event train?
%
% There is no designated respiratory baseline and no assumption that breathing is
% rhythmic, stationary, or describable by one frequency. Every cell is compared
% only against itself.
%
% ===================== WHY NOT PHASE, WHY NOT COHERENCE =====================
% COHERENCE needs breathing to occupy a stable frequency band; it weakens when the
% rate drifts, even for a cell that fires at the same point of every breath.
% PHASE (foot=0, peak=pi, next foot=2pi) is robust to rate drift but is NOT
% uniformly occupied in time: measured here, inspiration is ~15% of the cycle, so
% half the phase axis holds ~15% of the time and the other half ~85%. An
% unmodulated cell then looks strongly modulated until an ECDF correction is
% applied. Worse, under irregular breathing or a long pause, normalised phase
% erases the actual latency and the actual pause duration.
% ABSOLUTE TIME has neither problem: every trigger contributes to every time bin,
% so the null is flat with no correction of any kind.
% Both are kept as secondary descriptive measures; neither is the classifier.
%
% ========================= THE STATISTIC ====================================
% For the observed PETH r(t), with rs = lightly smoothed r and rbar = mean(r):
%
%     T_exc = max over the search window of [ rs(t) - rbar ]
%     T_sup = max over the search window of [ rbar - rs(t) ]
%
% Every shuffle recomputes the SAME quantity, including its own rbar and its own
% max. Because the maximum is re-taken under the null, the search across candidate
% latencies is paid for automatically -- no separate correction is needed.
%
%     M_exc = ( T_exc,obs - mean(T_exc,shuffle) ) / std(T_exc,shuffle)
%     p_exc = ( 1 + #{ T_exc,s >= T_exc,obs } ) / ( 1 + nShuffle )
%
% and the same, separately, for suppression. Excitation and suppression are never
% silently combined; a signed summary is provided but both originals are kept.
%
% Significance is called on the empirical p (rank-based, valid at any n).
% M is an EFFECT SIZE and is unreliable below ~40 events, where the null of a
% max-statistic is discrete and right-skewed and a z-score misdescribes it.
%
% ==================== QUIET BUT PRECISELY LOCKED CELLS ======================
% An all-breath PETH measures the average effect over all inspirations. A cell can
% be recruited on few breaths yet fire at a very reproducible latency when it does.
% So timing precision is measured separately, from event latencies rather than
% from the rate:
%   every event is assigned to its IMMEDIATELY PRECEDING accepted onset (so one
%   event is never counted for several neighbouring breaths), giving a latency dt;
%   recruitment = accepted breaths with an event near the preferred latency, over
%                 ALL accepted breaths;
%   precision   = how tightly those latencies concentrate.
% Precision is shuffle-tested with the preferred latency RE-ESTIMATED inside every
% shuffle, otherwise it is circular and inflated.
%
% ============================ IMPLEMENTATION ================================
% * Exposure is counted in FRAMES, not triggers. At 30 fps a frame lands every
%   33.3 ms, so a 50 ms bin holds one frame sometimes and two others; dividing by
%   trigger count leaves a sawtooth that is pure aliasing. Recordings here run
%   30-47 fps, so exposure accrues per recording at its own rate.
% * Only boundary-safe triggers are used -- those whose whole window lies inside
%   the recording. For them the circular FFT cross-correlation IS the linear one,
%   so there is no wraparound in the observed PETH while the shuffle keeps the
%   circular shift it needs. One line, exact, and the full correlogram is still
%   computed once per observation so each shuffle is only an index shift.
% * Triggers are screened by breath amplitude (see TRIGGER QC below).
% * Zero-event observations are RETAINED: they contribute duration, triggers and
%   exposure with zero counts. Dropping them inflates the cell's rate and its
%   recruitment denominator.
% * Bin width is FIXED and smoothing is specified in SECONDS, so latency is
%   comparable across datasets whose breath rates differ. Only the window EXTENT
%   is derived from the measured cycle.
%
% Input : cell_pool.mat
% Output: <dataset>\analysis_260727\breath_time\
%           breath_time_peth_data.mat   full observed + null distributions
%           breath_time_peth_cells.csv  one row per cell
% Figures: breath_time_peth_percell_260727.m, breath_time_modulation_scatter_260727.m
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir), 'coh_ca_breath'));

cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
% --- WHAT IS ZERO ---
TRIGGER = 'onset';     % 'onset' | 'peak'   -- run it BOTH ways, they answer
                       % different questions.
                       %
                       % 'onset' zero = inspiration onset. This is foot_idx, from
                       %         breath_insp_start_pc1.mat -> insp_start_idx: the foot
                       %         of the breath PC1, where the rise begins.
                       % 'peak'  zero = the inspiratory peak. This is peak_idx, from
                       %         breath_peak_pc1.mat -> insp_onset_idx. NOTE that
                       %         upstream field name says "onset" but it is the PEAK;
                       %         the two landmarks are 267 ms apart on this data.
                       %
                       % WHY BOTH. Onset and peak are separated by 267 ms with CV 0.15
                       % (SD ~40 ms). A cell locked to the peak is smeared by that
                       % jitter when aligned on the onset, and vice versa. Whichever
                       % alignment gives the larger M_exc and the smaller latency MAD
                       % is the landmark that cell is actually locked to -- which is a
                       % result, not a robustness check.
                       % Outputs go to breath_time\<TRIGGER>\ so both coexist.
if nargin >= 1 && ~isempty(triggerIn), TRIGGER = char(triggerIn); end
% --- window: extent from the data, resolution fixed ---
winPreFrac   = 0.25;   % window start = -winPreFrac  x median breath cycle
winPostFrac  = 1.25;   % window end   = +winPostFrac x median breath cycle
searchFrac   = 1.00;   % peak/trough searched over 0 .. searchFrac x median cycle
binWidth_s   = 0.050;  % FIXED. Must be >= the slowest frame period (33.3 ms at 30 fps).
                       %   Deliberately not derived from the cycle: an adaptive bin
                       %   made latency incomparable between datasets and quantised
                       %   it so coarsely that the IQR came out one bin wide.
smoothWidth_s= 0.150;  % smoothing in SECONDS, not bins, for the same reason

% --- trigger QC (recomputed here so the analysis can never drift from the diagnostic) ---
ampFrac      = 0.20;   % a cycle whose peak-minus-foot breath amplitude is below this
                       %   fraction of its RECORDING's median is a failed detection,
                       %   not a breath. Amplitude only -- a long respiratory pause is
                       %   never rejected for being long.

% --- inclusion ---
minRate_perMin = 1;    % pooled events per minute
minEvents      = 20;   % absolute floor as well: on short recordings a rate criterion
                       %   alone admitted 4-event cells and nothing was detectable

% --- precision / recruitment ---
relWin_s     = 0.100;  % +/- tolerance around the preferred latency

% --- suppression gate ---
% T_sup <= rbar by construction (rate cannot go below zero). When most smoothed bins
% are empty, nearly every shuffle hits that ceiling, sigma collapses and the z-score
% is noise. Require enough expected events per smoothing window for the statistic to
% have any dynamic range at all; otherwise report sup_tested = false rather than a
% column that looks meaningful and is not.
minExpectedPerSmoothWin = 5;

% --- nulls ---
% The smallest p a shuffle can return is 1/(nShuffle+1); BH needs p <= q/nTests at
% rank 1, so the shuffle count must scale with the number of tests or the top cells
% cannot be called at all.
shufflePerTest = 40;
nShuffleMin    = 1000;
nShuffleMax    = 20000;
fdr_q          = 0.05;
rngSeed        = 260727;
doSave         = true;
% =====================================================================

outDir = fullfile(cfg.outRoot, 'breath_time', lower(TRIGGER));
if ~isfolder(outDir), mkdir(outDir); end

rng(rngSeed);
fprintf('\n============ breath_time_peth_260727 [%s-triggered] ============\n', upper(TRIGGER));
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
pool = ensure_pool_260727();   % builds cell_link + cell_pool if they are missing
rec = pool.rec; obs = pool.obs; cells = pool.cells;

%% ---- 1. window from the measured breath cycle; resolution fixed ----
cycAll = [];
for k = 1:numel(rec)
    if ~rec(k).usable, continue; end
    f = sort(rec(k).foot_idx(:));
    if numel(f) > 1, cycAll = [cycAll; diff(f)/rec(k).fps]; end %#ok<AGROW>
end
medCyc = median(cycAll);
preLag_s = winPreFrac*medCyc;  postLag_s = winPostFrac*medCyc;
edges = -preLag_s : binWidth_s : postLag_s;
ctrs  = edges(1:end-1) + binWidth_s/2;
nB    = numel(ctrs);
srchM = ctrs >= 0 & ctrs <= searchFrac*medCyc;
smoothBins = max(1, round(smoothWidth_s/binWidth_s));
fprintf('median breath cycle %.3f s (%d cycles)\n', medCyc, numel(cycAll));
fprintf('window %.2f..%.2f s | %.0f ms bins (%d) | smoothing %.0f ms (%d bins) | search 0..%.2f s\n', ...
        edges(1), edges(end), 1000*binWidth_s, nB, 1000*smoothWidth_s, smoothBins, searchFrac*medCyc);

%% ---- 2. TRIGGER QC + boundary-safe trigger set, per recording ----
TRIG = cell(numel(rec),1);   % accepted, boundary-safe onset frames
qcTot = 0; qcRej = 0; bndDrop = 0;
for k = 1:numel(rec)
    if ~rec(k).usable, continue; end
    fs = rec(k).fps; T = rec(k).T;
    f = sort(rec(k).foot_idx(:));  p = sort(rec(k).peak_idx(:));
    if numel(f) < 2, continue; end
    bw  = rec(k).bw(:);
    bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
    amp = nan(numel(f)-1,1);  pkFr = nan(numel(f)-1,1);
    for i = 1:numel(f)-1
        pk = p(p > f(i) & p < f(i+1));
        if ~isempty(pk), amp(i) = bwz(pk(1)) - bwz(f(i)); pkFr(i) = pk(1); end
    end
    ok = amp > ampFrac*median(amp,'omitnan');        % amplitude ONLY; long pauses survive
    qcTot = qcTot + nnz(~isnan(amp));  qcRej = qcRej + nnz(~isnan(amp) & ~ok);
    % zero is either the foot of the accepted cycle (onset) or its peak
    if strcmpi(TRIGGER,'peak'), cand = pkFr(ok); else, cand = f(ok); end
    cand = cand(~isnan(cand));
    loFr = ceil(edges(1)*fs);  hiFr = floor(edges(end)*fs);
    safe = cand + loFr >= 1 & cand + hiFr <= T;       % whole window inside the recording
    bndDrop = bndDrop + nnz(~safe);
    TRIG{k} = cand(safe);
end
fprintf('trigger QC: %d of %d cycles rejected on amplitude (%.1f%%) | %d more dropped at recording edges\n', ...
        qcRej, qcTot, 100*qcRej/max(qcTot,1), bndDrop);

%% ---- 3. per-observation correlogram + exposure (ZERO-EVENT OBSERVATIONS KEPT) ----
fprintf('building cross-correlograms for %d observations...\n', numel(obs));
OB = struct('cell',{},'recName',{},'ccf',{},'lagIdx',{},'binOf',{},'expSec',{}, ...
            'T',{},'fps',{},'trig',{},'nTrig',{},'evFr',{},'nEv',{},'footShift',{});
for i = 1:numel(obs)
    o = obs(i);
    if ~o.usable || isnan(o.cell_id), continue; end
    k = o.rec; if isnan(k) || ~rec(k).usable, continue; end
    f = TRIG{k};
    if numel(f) < 2, continue; end                    % no usable triggers, not "no events"
    T = rec(k).T; fs = rec(k).fps;
    ev = full(double(o.spikes(:)));
    if numel(ev) < T, ev(end+1:T,1) = 0; end
    ev = ev(1:T);
    % NOTE: no `sum(ev)==0` test. A cell silent in one recording must still
    % contribute that recording's duration, triggers and exposure, or its rate and
    % its recruitment denominator are inflated by however many recordings it was
    % quiet in.

    trig = zeros(T,1); trig(f) = 1;
    ccf  = real(ifft(conj(fft(trig)) .* fft(ev)));    % ccf(m+1) = sum_t trig(t)*ev(t+m)
    % With only boundary-safe triggers in `trig`, f+m never leaves [1,T] across the
    % whole window, so this circular correlation IS the linear one -- no wraparound.

    loFr = ceil(edges(1)*fs);  hiFr = floor(edges(end)*fs);
    m  = (loFr:hiFr)';
    bo = discretize(m/fs, edges);
    keep = ~isnan(bo);  m = m(keep);  bo = bo(keep);

    e = struct('cell',o.cell_id, 'recName',rec(k).name, 'ccf',ccf, ...
               'lagIdx',mod(m,T)+1, 'binOf',bo, ...
               'expSec',accumarray(bo,1,[nB 1])*numel(f)/fs, ...
               'T',T, 'fps',fs, 'trig',f, 'nTrig',numel(f), ...
               'evFr',find(ev>0), 'nEv',nnz(ev>0), 'footShift',f-f(1));
    OB(end+1) = e; %#ok<AGROW>
end
obsCell = [OB.cell];
fprintf('  %d observations usable (%d of them with zero events, retained)\n', ...
        numel(OB), nnz([OB.nEv]==0));

%% ---- 4. how many cells will be tested -> how many shuffles ----
nCell = numel(cells);
nTestPre = 0;
for c = 1:nCell
    ii = find(obsCell == c);
    if isempty(ii), continue; end
    ev = sum([OB(ii).nEv]);
    du = sum(arrayfun(@(o) o.T/o.fps, OB(ii)));
    if ev >= minEvents && ev/max(du,eps)*60 >= minRate_perMin, nTestPre = nTestPre + 1; end
end
nShuffle = min(max(shufflePerTest*nTestPre, nShuffleMin), nShuffleMax);
pFloor = 1/(nShuffle+1);  bhThr = fdr_q/max(nTestPre,1);
fprintf('%d cells to test -> %d shuffles (p floor %.2e vs BH rank-1 threshold %.2e)\n', ...
        nTestPre, nShuffle, pFloor, bhThr);
if pFloor > bhThr
    warning('Shuffle resolution %.2e is coarser than the BH rank-1 threshold %.2e.', pFloor, bhThr);
end

%% ---- 5. per-cell statistics ----
R = init_results(nCell);
% Latencies are now signed distances to the NEAREST trigger, so they live in
% roughly [-halfCyc, +halfCyc] -- anything further away has a nearer trigger. The
% precision search runs over that same span, not over [0, cycle].
halfCyc  = medCyc/2;
searchLo = -halfCyc;  searchHi = halfCyc;
tt = 0;
for c = 1:nCell
    ii = find(obsCell == c);
    S = R(c);  S.cell_id = c;
    if isempty(ii), S.why_not_tested = "no usable observation"; R(c) = S; continue; end

    S.n_obs      = numel(ii);
    S.n_events   = sum([OB(ii).nEv]);
    S.duration_s = sum(arrayfun(@(o) o.T/o.fps, OB(ii)));
    S.n_accepted_breaths = sum([OB(ii).nTrig]);
    S.rate_per_min = S.n_events/max(S.duration_s,eps)*60;
    S.rec_names  = unique(string({OB(ii).recName}'),'stable');
    S.expSec     = sum(cat(2, OB(ii).expSec), 2)';

    if S.n_events < minEvents || S.rate_per_min < minRate_perMin
        S.why_not_tested = sprintf("%d events (need %d), %.2f/min (need %.2f)", ...
                                    S.n_events, minEvents, S.rate_per_min, minRate_perMin);
        S.peth = peth_counts(OB(ii), nB, 0) ./ max(S.expSec,eps);
        R(c) = S; continue;
    end
    S.tested = true;  tt = tt + 1;

    % ---------- observed ----------
    cnt = peth_counts(OB(ii), nB, 0);
    S.peth = cnt ./ max(S.expSec,eps);
    [S.T_exc, S.T_sup, S.preferred_latency_s, S.suppression_latency_s, ...
     S.mean_peth_rate, S.peak_rate, S.trough_rate] = peth_stats(S.peth, ctrs, srchM, smoothBins);
    S.peak_over_mean = S.peak_rate / max(S.mean_peth_rate,eps);

    % The PETH peak is found in [0, cycle]; a nearest-trigger latency cannot exceed
    % half a cycle, so a preferred latency in the second half of the cycle is the
    % SAME moment expressed as a negative latency to the NEXT breath. Wrap it, or
    % recruitment would look for events in a window no event can ever fall in.
    S.pref_latency_nearest_s = S.preferred_latency_s;
    if S.pref_latency_nearest_s >  halfCyc, S.pref_latency_nearest_s = S.pref_latency_nearest_s - medCyc; end
    if S.pref_latency_nearest_s < -halfCyc, S.pref_latency_nearest_s = S.pref_latency_nearest_s + medCyc; end

    dtObs = event_latencies(OB(ii), 0, halfCyc);
    S.n_assigned_events = numel(dtObs);
    S.dt_obs = single(dtObs);          % kept so the per-cell QC can plot the actual
                                       % latency distribution, not just its summary
    [S.precision_fraction, S.precision_tau_s] = max_window_frac(dtObs, relWin_s, searchLo, searchHi);
    S.precision_at_pref = mean(abs(dtObs - S.pref_latency_nearest_s) <= relWin_s);
    S.latency_median_s  = median(dtObs);
    S.latency_mad_s     = median(abs(dtObs - S.latency_median_s));
    [S.recruitment, S.n_breaths_recruited] = recruitment_of(OB(ii), 0, S.pref_latency_nearest_s, relWin_s, halfCyc);

    % suppression gate: expected events per smoothing window
    expPerWin = S.n_events * (smoothWidth_s/(edges(end)-edges(1)));
    S.sup_expected_per_window = expPerWin;
    S.sup_tested = expPerWin >= minExpectedPerSmoothWin;

    % ---------- nulls ----------
    nullExc = nan(nShuffle,1); nullSup = nan(nShuffle,1); nullPre = nan(nShuffle,1);
    pethAcc = zeros(1,nB); pethAcc2 = zeros(1,nB);
    for s = 1:nShuffle
        d = arrayfun(@(o) randi(o.T)-1, OB(ii));       % uniform circular shift per observation
        cs = peth_counts(OB(ii), nB, d);
        rs = cs ./ max(S.expSec,eps);
        [nullExc(s), nullSup(s)] = peth_stats(rs, ctrs, srchM, smoothBins);
        dts = event_latencies(OB(ii), d, halfCyc);
        nullPre(s) = max_window_frac(dts, relWin_s, searchLo, searchHi);  % tau RE-ESTIMATED
        pethAcc = pethAcc + rs;  pethAcc2 = pethAcc2 + rs.^2;
    end
    S.peth_shuffle_mean = pethAcc/nShuffle;
    S.peth_shuffle_sd   = sqrt(max(pethAcc2/nShuffle - S.peth_shuffle_mean.^2, 0));
    S.null_exc = single(nullExc); S.null_sup = single(nullSup); S.null_pre = single(nullPre);

    S.mod_exc_z = zscore_against(S.T_exc, nullExc);
    S.p_exc     = (1 + nnz(nullExc >= S.T_exc))/(1 + nShuffle);
    if S.sup_tested
        S.mod_sup_z = zscore_against(S.T_sup, nullSup);
        S.p_sup     = (1 + nnz(nullSup >= S.T_sup))/(1 + nShuffle);
    end
    S.precision_z = zscore_against(S.precision_fraction, nullPre);
    S.p_precision = (1 + nnz(nullPre >= S.precision_fraction))/(1 + nShuffle);

    R(c) = S;
    if mod(tt,10) == 0, fprintf('  %d/%d cells\n', tt, nTestPre); end
end

%% ---- 6. FDR, per family ----
tested = [R.tested];  ti = find(tested);
[q_exc, sig_exc] = deal(nan(1,nCell), false(1,nCell));
[q_sup, sig_sup] = deal(nan(1,nCell), false(1,nCell));
[q_pre, sig_pre] = deal(nan(1,nCell), false(1,nCell));
if ~isempty(ti)
    [sig_exc(ti), q_exc(ti)] = bh_fdr([R(ti).p_exc], fdr_q);
    si = ti([R(ti).sup_tested]);
    if ~isempty(si), [sig_sup(si), q_sup(si)] = bh_fdr([R(si).p_sup], fdr_q); end
    [sig_pre(ti), q_pre(ti)] = bh_fdr([R(ti).p_precision], fdr_q);
end
% signed summary; the two originals are always kept
resp_sign = zeros(1,nCell); mod_signed = nan(1,nCell);
for c = ti
    me = R(c).mod_exc_z; ms = R(c).mod_sup_z;
    if isnan(ms) || me >= ms, mod_signed(c) = me; resp_sign(c) = 1;
    else,                     mod_signed(c) = -ms; resp_sign(c) = -1; end
end

fprintf('\n---- results ----\n');
fprintf('  tested                       : %d cells\n', numel(ti));
fprintf('  EXCITATION significant       : %d  (BH q<%.2f)\n', nnz(sig_exc), fdr_q);
fprintf('  SUPPRESSION tested / signif  : %d / %d  (rest lacked the rate for the statistic to work)\n', ...
        nnz([R(ti).sup_tested]), nnz(sig_sup));
fprintf('  PRECISION significant        : %d\n', nnz(sig_pre));
if any(sig_exc)
    fprintf('  M_exc over significant : median %.2f  max %.2f\n', ...
            median([R(sig_exc).mod_exc_z]), max([R(sig_exc).mod_exc_z]));
    fprintf('  preferred latency      : median %.0f ms (IQR %.0f-%.0f)\n', ...
            1000*median([R(sig_exc).preferred_latency_s]), ...
            1000*prctile([R(sig_exc).preferred_latency_s],25), ...
            1000*prctile([R(sig_exc).preferred_latency_s],75));
    fprintf('  recruitment            : median %.3f | precision %.3f | latency MAD %.0f ms\n', ...
            median([R(sig_exc).recruitment]), median([R(sig_exc).precision_fraction]), ...
            1000*median([R(sig_exc).latency_mad_s]));
end
if numel(ti) > 3
    me = [R(ti).mod_exc_z]'; pz = [R(ti).precision_z]';
    ok = isfinite(me) & isfinite(pz);
    if nnz(ok) > 3
        rr = corr(me(ok), pz(ok));
        fprintf('  corr(M_exc, precision_z) = %.2f  %s\n', rr, ...
            ternary(abs(rr)>0.8, '-> nearly one axis; prefer the recruitment-vs-MAD plot', ...
                                 '-> genuinely two axes, the M_exc vs P_z scatter is informative'));
    end
end

%% ---- 7. save ----
if doSave
    rep = find([R.n_events] > 0);
    T = table([R(rep).cell_id]', [R(rep).n_obs]', [R(rep).n_events]', [R(rep).duration_s]', ...
              [R(rep).n_accepted_breaths]', [R(rep).rate_per_min]', [R(rep).tested]', ...
              [R(rep).mean_peth_rate]', [R(rep).peak_rate]', [R(rep).peak_over_mean]', ...
              1000*[R(rep).preferred_latency_s]', ...
              [R(rep).T_exc]', [R(rep).mod_exc_z]', [R(rep).p_exc]', q_exc(rep)', sig_exc(rep)', ...
              [R(rep).trough_rate]', 1000*[R(rep).suppression_latency_s]', ...
              [R(rep).T_sup]', [R(rep).mod_sup_z]', [R(rep).p_sup]', q_sup(rep)', sig_sup(rep)', ...
              [R(rep).sup_tested]', ...
              [R(rep).recruitment]', [R(rep).n_breaths_recruited]', ...
              [R(rep).precision_fraction]', [R(rep).precision_at_pref]', ...
              [R(rep).precision_z]', [R(rep).p_precision]', q_pre(rep)', sig_pre(rep)', ...
              1000*[R(rep).latency_median_s]', 1000*[R(rep).latency_mad_s]', ...
              1000*[R(rep).pref_latency_nearest_s]', ...
              [R(rep).n_assigned_events]', mod_signed(rep)', resp_sign(rep)', ...
              arrayfun(@(s) strjoin(cellstr(s.rec_names),'|'), R(rep),'UniformOutput',false)', ...
              [R(rep).why_not_tested]', ...
        'VariableNames', {'cell_id','n_obs','n_events','duration_s','n_accepted_breaths', ...
        'event_rate_per_min','tested','mean_peth_rate','peak_rate','peak_over_mean', ...
        'preferred_latency_ms','T_exc','mod_exc_z','p_exc','q_exc','sig_exc', ...
        'trough_rate','suppression_latency_ms','T_sup','mod_sup_z','p_sup','q_sup','sig_sup','sup_tested', ...
        'recruitment','n_breaths_recruited','precision_fraction','precision_at_pref', ...
        'precision_z','p_precision','q_precision','sig_precision', ...
        'latency_median_ms','latency_mad_ms','pref_latency_nearest_ms','n_assigned_events', ...
        'mod_signed_z','response_sign','rec_names','why_not_tested'});
    writetable(T, fullfile(outDir,'breath_time_peth_cells.csv'));

    params = struct('TRIGGER',TRIGGER,'winPreFrac',winPreFrac,'winPostFrac',winPostFrac,'searchFrac',searchFrac, ...
        'binWidth_s',binWidth_s,'smoothWidth_s',smoothWidth_s,'smoothBins',smoothBins, ...
        'ampFrac',ampFrac,'minRate_perMin',minRate_perMin,'minEvents',minEvents, ...
        'relWin_s',relWin_s,'minExpectedPerSmoothWin',minExpectedPerSmoothWin, ...
        'nShuffle',nShuffle,'fdr_q',fdr_q,'rngSeed',rngSeed,'medCycle_s',medCyc, ...
        'p_floor',pFloor); %#ok<NASGU>
    save(fullfile(outDir,'breath_time_peth_data.mat'), 'R','ctrs','edges', ...
         'q_exc','sig_exc','q_sup','sig_sup','q_pre','sig_pre','mod_signed','resp_sign', ...
         'params','cfg','-v7.3');
    fprintf('\nSaved breath_time_peth_data.mat + _cells.csv to\n  %s\n', outDir);
    fprintf('Next: breath_time_peth_percell_260727.m, breath_time_modulation_scatter_260727.m\n');
end
end

%% ========================= helpers =========================
function R = init_results(n)
f = {'cell_id','n_obs','n_events','duration_s','n_accepted_breaths','rate_per_min', ...
     'tested','why_not_tested','peth','expSec','peth_shuffle_mean','peth_shuffle_sd', ...
     'mean_peth_rate','peak_rate','trough_rate','preferred_latency_s','suppression_latency_s', ...
     'peak_over_mean','T_exc','T_sup','mod_exc_z','p_exc','mod_sup_z','p_sup', ...
     'sup_tested','sup_expected_per_window', ...
     'recruitment','n_breaths_recruited','precision_fraction','precision_at_pref', ...
     'precision_tau_s','precision_z','p_precision','latency_median_s','latency_mad_s', ...
     'pref_latency_nearest_s', ...
     'n_assigned_events','dt_obs','null_exc','null_sup','null_pre','rec_names'};
proto = struct();
for k = 1:numel(f), proto.(f{k}) = NaN; end
proto.tested = false; proto.sup_tested = false;
proto.why_not_tested = "";  proto.rec_names = strings(0,1);
proto.peth = []; proto.expSec = []; proto.peth_shuffle_mean = []; proto.peth_shuffle_sd = [];
proto.null_exc = single([]); proto.null_sup = single([]); proto.null_pre = single([]);
proto.dt_obs = single([]);
proto.n_obs = 0; proto.n_events = 0; proto.duration_s = 0; proto.n_accepted_breaths = 0;
proto.n_breaths_recruited = 0; proto.n_assigned_events = 0;
R = repmat(proto, 1, n);
end

function cnt = peth_counts(OBs, nB, d)
% Sum the correlogram over PETH bins. A circular shift of the event train by d
% frames is EXACTLY a circular index shift of the precomputed correlogram, so no
% recounting is ever done. d is scalar 0 (observed) or one shift per observation.
cnt = zeros(1,nB);
if isscalar(d), d = repmat(d,1,numel(OBs)); end
for j = 1:numel(OBs)
    o = OBs(j);
    idx = mod(o.lagIdx - 1 + d(j), o.T) + 1;
    cnt = cnt + accumarray(o.binOf, o.ccf(idx), [nB 1])';
end
end

function [Texc, Tsup, tPref, tSup, rbar, pk, tr] = peth_stats(rate, ctrs, srchM, smoothBins)
% Deviation of the lightly smoothed PETH from its OWN mean, maximised over the
% search window -- separately upward and downward. The same function runs on the
% observed and on every shuffle, so the maximisation is paid for by the null.
rbar = mean(rate,'omitnan');
rs = rate; if smoothBins > 1, rs = movmean(rate, smoothBins); end
d = rs - rbar;  d(~srchM) = -Inf;  [Texc, ie] = max(d);
u = rbar - rs;  u(~srchM) = -Inf;  [Tsup, is] = max(u);
if nargout > 2
    tPref = ctrs(ie);  tSup = ctrs(is);
    pk = rs(ie);       tr = rs(is);
end
end

function [dt, key] = event_latencies(OBs, d, maxAbsLag_s)
% Every event assigned to its NEAREST accepted trigger, giving a SIGNED latency.
%
% WHY NEAREST, NOT PRECEDING. With preceding-assignment an event 80 ms BEFORE an
% onset is attributed to the previous breath at +2.05 s, which is both wrong and
% destroys the timing precision of any pre-inspiratory cell -- exactly the cells a
% breathing circuit cares about. Nearest-assignment calls it -80 ms, which is what
% it is. Each event still belongs to exactly ONE breath, so nothing is double
% counted.
%
% Signed dt is therefore bounded by roughly +/- half the local inter-trigger
% interval: an event further away than that necessarily has a nearer trigger.
% `key` uniquely identifies the (observation, breath) an event was assigned to, so
% recruitment can count DISTINCT breaths without re-deriving any of this.
if isscalar(d), d = repmat(d,1,numel(OBs)); end
dt = []; key = [];
for j = 1:numel(OBs)
    o = OBs(j);
    if isempty(o.evFr), continue; end
    e = mod(o.evFr - d(j) - 1, o.T) + 1;             % shifted event frames
    f = o.trig(:);  nf = numel(f);
    b = discretize(e, [f; o.T+1]);                   % index of the PRECEDING trigger
    b(isnan(b)) = 1;                                 % before the first trigger -> compare to it
    ip = b;  in = min(b+1, nf);                      % preceding and following candidates
    dp = e - f(ip);                                  % >= 0 (except before the first trigger)
    dn = e - f(in);                                  % <= 0
    useNext = abs(dn) < abs(dp);
    v = dp;  v(useNext) = dn(useNext);
    bi = ip; bi(useNext) = in(useNext);
    v = v / o.fps;
    keep = abs(v) <= maxAbsLag_s;
    dt  = [dt;  v(keep)]; %#ok<AGROW>
    key = [key; j*1e7 + bi(keep)]; %#ok<AGROW>
end
end

function [frac, tau] = max_window_frac(dt, w, lo, hi)
% Largest fraction of latencies falling inside ANY window of half-width w whose
% centre lies in [lo,hi]. Used identically for the observed and for every shuffle,
% so the observed is not penalised for the null being allowed to optimise.
frac = NaN; tau = NaN;
if isempty(dt), return; end
x = sort(dt(:));  n = numel(x);
best = 0; bt = NaN;
for i = 1:n
    c = x(i) + w;                                    % window [x(i), x(i)+2w], centre c
    if c < lo || c > hi, continue; end
    k = sum(x >= x(i) & x <= x(i)+2*w);
    if k > best, best = k; bt = c; end
end
if best == 0                                          % no admissible centre: fall back
    frac = 0; tau = NaN; return;
end
frac = best/n; tau = bt;
end

function [recr, nRec] = recruitment_of(OBs, d, tPref, w, maxAbsLag_s)
% Accepted breaths with an event within +/-w of the preferred latency, over ALL
% accepted breaths -- including breaths in recordings where the cell was silent
% (those contribute triggers with no events, which is the point of keeping
% zero-event observations). Uses the same NEAREST-trigger assignment as the
% latencies, so recruitment and precision can never disagree about which breath an
% event belonged to.
[dt, key] = event_latencies(OBs, d, maxAbsLag_s);
nTot = sum(arrayfun(@(o) numel(o.trig), OBs));
near = abs(dt - tPref) <= w;
nRec = numel(unique(key(near)));                 % DISTINCT breaths, not events
recr = nRec / max(nTot,1);
end

function z = zscore_against(obsVal, nullVec)
nullVec = nullVec(isfinite(nullVec));
if isempty(nullVec) || ~isfinite(obsVal), z = NaN; return; end
s = std(nullVec);
if s <= 0, z = NaN; return; end
z = (obsVal - mean(nullVec))/s;
end

function [sig, q] = bh_fdr(p, level)
p = p(:); n = numel(p);
[ps, ord] = sort(p);
qs = min(1, ps .* n ./ (1:n)');
for i = n-1:-1:1, qs(i) = min(qs(i), qs(i+1)); end
q = nan(n,1); q(ord) = qs;
sig = reshape(q <= level, 1, []); q = reshape(q, 1, []);
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
