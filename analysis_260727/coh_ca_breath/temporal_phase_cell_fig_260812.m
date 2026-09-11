function [fig, proj, stats] = temporal_phase_cell_fig_260812(OBS, P)
%TEMPORAL_PHASE_CELL_FIG_260812  Per-CELL temporal-phase summary figure.
%
%   [fig, proj, stats] = temporal_phase_cell_fig_260812(OBS, P)
%
% OBS : struct array, ONE ENTRY PER RECORDING the cell appears in, with fields
%         folder   recording folder (needs *_ch1_dFF.mat + breath_peak_pc1.mat)
%         roi      ROI index within THAT recording's dFF matrix
%         recName  label for the title
%         group    genotype/ventral group ('' = sniff it from the path)
%         recDate  MMDD ('' outside the archive layout)
%       Pass a 1x1 OBS for a single recording -- pooling one thing is a no-op, so
%       the GUI and the batch run through exactly this code.
%
% THE CELL IS THE UNIT. Everything that can be pooled across the cell's
% recordings is pooled: the two spike histograms, the three triggered averages,
% the per-cycle heatmap and the Rayleigh log Z. Events and breath cycles from
% every recording are concatenated, each recording keeping ITS OWN occupancy
% weights and its own trigger set, exactly as the ventral analysis pools phase
% per cell. n_eff therefore accumulates across recordings, and the histogram null
% shifts each recording's spike train independently before re-pooling.
%
% Two things cannot be pooled and come from ONE recording, named in the title:
% the wide dF/F trace and the avg-projection crop. P.traceObs picks which (1 by
% default); the trace window is P.trace_xlim_sp, or 30 s from the middle of that
% recording when that is empty.
%
% NO ca-lag anywhere. Event frames are used raw.
%
% P : parameter struct -- see the USER-EDITABLE block of
%     temporal_phase_cell_gui_260812.m. P.statsOnly returns `stats` without
%     drawing, so a batch can select on the very number the figure would print.
%
% Runqi Zhang / 2026-08-12

lightDff = 0.30*P.dffColor + 0.70;
dffColor = P.dffColor; onsetCol = P.onsetCol; peakCol = P.peakCol;
ylim_dff = []; if isfield(P,'ylim_dff'), ylim_dff = P.ylim_dff; end
ylim_epc = []; if isfield(P,'ylim_epc'), ylim_epc = P.ylim_epc; end
statsOnly = isfield(P,'statsOnly') && P.statsOnly;
fig = []; stats = struct(); proj = struct('have_proj',false,'crop_img',[], ...
                                          'bnd_crop',{{}},'barLen_pr',NaN,'roi',NaN);

%% ===================== LOAD EVERY RECORDING =====================
R = struct('fps',{},'T',{},'dff',{},'bw',{},'peak',{},'foot',{},'spk',{}, ...
           'phi',{},'IBI',{},'f_pk',{},'px_um',{},'sm',{},'recName',{},'folder',{},'roi',{});
for i = 1:numel(OBS)
    R(end+1) = load_obs_local(OBS(i), P); %#ok<AGROW>
end
assert(~isempty(R), 'no usable recording for this cell');

nRec = numel(R);
% Pooled IBI: the median of every inter-onset interval from every recording, not
% the mean of per-recording medians, so a long recording carries more weight.
allIBI = [];
for i = 1:nRec, allIBI = [allIBI; diff(R(i).foot)/R(i).fps]; end %#ok<AGROW>
IBI = median(allIBI);
if ~isfinite(IBI) || IBI <= 0, IBI = median([R.IBI]); end
fpsRef = median([R.fps]);
f_pk   = median([R.f_pk]);

%% ===================== TRIGGER WINDOW =====================
if isfield(P,'trigWin_sec') && ~isempty(P.trigWin_sec) && P.trigWin_sec > 0
    trigWinTot = P.trigWin_sec;
else
    nIBI = 2; if isfield(P,'trigWinIBI') && P.trigWinIBI > 0, nIBI = P.trigWinIBI; end
    trigWinTot = nIBI * IBI;
end
winSec = trigWinTot/2;                                   % half-window, SECONDS
% Common time axis. Recordings may differ in fps, so every epoch is extracted in
% its own frames and then interpolated onto this one grid before pooling.
nHalf  = max(1, round(winSec*fpsRef));
tau    = (-nHalf:nHalf)/fpsRef;

% median inspiration-onset-to-peak interval, pooled, for the breath-strip marks
% P.breathObs restricts the BREATH STRIP -- and only the strip -- to a subset of
% the cell's recordings. [] = all of them, which is the default and what every
% figure to date used.
%
% WHY THIS EXISTS. The strip is a mean breath waveform, so pooling recordings
% that breathe at different rates averages one animal's cycle against another's
% and the strip washes out either side of t=0. Nothing else in the figure has
% that problem: the dF/F averages, the heatmap, the histograms and the Rayleigh
% all pool correctly because each recording contributes its own triggers and its
% own occupancy weights. So this narrows the strip alone and leaves the cell
% POOLED everywhere else -- it is a display fix for one panel, not a change to
% which data the cell is made of.
%
% The onset mark moves with it: dOnPk is the median onset-to-peak interval, and a
% mark derived from both recordings sitting on a waveform drawn from one of them
% would point at the wrong place on that waveform.
breathObs = [];
if isfield(P,'breathObs') && ~isempty(P.breathObs)
    breathObs = P.breathObs(:)';
    breathObs = breathObs(breathObs >= 1 & breathObs <= nRec);
end
useB = true(1,nRec);
if ~isempty(breathObs), useB = false(1,nRec); useB(breathObs) = true; end

dOnPk = NaN; dAll = [];
for i = 1:nRec
    if ~useB(i), continue; end
    if isempty(R(i).peak) || isempty(R(i).foot), continue; end
    b = discretize(R(i).peak, [R(i).foot(:); inf]); ok = ~isnan(b);
    if any(ok), dAll = [dAll; (R(i).peak(ok) - R(i).foot(b(ok)))/R(i).fps]; end %#ok<AGROW>
end
if ~isempty(dAll), dOnPk = median(dAll); end

%% ===================== POOLED TRIGGERED EPOCHS =====================
E = []; Ep = []; Epb = []; Esp = []; dtNear = []; nOnAll = 0; nPkAll = 0;
% onsRow: THIS cycle's own inspiration onset, in seconds before its peak, one
% entry per row of Ep. Negative because t = 0 is the peak. Distinct from dOnPk,
% which is the median over all cycles and is a single line on the strip -- these
% are per-cycle and are what the heatmap overlay draws. Built from the same
% boundary-filtered pk vector, in the same order, so it indexes Ep row-for-row.
onsRow = [];
Ez = []; Epz = [];
for i = 1:nRec
    r = R(i);  w = max(1, round(winSec*r.fps));
    on = r.foot(r.foot-w>=1 & r.foot+w<=r.T);
    pk = r.peak(r.peak-w>=1 & r.peak+w<=r.T);
    sp = find(r.spk>0);  sp = sp(sp-w>=1 & sp+w<=r.T);
    tk = (-w:w)/r.fps;
    E   = [E;   regrid_local(epochs_local(r.dff, on, w), tk, tau)];   %#ok<AGROW>
    Ep  = [Ep;  regrid_local(epochs_local(r.dff, pk, w), tk, tau)];   %#ok<AGROW>
    % Same epochs from the Z-SCORED trace. z-scoring is done WITHIN each
    % recording, before any pooling: dF/F amplitude varies by an order of
    % magnitude between cells, so a population mean of raw dF/F is set by
    % whichever cells are brightest rather than by how they are modulated.
    dz = (r.dff - mean(r.dff,'omitnan')) / max(std(r.dff,'omitnan'), eps);
    Ez  = [Ez;  regrid_local(epochs_local(dz, on, w), tk, tau)];      %#ok<AGROW>
    Epz = [Epz; regrid_local(epochs_local(dz, pk, w), tk, tau)];      %#ok<AGROW>
    % breath is in arbitrary per-video units, so z-score each recording before
    % pooling or the loudest video sets the shape of the average
    bwz = (r.bw - mean(r.bw)) / max(std(r.bw), eps);
    if useB(i)   % P.breathObs -- the strip may be drawn from a subset of recordings
        Epb = [Epb; regrid_local(epochs_local(bwz, pk, w), tk, tau)];  %#ok<AGROW>
    end
    Esp = [Esp; regrid_local(epochs_local(r.dff, sp, w), tk, tau)];   %#ok<AGROW>
    if numel(pk) > 1
        d = nan(numel(pk),1);
        for k = 1:numel(pk)
            o = pk; o(k) = [];
            [~,mi] = min(abs(o-pk(k))); d(k) = (o(mi)-pk(k))/r.fps;
        end
        dtNear = [dtNear; d]; %#ok<AGROW>
    elseif ~isempty(pk)
        dtNear = [dtNear; NaN]; %#ok<AGROW>
    end
    % Per-cycle onset, paired the same way the median on the strip is: each peak
    % takes the onset that PRECEDES it. A peak with no preceding onset (before the
    % first detected foot) gets NaN rather than being pushed onto a later onset,
    % which would draw a line at a time inspiration did not happen.
    if ~isempty(pk)
        oo = nan(numel(pk),1);
        bb = discretize(pk, [r.foot(:); inf]);
        okp = ~isnan(bb);
        oo(okp) = -(pk(okp) - r.foot(bb(okp)))/r.fps;
        onsRow = [onsRow; oo]; %#ok<AGROW>
    end
    nOnAll = nOnAll + numel(on);  nPkAll = nPkAll + numel(pk);
end
mu = mean(E,1); sd = std(E,0,1);
if ~isempty(Ep)
    mu_p = mean(Ep,1); sd_p = std(Ep,0,1); mu_pb = mean(Epb,1);
    switch lower(P.sortMode)
        case 'postmean', [~,si] = sort(mean(Ep(:,tau>=0),2),'descend');
        case 'dt',       [~,si] = sort(dtNear,'ascend','MissingPlacement','last');
        otherwise,       si = (1:size(Ep,1))';
    end
    Es = Ep(si,:); cl = prctile(Es(:), [0.5 99.5]);
    % Reorder the per-cycle onsets with the rows, or the overlay would draw each
    % cycle's onset against a different cycle's data the moment sortMode is not
    % 'none'.
    if numel(onsRow) == size(Ep,1), onsS = onsRow(si); else, onsS = []; end
else
    mu_p = []; sd_p = []; mu_pb = []; Es = zeros(0,numel(tau)); cl = [0 1]; onsS = [];
end
if ~isempty(Esp), mu_sp = mean(Esp,1); sd_sp = std(Esp,0,1); else, mu_sp = []; sd_sp = []; end

%% ===================== POOLED RAYLEIGH (all cycles, no lag) =====================
nBR = 36; if isfield(P,'rayPhaseBins'), nBR = P.rayPhaseBins; end
% Production uses NO lag (0). The knob exists only so this statistic can be
% reproduced against an archive that was computed WITH a lag -- validation, not
% analysis. It is deliberately separate from P.ca_lag_sec.
rayLag = 0; if isfield(P,'rayLagSec') && ~isempty(P.rayLagSec), rayLag = P.rayLagSec; end
aPool = []; wPool = []; nCycles = 0; nSpkTot = 0;
% Phase-binned profiles, pooled as SUMS so the pooling is frame-weighted: a long
% recording contributes proportionally more, which is what you want when the same
% cell was imaged for 3000 frames once and 6000 another time. These are the
% cluster features -- in PHASE, not seconds, because IBI spans 0.47-2.53 s across
% sessions and a time-domain waveform would mostly encode breath rate.
edgesP  = linspace(0, 2*pi, nBR+1);
occSum  = zeros(nBR,1);      % frames per phase bin, all recordings
dffSum  = zeros(nBR,1);      % sum of dF/F in each phase bin
evSum   = zeros(nBR,1);      % event count in each phase bin
for i = 1:nRec
    r = R(i);
    phiW  = mod(r.phi(:), 2*pi);
    valid = ~isnan(phiW);
    bidx = nan(r.T,1);
    if nnz(valid) >= 10
        bidx(valid) = min(discretize(phiW(valid), edgesP), nBR);
        occ = accumarray(bidx(valid), 1, [nBR 1]);        % FRAMES per bin, THIS rec
        ev = find(r.spk>0);
        if rayLag ~= 0
            ev = ev - round(rayLag*r.fps);  ev = ev(ev>=1 & ev<=r.T);
        end
        ev = ev(~isnan(bidx(ev)));
        % each recording contributes its events with ITS OWN occupancy weights
        aPool = [aPool; phiW(ev)];                        %#ok<AGROW>
        wPool = [wPool; 1./max(occ(bidx(ev)),1)];         %#ok<AGROW>

        occSum = occSum + occ / r.fps * r.fps;            % frames
        dffSum = dffSum + accumarray(bidx(valid), r.dff(valid), [nBR 1]);
        if ~isempty(ev), evSum = evSum + accumarray(bidx(ev), 1, [nBR 1]); end
    end
    pv = r.phi(~isnan(r.phi));
    if ~isempty(pv), nCycles = nCycles + (max(pv)-min(pv))/(2*pi); end
    nSpkTot = nSpkTot + nnz(r.spk>0);
end
%% ---- PEAK-REFERENCED phase profiles ----
% A SECOND phase convention, independent of the landmark one above. Here the
% cycle is defined peak-to-peak and the inspiratory PEAK sits at pi; inspiration
% onset is not used at all. That matters because onset is the softer landmark --
% it is the foot of a slow rise, so its detection is noisier than the peak's --
% and any analysis keyed on it inherits that jitter. Peak-to-peak also makes the
% phase a pure cycle coordinate rather than a two-landmark piecewise warp.
%
% dF/F is averaged FRAME-WEIGHTED per bin. A mean needs no occupancy divisor (it
% already divides by the frame count); only event RATES do.
%
% Two versions are returned: raw dF/F, and dF/F z-scored WITHIN each recording
% before pooling. The z-scored one is what makes a population average meaningful
% across cells whose absolute dF/F differs by an order of magnitude -- otherwise
% one bright cell sets the shape of the mean.
% EACH CYCLE IS RESAMPLED ONTO THE PHASE GRID, not binned into it. Binning frames
% breaks for fast breathers: at IBI 0.47 s and 30 Hz a cycle is ~14 frames, so 36
% phase bins are 2.5x oversampled -- most bins get 0 or 1 frame and a single noisy
% sample becomes a whole "bin mean". That produced a -0.19 dF/F spike in the
% Vglut2 population average from cells whose cycles are 14 frames long. Slow
% breathers (IBI 2.6 s, ~78 frames) were unaffected, which is exactly why it
% looked like one genotype misbehaving rather than a method problem.
% Interpolating within each cycle gives every cycle the same weight and leaves no
% empty bins, whatever the frame rate or the breathing rate.
gridPk   = (edgesP(1:end-1) + edgesP(2:end))/2;    % bin centres, 0..2pi
dffPkAcc = zeros(1,nBR); dffPkZAcc = zeros(1,nBR); nCycPk = 0;
onsPhPk  = [];    % where inspiration ONSET falls in this peak-referenced frame
for i = 1:nRec
    r = R(i);
    pk = r.peak(:);
    if numel(pk) < 2, continue; end
    z = (r.dff - mean(r.dff,'omitnan')) / max(std(r.dff,'omitnan'), eps);
    ft = r.foot(:);
    for k = 1:numel(pk)-1
        a = pk(k); b = pk(k+1);
        if b <= a || a < 1 || b > r.T, continue; end
        idx = (a:b).';
        % phase of each frame in this cycle, peak at pi
        u = mod( 2*pi*(idx - a)/(b - a) + pi, 2*pi );
        [us, o] = sort(u);
        if numel(us) < 4, continue; end
        [us, iu] = unique(us);
        dffPkAcc  = dffPkAcc  + interp1(us, r.dff(idx(o(iu))), gridPk, 'linear', 'extrap');
        dffPkZAcc = dffPkZAcc + interp1(us, z(idx(o(iu))),     gridPk, 'linear', 'extrap');
        nCycPk = nCycPk + 1;
        % onsets falling inside this cycle, on the same phase axis. This is a
        % property of the BREATH (the duty cycle), not of the cell -- it says
        % where inspiration begins relative to the peak-to-peak cycle.
        fi = ft(ft >= a & ft < b);
        if ~isempty(fi)
            onsPhPk = [onsPhPk; mod(2*pi*(fi - a)/(b - a) + pi, 2*pi)]; %#ok<AGROW>
        end
    end
end
% circular MEDIAN of the onset phase: the direction minimising total angular
% distance, on a 0.5 deg grid. Median not mean, because a handful of mis-detected
% onsets land half a cycle away and would drag a mean with them.
onsetPhasePk = NaN;
if numel(onsPhPk) >= 2
    gg = linspace(-pi, pi, 721); gg(end) = [];
    dd = abs(angle(exp(1i*(onsPhPk(:) - gg))));
    [~, kbest] = min(sum(dd, 1));
    onsetPhasePk = mod(gg(kbest), 2*pi);
end
if nCycPk > 0
    dffPhasePk  = dffPkAcc  / nCycPk;
    dffZPhasePk = dffPkZAcc / nCycPk;
else
    dffPhasePk  = nan(1,nBR);  dffZPhasePk = nan(1,nBR);
end
occPkSum = nCycPk;                      % cycles contributing, not frames

% dF/F: frame-weighted MEAN per phase bin (a mean is not occupancy-biased).
% Event RATE: events per SECOND in each bin, which DOES need the occupancy
% divisor -- raw counts would just retrace how long the cycle spends there.
dffPhase  = dffSum ./ max(occSum, eps);
ratePhase = evSum  ./ max(occSum/fpsRef, eps);
dffPhase(occSum == 0)  = NaN;
ratePhase(occSum == 0) = NaN;
logZ = NaN; nEff = 0; Rbar = 0; thPref = NaN;
if ~isempty(aPool)
    [thPref, Rbar, nEff] = wresultant_local(aPool, wPool);
    logZ = log(max(nEff*Rbar^2, eps));
end

% FIXED-WINDOW peak-triggered average, for population features. Deliberately in
% SECONDS on a window that does not depend on breath rate: GCaMP kinetics are a
% property of the indicator, so warping the axis by each cell's IBI would stretch
% the same rise and decay differently for a fast and a slow breather. This is
% separate from the display window (trigWin), which does scale with IBI.
featWin = 1.0; if isfield(P,'featWin_sec') && ~isempty(P.featWin_sec), featWin = P.featWin_sec; end
% FIXED sampling rate for the feature axis -- never this cell's own fps. Building
% tauF on fpsRef made the vector length depend on the recording (30 Hz -> 61
% samples, 42 Hz -> 85), so cells at different rates could not be stacked into
% one matrix and were silently lost by the caller. Only featPeakDff uses this
% axis; every drawn panel uses `tau`, which is unchanged.
featFs = 30; if isfield(P,'featFs') && ~isempty(P.featFs), featFs = P.featFs; end
nF   = max(1, round(featWin*featFs));
tauF = (-nF:nF)/featFs;
% Both triggers, on the same fixed axis, so the two can be compared directly.
EF = []; EO = [];
for i = 1:nRec
    r  = R(i);
    % Cut WIDER than the target, then interpolate inwards, so the native span
    % always brackets tauF and interp1 never extrapolates at the edges.
    wF = max(1, ceil(featWin*r.fps) + 1);
    tk = (-wF:wF)/r.fps;
    pk = r.peak(r.peak-wF>=1 & r.peak+wF<=r.T);
    on = r.foot(r.foot-wF>=1 & r.foot+wF<=r.T);
    if ~isempty(pk), EF = [EF; regrid_local(epochs_local(r.dff, pk, wF), tk, tauF)]; end %#ok<AGROW>
    if ~isempty(on), EO = [EO; regrid_local(epochs_local(r.dff, on, wF), tk, tauF)]; end %#ok<AGROW>
end
if isempty(EF)
    featPeakDff = nan(1,numel(tauF)); nFeatCyc = 0;
else
    featPeakDff = mean(EF,1,'omitnan'); nFeatCyc = size(EF,1);
end
if isempty(EO)
    featOnsetDff = nan(1,numel(tauF)); nFeatCycOn = 0;
else
    featOnsetDff = mean(EO,1,'omitnan'); nFeatCycOn = size(EO,1);
end

% Stage position for the spatial view. Same cell, so any recording's xy will do;
% the FIRST is used. z is deliberately NOT exported.
xy = [NaN NaN];
if isfield(R(1).sm,'motorPosition') && numel(R(1).sm.motorPosition) >= 2
    xy = R(1).sm.motorPosition(1:2);
end
durTot = sum(arrayfun(@(r) r.T/r.fps, R));

%% ===================== POOLED SPIKE HISTOGRAMS (computed) =====================
% Computed BEFORE the statsOnly return so the permutation p is available to a
% caller that never draws -- selecting cells on p requires it, and recomputing it
% in a second place would let the number drift from the one on the figure.
% Test window is 1 IBI; the plot is doubled to 2 IBI so neighbouring cycles show.
winSecH  = IBI;  testHalf = IBI/2;
nbf = 2; if isfield(P,'histBinFrames'), nbf = max(1,round(P.histBinFrames)); end
binW = nbf/fpsRef;
Mb   = floor(winSecH/binW);
ctrsC  = (-Mb:Mb)*binW;
edgesC = ((-Mb-0.5):(Mb+0.5))*binW;
nSh = 0; if isfield(P,'nShuffle'), nSh = P.nShuffle; end
shiftMinCyc = 3; if isfield(P,'shiftMinCyc'), shiftMinCyc = P.shiftMinCyc; end

trigSel = {'foot','peak'};
trigNm  = {'insp onset','insp peak'};
trigCol = {onsetCol, peakCol};
pHistAll = nan(1,2);
HIST = struct('spkH',{},'null_lo',{},'null_hi',{},'null_mu',{},'nTrig',{});

for q = 1:2
    cnt = zeros(1,numel(ctrsC)); nTrig = 0;
    nullH = zeros(nSh, numel(ctrsC)); okNull = nSh > 0;
    for i = 1:nRec
        r = R(i); wH = max(1,round(winSecH*r.fps));
        trg = r.(trigSel{q}); trg = trg(trg-wH>=1 & trg+wH<=r.T);
        if isempty(trg), continue; end
        ev = find(r.spk>0);
        cnt   = cnt + trig_hist_local(ev, trg, wH, edgesC, r.fps);
        nTrig = nTrig + numel(trg);
        if okNull && ~isempty(ev)
            % each recording's train is shifted INDEPENDENTLY, then re-pooled
            minShift = max([wH+1, round(shiftMinCyc/max(r.f_pk,eps)*r.fps), round(r.fps)]);
            if r.T - minShift <= minShift, okNull = false; continue; end
            for s = 1:nSh
                sh = randi([minShift, r.T-minShift]);
                es = mod(ev - 1 + sh, r.T) + 1;
                nullH(s,:) = nullH(s,:) + trig_hist_local(es, trg, wH, edgesC, r.fps);
            end
        end
    end
    spkH = []; null_lo = []; null_hi = []; null_mu = [];
    if nTrig > 0
        spkH = 100*cnt/nTrig;
        if okNull
            nullH = 100*nullH/nTrig;
            null_lo = prctile(nullH,2.5,1); null_hi = prctile(nullH,97.5,1);
            null_mu = mean(nullH,1);
            tm = abs(ctrsC) <= testHalf;
            if any(tm)
                sObs  = max(abs(spkH(tm) - null_mu(tm)));
                sNull = max(abs(nullH(:,tm) - null_mu(tm)), [], 2);
                pHistAll(q) = (1 + nnz(sNull >= sObs))/(1 + numel(sNull));
            end
        end
    end
    HIST(q) = struct('spkH',spkH,'null_lo',null_lo,'null_hi',null_hi, ...
                     'null_mu',null_mu,'nTrig',nTrig);
end

stats = struct('logZ',logZ, 'nEff',nEff, 'Rbar',Rbar, 'nEvents',numel(aPool), ...
               'nSpikes',nSpkTot, 'nCycles',nCycles, 'IBI',IBI, 'fps',fpsRef, ...
               'f_pk',f_pk, 'nRec',nRec, 'prefPhase',thPref, ...
               'phaseCtrs',(edgesP(1:end-1)+edgesP(2:end))/2, ...
               'dffPhase',dffPhase(:).', 'ratePhase',ratePhase(:).', ...
               'occPhase',occSum(:).', 'xy',xy(:).', 'durTot',durTot, ...
               'rateHz',nSpkTot/max(durTot,eps), ...
               'featTau',tauF, 'featPeakDff',featPeakDff, 'nFeatCyc',nFeatCyc, ...
               'featOnsetDff',featOnsetDff, 'nFeatCycOn',nFeatCycOn, ...
               'pOnset',pHistAll(1), 'pPeak',pHistAll(2), 'nShuffleUsed',nSh, ...
               ... % Triggered curves on the cell's OWN +/-1 IBI window, returned so
               ... % callers can pool across cells. tau is in seconds and spans
               ... % +/-trigWinIBI/2 IBI, so tau/IBI is the normalised cycle axis --
               ... % which is the only fair way to average cells whose IBI spans
               ... % 0.43-2.83 s. featTau/featPeakDff are NOT this: they are a fixed
               ... % 30 Hz, +/-1 s window for the population feature matrix.
               'tau',tau, 'muPeak',mu_p, 'muOnset',mu, ...
               'muPeakZ',mean(Epz,1), 'muOnsetZ',mean(Ez,1), ...
               ... % median inspiration onset->peak interval for this cell, in
               ... % SECONDS. On an onset-triggered axis this is where the peak
               ... % falls; it is a property of the breath, not of the cell.
               'dOnPk',dOnPk, ...
               'histCtrs',ctrsC, 'histOnset',HIST(1).spkH, 'histPeak',HIST(2).spkH, ...
               'histTestHalf',testHalf, 'histWinHalf',winSecH, ...
               ... % RAW pooled event phases and their occupancy weights, returned so
               ... % callers can compute the circular median and IQR from the events
               ... % themselves rather than from the 36-bin profile. Binned at 10 deg
               ... % a median is quantised to 10 deg, which is coarse next to the
               ... % between-genotype differences being measured (Vglut2 155 vs Sert
               ... % 212 deg). Each weight is 1/frames-in-that-bin from the event's
               ... % OWN recording, so pooling across recordings is concatenation.
               'evPhase',aPool(:).', 'evWeight',wPool(:).', ...
               ... % peak-to-peak phase profiles, insp PEAK at pi, onset unused
               'dffPhasePk',dffPhasePk(:).', 'dffZPhasePk',dffZPhasePk(:).', ...
               'occPhasePk',occPkSum(:).', 'onsetPhasePk',onsetPhasePk);
if statsOnly, return; end

%% ===================== TRACE / PROJECTION RECORDING =====================
% The trace and the crop must come from ONE recording, and they come from the
% same one so the picture is of the data being shown. Which one: the most
% REPRESENTATIVE, i.e. the crop with the highest mean normalised cross-correlation
% to the cell's other crops. That beats "first" (arbitrary) and beats "sharpest"
% (a high-zoom recording where the cell drifted out of focus, or where the mask
% wandered onto a neighbour, would still win on zoom but lose on NCC).
iT = pick_repr_obs_local(R, P);
if isfield(P,'traceObs') && ~isempty(P.traceObs) && P.traceObs>=1 && P.traceObs<=nRec
    iT = round(P.traceObs);                       % explicit override (the GUI)
end
rT = R(iT);
tT = (0:rT.T-1)'/rT.fps;
if isfield(P,'trace_xlim_sp') && ~isempty(P.trace_xlim_sp)
    wsp = [max(min(P.trace_xlim_sp),tT(1)), min(max(P.trace_xlim_sp),tT(end))];
else
    % 30 s from the middle of that recording
    mid = (tT(1)+tT(end))/2;  halfW = min(15, (tT(end)-tT(1))/2);
    wsp = [mid-halfW, mid+halfW];
end
mwsp = tT>=wsp(1) & tT<=wsp(2);

%% ===================== POOLED SPECTRA =====================
% Every recording contributes a spectrum. mtspectrumc's frequency grid depends on
% record length, and the recordings differ in duration (sometimes in fps), so each
% is interpolated onto one common grid and the LINEAR power is averaged before
% going to dB -- averaging dB would be a geometric mean and would understate the
% peaks. Individual recordings are drawn faintly behind the mean, so the panel
% shows how many went in and how much they disagree.
% The common grid INHERITS the finest native resolution rather than inventing
% one: Chronux's spacing is Fs/nfft (~0.004 Hz for a 6000-frame record), and the
% x-axis is logarithmic, so a coarse regrid would thin out exactly the decade
% where the breath peak lives.
specDeriv = true;
if isfield(P,'dffSpecDeriv') && ~isempty(P.dffSpecDeriv), specDeriv = logical(P.dffSpecDeriv); end
sBw = cell(nRec,1); fBw = cell(nRec,1);
sDd = cell(nRec,1); fDd = cell(nRec,1);
for i = 1:nRec
    r = R(i);
    pB.Fs=r.fps; pB.tapers=[P.TW_spec,2*P.TW_spec-1]; pB.pad=0;
    pB.fpass=[P.fmin,min(P.fmax,r.fps/2)]; pB.err=0;
    [s1,f1] = mtspectrumc(r.bw, pB);            sBw{i}=s1(:); fBw{i}=f1(:);
    % dF/F spectrum on the DERIVATIVE by default: dF/F is a slow, heavily
    % low-passed signal whose 1/f skirt buries the breathing band, and
    % differentiating flattens it. P.dffSpecDeriv = false gives the plain dF/F
    % spectrum instead. Default TRUE so the batch is unchanged.
    if specDeriv, xDd = diff(r.dff)*r.fps; else, xDd = r.dff; end
    [s2,f2] = mtspectrumc(xDd, pB);             sDd{i}=s2(:); fDd{i}=f2(:);
end
dfMin = inf; fHi = inf;
for i = 1:nRec
    if numel(fBw{i}) > 1, dfMin = min(dfMin, median(diff(fBw{i}))); end
    fHi = min(fHi, max(fBw{i}));
end
if ~isfinite(dfMin) || dfMin <= 0, dfMin = 0.005; end
fgrid  = (P.fmin : dfMin : min(fHi, P.fmax)).';
SbwAll = nan(numel(fgrid), nRec);
SddAll = nan(numel(fgrid), nRec);
for i = 1:nRec
    SbwAll(:,i) = interp1(fBw{i}, sBw{i}, fgrid, 'linear', NaN);
    SddAll(:,i) = interp1(fDd{i}, sDd{i}, fgrid, 'linear', NaN);
end
Sbw = mean(SbwAll, 2, 'omitnan');
Sdd = mean(SddAll, 2, 'omitnan');

%% ===================== AVG PROJECTION (trace recording) =====================
px_um = rT.px_um;
proj  = roi_crop_local(rT.folder, rT.roi, px_um, P);

%% ===================== ANNOTATION =====================
genoStr = OBS(iT).group;
if isempty(genoStr), genoStr = detect_genotype_local(rT.folder); end
gPath = detect_genotype_local(rT.folder);
if ~isempty(gPath) && ~isempty(genoStr) && ~strcmpi(gPath, genoStr)
    genoStr = sprintf('%s(%s)', genoStr, gPath);
end
if isempty(genoStr), genoStr = '(genotype ?)'; end
% Genotype + session + acquisition clock time -- "ChAT 260810 19:30". Without it
% the annotation names the genotype but not WHICH recording, and a GUI-rendered
% figure carries no session identity at all (the batch only gets it from the
% filename). The time is the TRACE recording's, like the stage xyz below.
genoTok = genoStr;      % the genotype alone, before anything is appended

% ROW 1 is identity: genotype, which cell, and how many recordings are pooled
% into it. n on the face of the annotation means a pooled cell declares itself
% without anyone having to read the suptitle or count the rec menu.
if isfield(P,'cellId') && ~isempty(P.cellId) && ~isnan(P.cellId)
    genoStr = sprintf('%s cell#%d (n=%d)', genoStr, P.cellId, nRec);
else
    genoStr = sprintf('%s (n=%d)', genoStr, nRec);
end

% ROW 2 is the acquisition date and time of the TRACE recording, on its own line
% rather than trailing the genotype. Left EMPTY when unavailable and then dropped
% from the list entirely -- a blank row would push everything down by one and
% leave a gap that reads as a missing value rather than an absent one.
acqStr = '';
% Default FALSE so the batch annotation is unchanged; the GUI sets it true.
showAcq = false;
if isfield(P,'showAcqTime') && ~isempty(P.showAcqTime), showAcq = logical(P.showAcqTime); end
if showAcq
    acqDT = acq_datetime_local(rT.folder);
    if ~isnat(acqDT)
        acqStr = sprintf('%s %s', datestr(acqDT,'yymmdd'), datestr(acqDT,'HH:MM'));
    elseif isfield(OBS,'recDate') && ~isempty(OBS(iT).recDate)
        acqStr = OBS(iT).recDate;      % MMDD only, no year/time
    end
end

sm = rT.sm;
% Stage xyz is the TRACE recording's (a cell sits at one place in each FOV, and
% the recordings can be at different depths). The duration after it is POOLED:
% total imaging time this cell contributed across every recording.
totDur = sum(arrayfun(@(r) r.T/r.fps, R));
if isfield(sm,'motorPosition') && numel(sm.motorPosition) >= 3
    % STAGE AXES ARE NOT ANATOMICAL AXES ON THE NEW MOUNT.
    % Everything acquired from 2026-07-21 onward sits on a mount rotated 90 deg,
    % where  -stage x = rostral  and  stage y = lateral. Printing raw stage x/y
    % therefore labels those sessions with the axes swapped. The mapping is keyed
    % on SESSION, never on genotype: Vglut2/0224 and Vglut2/1124 predate the
    % cutover and must not be moved.
    %   rotated : [lateral; rostral] = [0 1; -1 0] * [stage x; stage y]
    %   older   : identity
    % motorPosition is the only authority -- filenames disagree in these sessions.
    gRot = o_group_local(OBS(iT), rT.folder);
    isRot = is_rotated_session_local(gRot, OBS(iT).recDate);
    mp = sm.motorPosition(1:3);
    if isRot, latv = mp(2);  rosv = -mp(1);
    else,     latv = mp(1);  rosv =  mp(2); end
    % z: motor z grows NEGATIVE with depth (name z250 -> motor -245), so a
    % POSITIVE motor z is above the pia, which no cell can be. Those are clamped
    % to 0 rather than printed, since "z5" invites reading 5 um of depth when the
    % sign says the opposite. Clamping is DISPLAY ONLY -- mp(3) is untouched and
    % nothing downstream sees a modified value.
    % Note the clamp treats motor 0 as pia. That is an approximation: the stage
    % zero is set by hand at the surface, so a small positive reading is usually
    % that setting error rather than a cell genuinely above the pia. It is the
    % right thing to print either way; it is not a depth measurement.
    zv = min(mp(3), 0);
    % x is printed UNSIGNED: it is distance from the midline, 0 = medial, and the
    % two hemispheres are mirror images of one another. A signed x would make the
    % same anatomical position read as +1400 on one side and -1400 on the other,
    % which is a hemisphere label masquerading as a coordinate. y keeps its sign
    % because anterior and posterior are genuinely different directions along one
    % axis, not a mirror pair.
    xyzStr = sprintf('x%.0f y%.0f z%.0f, %.0fs', abs(latv), rosv, zv, totDur);
else
    xyzStr = sprintf('xyz n/a, %.0fs', totDur);
end
mwStr = '';
if isfield(sm,'laserPower_pct') && isscalar(sm.laserPower_pct) && isfinite(sm.laserPower_pct)
    acqDT = acq_date_local(rT.folder, sm);
    if ~isnat(acqDT)
        try
            mwStr = sprintf(', %.0fmW', laser_power_calibration(sm.laserPower_pct, acqDT));
        catch
            mwStr = sprintf(', %g%% (no mW)', sm.laserPower_pct);
        end
    else
        mwStr = sprintf(', %g%% (date?)', sm.laserPower_pct);
    end
end
% \mum renders as the mu glyph -- the annotation text is drawn with the tex
% interpreter, so this survives into the PDF as a real character.
if isfinite(px_um), resStr = sprintf('%.2f \\mum/px, %gHz%s', px_um, round(rT.fps), mwStr);
else,               resStr = sprintf('? \\mum/px, %gHz%s', round(rT.fps), mwStr); end

%% ===================== FIGURE =====================
% The grid is 3x finer than a panel so the trace row can be 2/3 of a whole
% panel-block: a normal panel is 6 row-units, the trace is 8 (was 12), the breath
% strip is 3 and the heatmap 9. 26 rows x 3 cols, no holes.
%
%   rows  1- 8   trace (spans all three columns)
%   rows  9-11   proj+annotation   PSD              spike-trig avg
%   rows 12-13   (strip floats, 2)  |                |
%   rows 14-26   heatmap (13 units) onset-trig avg   peak-trig avg
%                                   onset hist       peak hist
%
% Column 1 gives the projection 3 units at the TOP; the breath strip takes 2, and
% everything left over goes to the heatmap (13). sync_breath_strip is told both
% spans, so the strip is pinned to 2 grid rows instead of being a fixed fraction
% of a heatmap that keeps growing.
% P.guiLayout = true is the interactive two-column figure: the spike-triggered
% average and BOTH spike histograms are dropped, and column 2 carries the two
% triggered averages -- onset above, peak below.
% Everything the removed panels measured is still computed and returned in stats
% -- only the drawing goes. Default FALSE = the three-column batch layout.
guiLayout = false;
if isfield(P,'guiLayout') && ~isempty(P.guiLayout), guiLayout = logical(P.guiLayout); end

if guiLayout
    %   rows  1- 8   trace (spans both columns)
    %   rows  9-12   proj+annotation   PSD
    %   rows 13-14   (strip floats over the heatmap)
    %   rows 15-26   heatmap (12)      onset-trig avg (rows 15-20)
    %                                  peak-trig avg  (rows 21-26)
    %
    % 2026-08-17: BOTH histograms are gone and column 2 is now the two triggered
    % averages, onset above peak. The histograms are still computed -- their
    % permutation p is what the "onset **, peak **" line in the annotation
    % reports -- only the drawing goes, exactly as with the panels dropped above.
    %
    % The heatmap gave one grid row to the projection: 13 -> 12 units (-7.7%,
    % the nearest whole row to the 10% asked for) and the projection 3 -> 4
    % (+33%). Column 1 still totals 18 rows, so nothing else moves.
    nRows = 26; nCols = 2;   span_pan = [6 1];
    tile_trace = 1;   span_trace  = [8 2];
    tile_proj  = 17;  span_proj   = [4 1];
    tile_psd   = 18;
    tile_pthm  = 29;  span_pthm   = [12 1];
    tile_ptavg  = 30;                                        % onset-triggered
    tile_pktavg = 42;                                        % peak-triggered
    tile_stavg = [];  tile_onhist = [];  tile_pkhist = [];   % dropped
else
    nRows = 26; nCols = 3;   span_pan = [6 1];
    tile_trace = 1;   span_trace  = [8 3];
    tile_proj  = 25;  span_proj   = [3 1];
    tile_psd   = 26;  tile_stavg  = 27;
    tile_ptavg = 44;  tile_pktavg = 45;
    tile_pthm  = 40;  span_pthm   = [13 1];
    tile_onhist= 62;  tile_pkhist = 63;
end
nUnitHm    = span_pthm(1);   nUnitStrip = 2;

% Width is 2/3 of what it was (0.9 -> 0.6). The trace spans the same three
% columns, so narrowing the figure narrows the columns and the trace together and
% they stay aligned by construction.
% Two columns at 2/3 the width keeps each column -- and therefore every bottom
% panel -- the same physical size it had in the three-column layout. Narrowing the
% figure rather than widening the columns is what makes the trace row, which spans
% all columns, shrink to match the panels underneath it.
figW = 0.6;  if guiLayout, figW = 0.6 * 2/3; end
fig = figure('Color','w','Units','normalized','Position',[0.03 0.03 figW 0.8], ...
             'Name',sprintf('cell summary (%d recording(s))', nRec));
tl = tiledlayout(fig,nRows,nCols,'TileSpacing','compact','Padding','compact');
cellTag = '';
if isfield(P,'cellId') && ~isempty(P.cellId) && ~isnan(P.cellId)
    cellTag = sprintf('cell %d   |   ', P.cellId);
end
recTag = rT.recName;
if nRec > 1, recTag = sprintf('%s   (+%d more recording(s), pooled)', rT.recName, nRec-1); end
% The GUI figure is 2/3 as wide, so the title has to come down with it or it is
% clipped at both ends. Shorter separators and a smaller font, not a smaller
% amount of information.
if guiLayout
    ttlLine1 = sprintf('%s%s  ROI%d', cellTag, recTag, rT.roi);
    ttlLine2 = sprintf('%d rec pooled | IBI %.2f s, %gHz | trace %.1f-%.1f s | win %.2f s', ...
                       nRec, IBI, round(fpsRef), wsp(1), wsp(2), 2*winSec);
    ttlSize  = 9;
else
    ttlLine1 = sprintf('%s%s   ROI%d', cellTag, recTag, rT.roi);
    ttlLine2 = sprintf(['%d recording(s) pooled  |  IBI %.2f s, %gHz  |  trace %.1f-%.1f s  ' ...
                        '|  trig win %.2f s (%.1f x IBI)'], nRec, IBI, round(fpsRef), ...
                        wsp(1), wsp(2), 2*winSec, 2*winSec/IBI);
    ttlSize  = get(groot,'defaultAxesFontSize') + 2;
end
title(tl, {ttlLine1, ttlLine2}, 'Interpreter','none','FontWeight','bold', ...
      'FontSize', ttlSize);

% ---- wide trace (one recording) ----
% x is time WITHIN THE SHOWN WINDOW, so a 30 s clip always reads 0-30 whatever
% part of the recording it was cut from. The absolute window is in the figure
% title ("trace 84.5-114.5 s"), so the provenance is not lost.
ax1 = nexttile(tl,tile_trace,span_trace);
tRel = tT - wsp(1);
% Breath colour in the wide trace only. Grey elsewhere (PSD, cycle strip) is a
% "reference signal" convention; here P.traceBreathCol overrides it, and the axis
% colour follows the line so the two cannot disagree.
bwCol = [0.6 0.6 0.6];
if isfield(P,'traceBreathCol') && ~isempty(P.traceBreathCol)
    bwCol = P.traceBreathCol;
end
if isfield(P,'traceBreathCol') && ~isempty(P.traceBreathCol)
    % SINGLE-AXIS path. yyaxis always draws the right-hand side ABOVE the left,
    % whatever order things were created in, so a breath trace on the right can
    % never sit under the dF/F -- uistack cannot cross the two sides either. With
    % a black breath line that matters: it covers the data it is context for.
    % So the breath is rescaled onto the dF/F axis and both live on ONE axis,
    % where z-order is just child order. The right axis is kept for its LABEL
    % only (it already had no ticks), so the panel looks unchanged.
    yyaxis(ax1,'right'); set(ax1,'YColor',bwCol,'YTick',[]); ylabel(ax1,'breath');
    yyaxis(ax1,'left');  hold(ax1,'on');
    plot(ax1, tRel, rT.dff, '-','Color',dffColor,'LineWidth',0.8);
    set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
    if nnz(mwsp)>2, set_ylim_win(ax1, rT.dff, mwsp); end
    dl = ylim(ax1);
    wbw = rT.bw(mwsp); if isempty(wbw), wbw = rT.bw; end
    lo = min(wbw); hi = max(wbw);
    bwS = (rT.bw - lo) / max(hi - lo, eps);          % 0..1 across the window
    bwS = dl(1) + bwS * (dl(2) - dl(1));
    hBw = plot(ax1, tRel, bwS, '-','Color',bwCol,'LineWidth',0.6);
    uistack(hBw,'bottom');
    ylim(ax1, dl);                                    % the breath must not rescale it
    hold(ax1,'off');
else
    yyaxis(ax1,'right'); plot(ax1, tRel, rT.bw, '-','Color',bwCol,'LineWidth',0.6);
    set(ax1,'YColor',bwCol,'YTick',[]); ylabel(ax1,'breath');
    if nnz(mwsp)>2, set_ylim_win(ax1, rT.bw, mwsp); end
    yyaxis(ax1,'left'); plot(ax1, tRel, rT.dff, '-','Color',dffColor,'LineWidth',0.8);
    set(ax1,'YColor','k'); ylabel(ax1,'\DeltaF/F');
    if nnz(mwsp)>2, set_ylim_win(ax1, rT.dff, mwsp); end
end
xlim(ax1,[0 wsp(2)-wsp(1)]); xlabel(ax1,'Time (s)'); box(ax1,'off');

% ---- PSD ----
ax2 = nexttile(tl,tile_psd,span_pan); hold(ax2,'on');
if nRec > 1                                   % individual recordings, faint
    plot(ax2, fgrid, 10*log10(SbwAll), '-', 'Color',[0.5 0.5 0.5 0.25], 'LineWidth',0.5);
    plot(ax2, fgrid, 10*log10(SddAll), '-', 'Color',[dffColor 0.25],    'LineWidth',0.5);
end
plot(ax2, fgrid, 10*log10(Sbw), 'Color',[0.5 0.5 0.5], 'LineWidth',1.2);
plot(ax2, fgrid, 10*log10(Sdd), 'Color',dffColor,      'LineWidth',1.2);
set(ax2,'XScale','log'); xlim(ax2,[P.fmin P.fmax]); xticks(ax2,[0.1 0.3 1 3 10]);
set(ax2,'XMinorTick','off'); xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'power (dB)');
% No title: the axes already say 'power (dB)' against 'Frequency (Hz)', and the
% recording count now lives on row 1 of the annotation.

% ---- triggered averages (pooled) ----
% P.trialOverlay = true draws every single trial faintly instead of a +/- SD
% ribbon. Default FALSE so the batch keeps the ribbon it was rendered with.
trialOv = false;
if isfield(P,'trialOverlay') && ~isempty(P.trialOverlay), trialOv = logical(P.trialOverlay); end
if trialOv, sdTag = ''; else, sdTag = '\pm SD, '; end

% ---- onset-triggered ----
if ~isempty(tile_ptavg)
ax5 = nexttile(tl,tile_ptavg,span_pan); hold(ax5,'on');
if ~isempty(mu)
    if trialOv
        overlay_trials_local(ax5, tau, E, mu, dffColor, ylim_dff);
    else
        fill(ax5,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],lightDff,'EdgeColor','none');
        plot(ax5, tau, mu, '-','Color',dffColor,'LineWidth',1.5);
        xlim(ax5,[tau(1) tau(end)]);
    end
end
xline(ax5,0,'Color',onsetCol,'LineWidth',0.8);
ylabel(ax5,'\DeltaF/F'); xlabel(ax5,'time from insp onset (s)'); grid(ax5,'on');
if ~isempty(ylim_dff), ylim(ax5, ylim_dff); end
title(ax5, 'onset-triggered');
end

% ---- spike-triggered ----
% P.spkWinIBI gives THIS panel its own width, in IBI, independent of the shared
% trigger window: spkWinIBI = 1 means +/-0.5 IBI. The spike-triggered average is
% the only panel not locked to the breath cycle, so a full 2-IBI span mostly shows
% neighbouring cycles rather than the transient. Empty = share the common window,
% which is what the batch does.
tauS = tau; EspS = Esp; muS = mu_sp; sdS = sd_sp;
if isfield(P,'spkWinIBI') && ~isempty(P.spkWinIBI) && P.spkWinIBI > 0
    selS = abs(tau) <= 0.5*P.spkWinIBI*IBI;
    % ~isempty(mu_sp): a cell with NO detected events has Esp = [] and therefore
    % mu_sp = [] (0x0), while selS is 1xnumel(tau). Trimming then indexes column
    % 1 of an empty matrix and throws. The panel is not drawn for such a cell
    % anyway (and P.guiLayout drops it outright), so skip the trim.
    if nnz(selS) > 2 && ~isempty(mu_sp) && size(mu_sp,2) == numel(tau)
        tauS = tau(selS);  muS = mu_sp(:,selS);  sdS = sd_sp(:,selS);
        if ~isempty(Esp), EspS = Esp(:,selS); end
    end
end
if ~isempty(tile_stavg)
ax_sp = nexttile(tl,tile_stavg,span_pan); hold(ax_sp,'on');
if ~isempty(muS)
    if trialOv
        overlay_trials_local(ax_sp, tauS, EspS, muS, dffColor, ylim_dff);
    else
        fill(ax_sp,[tauS fliplr(tauS)],[muS+sdS fliplr(muS-sdS)],lightDff,'EdgeColor','none');
        plot(ax_sp, tauS, muS, 'Color',dffColor,'LineWidth',1.5);
        xlim(ax_sp,[tauS(1) tauS(end)]);
    end
    xline(ax_sp,0,'k--','LineWidth',0.8);
end
xlabel(ax_sp,'time from spike (s)'); ylabel(ax_sp,'\DeltaF/F'); grid(ax_sp,'on');
if ~isempty(ylim_dff), ylim(ax_sp, ylim_dff); end
title(ax_sp, sprintf('spike-triggered average (%sn=%d spikes)', sdTag, size(Esp,1)));
end

% ---- peak-triggered ----
ax_pk = nexttile(tl,tile_pktavg,span_pan); hold(ax_pk,'on');
if ~isempty(mu_p)
    if trialOv
        overlay_trials_local(ax_pk, tau, Ep, mu_p, dffColor, ylim_dff);
    else
        fill(ax_pk,[tau fliplr(tau)],[mu_p+sd_p fliplr(mu_p-sd_p)],lightDff,'EdgeColor','none');
        plot(ax_pk, tau, mu_p, 'Color',dffColor,'LineWidth',1.5);
        xlim(ax_pk,[tau(1) tau(end)]);
    end
end
% ONE trigger line per panel, in that panel's own colour: blue = peak here, red =
% onset on the panel above. The mean-onset landmark that used to be drawn here in
% red as well was there only because the onset-triggered panel had been dropped
% from the GUI layout; that panel is back (2026-08-17), so a second coloured line
% on this axes now says nothing the neighbouring panel does not say better, and
% two trigger colours on one panel invite reading the wrong one as t=0.
xline(ax_pk,0,'Color',peakCol,'LineWidth',0.8);
xlabel(ax_pk,'time from insp peak (s)'); ylabel(ax_pk,'\DeltaF/F'); grid(ax_pk,'on');
if ~isempty(ylim_dff), ylim(ax_pk, ylim_dff); end
title(ax_pk, 'peak-triggered');

% ---- per-cycle heatmap (pooled cycles) ----
ax4 = nexttile(tl,tile_pthm,span_pthm);
if isempty(Es)
    axis(ax4,'off'); text(ax4,0.5,0.5,'(no breath peaks)','Horizontal','center');
else
    imagesc(ax4, tau, 1:size(Es,1), Es); axis(ax4,'tight');
    colormap(ax4, gray(256)); caxis(ax4, cl); set(ax4,'YDir','reverse');
    xlim(ax4,[tau(1) tau(end)]);
    % ---- per-cycle inspiration onset ----
    % One point per row, at THAT breath's own onset, joined into a line so the
    % cycle-to-cycle scatter is visible as the wobble of the line. This is the
    % per-cycle counterpart of the single median mark on the strip above: where
    % the strip says "onset is typically here", this says where it actually was,
    % breath by breath. Rows with no preceding onset are NaN and simply break the
    % line rather than being drawn at a wrong time.
    % Defaults to onsetCol, so this line and every other onset mark on the figure
    % are ONE colour set in ONE place and cannot drift apart. P.onsetRowCol still
    % overrides if the per-cycle line ever needs to differ again. Clipping to the
    % axes is left to MATLAB -- an onset earlier than the window edge should
    % disappear, not be pinned to the edge.
    % P.onsetRowStyle = 'raster' draws one short tick per row INSTEAD of the
    % joined line: same x positions, same colour, but nothing is drawn between
    % consecutive rows. A joined line interpolates across rows, which invents a
    % continuous "onset trajectory" the data does not have -- rows are separate
    % breaths and, under a sort, not even in time order. The raster makes each
    % row's onset its own mark and lets the cloud speak for itself. 'line' is
    % still the default so the batch and the 265 archived figures are untouched.
    if ~isempty(onsS) && any(isfinite(onsS))
        onsCol = onsetCol;
        if isfield(P,'onsetRowCol') && ~isempty(P.onsetRowCol), onsCol = P.onsetRowCol; end
        onsSty = 'line';
        if isfield(P,'onsetRowStyle') && ~isempty(P.onsetRowStyle), onsSty = lower(P.onsetRowStyle); end
        hold(ax4,'on');
        switch onsSty
            case 'raster'
                % Tick height in ROWS (P.onsetRowTickFrac, default 0.5 = half a
                % row, since row i spans i-0.5 .. i+0.5). 1 gives touching ticks,
                % which at 166 rows in a 600 px panel is back to a solid column;
                % anything below ~0.4 goes sub-pixel and starts dropping out.
                % Drawn as ONE NaN-separated line, so thousands of rows stay a
                % single graphics object.
                tickH = 0.5;
                if isfield(P,'onsetRowTickFrac') && ~isempty(P.onsetRowTickFrac)
                    tickH = P.onsetRowTickFrac;
                end
                rr = (1:numel(onsS))';
                ok = isfinite(onsS(:));
                xr = [onsS(ok)'; onsS(ok)'; nan(1,nnz(ok))];
                yr = [rr(ok)'-tickH/2; rr(ok)'+tickH/2; nan(1,nnz(ok))];
                plot(ax4, xr(:), yr(:), '-', 'Color', onsCol, 'LineWidth', 0.8);
            otherwise
                plot(ax4, onsS(:), (1:numel(onsS))', '-', 'Color', onsCol, 'LineWidth', 0.8);
        end
        hold(ax4,'off');
    end
    xlabel(ax4,'time from insp peak (s)');
    % The row order is whatever P.sortMode says, so the axis has to name it.
    % A bare 'breath #' on duration-sorted rows reads as chronological and is
    % simply wrong. Under 'none' the rows ARE chronological, but for a pooled
    % multi-recording cell that means chronological WITHIN each recording, with
    % the recordings stacked in registry order -- there is no common clock
    % across separate acquisitions.
    switch lower(P.sortMode)
        case 'postmean', ylabel(ax4,'breath (sorted by post-peak \DeltaF/F)');
        case 'dt',       ylabel(ax4,'breath (sorted by cycle duration)');
        otherwise,       ylabel(ax4,'breath # (chronological)');
    end
    cb=colorbar(ax4); cb.Label.String='\DeltaF/F';
end

%% ---- pooled spike histograms: DRAW the precomputed HIST ----
tileHst = {tile_onhist, tile_pkhist};
for q = 1:2
    if isempty(tileHst{q}), continue; end     % onset histogram dropped in the GUI layout
    axq = nexttile(tl, tileHst{q}, span_pan);
    H = HIST(q);
    if H.nTrig == 0
        axis(axq,'off');
        text(axq,0.5,0.5,sprintf('(no %s triggers)',trigNm{q}),'Horizontal','center');
        continue;
    end
    hold(axq,'on');
    if ~isempty(H.null_hi)
        fill(axq,[ctrsC fliplr(ctrsC)],[H.null_hi fliplr(H.null_lo)],[0.6 0.6 0.6], ...
             'FaceAlpha',0.30,'EdgeColor','none','HandleVisibility','off');
        plot(axq, ctrsC, H.null_mu, '--','Color',[0.4 0.4 0.4],'LineWidth',0.8);
    end
    bar(axq, ctrsC, H.spkH, 1, 'FaceColor','k','FaceAlpha',0.85,'EdgeColor','none');
    set(axq,'YColor','k');
    if ~isempty(ylim_epc), ylim(axq, ylim_epc);
    else,                  ylim(axq, [0 max([H.spkH, H.null_hi, 1])*1.10]); end
    ylabel(axq,'spk/cyc %');
    xline(axq, 0, 'Color',trigCol{q},'LineWidth',1);
    % Mean inspiration onset, the same red landmark as on the peak-triggered
    % average and the breath strip. Only meaningful on the PEAK-triggered
    % histogram (q=2) -- on the onset-triggered one, onset already IS zero.
    if guiLayout && q == 2 && isfinite(dOnPk) && abs(dOnPk) <= winSecH
        xline(axq, -dOnPk, '-', 'Color', onsetCol, 'LineWidth', 1);
    end
    xline(axq, -testHalf, ':', 'Color',[0.45 0.45 0.45],'LineWidth',0.8);
    xline(axq,  testHalf, ':', 'Color',[0.45 0.45 0.45],'LineWidth',0.8);
    xlim(axq,[-winSecH winSecH]);
    xlabel(axq, sprintf('time from %s (s)', trigNm{q}));
    title(axq, sprintf('%d spk, %.1f cyc, %.2f spk/cyc', nSpkTot, nCycles, ...
          nSpkTot/max(nCycles,eps)), 'Interpreter','none','FontSize',9);
    box(axq,'on');
end

%% ---- annotation + projection panel ----
% P.showRayleigh = false drops the Rayleigh log Z line. logZ is the ONLY
% Rayleigh-derived quantity this figure ever drew -- Rbar, nEff and prefPhase are
% computed and returned in stats but never plotted, and the polar plot and phase
% histogram were removed earlier -- so this one flag removes all of it. logZ keeps
% being computed and returned either way, because cell SELECTION (logZ >= 1/2/3)
% and the population scripts depend on it.
% Default TRUE so the batch is unchanged; the GUI sets it false.
showRay = true;
if isfield(P,'showRayleigh') && ~isempty(P.showRayleigh), showRay = logical(P.showRayleigh); end

% P.boldGenoOnly = true makes the GENOTYPE the only bold text on the whole figure:
% every panel title, axis label and the suptitle drop to normal weight, and the
% genotype token alone is bolded with a TeX \bf group so the date and time beside
% it stay normal. Default FALSE = the batch's existing weighting.
boldGenoOnly = false;
if isfield(P,'boldGenoOnly') && ~isempty(P.boldGenoOnly)
    boldGenoOnly = logical(P.boldGenoOnly);
end
if boldGenoOnly && ~isempty(genoTok)
    % escape nothing: genotype tokens are plain words, at most "IO(Vglut2)"
    genoStr = regexprep(genoStr, ['^' regexptranslate('escape',genoTok)], ...
                        ['{\\bf ' genoTok '}'], 'once');
end

% The onset/peak permutation stars were removed from the annotation 2026-08-17
% (RZ). pHistAll is still computed and still returned in stats.
annotLines = {genoStr};
if ~isempty(acqStr), annotLines{end+1} = acqStr; end
annotLines{end+1} = xyzStr;
annotLines{end+1} = resStr;
if showRay
    if isfinite(logZ), statStr = sprintf('Rayleigh Log(Z)=%.2f', logZ);
    else,              statStr = sprintf('Rayleigh Log(Z) n/a (%d events)', numel(aPool)); end
    annotLines{end+1} = statStr;
end

% The picture and the annotation are TWO separately positioned axes, not one.
% Sharing an axes made the text spacing a fraction of the crop's pixel height, so
% shrinking this tile collapsed the annotation lines on top of each other. Now the
% image is top-aligned in the tile at its natural aspect and the text sits beside
% it with spacing in normalized units, independent of the crop.
axTile = nexttile(tl,tile_proj,span_proj);
axis(axTile,'off'); drawnow;
tp = get(axTile,'OuterPosition');    % the tile box in figure units
delete(axTile);

% Left-align the picture with the HEATMAP'S PLOT BOX, not with the tile edge. The
% tile edge sits left of the heatmap's y-label and tick labels, so aligning to it
% pushes the crop out past the column. The breath strip is positioned from the
% same heatmap Position, so picture, strip and heatmap all share one left edge.
if isgraphics(ax4)
    p4 = get(ax4,'Position');
    if p4(1) > tp(1) && p4(1) < tp(1)+tp(3)
        tp(3) = (tp(1)+tp(3)) - p4(1);
        tp(1) = p4(1);
    end
end

if proj.have_proj, [Hc,Wc] = size(proj.crop_img); else, Hc = 100; Wc = 100; end
% Square pixels: a crop that is Wc/Hc wide in pixels must be (Wc/Hc)*(figH/figW)
% wide in normalized figure units, because normalized x and y span different
% numbers of screen pixels.
gapX = 0.012;
fpx  = getpixelposition(fig);
aspF = fpx(4)/max(fpx(3),eps);
imgH = tp(4);
imgW = imgH * (Wc/Hc) * aspF;
if imgW > 0.42*tp(3)                            % never let it eat the text room
    imgW = 0.42*tp(3);  imgH = imgW / ((Wc/Hc)*aspF);
end
ax_pr = axes(fig,'Position',[tp(1), tp(2)+tp(4)-imgH, imgW, imgH]);
hold(ax_pr,'on');
if proj.have_proj
    imagesc(ax_pr, proj.crop_img); colormap(ax_pr, gray(256)); caxis(ax_pr,[0 1]);
    % P.outlineImgPx (image pixels) scales the outline to the displayed crop, so a
    % 195 px crop and a 400 px crop get the same apparent thickness. Empty = the
    % old fixed-points behaviour, which is what the batch uses.
    % The panel magnification is imgH/Hc, so a SHORT crop is blown up more and a
    % fixed-points outline comes out relatively thinner on it -- that is why the
    % outline looks different from cell to cell. Scaling by image pixels removes
    % that. The max() is a visibility floor: this panel often shows the crop at
    % less than 1:1, where a pure image-pixel rule computes a sub-pixel line.
    % Outline thickness: P.outlineFrac as a FRACTION of the crop width takes
    % priority over P.outlineImgPx (absolute image pixels). With the crop
    % resampled to a fixed pixel size the two are equivalent, but the fraction
    % keeps its proportions if crop_um or cropPxUm change.
    wantPx = [];
    if isfield(P,'outlineFrac') && ~isempty(P.outlineFrac) && P.outlineFrac > 0
        wantPx = P.outlineFrac * size(proj.crop_img,2);
    elseif isfield(P,'outlineImgPx') && ~isempty(P.outlineImgPx)
        wantPx = P.outlineImgPx;
    end
    lwOut = P.outlineLW;
    if ~isempty(wantPx)
        drawnow;   % Position is only meaningful once the axes is laid out
        lwOut = max(P.outlineLW, ...
                    outline_lw_local(ax_pr, size(proj.crop_img,2), wantPx, P.outlineLW));
    end
    % P.outlineAlpha < 1 lets the cell show through the outline. MATLAB takes a
    % 4-element Color as [r g b alpha] on line objects.
    outCol = [1 1 0];
    if isfield(P,'outlineAlpha') && ~isempty(P.outlineAlpha) && P.outlineAlpha < 1
        outCol = [1 1 0 max(0, min(1, P.outlineAlpha))];
    end
    for k = 1:numel(proj.bnd_crop)
        plot(ax_pr, proj.bnd_crop{k}(:,2), proj.bnd_crop{k}(:,1), '-','Color',outCol, ...
             'LineWidth',lwOut);
    end
    sbCorner = 'lower left';
    if isfield(P,'scalebarCorner') && ~isempty(P.scalebarCorner)
        sbCorner = P.scalebarCorner;
    end
    sbThick = 0.02;
    if isfield(P,'scalebarThickFrac') && ~isempty(P.scalebarThickFrac)
        sbThick = P.scalebarThickFrac;
    end
    rectangle(ax_pr,'Position', ...
              scalebar_pos_260815(Wc, Hc, proj.barLen_pr, sbCorner, sbThick), ...
              'FaceColor','w','EdgeColor','none');
    xlim(ax_pr,[0.5 Wc+0.5]); ylim(ax_pr,[0.5 Hc+0.5]);
else
    text(ax_pr, 0.5, 0.5, '(no mask / avg proj)','Horizontal','center','FontSize',8);
end
set(ax_pr,'YDir','reverse'); axis(ax_pr,'off'); hold(ax_pr,'off');

ax_tx = axes(fig,'Position',[tp(1)+imgW+gapX, tp(2), max(tp(3)-imgW-gapX,0.05), tp(4)]);
axis(ax_tx,'off'); xlim(ax_tx,[0 1]); ylim(ax_tx,[0 1]);
nL = numel(annotLines);
yN = linspace(0.92, 0.08, nL);                  % even spacing, top to bottom
for k = 1:nL
    text(ax_tx, 0, yN(k), annotLines{k}, 'FontSize',9, 'Interpreter','tex', ...
         'HorizontalAlignment','left','VerticalAlignment','middle', ...
         'FontWeight', ternary_local(k==1 && ~boldGenoOnly,'bold','normal'));
end

%% ---- cycle-average breath strip, floated on the heatmap ----
% n is the strip's OWN cycle count (size(Epb,1)), not the heatmap's: under
% P.breathObs the two differ, and printing the heatmap's count over a waveform
% built from fewer cycles would misstate what the strip averages. When the strip
% is restricted the title says which recordings it came from.
stripN = size(Epb,1);
if isempty(breathObs) || numel(breathObs) == nRec
    strip_ttl = sprintf('peak-triggered, n=%d cycles', stripN);
else
    strip_ttl = sprintf('peak-triggered, n=%d cycles (rec %s of %d)', ...
                        stripN, strjoin(cellstr(num2str(breathObs(:))).', ','), nRec);
end
if ~isempty(Epb) && ~isempty(mu_pb) && ~isempty(Es)
    drawnow;
    ax_b = axes(fig);
    plot(ax_b, tau, mu_pb, '-','Color',[0.35 0.35 0.35],'LineWidth',1.2);
    xlim(ax_b,[tau(1) tau(end)]);
    padb = 0.08*max(range_local(mu_pb), eps);
    ylim(ax_b,[min(mu_pb)-padb, max(mu_pb)+padb]);
    xline(ax_b, 0, '-', 'Color', peakCol, 'LineWidth', 1);
    if isfinite(dOnPk)
        xline(ax_b, -dOnPk, '-', 'Color', onsetCol, 'LineWidth', 1);
        % The NEXT cycle's onset (at IBI - dOnPk) is drawn only in the batch
        % layout. It marks where the following breath begins, which is useful when
        % the window spans two cycles and clutter otherwise.
        if ~guiLayout && (IBI - dOnPk) <= tau(end)
            xline(ax_b, IBI - dOnPk, '-', 'Color', onsetCol, 'LineWidth', 1);
        end
    end
    set(ax_b,'XTick',[],'YTick',[]); box(ax_b,'off');
    ylabel(ax_b,'breath','FontSize',8);
    title(ax_b, strip_ttl, 'FontSize',9);
    sync_breath_strip(ax_b, ax4, nUnitHm, nUnitStrip);
    set(fig,'SizeChangedFcn', @(s,e) sync_breath_strip(ax_b, ax4, nUnitHm, nUnitStrip));
end

%% ---- font weight sweep: the genotype is the only bold thing ----
% Done LAST, after every panel exists, so nothing added earlier keeps MATLAB's
% default bold title. The genotype stays bold because its weight comes from a TeX
% \bf group inside the string, which FontWeight does not override.
if boldGenoOnly
    set(findall(fig,'Type','text'), 'FontWeight','normal');
    for ax = findall(fig,'Type','axes').'
        set([get(ax,'Title'), get(ax,'XLabel'), get(ax,'YLabel')], 'FontWeight','normal');
    end
    for ax = findall(fig,'Type','colorbar').'
        set(get(ax,'Label'), 'FontWeight','normal');
    end
    if isgraphics(tl) && ~isempty(get(tl,'Title'))
        set(get(tl,'Title'), 'FontWeight','normal');
    end
end
end

% =======================================================================
function r = load_obs_local(o, P)
%LOAD_OBS_LOCAL  One recording's aligned traces, events and breath phase.
df = dir(fullfile(o.folder,'*_ch1_dFF.mat'));
bp = dir(fullfile(o.folder,'breath_peak_pc1.mat'));
ip = dir(fullfile(o.folder,'breath_insp_start_pc1.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', o.folder);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', o.folder);
[fps, sm] = detect_session_fps(o.folder, P.fallback_fps);
D  = load(fullfile(df(1).folder, df(1).name),'dFF');
BP = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(D.dFF);
assert(o.roi>=1 && o.roi<=size(dff_all,2),'ROI %d out of range (1..%d) in %s', ...
       o.roi, size(dff_all,2), o.folder);

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
% Vglut2/1124 used a rising-edge trigger, so breath and events are one frame
% late. GATE ON THE GENOTYPE FOLDER, not on the ventral group label: the IO sites
% of that session live at <root>\Vglut2\1124\IO\ and carry the group name 'IO',
% so testing o.group skipped the fix for exactly those recordings. That was worth
% ~1 frame of phase and moved 5 cells across a logZ cut.
gFix = o.group;
if isempty(gFix) || strcmpi(gFix,'IO'), gFix = detect_genotype_local(o.folder); end
if strcmpi(gFix,'Vglut2') && strcmp(o.recDate,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end

spk = [];
sp_file = fullfile(o.folder,'ca_spike_data.mat');
if isfile(sp_file)
    CA = load(sp_file,'roi_spikes');
    if isfield(CA,'roi_spikes') && o.roi <= numel(CA.roi_spikes)
        spk = double(CA.roi_spikes(o.roi).spike_train(:));
    end
end
T = min([size(dff_all,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
dff = dff_all(1:T, o.roi); bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); else, ev_foot = zeros(T,1); end
if isempty(spk), spk = zeros(T,1); end
if numel(spk) < T, spk(end+1:T) = 0; end
spk = spk(1:T);

peak = find(ev>0); foot = find(ev_foot>0);
phi  = piecewise_phase_local(peak, foot, T);

pB2.Fs=fps; pB2.tapers=[P.TW_spec,2*P.TW_spec-1]; pB2.pad=0;
pB2.fpass=[P.fmin,min(P.fmax,fps/2)]; pB2.err=0;
[S2,f2] = mtspectrumc(bw, pB2); S2=S2(:); f2=f2(:);
mm = f2>=P.f_breath_search(1) & f2<=P.f_breath_search(2);
[~,rl] = max(S2(mm)); f_pk = f2(find(mm,1)+rl-1);

IBI = NaN;
if numel(foot) >= 2, IBI = median(diff(foot))/fps; end
if ~isfinite(IBI) || IBI <= 0, IBI = 1/max(f_pk,eps); end

r = struct('fps',fps,'T',T,'dff',dff,'bw',bw,'peak',peak,'foot',foot,'spk',spk, ...
           'phi',phi,'IBI',IBI,'f_pk',f_pk,'px_um',px_um,'sm',sm, ...
           'recName',o.recName,'folder',o.folder,'roi',o.roi);
end

function pr = roi_crop_local(folderPath, roi, px_um, P)
%ROI_CROP_LOCAL  Avg-projection crop around one ROI, plus its outline+scalebar.
%  Square box centred on the ROI centroid, half-width = ROI radius + P.pad_um, so
%  every recording of a cell covers about the same PHYSICAL extent regardless of
%  zoom -- which is what lets the crops be compared to each other.
pr = struct('have_proj',false,'crop_img',[],'bnd_crop',{{}},'barLen_pr',NaN,'roi',roi);
if ~isfinite(px_um), return; end
sam = dir(fullfile(folderPath,'*_cpSAM_output.mat'));
if isempty(sam), return; end
SAM = load(fullfile(sam(1).folder, sam(1).name), 'maskL');
if ~isfield(SAM,'maskL'), return; end
maskL  = SAM.maskL;
avgimg = read_avgproj_local(folderPath, size(maskL));
mask   = (maskL == roi);
if ~any(mask(:)) || isempty(avgimg), return; end
[Himg,Wimg] = size(avgimg);
[yy,xx] = find(mask); cx = mean(xx); cy = mean(yy); %#ok<ASGLU>
b0 = bwboundaries(mask,'noholes'); ball = cat(1, b0{:});

% CROP SIZE
%   P.crop_um : total crop WIDTH in microns -- a fixed physical field of view,
%               independent of how big the ROI is. Every cell then gets the same
%               scale, which is what makes crops comparable between cells and
%               keeps the ROI outline a consistent thickness.
%   P.pad_um  : legacy -- a MARGIN added to the ROI's own radius, so the crop was
%               2*(radius + pad) and the number matched no visible dimension.
%               Kept as the default so the batch is unchanged.
if isfield(P,'crop_um') && ~isempty(P.crop_um) && P.crop_um > 0
    half_px = 0.5 * P.crop_um / px_um;
else
    half_px = max(hypot(ball(:,2)-cx, ball(:,1)-cy)) + P.pad_um/px_um;
end
half_px = max(4, round(half_px));

% CENTRING
%   The window is built around the centroid and NOT clamped to the image. Clamping
%   (max(1,..)/min(W,..)) silently cut one side for any ROI near the FOV edge, so
%   the cell slid off centre -- into the corner for a corner ROI -- and the crop
%   was no longer the requested size either. Out-of-image area is filled instead,
%   so the ROI is always dead centre and the scale is always what was asked for.
cxr = round(cx); cyr = round(cy);
x0 = cxr - half_px; x1 = cxr + half_px;
y0 = cyr - half_px; y1 = cyr + half_px;
craw = nan(y1-y0+1, x1-x0+1);
xs = max(x0,1):min(x1,Wimg);
ys = max(y0,1):min(y1,Himg);
if isempty(xs) || isempty(ys), return; end
craw(ys-y0+1, xs-x0+1) = double(avgimg(ys, xs));

loi = prctile(craw(:),P.clip_pct(1)); hii = prctile(craw(:),P.clip_pct(2));
cimg = min(max((craw-loi)/max(hii-loi,eps),0),1) .^ P.gamma_val;
cimg(isnan(cimg)) = 0;        % outside the FOV -> black, not a broken pixel
bd = bwboundaries(imdilate(mask,strel('square',3)),'noholes');
for k = 1:numel(bd), bd{k} = [bd{k}(:,1)-y0+1, bd{k}(:,2)-x0+1]; end

% RESAMPLE TO A UNIFORM PIXEL SIZE
%   P.cropPxUm fixes the crop's microns-per-pixel, so every cell -- whatever zoom
%   it was acquired at -- comes out the same number of pixels for a given crop_um.
%   That is what makes the scale bar and the ROI outline a CONSTANT width: both
%   are specified in image pixels, and image pixels now mean a fixed physical
%   size. Without this, native px_um ranges ~0.4-1.8 here, so the same settings
%   drew a bar and an outline that changed thickness from cell to cell.
%   Stretch first, resample after: interpolating raw counts and then taking
%   percentiles would let the interpolation shift the black and white points.
%   Empty = keep native pixels, which is what the batch does.
pxOut = px_um;
if isfield(P,'cropPxUm') && ~isempty(P.cropPxUm) && P.cropPxUm > 0
    s = px_um / P.cropPxUm;
    % Resize to an EXPLICIT output size when the crop is a known physical width,
    % so every cell lands on exactly the same pixel grid. Scaling by the factor
    % alone leaves +/-2 px of slop, because half_px was rounded to whole native
    % pixels before scaling (measured 141-145 px for a 70 um crop).
    nOut = [];
    if isfield(P,'crop_um') && ~isempty(P.crop_um) && P.crop_um > 0
        nOut = max(8, round(P.crop_um / P.cropPxUm));
        s = nOut / size(cimg,1);              % actual scale used, for the outline
    end
    if isfinite(s) && s > 0 && abs(s-1) > 1e-6
        % bicubic: sharper than bilinear on somata. It overshoots at edges, which
        % is why the result is clamped back into [0 1] -- otherwise a bright soma
        % rim would ring past white and past black.
        if isempty(nOut), cimg = imresize(cimg, s, 'bicubic');
        else,             cimg = imresize(cimg, [nOut nOut], 'bicubic'); end
        cimg = min(max(cimg, 0), 1);
        % imresize maps pixel CENTRES, so a coordinate transforms as
        % (c - 0.5)*s + 0.5, not c*s -- using c*s drifts the outline half a
        % pixel off the cell, which is visible at these crop sizes.
        for k = 1:numel(bd), bd{k} = (bd{k} - 0.5)*s + 0.5; end
        pxOut = P.cropPxUm;
    end
end

pr = struct('have_proj',true,'crop_img',cimg,'bnd_crop',{bd}, ...
            'barLen_pr',max(1,round(P.scalebar_um/pxOut)),'roi',roi,'px_um',pxOut);
end

function iBest = pick_repr_obs_local(R, P)
%PICK_REPR_OBS_LOCAL  Index of the most representative recording of this cell:
%  the ROI crop with the highest MEAN NCC against the cell's other crops.
%  Crops are resampled to a common pixel size first -- they already cover the same
%  physical extent (radius + pad_um), so the resample only removes the zoom
%  difference, and correlation is then at zero lag because every crop is centred
%  on its own ROI centroid. Recordings with no usable crop score -Inf; if none
%  has a crop, fall back to the first.
iBest = 1;
n = numel(R);
if n < 2, return; end
SZ = 96;
C = cell(n,1); ok = false(n,1);
for i = 1:n
    pr = roi_crop_local(R(i).folder, R(i).roi, R(i).px_um, P);
    if ~pr.have_proj || any(size(pr.crop_img) < 4), continue; end
    img = imresize(pr.crop_img, [SZ SZ]);
    img = img - mean(img(:));
    s = std(img(:));
    if ~isfinite(s) || s <= 0, continue; end
    C{i} = img / s;  ok(i) = true;
end
if nnz(ok) < 2, if any(ok), iBest = find(ok,1); end, return; end
idx = find(ok);
score = -inf(n,1);
for a = 1:numel(idx)
    v = zeros(numel(idx)-1,1); m = 0;
    for b = 1:numel(idx)
        if a == b, continue; end
        m = m + 1;
        v(m) = sum(C{idx(a)}(:) .* C{idx(b)}(:)) / (SZ*SZ);   % zero-lag NCC
    end
    score(idx(a)) = mean(v);
end
[~, iBest] = max(score);
end

function Eo = epochs_local(x, idx, w)
%EPOCHS_LOCAL  [numel(idx) x 2w+1] matrix of x around each idx.
Eo = zeros(numel(idx), 2*w+1);
for k = 1:numel(idx), Eo(k,:) = x(idx(k)-w : idx(k)+w); end
end

function Eg = regrid_local(Ein, tk, tau)
%REGRID_LOCAL  Put one recording's epochs on the pooled time axis. Recordings can
%  differ in fps, so epochs are cut in native frames then interpolated here.
if isempty(Ein), Eg = zeros(0,numel(tau)); return; end
if numel(tk)==numel(tau) && max(abs(tk-tau)) < 1e-9, Eg = Ein; return; end
Eg = interp1(tk(:), Ein.', tau(:), 'linear', NaN).';
end

function h = trig_hist_local(evIdx, trigIdx, winH, edgesC, fps)
%TRIG_HIST_LOCAL  Counts of evIdx relative to every trigger, |dt| <= winH frames.
h = zeros(1, numel(edgesC)-1);
if isempty(evIdx) || isempty(trigIdx), return; end
d = double(evIdx(:)).' - double(trigIdx(:));
d = d(abs(d) <= winH);
if isempty(d), return; end
h = histcounts(d(:)/fps, edgesC);
end

function lw = outline_lw_local(ax, cropWpx, wantImgPx, fallbackLW)
%OUTLINE_LW_LOCAL  LineWidth (points) that draws wantImgPx IMAGE pixels thick.
%
%  MATLAB LineWidth is in POINTS -- fixed on the page -- but the crop is scaled to
%  fit whatever axes it lands in. So one setting gives a different apparent
%  thickness for every crop size, and different again between the figure panel and
%  the full-resolution standalone PNG: measured 3.0 px vs 8.0 px for the same
%  195x195 crop at outlineLW = 0.8.
%
%  Converting through the axes' on-page width makes the outline a constant
%  fraction of the IMAGE, which is what actually looks the same. Export resolution
%  cancels out, because exportgraphics scales the axes and the line together.
    lw = fallbackLW;
    if cropWpx <= 0, return; end
    try
        old = get(ax,'Units');  set(ax,'Units','points');
        p   = get(ax,'Position'); set(ax,'Units',old);
        if p(3) > 0, lw = max(0.15, wantImgPx * p(3) / cropWpx); end
    catch
        % keep the fallback
    end
end

% -----------------------------------------------------------------------
function overlay_trials_local(ax, tau, M, mu, col, ylim_dff)
%OVERLAY_TRIALS_LOCAL  Every single-trial trace, faint, with the mean on top.
%  Replaces the +/- SD ribbon. The ribbon compressed the spread into two numbers
%  and hid how many trials actually carry the response -- a handful of large
%  events and a hundred flat ones give the same band as a uniform modest response.
%
%  Y-LIMITS: single trials are far noisier than their mean, so a min/max scale
%  gets set by one outlier trial and flattens everything. When no manual limit is
%  given, clip to the 1-99th percentile of the trial values so the bulk of the
%  traces and the mean are both legible.
    if isempty(M), return; end
    % Single trials only -- the mean line was removed 2026-08-17 (RZ). `mu` is
    % still passed in so putting it back is one line, not a re-plumb.
    plot(ax, tau, M.', '-', 'Color',[col 0.20], 'LineWidth',0.30);
    xlim(ax, [tau(1) tau(end)]);
    if isempty(ylim_dff)
        q = prctile(M(:), [1 99]);
        if all(isfinite(q)) && q(2) > q(1)
            pad = 0.08*(q(2)-q(1));
            ylim(ax, [q(1)-pad, q(2)+pad]);
        end
    end
end

% -----------------------------------------------------------------------
function g = o_group_local(o, folderPath)
%O_GROUP_LOCAL  Genotype folder for a recording, not the ventral group label.
%  The mount is a property of the SESSION, so it must be keyed on the genotype
%  directory. The group label cannot be used: the IO sites of Vglut2/1124 carry
%  group 'IO', and keying on that would apply the wrong mount to them -- the same
%  trap that made the 1124 one-frame breath fix miss those recordings.
    g = '';
    if isstruct(o) && isfield(o,'group'), g = o.group; end
    if isempty(g) || strcmpi(g,'IO'), g = detect_genotype_local(folderPath); end
end

% -----------------------------------------------------------------------
function tf = is_rotated_session_local(genotype, dateStr)
%IS_ROTATED_SESSION_LOCAL  Was this session mounted rotated 90 deg?
%  THE RULE IS TEMPORAL: everything acquired from 2026-07-21 onward is on the new
%  mount. The list is explicit rather than derived, because the MMDD folder names
%  cannot be ordered without a year -- 1124 is a PREVIOUS year and is NOT rotated,
%  while 0721 of 2026 is -- and the raw tifs that carry the true date live on
%  drives that are often not mounted.
%
%  Do NOT try to infer the mount from where the stage coordinates land. Vgat/0730
%  sits squarely inside the old-mount cloud and is still rotated.
%
%  ADD ANY NEW SESSION HERE. Assume new data is rotated.
    ROT = { 'Sert','0721'
            'Vglut2','0728'
            'Vgat','0730'
            'Sst','0806'      % 2026-08-06, confirmed from the ScanImage epoch
            'Sst','0807'      % 2026-08-07, confirmed
            'Vglut2','0810'   % 2026-08-10, confirmed
            'Vglut2','0824' };% 2026-08-24, vagotomised session (temporal rule)
    tf = any(strcmpi(ROT(:,1), genotype) & strcmpi(ROT(:,2), dateStr));
end

% -----------------------------------------------------------------------
function dt = acq_datetime_local(folderPath)
%ACQ_DATETIME_LOCAL  Acquisition date/time from the ScanImage header.
%  ScanImage writes   epoch = [YYYY MM DD HH MM SS.sss]   into the
%  ImageDescription of the FIRST frame. Only the RAW tif carries it -- every
%  processed tif in this pipeline has its metadata stripped -- so this has to
%  reach the raw file through _meta.mat/source_tif, which points at the
%  acquisition drive (C:\260810_...\phys\...), not at the archive copy.
%
%  Returns NaT when the raw file cannot be reached, which is NORMAL: the archive
%  is self-contained and the acquisition drive is often absent. Callers must fall
%  back to the MMDD from the folder, never error.
%
%  Cached per folder: imfinfo on a multi-GB tif is cheap (first IFD only) but the
%  batch calls this 269 times.
    persistent CACHE
    if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','any'); end
    key = lower(folderPath);
    if isKey(CACHE, key), dt = CACHE(key); return; end

    dt = NaT;
    try
        src = '';
        dm = dir(fullfile(folderPath,'*_meta.mat'));
        if ~isempty(dm)
            M = load(fullfile(folderPath, dm(1).name), 'source_tif');
            if isfield(M,'source_tif') && isfile(M.source_tif), src = M.source_tif; end
        end
        if isempty(src)
            % a raw tif sitting in the folder itself -- exclude processed ones,
            % which have no metadata and would silently give nothing
            dd = dir(fullfile(folderPath,'*.tif'));
            if ~isempty(dd)
                bad = contains({dd.name}, {'_preproc','_MC','_AVG','_ROI','minusDark'});
                dd = dd(~bad);
                if ~isempty(dd), src = fullfile(folderPath, dd(1).name); end
            end
        end
        if ~isempty(src)
            info = imfinfo(src);
            tok = regexp(info(1).ImageDescription, 'epoch\s*=\s*\[([^\]]*)\]', ...
                         'tokens','once');
            if ~isempty(tok)
                v = sscanf(tok{1}, '%f');
                if numel(v) >= 6
                    dt = datetime(v(1), v(2), v(3), v(4), v(5), v(6));
                end
            end
        end
    catch
        dt = NaT;                 % unreadable header is not an error here
    end
    CACHE(key) = dt;
end

% -----------------------------------------------------------------------
function [th, Rbar, nEff] = wresultant_local(a, w)
%WRESULTANT_LOCAL  Weighted circular resultant with Kish's effective sample size.
%  Z = nEff*Rbar^2 with nEff = (sum w)^2/sum(w^2): using raw n with weighted Rbar
%  breaks the Rayleigh, because a few events in a sparse inspiratory bin carry
%  huge leverage. Pooling across recordings just concatenates a and w, so nEff
%  accumulates exactly as the ventral analysis intends.
S1 = sum(w); S2 = sum(w.^2);
if S1 <= 0, th = NaN; Rbar = 0; nEff = 0; return; end
nEff = S1^2 / max(S2, eps);
v    = sum(w(:) .* exp(1i*a(:))) / S1;
th   = angle(v);
Rbar = min(abs(v), 1);
end

function s = fmt_pval_local(p, nSh)
%FMT_PVAL_LOCAL  Stars: * <0.05, ** <0.01, *** <0.001. At nShuffle=1200 the p
%  floors at 1/1201 = 8.3e-4, so *** is reachable; below ~1000 shuffles it is not
%  and the star is capped rather than claiming a threshold that cannot resolve.
if ~isfinite(p),  s = 'n/a';
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
else,             s = 'n.s.';
end
if nSh > 0 && p < 0.001 && 1/(1+nSh) >= 0.001, s = '**'; end
end

function out = ternary_local(cond, a, b)
if cond, out = a; else, out = b; end
end

function r = range_local(y)
r = max(y(:)) - min(y(:));
end

function set_ylim_win(ax, y, mwsp)
lo = min(y(mwsp)); hi = max(y(mwsp)); pad = 0.05*max(hi-lo, eps);
ylim(ax, [lo-pad hi+pad]);
end

function sync_breath_strip(ax_b, ax4, nUnitHm, nUnitStrip)
%SYNC_BREATH_STRIP  Park the breath strip on top of the heatmap, same x and width.
%  tiledlayout aligns tile OUTER boxes but lets each axes size its own inner box,
%  so a tile would not line up with an axes carrying a colorbar.
%  The strip's height is nUnitStrip GRID ROWS, derived from the heatmap's own
%  height and its span -- not a fixed fraction of the heatmap, so growing the
%  heatmap does not stretch the strip with it.
if ~isgraphics(ax_b) || ~isgraphics(ax4), return; end
if nargin < 3 || isempty(nUnitHm),    nUnitHm    = 9; end
if nargin < 4 || isempty(nUnitStrip), nUnitStrip = 3; end
p    = get(ax4,'Position');
unit = p(4)/nUnitHm;                      % one grid row, in figure units
h    = nUnitStrip*unit;
gap  = 0.28*unit;
set(ax_b,'Units','normalized','Position',[p(1), p(2)+p(4)+gap, p(3), h-gap]);
end

function avgimg = read_avgproj_local(folderPath, szWant)
%READ_AVGPROJ_LOCAL  Average projection matching the cpSAM mask.
%  *_AVG_for_CP.tif FIRST: it is the frame Cellpose segmented, so it matches
%  maskL by construction, it exists in 111/111 archived recordings, and it is one
%  frame rather than a 1.5 GB stack. Only 61/111 still have the movie, which is
%  why preferring the movie silently dropped this panel for the rest.
avgimg = [];
if nargin < 2, szWant = []; end
cand = [ dir(fullfile(folderPath,'*_MC_MC_AVG_for_CP.tif'))
         dir(fullfile(folderPath,'AVG_*_ch1_preproc_MC_MC.tif'))
         dir(fullfile(folderPath,'AVG_*_MC_MC.tif')) ];
mv = dir(fullfile(folderPath,'*_ch1_preproc_MC_MC.tif'));
mv = mv(~contains({mv.name},'AVG','IgnoreCase',true));
cand = [cand; mv];
for k = 1:numel(cand)
    try
        V = tiffreadVolume(fullfile(cand(k).folder, cand(k).name));
    catch
        continue;
    end
    img = mean(double(V),3);
    if isempty(szWant) || isequal(size(img), szWant), avgimg = img; return; end
end
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
%PIECEWISE_PHASE_LOCAL  Insp onset = 0, peak = pi, next onset = 2pi, linear in
%  TIME between consecutive landmarks. Strict alternation enforced; NaN outside
%  the first/last landmark. Phase accumulates past 2pi (wrap with mod).
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

function dt = acq_date_local(folderPath, sm)
%ACQ_DATE_LOCAL  Acquisition date, for choosing the laser-power table.
dt = NaT;
tok = regexp(folderPath, '(?<!\d)(2\d)(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?!\d)', 'tokens');
for i = 1:numel(tok)
    cand = datetime(2000+str2double(tok{i}{1}), str2double(tok{i}{2}), str2double(tok{i}{3}));
    if cand >= datetime(2020,1,1) && cand <= datetime(2030,1,1), dt = cand; return; end
end
% archive layout is <root>\<Genotype>\<MMDD>\... -- MMDD carries no year, so take
% the most recent occurrence that is not in the future ('0124' = this January,
% '1124' = last November). Only the laser table depends on it.
parts = strsplit(regexprep(folderPath,'[\\/]+$',''), {'\','/'});
for i = 1:numel(parts)
    if ~isempty(regexp(parts{i}, '^\d{4}$', 'once'))
        mo = str2double(parts{i}(1:2)); dy = str2double(parts{i}(3:4));
        if mo >= 1 && mo <= 12 && dy >= 1 && dy <= 31
            today0 = datetime('today');
            cand = datetime(year(today0), mo, dy);
            if cand > today0, cand = datetime(year(today0)-1, mo, dy); end
            dt = cand; return;
        end
    end
end
tf = dir(fullfile(folderPath,'*.tif'));
tf = tf(~contains({tf.name},{'_preproc','_MC','_AVG'},'IgnoreCase',true));
if ~isempty(tf), dt = datetime(tf(1).datenum,'ConvertFrom','datenum'); return; end
if isfield(sm,'source_tif') && ~isempty(sm.source_tif) && isfile(sm.source_tif)
    f = dir(sm.source_tif);
    if ~isempty(f), dt = datetime(f.datenum,'ConvertFrom','datenum'); end
end
end

function g = detect_genotype_local(fp)
%DETECT_GENOTYPE_LOCAL  Genotype token from anywhere in the path. IO needs a
%  delimiter on both sides or it matches inside ordinary words.
g = '';
toks = {'Vglut2','Vgat','ChAT','Sert','Sst'};
for i = 1:numel(toks)
    if ~isempty(regexpi(fp, toks{i}, 'once')), g = toks{i}; return; end
end
if ~isempty(regexpi(fp, '(^|[\\/_\-\s])IO([\\/_\-\s]|$)', 'once')), g = 'IO'; end
end
