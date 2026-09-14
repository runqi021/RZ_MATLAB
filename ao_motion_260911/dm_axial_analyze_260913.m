function R = dm_axial_analyze_260913(sessDir, opts)
%DM_AXIAL_ANALYZE_260913  Offline analysis of a dm_axial_motion_acq_260913 session.
%
%   R = dm_axial_analyze_260913('D:\Data\RZ\ao_axial_260913\axial_260913_2130')
%
% Answers the go/no-go question: HOW FAR DOES THE BRAIN MOVE IN Z PER BREATH?
% Nothing here touches the rig. Everything is recomputable from the logged raw
% frames, so gates and parameters can be retuned without re-imaging.
%
% WHAT IT DOES, IN ORDER
%   1  Build the reference volume from block B, BREATH-GATED (see below).
%   2  Calibrate the coefficient axis in um against the block A stage stack --
%      this is the tissue repeat of dm_zstack_demo_260911, for free, from data
%      that had to be taken anyway.
%   3  Run the gradient estimator over block C, in chunks, -> dz(t), dx(t), dy(t).
%   4  Take breath phase FROM THE LATERAL SHIFTS THEMSELVES, fit dz(phase),
%      report peak-to-peak axial excursion against a circular-shift null.
%   5  Check block D: recover a defocus staircase the estimator was told about.
%
% WHY BREATH PHASE COMES OUT OF THE IMAGE AND NOT A CAMERA
% This machine has no conda and no ffmpeg, so the Basler and thermal legs do
% not run here. It turns out not to matter: the estimator already solves for
% dx, dy per strip at ~240 Hz, and lateral breathing motion is large and
% obvious. Phase from lateral, displacement from axial -- two independent
% quantities out of the same frames, no extra hardware, nothing to synchronise.
%
% WHY THE REFERENCE IS GATED -- AND HOW MUCH IT ACTUALLY MATTERS
% Block B was acquired while the animal breathed, so each plane is the true
% structure convolved with the breathing displacement.
%
% MEASURED 2026-09-13 (test_axial_estimator_260913, tests 6 and 6b): this does
% NOT bias the recovered slope. 4 um of reference blur gives slope 1.004; a
% reference averaged from 2 frames per plane gives 0.995. Smoothing is a
% deterministic linear operation, so least squares stays consistent, and the
% fit averages over thousands of pixels so reference noise washes out too.
% Blur costs a constant offset and some sensitivity, not a gain error.
%
% Gating is therefore DEFENSIVE and cheap, not load-bearing. Set opts.gate =
% false and the numbers should barely move -- if they move a lot, something
% else is wrong and that is worth knowing. The gate is on the MODE of the
% lateral displacement, not the median: a breathing trace dwells at
% end-expiration, so the density peak is end-expiration and the median is not.
%
% Runqi Zhang / 2026-09-13.  NOT YET RUN ON REAL DATA.

if nargin < 2, opts = struct(); end
opts = fill_defaults(opts);

assert(exist(sessDir,'dir') == 7, 'no such session folder: %s', sessDir);
M = load(fullfile(sessDir,'meta.mat'));
fprintf('\n=== %s ===\n', sessDir);
fprintf('channel %d | baseline %s | nominal gain %.3f um/coeff\n', ...
        M.channel, M.baselineVar, M.gainUmPerCoeff);

R = struct('sessDir', sessDir, 'meta', M, 'opts', opts);
outDir = fullfile(sessDir, 'analysis');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ===================== 1. REFERENCE VOLUME (block B) =====================
fprintf('\n--- 1. reference volume from block B (gated) ---\n');
assert(isfield(M,'B'), 'session has no block B -- there is no reference to register against');

Bz = M.B.zNom(M.B.got);
Bf = M.B.files(M.B.got);
K  = numel(Bz);
refVol = [];  keptFrac = nan(1,K);
for k = 1:K
    F = si_read_all(Bf{k}, M.channel);
    [mu, kept] = gate_mean(F, opts);
    keptFrac(k) = kept;
    if isempty(refVol), refVol = zeros([size(mu) K], 'single'); end
    refVol(:,:,k) = mu; %#ok<AGROW>
end
fprintf('  %d planes, %.1f to %.1f um nominal, gate kept %.0f%% of frames (min %.0f%%)\n', ...
        K, Bz(1), Bz(end), 100*mean(keptFrac), 100*min(keptFrac));
if min(keptFrac) * M.B.nFrames < 4
    warning('dm_axial:thinGate', ...
        ['the tightest plane kept only %.0f frames. Synthetic testing says even ' ...
         '2 frames/plane holds the slope to 1%%, so this is unlikely to bias ' ...
         'anything -- but it does cost sensitivity. Widen opts.gateWidth or set ' ...
         'opts.gate = false.'], min(keptFrac)*M.B.nFrames);
end
R.refVol = refVol;  R.refZnom = Bz;  R.keptFrac = keptFrac;

%% ============ 2. CALIBRATE THE COEFFICIENT AXIS IN TISSUE (block A) ============
% dm_zstack_demo_260911 measured ~7 um per um of coefficient on a STATIC TEST
% SAMPLE. Tissue scatters and the PSF degrades with depth, so the gain has to be
% confirmed in the real medium before any um number here means anything.
umPerCoeff = M.gainUmPerCoeff;
if isfield(M,'A') && opts.doCalib
    fprintf('\n--- 2. tissue gain: DM planes registered into the STAGE stack ---\n');
    Az = M.A.zMeas;                          % encoder-derived um, the ruler
    stageVol = [];
    for k = 1:numel(M.A.files)
        F = si_read_all(M.A.files{k}, M.channel);
        mu = gate_mean(F, opts);
        if isempty(stageVol), stageVol = zeros([size(mu) numel(M.A.files)],'single'); end
        stageVol(:,:,k) = mu; %#ok<AGROW>
    end
    Ecal = dm_axial_estimator_260913('prep', stageVol, Az, opts.est);

    matched = nan(1,K);  rpk = nan(1,K);
    for k = 1:K
        [matched(k), rpk(k)] = dm_axial_estimator_260913('coarse', Ecal, refVol(:,:,k));
    end
    cUsed = M.B.coeff(M.B.got);
    ok = isfinite(matched) & isfinite(cUsed);
    if nnz(ok) >= 3
        p = polyfit(cUsed(ok), matched(ok), 1);
        umPerCoeff = p(1);
        resid = matched(ok) - polyval(p, cUsed(ok));
        fprintf('  matched_z = %.3f * coeff %+.2f\n', p(1), p(2));
        fprintf('  TISSUE GAIN %.2f um per um of coefficient  (bench was 6.66-7.04)\n', p(1));
        fprintf('  residual RMS %.2f um over %d planes, median peak r %.3f\n', ...
                sqrt(mean(resid.^2)), nnz(ok), median(rpk(ok)));
        R.calib = struct('coeff',cUsed,'matched',matched,'rpk',rpk,'fit',p, ...
                         'residRms',sqrt(mean(resid.^2)));
    else
        warning('dm_axial:calibFail','too few matched planes to fit a tissue gain');
    end
else
    fprintf('\n--- 2. skipped (no block A, or opts.doCalib = false) ---\n');
    fprintf('  using the BENCH gain %.3f um/coeff. Every um below inherits its error.\n', umPerCoeff);
end
R.umPerCoeff = umPerCoeff;

% The reference z axis is what the estimator reports in. Rescale nominal um
% (which was built with the bench gain) onto the measured tissue gain.
refZ = Bz * (umPerCoeff / M.gainUmPerCoeff);
R.refZ = refZ;

%% ===================== 3. ESTIMATOR OVER THE TIME SERIES =====================
fprintf('\n--- 3. estimating dz over block C ---\n');
assert(isfield(M,'C'), 'session has no block C -- there is no measurement to analyse');

E = dm_axial_estimator_260913('prep', refVol, refZ, opts.est);
R.E = rmfield(E, {'ref','Rz','M','idx','colN','coarseRef','regRef','win'});  % keep it small

C = run_series(M.C.file, M.channel, E, opts);
R.C = C;
S  = E.S;
fps = pick_fps(C.t, M.fpsHeader);
fprintf('  %d frames, %.2f Hz (from %s), %d strips -> %.1f Hz effective\n', ...
        C.nFrames, fps, C.fpsSource, S, fps*S);
fprintf('  coarse re-lock fired on %d of %d frames\n', nnz(C.coarseUsed), C.nFrames);

% Strip time base. A frame is not a snapshot: strip s is acquired at
% (s-0.5)/S of the way through the frame period, and ignoring that shears
% every phase estimate by up to a whole frame.
tStrip = reshape(C.t(:).' + ((1:S)'-0.5)/S / fps, [], 1);
dz = C.dzSerial;  dx = C.dxSerial;  dy = C.dySerial;
good = isfinite(dz) & isfinite(tStrip);
if isfield(C,'residSerial')
    rmax = qtile(C.residSerial(isfinite(C.residSerial)), opts.residReject);
    good = good & (C.residSerial <= rmax);
end
fprintf('  %d of %d strip estimates usable (%.1f%%)\n', nnz(good), numel(dz), 100*nnz(good)/numel(dz));

%% ============ 4. THE MEASUREMENT: how big is dz, and how fast ============
% PHASE-FREE ON PURPOSE. The go/no-go needs exactly two things about dz(t) --
% how large it is, and what frequencies it contains. Neither requires a breath
% sensor, a phase model, or any assumption that the motion is periodic.
%
% A phase-binned peak-to-peak silently assumes periodicity and discards
% everything not locked to it: sniff bouts, apnea, the cardiac component, and
% drift. Those all move the focus just as much. Section 4c below still reports
% the phase description because it is informative, but the numbers that decide
% whether the DM has enough stroke are these.
fprintf('\n--- 4. dz amplitude and spectrum (no phase model) ---\n');
dzd0 = detrend_slow(dz, tStrip, good, opts.driftHz);
v = dzd0(good);
fsStrip = (numel(tStrip)-1) / range(tStrip(isfinite(tStrip)));

pctl = [1 5 50 95 99];
qv = arrayfun(@(p) qtile(v, p/100), pctl);
fprintf('    RMS                         : %6.2f um\n', sqrt(mean(v.^2)));
fprintf('    5-95%% spread               : %6.2f um\n', qv(4)-qv(2));
fprintf('    1-99%% spread               : %6.2f um\n', qv(5)-qv(1));
fprintf('    full range (min to max)     : %6.2f um\n', max(v)-min(v));

% Spectrum: where the motion actually lives decides how fast the loop must be.
nfft = 2^nextpow2(min(numel(v), 2^15));
vv = v(1:min(numel(v),nfft));
vv = (vv - mean(vv)) .* (0.5 - 0.5*cos(2*pi*(0:numel(vv)-1).'/(numel(vv)-1)));
P = abs(fft(vv, nfft)).^2;
fAx = (0:nfft-1).' * fsStrip / nfft;
h = fAx <= fsStrip/2;
Ph = P(h);  fh = fAx(h);
[~, ipk] = max(Ph(fh > 0.3));
fsub = fh(fh > 0.3);  fPeak = fsub(ipk);
cum = cumsum(Ph) / sum(Ph);
f90 = fh(find(cum >= 0.90, 1));
fprintf('    dominant frequency          : %6.2f Hz\n', fPeak);
fprintf('    90%% of power below         : %6.2f Hz\n', f90);

R.amplitude = struct('rms',sqrt(mean(v.^2)),'pctl',pctl,'q',qv, ...
                     'range',max(v)-min(v),'fPeak',fPeak,'f90',f90, ...
                     'spec',Ph,'specF',fh,'fsStrip',fsStrip);

%% ============ 4b. WHAT A REACTIVE LOOP WOULD ACTUALLY LEAVE ============
% SIMULATED on the measured dz(t), not predicted from a formula.
%
% The closed-form 2*sin(pi*f*tau) assumes a single sinusoid. Real breathing is
% not one: inspiration is sharp, so there are harmonics, and there are sniffs,
% pauses and drift. Those are exactly the components a delayed loop handles
% worst, so they must be in the test. Replaying the loop over the real trace
% costs nothing and answers the question honestly.
%
% The loop being simulated is the whole design: image in -> compare to the
% reference -> dz -> send -dz. Nothing external enters.
%
% VELOCITY EXTRAPOLATION. The error is pure DELAY, so extrapolate forward with
% the velocity already available from consecutive estimates:
%     send  -( e[k] + gx*(e[k] - e[k-1]) )
% Two previous estimates and a subtraction. This is NOT a breath model -- no
% periodicity is assumed, nothing external is used, and it keeps working during
% sniffs and apnea, which is where a phase model fails worst. It trades noise
% for lag, and there is noise budget to spend (estimator scatter ~0.07 um
% against um-scale motion).
fprintf('\n--- 4b. simulated reactive loop on the measured trace ---\n');
vClean = dzd0;
vClean(~good) = interp1(find(good), dzd0(good), find(~good), 'linear', 'extrap');
vClean = vClean(isfinite(vClean));
rms0 = sqrt(mean((vClean - mean(vClean)).^2));
fprintf('  uncorrected RMS %.3f um. Numbers below are %% of that left behind.\n', rms0);
fprintf('  (a loop cannot do better than the estimator noise floor, ~0.1 um)\n\n');
fprintf('%10s %8s %10s %10s %10s\n','loop rate','lines','no extrap','gx=0.5','gx=1.0');

rates = unique(round([fps, 60, fps*2, 125, fps*S]));
rates = rates(rates >= 5 & rates <= fsStrip/2);
loopTab = nan(numel(rates), 3);
for i = 1:numel(rates)
    for j = 1:3
        gx = [0 0.5 1.0];
        loopTab(i,j) = sim_loop(vClean, fsStrip, rates(i), gx(j)) / rms0;
    end
    linesAt = round(fps * 512 / rates(i));   % 512 lines at the measured fps
    fprintf('%8.0f Hz %8d %9.0f%% %9.0f%% %9.0f%%%s\n', rates(i), linesAt, ...
            100*loopTab(i,1), 100*loopTab(i,2), 100*loopTab(i,3), ...
            ternary(loopTab(i,1) >= 1, '   <-- plain loop AMPLIFIES', ''));
end
fprintf('\n  "lines" = frame height needed to reach that rate, at the current\n');
fprintf('  line rate and 512 px/line. Halving lines halves the FOV height.\n');
R.loop = struct('rates',rates,'residFrac',loopTab,'rms0',rms0,'gx',[0 0.5 1.0]);

%% ===================== 4c. BREATH PHASE (description only) =====================
fprintf('\n--- 4c. breath phase from the lateral shifts (descriptive) ---\n');
lat = lateral_pc1(dx, dy, good);
[phi, fBreath, latBP] = breath_phase(lat, tStrip, good, opts);
fprintf('  breathing at %.2f Hz (peak of the lateral spectrum, %.1f-%.1f Hz band)\n', ...
        fBreath, opts.breathBand(1), opts.breathBand(2));

dzd = detrend_slow(dz, tStrip, good, opts.driftHz);   % remove thermal drift only
[binC, binM, binSem, nBin] = phase_bin(phi(good), dzd(good), opts.nPhaseBin);
p2p = max(binM) - min(binM);

X = [ones(nnz(good),1), cos(phi(good)), sin(phi(good))];
bfit = X \ dzd(good);
amp  = hypot(bfit(2), bfit(3));
mu   = atan2(-bfit(3), bfit(2));

fprintf('\n  *** AXIAL EXCURSION PER BREATH ***\n');
fprintf('  peak-to-peak, %d-bin phase average : %6.2f um\n', opts.nPhaseBin, p2p);
fprintf('  cosine amplitude (half of p2p)     : %6.2f um  at phase %+.0f deg\n', amp, rad2deg(mu));
fprintf('  raw 5-95%% spread of dz            : %6.2f um\n', ...
        qtile(dzd(good),0.95) - qtile(dzd(good),0.05));

% --- circular-shift null ------------------------------------------------
% Is the modulation real, or would any signal with this autocorrelation give
% it? Shift dz against phase and refit. Shifting preserves the autocorrelation
% of dz exactly, which a random permutation does not.
nS = opts.nShuffle;
nullAmp = nan(1,nS);
dzg = dzd(good);  phig = phi(good);  ng = numel(dzg);
minShift = max(10, round(ng/20));
for i = 1:nS
    sh = minShift + randi(ng - 2*minShift);
    Xs = [ones(ng,1), cos(phig), sin(phig)];
    bs = Xs \ circshift(dzg, sh);
    nullAmp(i) = hypot(bs(2), bs(3));
end
pVal = (1 + nnz(nullAmp >= amp)) / (nS + 1);
fprintf('  null (circular shift, n=%d): %.3f +- %.3f um, p = %.4f\n', ...
        nS, mean(nullAmp), std(nullAmp), pVal);
if pVal >= 0.05
    fprintf('  >> NOT significant. Either the brain is not moving axially, or\n');
    fprintf('     the estimator is not sensitive enough. Block D tells you which.\n');
end

R.phase = struct('phi',phi,'fBreath',fBreath,'binC',binC,'binM',binM, ...
                 'binSem',binSem,'nBin',nBin,'p2p',p2p,'amp',amp,'mu',mu, ...
                 'nullAmp',nullAmp,'p',pVal,'lat',lat,'latBP',latBP);
R.dz = dz;  R.dzDetrended = dzd;  R.tStrip = tStrip;  R.good = good;

%% ============ 4b. IS IT REALLY AXIAL MOTION, OR A BREATH-LOCKED ARTIFACT? ============
% THE CONFOUND THIS SESSION EXISTS TO RULE OUT.
%
% A breath-locked signal is not necessarily breath-locked MOTION. This prep has
% already produced one that was not: deep Vgat, top SVD mode at r = 0.96 with
% chest breathing, SURVIVED motion correction, 4x stronger in background than in
% somata -> optical, not neural, and not displacement. Anything of that shape
% will project onto dR/dz and be reported here as axial motion. The gain term
% only absorbs GLOBAL brightness; a spatially structured intensity change walks
% straight through it.
%
% The discriminator: a real axial displacement of dz must change the image by
% dz * dR/dz -- a pattern we have MEASURED, not assumed. So project the observed
% breath-triggered frame difference onto dR/dz and ask how much of it that
% explains. Genuine motion is ~all of it. An artifact is ~none of it and leaves
% a large structured residual.
fprintf('\n--- 4b. axial motion, or a breath-locked optical artifact? ---\n');
try
    % Mean frame at the two phase extremes of the fitted dz(phase).
    phiHi = mu;  phiLo = wrapToPiLocal(mu + pi);
    dObs  = phase_extreme_diff(M.C.file, M.channel, phi, good, phiHi, phiLo, opts, E);

    kMid = E.kMid;
    g0 = double(E.Rz(:,:,kMid));
    d0 = double(dObs);
    sel = isfinite(g0) & isfinite(d0);
    % Least-squares projection: how many um of real defocus does the observed
    % difference correspond to, and how much of its variance does that explain?
    a  = (g0(sel).' * d0(sel)) / (g0(sel).' * g0(sel));
    ve = 1 - sum((d0(sel) - a*g0(sel)).^2) / sum(d0(sel).^2);
    rr = corr_local(g0(sel), d0(sel));

    fprintf('  breath-extreme difference projected onto dR/dz:\n');
    fprintf('    implied displacement        : %+.2f um\n', a);
    fprintf('    correlation with dR/dz      : %+.3f\n', rr);
    fprintf('    variance explained by dR/dz : %.1f%%\n', 100*ve);
    if ve > 0.5
        fprintf('  >> CONSISTENT WITH REAL AXIAL MOTION. The image changes the way a\n');
        fprintf('     focus shift makes it change.\n');
    elseif ve > 0.2
        fprintf('  >> MIXED. Some of the breath-locked change is axial, some is not.\n');
        fprintf('     Treat the p2p above as an UPPER BOUND.\n');
    else
        fprintf('  >> NOT AXIAL MOTION. The breath-locked change does not look like a\n');
        fprintf('     focus shift. This is the optical-artifact signature -- do NOT\n');
        fprintf('     quote the p2p above as a displacement.\n');
    end
    R.signature = struct('impliedUm',a,'corr',rr,'varExp',ve,'diff',dObs);
catch ME
    warning('dm_axial:signature','signature test failed: %s', ME.message);
end

%% ===================== 5. POSITIVE CONTROL (block D) =====================
if isfield(M,'D') && opts.doInject
    fprintf('\n--- 5. positive control: recovering the injected staircase ---\n');
    D = run_series(M.D.file, M.channel, E, opts);
    fpsD  = pick_fps(D.t, M.fpsHeader);
    tD    = reshape(D.t(:).' + ((1:S)'-0.5)/S / fpsD, [], 1);
    dzD   = D.dzSerial;

    cmd = M.D.stepsUm(:).' * (umPerCoeff / M.gainUmPerCoeff);  % in measured um
    tS  = M.D.tSendSec(:).';
    nSt = numel(cmd);
    meas = nan(1,nSt);  sdm = nan(1,nSt);
    for k = 1:nSt
        t1 = tS(k) + opts.injSettleS;
        if k < nSt, t2 = tS(k+1); else, t2 = tS(k) + M.D.holdS; end
        sel = tD >= t1 & tD < t2 & isfinite(dzD);
        if nnz(sel) > 5
            meas(k) = median(dzD(sel));
            sdm(k)  = std(dzD(sel));
        end
    end
    ok = isfinite(meas);
    % Referenced to the 0 um steps so a constant offset does not enter.
    base = median(meas(ok & cmd == 0));
    if ~isfinite(base), base = 0; end
    meas = meas - base;

    fprintf('%12s %12s %12s\n','commanded um','measured um','sd um');
    for k = 1:nSt
        fprintf('%12.2f %12.2f %12.2f\n', cmd(k), meas(k), sdm(k));
    end
    if nnz(ok) >= 3
        pD = polyfit(cmd(ok), meas(ok), 1);
        rr = meas(ok) - polyval(pD, cmd(ok));
        fprintf('\n  measured = %.3f * commanded %+.2f um, residual RMS %.2f um\n', ...
                pD(1), pD(2), sqrt(mean(rr.^2)));
        fprintf('  SLOPE 1 = the estimator recovers what it was told. ');
        if abs(pD(1)-1) > 0.25
            fprintf('THIS ONE DOES NOT.\n');
            fprintf(['  On synthetic data the estimator holds slope to 1%% even with a\n' ...
                     '  blurred or 2-frame reference, so a bad slope HERE is not the gate.\n' ...
                     '  Look instead at: the tissue gain from block A (is umPerCoeff right?),\n' ...
                     '  whether the steps left the linear zone (opts.est.trustUm), and\n' ...
                     '  whether the sample moved between blocks B and D.\n']);
        else
            fprintf('This one does.\n');
        end
        R.inject = struct('cmd',cmd,'meas',meas,'sd',sdm,'fit',pD, ...
                          'residRms',sqrt(mean(rr.^2)));
    end
else
    fprintf('\n--- 5. skipped (no block D) ---\n');
    fprintf('  WITHOUT THE POSITIVE CONTROL a null result in step 4 is uninterpretable.\n');
end

%% ===================== 6. FIGURE + SAVE =====================
make_figure(R, outDir, opts);
save(fullfile(outDir,'axial_analysis.mat'), '-struct', 'R', '-v7.3');
fprintf('\nsaved %s\n', fullfile(outDir,'axial_analysis.mat'));

fprintf('\n=========================== VERDICT ===========================\n');
% The stroke question is answered by the FULL excursion, not a phase average --
% the mirror has to reach wherever the tissue actually goes, including the
% sniffs and the outliers that a phase bin averages away.
exc = R.amplitude.q(5) - R.amplitude.q(1);      % 1-99 percentile
fprintf('axial excursion (1-99%%)   : %.2f um   [RMS %.2f, full range %.2f]\n', ...
        exc, R.amplitude.rms, R.amplitude.range);
fprintf('motion lives below         : %.1f Hz (peak %.1f Hz)\n', R.amplitude.f90, R.amplitude.fPeak);
fprintf('DM authority (measured)    : -14 to +24 um of focus\n');
if exc < 8
    fprintf('=> STROKE COMFORTABLE. The DM has several times what is needed.\n');
elseif exc < 20
    fprintf('=> STROKE WORKABLE, but re-centre the operating point first.\n');
else
    fprintf('=> STROKE TIGHT. %.1f um against a -14/+24 envelope leaves little margin.\n', exc);
end

% The loop-rate question is separate and is usually the binding one.
tauF = 1.5/fps;  tauS = 1.5/(fps*S);
fprintf('\nloop rate needed (residual = 2*sin(pi*f*tau), >1 means it AMPLIFIES):\n');
fprintf('  per frame (%.0f Hz)  : %3.0f%% of the motion left at %.1f Hz\n', ...
        fps, 100*2*sin(pi*R.amplitude.f90*tauF), R.amplitude.f90);
fprintf('  per strip (%.0f Hz)  : %3.0f%% of the motion left at %.1f Hz\n', ...
        fps*S, 100*2*sin(pi*R.amplitude.f90*tauS), R.amplitude.f90);
fprintf('  (image in -> dz -> send -dz. No breath signal enters the loop.)\n');
fprintf('===============================================================\n');
end

%% ========================== SERIES RUNNER ==========================
function out = run_series(tifPath, channel, E, opts)
%RUN_SERIES  Stream a logged TIFF through the estimator in chunks.
assert(exist(tifPath,'file') == 2, 'missing time series file: %s', tifPath);
[t, nCh, chIdx, fpsHdr] = si_header(tifPath, channel); %#ok<ASGLU>

acc = struct('dz',[],'dx',[],'dy',[],'gain',[],'resid',[],'coarseUsed',[],'t',[]);
zPrev = NaN;

    function feed(F, tF)
        o = dm_axial_estimator_260913('run', E, F, zPrev);
        zPrev = o.zEnd;
        acc.dz    = [acc.dz,   o.dz];
        acc.dx    = [acc.dx,   o.dx];
        acc.dy    = [acc.dy,   o.dy];
        acc.gain  = [acc.gain, o.gain];
        acc.resid = [acc.resid, o.resid];
        acc.coarseUsed = [acc.coarseUsed, o.coarseUsed];
        acc.t = [acc.t, tF];
        fprintf('    %d frames\r', numel(acc.t));
    end

si_stream(tifPath, channel, opts.chunkFrames, @feed);
fprintf('\n');

out = acc;
out.nFrames    = numel(acc.t);
out.dzSerial   = acc.dz(:);
out.dxSerial   = acc.dx(:);
out.dySerial   = acc.dy(:);
out.residSerial= acc.resid(:);
if all(isfinite(acc.t)) && numel(acc.t) > 2 && range(acc.t) > 0
    out.fpsSource = 'TIFF frameTimestamps_sec';
else
    out.fpsSource = 'SI header scanFrameRate';
    out.t = (0:out.nFrames-1) / fpsHdr;
end
end

function fps = pick_fps(t, fpsHdr)
% NEVER hardcode or default the frame rate. Timestamps first, header second.
if numel(t) > 2 && all(isfinite(t)) && range(t) > 0
    fps = (numel(t)-1) / range(t);
else
    fps = fpsHdr;
end
assert(isfinite(fps) && fps > 0, 'could not establish a frame rate');
end

%% ========================== SI TIFF I/O ==========================
function [t, nCh, chIdx, fps] = si_header(path, channel)
%SI_HEADER  channelSave / frame rate out of the ScanImage Software tag.
tf = Tiff(path,'r');
c = onCleanup(@() close(tf));
sw = '';
try, sw = tf.getTag('Software'); catch, end
if isempty(sw)
    try, sw = tf.getTag('ImageDescription'); catch, end
end
chSave = parse_vec(sw, 'SI\.hChannels\.channelSave');
if isempty(chSave), chSave = channel; end
nCh   = numel(chSave);
chIdx = find(chSave == channel, 1);
assert(~isempty(chIdx), ...
    ['channel %d was not saved in %s (saved: %s).\n' ...
     'The pages interleave channels PER FRAME, so picking the wrong index ' ...
     'silently averages the OTHER label instead of erroring.'], ...
    channel, path, mat2str(chSave));
fps = parse_num(sw, 'SI\.hRoiManager\.scanFrameRate');
if isempty(fps) || ~isfinite(fps), fps = NaN; end
t = [];
end

function si_stream(path, channel, chunkFrames, feed)
%SI_STREAM  Walk a logged SI TIFF once, handing the caller chunks of frames.
%
% Pages are walked with nextDirectory. NEVER setDirectory(t,k) in a loop --
% it is quadratic in the page count and a 2 min series is thousands of pages.
[~, nCh, chIdx] = si_header(path, channel);
tf = Tiff(path,'r');
c  = onCleanup(@() close(tf));

buf = [];  tbuf = [];  nIn = 0;  page = 0;
while true
    page = page + 1;
    isWanted = (mod(page-1, nCh) + 1) == chIdx;
    if isWanted
        im = tf.read();
        ts = page_time(tf);
        if isempty(buf)
            buf = zeros([size(im) chunkFrames], 'single');
            tbuf = nan(1, chunkFrames);
        end
        nIn = nIn + 1;
        buf(:,:,nIn) = single(im); %#ok<AGROW>
        tbuf(nIn) = ts;
        if nIn == chunkFrames
            feed(buf, tbuf);
            nIn = 0;
        end
    end
    if tf.lastDirectory(), break; end
    tf.nextDirectory();
end
if nIn > 0
    feed(buf(:,:,1:nIn), tbuf(1:nIn));
end
end

function F = si_read_all(path, channel)
%SI_READ_ALL  Read every frame of one channel from a short logged TIFF.
assert(exist(path,'file') == 2, 'missing file: %s', path);
out = {};
    function grab(B, ~)
        out{end+1} = B; %#ok<AGROW>
    end
si_stream(path, channel, 64, @grab);
assert(~isempty(out), 'no frames read from %s', path);
F = cat(3, out{:});
end

function ts = page_time(tf)
% frameTimestamps_sec lives in the PER-PAGE ImageDescription, not the header.
ts = NaN;
try
    d = tf.getTag('ImageDescription');
    v = parse_num(d, 'frameTimestamps_sec');
    if ~isempty(v), ts = v; end
catch
end
end

function v = parse_num(s, key)
v = [];
if isempty(s), return; end
m = regexp(s, [key '\s*=\s*([-\d\.eE+]+)'], 'tokens', 'once');
if ~isempty(m), v = str2double(m{1}); end
end

function v = parse_vec(s, key)
v = [];
if isempty(s), return; end
m = regexp(s, [key '\s*=\s*\[([^\]]*)\]'], 'tokens', 'once');
if ~isempty(m)
    v = str2double(regexp(m{1}, '[-\d\.eE+]+', 'match'));
    v = v(isfinite(v));
    return
end
m = regexp(s, [key '\s*=\s*([-\d\.eE+]+)'], 'tokens', 'once');
if ~isempty(m), v = str2double(m{1}); end
end

%% ==================== BREATH-EXTREME FRAME DIFFERENCE ====================
function dObs = phase_extreme_diff(tifPath, channel, phi, good, phiHi, phiLo, opts, E)
%PHASE_EXTREME_DIFF  Mean frame at one breath phase minus mean frame at the other.
%
% Frame-level phase is the CIRCULAR mean of its strips' phases -- a plain mean
% of angles is wrong across the +-pi wrap and would scatter frames from one
% extreme into the other, diluting exactly the difference we are measuring.
S = E.S;
nStrip = numel(phi);
N = floor(nStrip / S);
P = reshape(phi(1:N*S), S, N);
G = reshape(good(1:N*S), S, N);
P(~G) = NaN;
phiFrame = angle(mean(exp(1i*P), 1, 'omitnan'));

w = opts.extremeHalfWidth;
selHi = abs(angdiff_local(phiFrame, phiHi)) <= w;
selLo = abs(angdiff_local(phiFrame, phiLo)) <= w;
assert(nnz(selHi) > 10 && nnz(selLo) > 10, ...
    'too few frames at the phase extremes (%d hi, %d lo)', nnz(selHi), nnz(selLo));

accHi = [];  accLo = [];  nHi = 0;  nLo = 0;  base = 0;
    function feed(F, ~)
        n = size(F,3);
        idx = base + (1:n);
        base = base + n;
        idx = idx(idx <= N);
        if isempty(idx), return; end
        F = F(:,:,1:numel(idx));
        h = selHi(idx);  l = selLo(idx);
        if isempty(accHi)
            accHi = zeros(size(F,1), size(F,2));
            accLo = accHi;
        end
        if any(h), accHi = accHi + sum(double(F(:,:,h)),3);  nHi = nHi + nnz(h); end
        if any(l), accLo = accLo + sum(double(F(:,:,l)),3);  nLo = nLo + nnz(l); end
    end
si_stream(tifPath, channel, opts.chunkFrames, @feed);

mHi = accHi / max(nHi,1);
mLo = accLo / max(nLo,1);
% Remove the global brightness difference: that is what the gain term already
% absorbs, and leaving it in would let a pure intensity fluctuation correlate
% with dR/dz wherever the reference is bright.
mHi = mHi - mean(mHi(:));
mLo = mLo - mean(mLo(:));
dObs = mHi - mLo;
end

function a = wrapToPiLocal(x)
a = mod(x + pi, 2*pi) - pi;
end

function d = angdiff_local(a, b)
d = mod(a - b + pi, 2*pi) - pi;
end

function r = corr_local(x, y)
x = double(x(:)) - mean(double(x(:)));
y = double(y(:)) - mean(double(y(:)));
r = (x.'*y) / (norm(x)*norm(y) + eps);
end

%% ========================== GATING ==========================
function [mu, keptFrac] = gate_mean(F, opts)
%GATE_MEAN  Average only the frames from one part of the breath cycle.
%
% The gate is on the MODE of the lateral displacement. A breathing trace
% dwells at end-expiration, so the density peak IS end-expiration; the median
% sits partway up the excursion and the mean is worse still. Picking the mode
% gives the sharpest, most repeatable reference plane.
n = size(F,3);
if n < 4 || ~opts.gate
    mu = mean(F,3);  keptFrac = 1;  return
end
ref = mean(F,3);
d = zeros(n,2);
for i = 1:n
    d(i,:) = dm_axial_estimator_260913('shift', ref, double(F(:,:,i)));
end
d = -d;                                  % correction -> displacement
d = d - mean(d,1);
[~, ~, V] = svd(d, 0);
proj = d * V(:,1);                       % 1-D lateral displacement

edges = linspace(min(proj), max(proj), max(4, round(sqrt(n))));
if range(proj) < 1e-6
    mu = mean(F,3);  keptFrac = 1;  return
end
cnt = histcounts(proj, edges);
[~, b] = max(cnt);
centre = (edges(b) + edges(b+1)) / 2;
half = opts.gateWidth * std(proj);
sel = abs(proj - centre) <= half;
if nnz(sel) < opts.gateMinFrames
    [~, ord] = sort(abs(proj - centre));
    sel = false(n,1);  sel(ord(1:min(opts.gateMinFrames,n))) = true;
end
mu = mean(F(:,:,sel), 3);
keptFrac = nnz(sel) / n;
end

%% ========================== PHASE ==========================
function lat = lateral_pc1(dx, dy, good)
d = [dx(:), dy(:)];
d(~good,:) = NaN;
mu = mean(d(good,:), 1);
d = d - mu;
[~,~,V] = svd(d(good,:), 0);
lat = d * V(:,1);
end

function [phi, fPk, xbp] = breath_phase(lat, t, good, opts)
%BREATH_PHASE  Bandpass the lateral trace, take the analytic-signal phase.
%
% hilbert() and the filter designers live in the Signal Processing Toolbox,
% which is not guaranteed on the rig MATLAB, so both are written out here.
x = lat(:);
x(~good) = 0;
x = x - mean(x(good));
n = numel(x);
fs = (numel(t)-1) / range(t(isfinite(t)));

X = fft(x);
f = (0:n-1)' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
keep = abs(f) >= opts.breathBand(1) & abs(f) <= opts.breathBand(2);
Xb = X;  Xb(~keep) = 0;
xbp = real(ifft(Xb));

pw = abs(X).^2;
band = f >= opts.breathBand(1) & f <= opts.breathBand(2);
[~, i] = max(pw .* band);
fPk = f(i);

phi = angle(analytic(xbp));
end

function z = analytic(x)
% Analytic signal by the FFT one-sided recipe -- hilbert() without the toolbox.
x = x(:);  n = numel(x);
X = fft(x);
h = zeros(n,1);
if mod(n,2) == 0
    h(1) = 1;  h(n/2+1) = 1;  h(2:n/2) = 2;
else
    h(1) = 1;  h(2:(n+1)/2) = 2;
end
z = ifft(X .* h);
end

function y = detrend_slow(x, t, good, fHz)
%DETREND_SLOW  Remove drift below fHz, keep everything at breath rate.
%
% Thermal drift wanders ~1 um over tens of minutes (measured 2026-09-11 as the
% scatter of the stack-to-stack offsets). That is a real displacement but it is
% not breathing, and leaving it in inflates every peak-to-peak number.
y = x(:);
n = numel(y);
fs = (numel(t)-1) / range(t(isfinite(t)));
yi = y;  yi(~good) = interp1(find(good), y(good), find(~good), 'linear', 'extrap');
Y = fft(yi);
f = (0:n-1)' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
Y(abs(f) < fHz) = 0;
y = real(ifft(Y));
y(~good) = NaN;
end

function [c, m, sem, nb] = phase_bin(phi, v, nBin)
edges = linspace(-pi, pi, nBin+1);
c = (edges(1:end-1) + edges(2:end)) / 2;
m = nan(1,nBin);  sem = nan(1,nBin);  nb = zeros(1,nBin);
for k = 1:nBin
    s = phi >= edges(k) & phi < edges(k+1) & isfinite(v);
    nb(k) = nnz(s);
    if nb(k) > 0
        m(k) = mean(v(s));
        sem(k) = std(v(s)) / sqrt(nb(k));
    end
end
end

%% ========================== FIGURE ==========================
function make_figure(R, outDir, opts)
% NVIDIA's OpenGL path drops ~1% of exportgraphics PNGs as blank 9,223-byte
% stubs on this machine. Force software rendering for the export.
try, opengl('software'); catch, end

hf = figure('Color','w','Position',[80 80 1400 900],'Visible','off');

% --- dz trace ---
ax1 = subplot(3,2,1:2);
t = R.tStrip;  g = R.good;
plot(ax1, t(g), R.dzDetrended(g), '-', 'Color',[0.15 0.35 0.75], 'LineWidth',0.5);
xlabel(ax1,'time (s)'); ylabel(ax1,'axial displacement (\mum)');
title(ax1, sprintf('%s   |   dz(t), drift removed above %.2f Hz', ...
      strrep(R.meta.stamp,'_','\_'), opts.driftHz));
xlim(ax1, [t(find(g,1)) min(t(find(g,1))+20, t(find(g,1,'last')))]);
grid(ax1,'on');

% --- lateral, for comparison ---
ax2 = subplot(3,2,3);
plot(ax2, t(g), R.phase.latBP(g), '-', 'Color',[0.6 0.3 0.1], 'LineWidth',0.5);
xlabel(ax2,'time (s)'); ylabel(ax2,'lateral PC1 (px)');
title(ax2, sprintf('lateral, %.2f Hz', R.phase.fBreath));
xlim(ax2, xlim(ax1)); grid(ax2,'on');

% --- dz vs phase ---
ax3 = subplot(3,2,4);
errorbar(ax3, rad2deg(R.phase.binC), R.phase.binM, R.phase.binSem, ...
         'o-','Color',[0.15 0.35 0.75],'MarkerFaceColor',[0.15 0.35 0.75],'LineWidth',1.2);
hold(ax3,'on');
pp = linspace(-pi,pi,200);
plot(ax3, rad2deg(pp), mean(R.phase.binM,'omitnan') + R.phase.amp*cos(pp - R.phase.mu), ...
     'k--','LineWidth',1);
xlabel(ax3,'breath phase (deg)'); ylabel(ax3,'dz (\mum)');
title(ax3, sprintf('p2p %.2f \\mum, amp %.2f \\mum, p = %.4f', ...
      R.phase.p2p, R.phase.amp, R.phase.p));
xlim(ax3,[-180 180]); grid(ax3,'on');

% --- null ---
ax4 = subplot(3,2,5);
histogram(ax4, R.phase.nullAmp, 30, 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none');
hold(ax4,'on');
yl = ylim(ax4);
plot(ax4, [R.phase.amp R.phase.amp], yl, 'r-','LineWidth',2);
xlabel(ax4,'cosine amplitude (\mum)'); ylabel(ax4,'count');
title(ax4,'circular-shift null');

% --- positive control ---
ax5 = subplot(3,2,6);
if isfield(R,'inject')
    plot(ax5, R.inject.cmd, R.inject.meas, 'o', 'MarkerFaceColor',[0.15 0.35 0.75], ...
         'MarkerEdgeColor','none','MarkerSize',8);
    hold(ax5,'on');
    xl = [min(R.inject.cmd) max(R.inject.cmd)];
    plot(ax5, xl, xl, 'k:','LineWidth',1);
    plot(ax5, xl, polyval(R.inject.fit, xl), 'r-','LineWidth',1.2);
    xlabel(ax5,'commanded (\mum)'); ylabel(ax5,'measured (\mum)');
    title(ax5, sprintf('injected staircase: slope %.3f, resid %.2f \\mum', ...
          R.inject.fit(1), R.inject.residRms));
    axis(ax5,'equal'); grid(ax5,'on');
else
    text(ax5, 0.5, 0.5, 'no block D', 'HorizontalAlignment','center');
    axis(ax5,'off');
end

png = fullfile(outDir, 'axial_motion.png');
% print(), NOT exportgraphics -- the rig runs R2017b (SI2018b's MATLAB) and
% exportgraphics is R2020a+.
print(hf, png, '-dpng', '-r200');
close(hf);

% An OpenGL stub comes back as a ~9 kB blank. Verify by SIZE, never by the
% file merely existing.
d = dir(png);
if isempty(d) || d.bytes < 20000
    warning('dm_axial:exportStub', ...
        'axial_motion.png is %d bytes -- almost certainly an OpenGL blank. Re-export.', ...
        max(0, sum([d.bytes])));
else
    fprintf('wrote %s (%.0f kB)\n', png, d.bytes/1024);
end
end

%% ========================== DEFAULTS ==========================
function o = fill_defaults(o)
d = struct( ...
    'gate',          true, ...
    'gateWidth',     0.5, ...    % half-width of the gate, in SD of the lateral projection
    'gateMinFrames', 8, ...
    'doCalib',       true, ...
    'doInject',      true, ...
    'chunkFrames',   200, ...
    'nPhaseBin',     24, ...     % matches the repo's phase-tuning convention
    'nShuffle',      500, ...
    'breathBand',    [1 10], ...
    'driftHz',       0.3, ...
    'residReject',   0.99, ...
    'injSettleS',    1.0, ...
    'extremeHalfWidth', pi/6, ...   % phase window for the 4b signature test
    'est',           struct());
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i}) || isempty(o.(f{i})), o.(f{i}) = d.(f{i}); end
end
end

function residRms = sim_loop(v, fsFine, rate, gx)
%SIM_LOOP  Replay a delayed reactive loop over a measured displacement trace.
%
% v      displacement, sampled at fsFine
% rate   loop rate (frames/s); the estimate is the MEAN over each loop period,
%        because a frame integrates over its own acquisition
% gx     velocity-extrapolation gain (0 = plain reactive)
%
% Timing, which is the whole point: the estimate for period k has its centroid
% at the MIDDLE of period k but cannot be applied until period k has finished.
% It is therefore applied during period k+1. That is one full loop period of
% delay, and it is what the correction has to overcome.
n = numel(v);
step = max(1, round(fsFine / rate));
nP = floor(n / step);
if nP < 4, residRms = NaN; return; end

vp = reshape(v(1:nP*step), step, nP);
e  = mean(vp, 1);                    % what the loop measures each period

c = zeros(1, nP);                    % correction applied during each period
for k = 3:nP
    pred = e(k-1) + gx * (e(k-1) - e(k-2));
    c(k) = pred;
end

resid = vp - repmat(c, step, 1);     % what is left of the motion
resid = resid(:, 3:end);
resid = resid(:) - mean(resid(:));
residRms = sqrt(mean(resid.^2));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function q = qtile(x, p)
x = x(isfinite(x));
if isempty(x), q = NaN; return; end
x = sort(double(x(:)));
n = numel(x);
i = max(1, min(n, round(p*(n-1)) + 1));
q = x(i);
end
