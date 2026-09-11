function ca_spike_deconv_bench_260722()
%CA_SPIKE_DECONV_BENCH_260722  Does calibrated deconvolution alone fix detection?
%
% Step 2 of the calibrated-deconvolution spike detector. Benchmarks three
% detection statistics against the supervised labels on all labelled sessions:
%
%   raw    -- dF/F itself, thresholded in robust noise units   (the control:
%             this is what calcium_spike_gui.m does today, minus the fact that
%             its threshold is in absolute dF/F rather than sigma)
%   mf     -- matched filter with the kernel measured in step 1
%             (ca_spike_kernel_fit_260722.m). Handles the 2-4 frame rise that
%             AR(1) cannot represent.
%   oasis  -- OASIS AR(1) at the per-group g measured in step 1.
%
% Scored two ways, per session and pooled:
%
%   recall       -- fraction of high-confidence labels (>= SNR_MIN sigma)
%                   recovered within TOL_FR frames
%   neg rate     -- the SAME detector run on the sign-flipped trace. Calcium is
%                   one-sided, so anything it finds there is noise or artifact.
%                   Counted only OUTSIDE a window after each labelled event,
%                   because the 20 s sliding-median baseline rebounds after a
%                   real transient and digs a genuine negative bowl -- without
%                   that mask the null is biased upward.
%
% The question this answers: at 95% recall, how much does the temporal model
% alone cut the false-positive rate? If it lands near zero, the spatial veto
% (step 3, requires reading the MC TIFFs) is optional polish rather than a
% requirement.
%
% Read-only w.r.t. session folders. Writes to <ROOT_DIR>\_spike_deconv_bench_260722\

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = 'D:\Ventral_surface_summary';
KERNEL_MAT = fullfile(ROOT_DIR, '_spike_kernel_260722', 'spike_kernel_fit.mat');
OUT_SUB    = '_spike_deconv_bench_260722';

SNR_MIN    = 4;                 % labels at/above this are the recall target
TOL_FR     = 4;                 % frames: detection-to-label match tolerance
MIN_DIST_S = 0.2;               % s: refractory for every detector (fair comparison)
THETA      = 1.5:0.25:8;        % threshold sweep, in robust sigma of each statistic
RECALL_TGT = 0.95;              % secondary readout: threshold reaching this recall
FAR_TGT    = [1.0 0.1];         % PRIMARY readout: recall at these false-alarm
                                % rates (detections/ROI-min on the flipped
                                % trace). Comparing methods at equal threshold
                                % is meaningless -- each statistic has its own
                                % noise scale -- so compare at equal FAR.
EXCL_POST  = 2.5;               % s after a label to ignore when counting negatives
EXCL_PRE   = 0.2;               % s before a label likewise

DIAG_THETA = 3.0;               % sigma: threshold used for the alignment self-check
DIAG_WIN   = 60;                % frames: search window for that check

KER_LEN_TAU       = 4;          % matched-filter kernel length, in tau
KER_LEN_TAU_SHORT = 1.5;        % ...and the short variant ('mfs')

ONLY_LABELLED_ROIS = true;      % Score ONLY ROIs carrying >= 1 label.
                                % ca_spike_data.mat cannot distinguish "inspected
                                % and silent" from "never inspected": both give
                                % ifSpike = false. Sessions like Vglut2/1124/IO
                                % have 190-222 ROIs but ~0.6 labels/ROI-min, so
                                % most ROIs were never curated. Counting them
                                % makes real events in uninspected ROIs look
                                % free -- they inflate no recall and trip no
                                % null -- which is how raw dF/F scored 117x the
                                % label rate while appearing to be clean.

DO_OASIS   = true;              % set false to skip the Python round-trip entirely
PYTHON_EXE = fullfile(getenv('USERPROFILE'), '.conda','envs','oasis','python.exe');
FALLBACK_FPS = 30;
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

outDir = fullfile(ROOT_DIR, OUT_SUB);
if ~isfolder(outDir), mkdir(outDir); end

assert(isfile(KERNEL_MAT), ['Missing step-1 output: %s\n' ...
    'Run ca_spike_kernel_fit_260722.m first.'], KERNEL_MAT);
Kf = load(KERNEL_MAT, 'K'); K = Kf.K;
fprintf('Kernel from step 1: %s\n', strjoin(fieldnames(K)', ', '));

d = dir(fullfile(ROOT_DIR, '**', 'ca_spike_data.mat'));
assert(~isempty(d), 'No ca_spike_data.mat under %s', ROOT_DIR);

METHODS  = {'raw','mf','mfs','oasis'};
if ~DO_OASIS, METHODS = METHODS(1:3); end
nM = numel(METHODS); nT = numel(THETA);

R = struct('session',{},'group',{},'nROI',{},'nROIscored',{},'nLab',{}, ...
           'nLabHi',{},'minutes',{},'hits',{},'hitsD',{},'nPos',{}, ...
           'nNeg',{},'negMin',{},'lag',{});
oasisDead = false;

for i = 1:numel(d)
    fo  = d(i).folder;
    rel = fo(numel(ROOT_DIR)+2:end);
    grp = strtok(rel, filesep);
    if strcmp(grp,'Vglut2_test'), grp = 'Vglut2'; end
    fprintf('\n[%2d/%2d] %s\n', i, numel(d), rel);

    dh = dir(fullfile(fo,'*_dFF.mat'));
    if isempty(dh), fprintf(2,'  no _dFF.mat -- skip\n'); continue; end
    L = load(fullfile(dh(1).folder,dh(1).name), 'dFF','F_roi');
    if ~isfield(L,'dFF'), fprintf(2,'  no dFF var -- skip\n'); continue; end
    dFF = L.dFF; [T,N] = size(dFF);

    S = load(fullfile(fo,'ca_spike_data.mat'),'roi_spikes');
    rs = S.roi_spikes;
    if numel(rs) ~= N, fprintf(2,'  ROI mismatch -- skip\n'); continue; end

    fps  = detect_session_fps(fo, FALLBACK_FPS);
    mdFr = max(1, round(MIN_DIST_S*fps));
    sig  = median(abs(diff(dFF,1,1)),1) * 1.4826 / sqrt(2);   % as in step 1

    if ~isfield(K, grp), fprintf(2,'  no kernel for %s -- skip\n', grp); continue; end
    kk  = K.(grp);
    ker = build_kernel(kk.rise_ms/1000, kk.tau_s, fps, KER_LEN_TAU);
    [~, kLag] = max(ker);  kLag = kLag - 1;    % onset -> peak, in frames
    % Short variant: a 4-tau kernel spans 2-3.6 s here, so in dense recordings
    % neighbouring events smear into one hump that no refractory can split.
    % The discriminative part is the rise plus early decay, so test that too.
    kerS = build_kernel(kk.rise_ms/1000, kk.tau_s, fps, KER_LEN_TAU_SHORT);
    [~, kLagS] = max(kerS); kLagS = kLagS - 1;
    gAR = exp(-(1/fps)/kk.tau_s);

    % ---- OASIS: one call for the session (positive), one for the mirror ----
    % Deconvolve dF/F, not raw F: its baseline is already ~0 from the sliding
    % median, so b=0 is a fair assumption and the sign-flipped null becomes the
    % same operation the other two methods get (just -dFF), instead of the
    % mirrored-fluorescence workaround.
    oaPos = []; oaNeg = []; snPos = []; snNeg = [];
    if DO_OASIS && ~oasisDead
        Fp = dFF;
        Fn = -dFF;
        try
            [oaPos, snPos] = run_oasis(Fp, gAR, PYTHON_EXE, repoRoot, sprintf('%03dp',i));
            [oaNeg, snNeg] = run_oasis(Fn, gAR, PYTHON_EXE, repoRoot, sprintf('%03dn',i));
            fprintf('  OASIS ok (g=%.3f)\n', gAR);
        catch ME
            fprintf(2,'  OASIS failed (%s) -- dropping OASIS for the rest of the run\n', ...
                ME.message);
            oasisDead = true; oaPos = []; oaNeg = []; snPos = []; snNeg = [];
        end
    end

    hits = zeros(nM,nT); hitsD = zeros(nM,nT);
    nPos = zeros(nM,nT); nNeg = zeros(nM,nT);
    nLab = 0; nLabHi = 0; okMin = 0; nScored = 0;
    lagAcc = repmat({[]}, nM, 1);

    for k = 1:N
        x   = dFF(:,k);
        idx = rs(k).spike_idx; idx = idx(idx>=1 & idx<=T);
        if ONLY_LABELLED_ROIS && isempty(idx), continue; end
        hi  = idx(x(idx)/sig(k) >= SNR_MIN);
        nLab = nLab + numel(idx); nLabHi = nLabHi + numel(hi);

        % frames where a negative detection would just be baseline rebound
        excl = false(T,1);
        for e = 1:numel(idx)
            a = max(1, idx(e)-round(EXCL_PRE*fps));
            b = min(T, idx(e)+round(EXCL_POST*fps));
            excl(a:b) = true;
        end
        okMin   = okMin + sum(~excl)/fps/60;
        nScored = nScored + 1;

        for m = 1:nM
            scl = [];   % explicit scale; [] = use robust_sigma of the statistic
            switch METHODS{m}
                case 'raw'
                    sp = x;                  sq = -x;                  lag = 0;
                case 'mf'
                    sp = matched(x, ker);    sq = matched(-x, ker);    lag = kLag;
                case 'mfs'
                    sp = matched(x, kerS);   sq = matched(-x, kerS);   lag = kLagS;
                case 'oasis'
                    if isempty(oaPos), continue; end
                    % The spike train marks the transient ONSET, so it needs the
                    % same onset->peak shift as the matched filter. And it is
                    % mostly zeros, which makes MAD degenerate (median = 0), so
                    % scale by the trace noise sn instead -- also the natural
                    % unit, since OASIS amplitudes are in dF/F.
                    sp = oaPos(:,k);         sq = oaNeg(:,k);          lag = kLag;
                    scl = [snPos(min(k,end)), snNeg(min(k,end))];
            end
            if isempty(scl)
                zp = sp / robust_sigma(sp);
                zq = sq / robust_sigma(sq);
            else
                zp = sp / max(scl(1), eps);
                zq = sq / max(scl(2), eps);
            end
            zn = zq;

            % --- alignment self-check ------------------------------------
            % Signed offset from each label to its nearest detection, over a
            % window much wider than TOL_FR. If the statistic is correctly
            % aligned this is ~0; a systematic offset means the lag correction
            % is wrong, which would masquerade as terrible recall.
            dk = findpeaks_idx(zp, DIAG_THETA, mdFr) + lag;
            dk = dk(dk>=1 & dk<=T);
            if ~isempty(hi) && ~isempty(dk)
                Dd = hi(:) - dk(:)';
                [~, wj] = min(abs(Dd), [], 2);
                off = -Dd(sub2ind(size(Dd), (1:numel(hi))', wj));
                lagAcc{m} = [lagAcc{m}; off(abs(off) <= DIAG_WIN)];
            end

            for ti = 1:nT
                pk = findpeaks_idx(zp, THETA(ti), mdFr) + lag;
                pk = pk(pk>=1 & pk<=T);
                nPos(m,ti) = nPos(m,ti) + numel(pk);
                if ~isempty(hi) && ~isempty(pk)
                    D = abs(hi(:) - pk(:)');
                    hits(m,ti)  = hits(m,ti)  + sum(min(D,[],2) <= TOL_FR); % labels found
                    hitsD(m,ti) = hitsD(m,ti) + sum(min(D,[],1) <= TOL_FR); % detections that landed on a label
                end
                nk = findpeaks_idx(zn, THETA(ti), mdFr) + lag;
                nk = nk(nk>=1 & nk<=T);
                nNeg(m,ti) = nNeg(m,ti) + sum(~excl(nk));
            end
        end
    end

    lagMed = cellfun(@(v) medOrNaN(v), lagAcc);
    R(end+1) = struct('session',rel, 'group',grp, 'nROI',N, 'nROIscored',nScored, ...
        'nLab',nLab, 'nLabHi',nLabHi, 'minutes',okMin, 'hits',hits, ...
        'hitsD',hitsD, 'nPos',nPos, 'nNeg',nNeg, ...
        'negMin',nNeg/max(okMin,eps), 'lag',lagMed); %#ok<AGROW>

    for m = 1:nM
        if all(nPos(m,:)==0)
            fprintf('  %-6s (no data)\n', METHODS{m}); continue;
        end
        rcv = hits(m,:)/max(nLabHi,1);
        prv = hitsD(m,:)./max(nPos(m,:),1);
        [f1, jf] = max(2*rcv.*prv./max(rcv+prv,eps));
        fprintf(['  %-6s bestF1=%.2f @th%.2f (P=%.2f R=%.2f)  det=%5.1f/min ' ...
                 '(lab %4.1f)  lag=%+.0f fr\n'], ...
            METHODS{m}, f1, THETA(jf), prv(jf), rcv(jf), ...
            nPos(m,jf)/max(okMin,eps), nLabHi/max(okMin,eps), lagMed(m));
    end
end

assert(~isempty(R), 'No sessions produced results.');

%% ------------------------------ pooled summary -----------------------------
Hp = zeros(nM,nT); Hd = zeros(nM,nT); Np = zeros(nM,nT); Pp = zeros(nM,nT);
nHi = 0; mins = 0; nSc = 0;
for i = 1:numel(R)
    Hp = Hp + R(i).hits;  Hd = Hd + R(i).hitsD;
    Np = Np + R(i).nNeg;  Pp = Pp + R(i).nPos;
    nHi = nHi + R(i).nLabHi; mins = mins + R(i).minutes; nSc = nSc + R(i).nROIscored;
end
recall = Hp / max(nHi,1);
precis = Hd ./ max(Pp,1);
negmin = Np / max(mins,eps);
f1all  = 2*recall.*precis ./ max(recall+precis, eps);

fprintf(['\n================= POOLED (%d sessions, %d scored ROIs, %d high-conf ' ...
         'labels, %.0f ROI-min) =================\n'], numel(R), nSc, nHi, mins);
fprintf('%-6s %8s %8s %8s %8s %10s\n','method','bestF1','th','precis','recall','det/min');
for m = 1:nM
    if all(Pp(m,:)==0), fprintf('%-6s %8s\n', METHODS{m}, '(no data)'); continue; end
    [f1, jf] = max(f1all(m,:));
    fprintf('%-6s %8.2f %8.2f %8.2f %8.2f %10.1f\n', METHODS{m}, f1, THETA(jf), ...
        precis(m,jf), recall(m,jf), Pp(m,jf)/max(mins,eps));
end
fprintf('(hand-labelled rate = %.1f events /ROI-min)\n', nHi/max(mins,eps));

fprintf('%-6s %12s %12s %12s %10s\n', 'method', ...
    sprintf('recall@%.1f',FAR_TGT(1)), sprintf('recall@%.1f',FAR_TGT(2)), ...
    'theta*(95%)','lag (fr)');
LAG   = [R.lag];                       % nM x nSessions
alive = false(1,nM);
for m = 1:nM
    alive(m) = ~(all(Np(m,:)==0) && all(Hp(m,:)==0));
    if ~alive(m)
        fprintf('%-6s %12s %12s %12s %10s\n', METHODS{m}, '--','--','--','(no data)');
        continue;
    end
    r1 = recall_at_far(THETA, recall(m,:), negmin(m,:), FAR_TGT(1));
    r2 = recall_at_far(THETA, recall(m,:), negmin(m,:), FAR_TGT(2));
    th = op_point(THETA, recall(m,:), negmin(m,:), RECALL_TGT);
    fprintf('%-6s %12.2f %12.2f %12.2f %10.1f\n', METHODS{m}, r1, r2, th, ...
        median(LAG(m,:), 'omitnan'));
end
fprintf(['\nrecall@X = fraction of >=%g-sigma labels recovered while the SAME detector\n' ...
         'fires at most X times per ROI-minute on the sign-flipped trace.\n' ...
         'lag = median frames from label to nearest detection at %.1f sigma; far from 0\n' ...
         'means a misaligned statistic, not a bad detector.\n'], SNR_MIN, DIAG_THETA);

fprintf('\n---- per group: best F1 (labelled ROIs only) ----\n');
groups = unique({R.group});
fprintf('%-8s %6s', 'group', 'nLab');
for m=1:nM, fprintf(' %14s', METHODS{m}); end, fprintf('\n');
for gi = 1:numel(groups)
    sel = strcmp({R.group}, groups{gi});
    Hg = zeros(nM,nT); Dg = zeros(nM,nT); Pg = zeros(nM,nT); hg = 0;
    for i = find(sel)
        Hg = Hg + R(i).hits;  Dg = Dg + R(i).hitsD;
        Pg = Pg + R(i).nPos;  hg = hg + R(i).nLabHi;
    end
    fprintf('%-8s %6d', groups{gi}, hg);
    for m = 1:nM
        if ~alive(m), fprintf(' %14s', '--'); continue; end
        rcg = Hg(m,:)/max(hg,1);
        prg = Dg(m,:)./max(Pg(m,:),1);
        fprintf(' %14.2f', max(2*rcg.*prg./max(rcg+prg,eps)));
    end
    fprintf('\n');
end

%% ------------------------------ figure -------------------------------------
f = figure('Color','w','Position',[60 60 1500 460]);
cols = lines(nM);
subplot(1,3,1); hold on; grid on; box on;
for m = 1:nM
    if all(Pp(m,:)==0), continue; end
    plot(recall(m,:), precis(m,:), '-o', 'Color', cols(m,:), 'LineWidth',1.6, 'MarkerSize',3);
end
xlabel(sprintf('recall of >=%g\\sigma labels', SNR_MIN)); ylabel('precision');
xlim([0 1]); ylim([0 1]);
title('precision-recall (labelled ROIs only)');
legend(METHODS,'Location','southwest','Box','off');

subplot(1,3,2); hold on; grid on; box on;
for m = 1:nM
    if all(Np(m,:)==0) && all(Hp(m,:)==0), continue; end
    plot(negmin(m,:), recall(m,:), '-o', 'Color', cols(m,:), 'LineWidth',1.6, 'MarkerSize',3);
end
yline(RECALL_TGT,'k--'); set(gca,'XScale','log');
xlabel('false det /ROI-min (sign-flipped null)');
ylabel(sprintf('recall of >=%g\\sigma labels', SNR_MIN));
title('trade-off vs the flipped-trace null');

subplot(1,3,3); hold on; grid on; box on;
for m = 1:nM
    if all(Pp(m,:)==0), continue; end
    plot(THETA, Pp(m,:)/max(mins,eps), '-', 'Color', cols(m,:), 'LineWidth',1.6);
end
yline(nHi/max(mins,eps), 'k--', 'hand-labelled rate');
set(gca,'YScale','log');
xlabel('threshold (robust \sigma)'); ylabel('detections /ROI-min');
title('detection rate vs labels');
sgtitle('Step 2: does the temporal model alone reject false events?');

saveas(f, fullfile(outDir,'deconv_bench_summary.png'));
save(fullfile(outDir,'deconv_bench.mat'), 'R','THETA','METHODS','SNR_MIN', ...
     'TOL_FR','MIN_DIST_S','EXCL_PRE','EXCL_POST','RECALL_TGT','-v7.3');
fprintf('\nSaved: %s\n', outDir);
end

%% =============================== LOCAL FUNCTIONS ===============================
function k = build_kernel(riseSec, tauSec, fps, lenTau)
% Unit-energy kernel  (1-exp(-t/tr)) .* exp(-t/td)  whose 10-90% rise matches
% the measured value. tr is found by grid search since the crossing has no
% closed form.
    t  = (0:round(lenTau*tauSec*fps))'/fps;
    trGrid = logspace(log10(1/fps/10), log10(max(riseSec,1/fps)*3), 200);
    best = trGrid(1); bestErr = inf;
    for tr = trGrid
        kk = (1-exp(-t/tr)) .* exp(-t/tauSec);
        [pk, ip] = max(kk);
        up = kk(1:ip);
        i10 = find(up <= 0.1*pk, 1, 'last');
        i90 = find(up <= 0.9*pk, 1, 'last');
        if isempty(i10) || isempty(i90), continue; end
        err = abs((i90-i10)/fps - riseSec);
        if err < bestErr, bestErr = err; best = tr; end
    end
    k = (1-exp(-t/best)) .* exp(-t/tauSec);
    k = k / norm(k);
end

function s = matched(x, k)
% Correlate with the kernel: s(t) = sum_tau k(tau+1)*x(t+tau), i.e. the score
% for a template STARTING at t. Peaks therefore land on the event onset, and
% the caller shifts by the kernel's onset-to-peak lag to reach the labelled peak.
%
% Do NOT use conv(...,'same') here: it recentres by floor(K/2), which for these
% kernels is 30-55 frames -- far beyond any sane match tolerance.
    K  = numel(k);
    yf = conv(x - median(x), flipud(k));      % 'full', length T+K-1
    s  = yf(K : K + numel(x) - 1);            % s(t) = yf(t+K-1)
end

function m = medOrNaN(v)
    if isempty(v), m = NaN; else, m = median(v); end
end

function sg = robust_sigma(x)
    sg = 1.4826 * median(abs(x - median(x)));
    if sg <= 0 || ~isfinite(sg), sg = std(x); end
    if sg <= 0 || ~isfinite(sg), sg = 1; end
end

function loc = findpeaks_idx(z, thr, minDistFr)
    [~, loc] = findpeaks(z, 'MinPeakHeight', thr, 'MinPeakDistance', minDistFr);
    loc = loc(:);
end

function [rc, th] = recall_at_far(THETA, recall, negrate, farTarget)
% Best recall achievable while the sign-flipped null fires at most farTarget
% times per ROI-minute. negrate falls with threshold, so the qualifying set is
% the high-threshold tail and its FIRST entry is the most permissive threshold
% that still meets the budget. Threshold-free by construction, so statistics
% with different noise scales can be compared directly.
    idx = find(negrate <= farTarget);
    if isempty(idx), rc = NaN; th = NaN; return; end
    j  = idx(1);
    rc = recall(j);
    th = THETA(j);
end

function [th, rc, ng] = op_point(THETA, recall, negrate, target)
% Lowest threshold whose recall still meets the target; if none does, the
% threshold with the best recall.
    ok = find(recall >= target);
    if isempty(ok), [~, j] = max(recall); else, j = ok(end); end
    th = THETA(j); rc = recall(j); ng = negrate(j);
end

function [sp, sn] = run_oasis(Y, g, pyExe, repoRoot, tag)
% OASIS AR(1) via a subprocess. MATLAB R2021b's pyenv only supports CPython
% 3.7-3.9 and the conda envs holding oasis-deconv here are 3.10, so
% helper.oasis_deconv_and_dff_AR1 (which goes through pyenv) cannot run on this
% machine. Shell out instead: -v7 .mat in, .mat out.
    script = fullfile(repoRoot, 'ca_oasis_run.py');
    assert(isfile(script), 'Missing %s', script);
    tmp = tempdir;
    fin  = fullfile(tmp, sprintf('oasis_in_%s.mat',  tag));
    fout = fullfile(tmp, sprintf('oasis_out_%s.mat', tag));
    save(fin, 'Y', 'g', '-v7');               % -v7: scipy.io cannot read v7.3
    cmd = sprintf('"%s" "%s" --in "%s" --out "%s" --g %.6f --lam 0', ...
                  pyExe, script, fin, fout, g);
    [st, msg] = system(cmd);
    if st ~= 0 || ~isfile(fout)
        delete_if(fin); delete_if(fout);
        error('ca_oasis_run failed (status %d): %s', st, strtrim(msg));
    end
    O  = load(fout, 'S', 'sn');
    sp = O.S;
    sn = O.sn(:)';                            % per-ROI noise, for scaling
    delete_if(fin); delete_if(fout);
end

function delete_if(f)
    if isfile(f)
        try
            delete(f);
        catch
        end
    end
end
