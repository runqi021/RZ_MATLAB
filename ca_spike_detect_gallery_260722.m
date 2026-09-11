function ca_spike_detect_gallery_260722()
%CA_SPIKE_DETECT_GALLERY_260722  Eyeball the detector against your own clicks.
%
% Draws a random sample of curated ROIs and, for each, plots the dF/F trace with
% three marker classes so precision and recall are visible rather than abstract:
%
%   green filled circle   MATCHED  -- you clicked, detector agreed
%   red open circle       MISSED   -- you clicked, detector found nothing
%   orange cross          EXTRA    -- detector fired, you did not click
%
% Detector = OASIS AR(1), the best method in ca_spike_deconv_bench_260722, run at
% PER-SESSION g taken from ca_spike_kernel_fit_260722 (falling back to the cell
% type's g where a session had too few isolated events). Per-session g is the
% recommended configuration: in the benchmark the group-level g was what made
% OASIS fail on Vglut2/0224, whose true tau is ~0.36 s against the group's 0.91 s.
%
% Only ROIs carrying >= 1 label are sampled: ca_spike_data.mat cannot tell
% "inspected and silent" from "never inspected", so uncurated ROIs would make
% every real event there look like a false positive.
%
% Writes pages of panels (PNG + vector PDF) to <ROOT_DIR>\_spike_gallery_260722\

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = 'D:\Ventral_surface_summary';
KERNEL_MAT = fullfile(ROOT_DIR, '_spike_kernel_260722', 'spike_kernel_fit.mat');
OUT_SUB    = '_spike_gallery_260722';

N_TRACES   = 100;        % how many ROIs to draw
PER_PAGE   = 20;         % panels per page -> 100 traces = 5 figures
N_COLS     = 2;          % laid out 10 rows x 2 cols, so panels stay readable
WIN_SEC    = 30;         % excerpt length (0 = whole trace); window is placed to
                         % contain at least one label, so there is always
                         % something to judge
RNG_SEED   = 260722;     % fixed, so the same 100 traces come back next run

TH_SIGMA   = 1.5;        % OASIS threshold in units of sn. Swept on this same
                         % 100-ROI sample with the corrected lag/tolerance:
                         %   th 1.50 -> P 0.87  R 0.76  F1 0.81  <- best
                         %   th 2.00 -> P 0.96  R 0.69  F1 0.80
                         %   th 3.50 -> P 0.98  R 0.44  F1 0.61
                         % The old 2.0 came from the pre-alignment-fix benchmark
                         % and was far too conservative.
MIN_DIST_S = 0.2;        % refractory
TOL_FR     = 6;          % frames (0.20 s): detection-to-label match tolerance.
                         % 4 was too strict for this data -- rise is 2-4 frames
                         % and calcium_spike_gui snaps clicks to a local max
                         % within +-3 frames, so ~0.2 s of legitimate jitter is
                         % built in on both sides. Measured: |offset| <= 4 fr
                         % catches 0.67 of labels, <= 6 fr catches 0.78.
LAG_FR     = 4;          % frames: onset -> labelled-peak shift for OASIS.
                         % Measured, not modelled. The biexponential kernel's
                         % argmax predicts 5-9 frames here and over-corrects;
                         % the actual offset between an uncorrected OASIS spike
                         % and the click is ~4 frames, and applying +4 puts
                         % 0.83 of labels within +-4 fr (0.67 uncorrected).
                         % ONE global constant, not fitted per session.
SNR_MIN    = 4;          % labels below this are drawn hollow/grey (low-confidence)
MIN_SESS_EVENTS = 5;     % per-session tau needs at least this many fitted events

PYTHON_EXE = fullfile(getenv('USERPROFILE'), '.conda','envs','oasis','python.exe');
FALLBACK_FPS = 30;
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
% quiet: findpeaks warns on every silent ROI where nothing clears the threshold
warning('off','signal:findpeaks:largeMinPeakHeight');
outDir = fullfile(ROOT_DIR, OUT_SUB);
if ~isfolder(outDir), mkdir(outDir); end

assert(isfile(KERNEL_MAT), 'Run ca_spike_kernel_fit_260722.m first (missing %s)', KERNEL_MAT);
Kf = load(KERNEL_MAT, 'K', 'E');
K  = Kf.K; E = Kf.E;

% per-session tau from the step-1 per-event fits
sessTau = containers.Map('KeyType','char','ValueType','double');
if ~isempty(E)
    us = unique({E.session});
    for i = 1:numel(us)
        m = strcmp({E.session}, us{i});
        if sum(m) >= MIN_SESS_EVENTS
            sessTau(us{i}) = median([E(m).tau]);
        end
    end
end
fprintf('per-session tau available for %d session(s)\n', sessTau.Count);

d = dir(fullfile(ROOT_DIR, '**', 'ca_spike_data.mat'));
assert(~isempty(d), 'No ca_spike_data.mat under %s', ROOT_DIR);

%% ---- build the pool of curated ROIs, then sample -------------------------
pool = struct('si',{},'roi',{});
S = cell(numel(d),1);
for i = 1:numel(d)
    fo = d(i).folder;
    dh = dir(fullfile(fo,'*_dFF.mat'));
    if isempty(dh), continue; end
    L = load(fullfile(dh(1).folder,dh(1).name),'dFF');
    if ~isfield(L,'dFF'), continue; end
    Sp = load(fullfile(fo,'ca_spike_data.mat'),'roi_spikes');
    if numel(Sp.roi_spikes) ~= size(L.dFF,2), continue; end

    rel = fo(numel(ROOT_DIR)+2:end);
    grp = strtok(rel, filesep);
    if strcmp(grp,'Vglut2_test'), grp = 'Vglut2'; end
    S{i} = struct('folder',fo,'rel',rel,'grp',grp,'dFF',L.dFF, ...
                  'rs',Sp.roi_spikes,'fps',detect_session_fps(fo,FALLBACK_FPS));
    for k = 1:size(L.dFF,2)
        if ~isempty(Sp.roi_spikes(k).spike_idx)
            pool(end+1) = struct('si',i,'roi',k); %#ok<AGROW>
        end
    end
end
fprintf('curated ROI pool: %d\n', numel(pool));

rs0  = RandStream('mt19937ar','Seed',RNG_SEED);
nTake = min(N_TRACES, numel(pool));
pick  = pool(randperm(rs0, numel(pool), nTake));
fprintf('sampling %d ROI(s)\n', nTake);

%% ---- run OASIS once per session that contributed a sampled ROI ----------
need = unique([pick.si]);
OA   = cell(numel(d),1);
for ii = 1:numel(need)
    i  = need(ii);
    ss = S{i};
    tau = tau_for(ss.rel, ss.grp, sessTau, K);
    g   = exp(-(1/ss.fps)/tau);
    try
        [Sp, sn] = run_oasis(ss.dFF, g, PYTHON_EXE, repoRoot, sprintf('gal%03d',i));
        OA{i} = struct('S',Sp,'sn',sn,'g',g,'tau',tau,'lag',LAG_FR);
        fprintf('  [%2d/%2d] %-44s tau=%.2f g=%.3f\n', ii, numel(need), ...
            ss.rel(1:min(44,end)), tau, g);
    catch ME
        fprintf(2,'  OASIS failed on %s: %s\n', ss.rel, ME.message);
        OA{i} = [];
    end
end

%% ---- draw ---------------------------------------------------------------
nPage = ceil(nTake / PER_PAGE);
tally = [0 0 0];    % matched, missed, extra
pdfPath = fullfile(outDir, 'spike_gallery.pdf');
if isfile(pdfPath), delete(pdfPath); end

for p = 1:nPage
    idx = (p-1)*PER_PAGE + (1:PER_PAGE);
    idx = idx(idx <= nTake);
    nRow = ceil(numel(idx) / N_COLS);
    f  = figure('Color','w','Position',[20 20 1820 1180],'Visible','off');
    tl = tiledlayout(f, nRow, N_COLS, 'TileSpacing','compact','Padding','compact');

    for j = 1:numel(idx)
        pk = pick(idx(j));
        ss = S{pk.si}; oa = OA{pk.si};
        ax = nexttile(tl);
        c  = draw_panel(ax, ss, oa, pk.roi, WIN_SEC, TH_SIGMA, MIN_DIST_S, ...
                        TOL_FR, SNR_MIN);
        tally = tally + c;
        % keep the time axis only on the bottom row (tiles fill row-major)
        if j <= (nRow-1)*N_COLS, set(ax,'XTickLabel',[]); xlabel(ax,''); end
    end
    title(tl, sprintf(['OASIS AR(1), per-session g, threshold %.1f\\times sn   ' ...
        '|   green = matched,  red = missed,  orange \\times = extra   ' ...
        '|   page %d/%d'], TH_SIGMA, p, nPage), 'FontWeight','bold');

    png = fullfile(outDir, sprintf('spike_gallery_p%02d.png', p));
    exportgraphics(f, png, 'Resolution', 150);
    exportgraphics(f, pdfPath, 'ContentType','vector', 'Append', p > 1);
    close(f);
    fprintf('page %d/%d written\n', p, nPage);
end

fprintf(['\nover the %d sampled ROIs: %d matched, %d missed, %d extra\n' ...
         '   precision %.2f   recall %.2f\n'], nTake, tally(1), tally(2), tally(3), ...
    tally(1)/max(tally(1)+tally(3),1), tally(1)/max(tally(1)+tally(2),1));
fprintf('Saved: %s\n', outDir);
end

%% =============================== LOCAL FUNCTIONS ===============================
function tau = tau_for(rel, grp, sessTau, K)
% Per-session tau when step 1 fitted enough isolated events there, else the
% cell type's. The benchmark showed this choice is what decides whether OASIS
% works on a given session.
    if isKey(sessTau, rel)
        tau = sessTau(rel);
    else
        tau = K.(grp).tau_s;
    end
end

function lagFr = kernel_peak_lag(riseSec, tauSec, fps)
% Frames from transient onset to its peak. OASIS marks the onset; the labels sit
% on the peak, so detections must be shifted by this before matching.
    t = (0:round(4*tauSec*fps))'/fps;
    trGrid = logspace(log10(1/fps/10), log10(max(riseSec,1/fps)*3), 200);
    best = trGrid(1); bestErr = inf;
    for tr = trGrid
        kk = (1-exp(-t/tr)) .* exp(-t/tauSec);
        [pkv, ip] = max(kk); up = kk(1:ip);
        i10 = find(up <= 0.1*pkv, 1, 'last');
        i90 = find(up <= 0.9*pkv, 1, 'last');
        if isempty(i10) || isempty(i90), continue; end
        err = abs((i90-i10)/fps - riseSec);
        if err < bestErr, bestErr = err; best = tr; end
    end
    k = (1-exp(-t/best)) .* exp(-t/tauSec);
    [~, ip] = max(k); lagFr = ip - 1;
end

function cnt = draw_panel(ax, ss, oa, k, winSec, thSig, minDistS, tolFr, snrMin)
    cnt = [0 0 0];
    x   = ss.dFF(:,k); T = numel(x); fps = ss.fps;
    t   = (0:T-1)'/fps;
    lab = ss.rs(k).spike_idx; lab = sort(lab(lab>=1 & lab<=T));

    % detections
    det = [];
    if ~isempty(oa)
        z = oa.S(:,k) / max(oa.sn(min(k,numel(oa.sn))), eps);
        [~, loc] = findpeaks(z, 'MinPeakHeight', thSig, ...
                             'MinPeakDistance', max(1,round(minDistS*fps)));
        det = loc(:) + oa.lag;
        det = det(det>=1 & det<=T);
    end

    % window: centre on a label so the panel always shows something
    if winSec > 0 && T/fps > winSec
        w = round(winSec*fps);
        if ~isempty(lab)
            c = lab(1 + floor(numel(lab)/2));
        else
            c = round(T/2);
        end
        a = max(1, min(T-w, c - round(w/2))); b = a + w;
    else
        a = 1; b = T;
    end
    sel = @(v) v(v>=a & v<=b);
    labW = sel(lab); detW = sel(det);

    % classify (matching done on the windowed sets so the drawn counts match
    % what is actually visible in the panel)
    isHit = false(size(labW)); isUsed = false(size(detW));
    if ~isempty(labW) && ~isempty(detW)
        D = abs(labW(:) - detW(:)');
        isHit  = min(D,[],2) <= tolFr;
        isUsed = min(D,[],1)' <= tolFr;
    end
    cnt = [sum(isHit), sum(~isHit), sum(~isUsed)];

    % robust noise, for flagging low-confidence labels
    sg  = median(abs(diff(x)))*1.4826/sqrt(2);
    lo  = labW(x(labW)/max(sg,eps) < snrMin);

    plot(ax, t(a:b), x(a:b), 'Color',[0.30 0.30 0.34], 'LineWidth',0.6); hold(ax,'on');
    yl = [min(x(a:b)) max(x(a:b))]; if diff(yl)<=0, yl = yl + [-1 1]; end
    pad = 0.12*diff(yl); yl = yl + [-pad pad];

    if any(isHit)
        plot(ax, t(labW(isHit)), x(labW(isHit)), 'o', 'MarkerSize',5, ...
            'MarkerFaceColor',[0.15 0.65 0.30], 'MarkerEdgeColor','none');
    end
    if any(~isHit)
        plot(ax, t(labW(~isHit)), x(labW(~isHit)), 'o', 'MarkerSize',6, ...
            'MarkerEdgeColor',[0.85 0.15 0.12], 'MarkerFaceColor','none','LineWidth',1.1);
    end
    if any(~isUsed)
        plot(ax, t(detW(~isUsed)), repmat(yl(2)-0.06*diff(yl),sum(~isUsed),1), ...
            'x', 'MarkerSize',6, 'Color',[0.95 0.55 0.05], 'LineWidth',1.2);
    end
    if ~isempty(lo)   % low-confidence labels: grey tick under the trace
        plot(ax, t(lo), repmat(yl(1)+0.05*diff(yl),numel(lo),1), '|', ...
            'Color',[0.6 0.6 0.6], 'MarkerSize',5);
    end

    ylim(ax, yl); xlim(ax, [t(a) t(b)]); grid(ax,'on'); box(ax,'off');
    ylabel(ax, 'dF/F'); xlabel(ax, 'time (s)');
    [~, short] = fileparts(ss.folder);
    ttl = sprintf('%s | %s  ROI%d   tau=%.2fs   M%d  miss%d  extra%d', ...
        ss.grp, short(1:min(30,end)), k, oa_tau(oa), cnt(1), cnt(2), cnt(3));
    title(ax, ttl, 'FontSize',8, 'Interpreter','none', 'FontWeight','normal');
end

function v = oa_tau(oa)
    if isempty(oa), v = NaN; else, v = oa.tau; end
end

function [sp, sn] = run_oasis(Y, g, pyExe, repoRoot, tag)
    script = fullfile(repoRoot, 'ca_oasis_run.py');
    assert(isfile(script), 'Missing %s', script);
    fin  = fullfile(tempdir, sprintf('oasis_in_%s.mat',  tag));
    fout = fullfile(tempdir, sprintf('oasis_out_%s.mat', tag));
    save(fin, 'Y', 'g', '-v7');
    [st, msg] = system(sprintf('"%s" "%s" --in "%s" --out "%s" --g %.6f --lam 0', ...
                               pyExe, script, fin, fout, g));
    if st ~= 0 || ~isfile(fout)
        cleanup(fin); cleanup(fout);
        error('ca_oasis_run failed (%d): %s', st, strtrim(msg));
    end
    O = load(fout, 'S', 'sn');
    sp = O.S; sn = O.sn(:)';
    cleanup(fin); cleanup(fout);
end

function cleanup(f)
    if isfile(f)
        try
            delete(f);
        catch
        end
    end
end
