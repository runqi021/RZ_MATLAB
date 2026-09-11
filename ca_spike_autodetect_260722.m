function ca_spike_autodetect_260722()
%CA_SPIKE_AUTODETECT_260722  Fully unsupervised calcium event detection.
%
% Uses NO hand labels anywhere -- not for the kernel, not for the threshold, not
% for validation. Everything is estimated from each trace itself, so all 868
% ROIs are usable rather than only the 362 that carry clicks.
%
% Per ROI:
%   1. sn      robust noise from the frame-to-frame difference
%   2. g       AR(1) coefficient from the autocovariance at lags >= 1 (lag 0 is
%              the only lag white noise touches, so excluding it makes this
%              noise-free). Reproduces the supervised kernel fit to ~0.1 s.
%   3. regime  rhythmic vs sparse, from spectral concentration + event rate
%   4a. SPARSE   -> OASIS AR(1) at that ROI's own g; event = transient onset
%   4b. RHYTHMIC -> band-pass around the ROI's own rhythm, segment cycles, walk
%                   back to the foot; event = inspiration onset, one per cycle
%   5. threshold chosen PER ROI by the sign-flipped null: calcium cannot go
%      down, so detections on -dF/F are false by construction. Take the most
%      permissive threshold whose downward count stays under TARGET_FDR of the
%      upward count. No labels, and it self-calibrates to each ROI's noise.
%
% Writes per session `ca_spike_auto.mat` and, at ROOT_DIR, a summary + gallery.

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR    = 'D:\Ventral_surface_summary';
OUT_SUB     = '_spike_autodetect_260722';

TARGET_FDR  = 0.05;      % allowed sign-flipped detections, as a fraction
TH_GRID     = 5:-0.1:1;  % searched high -> low; first value meeting the FDR wins
TH_FLOOR    = 1.0;       % never go below this many sn

RHYTHM_MIN  = 0.25;      % spectral concentration above this -> rhythmic regime
RATE_MIN    = 20;        % ...or events/min above this with rhythm > RHYTHM_MIN/2
ASYM_MIN    = 1.3;       % below this an ROI shows no calcium asymmetry -> skip

MIN_DIST_S  = 0.15;      % refractory for the sparse detector
FOOT_FRAC   = 0.20;      % walk back to baseline + this fraction of amplitude
GAL_N       = 100;       % ROIs drawn for the gallery
GAL_PER_PAGE= 20;
GAL_COLS    = 2;
GAL_WIN_SEC = 30;
RNG_SEED    = 260722;

PYTHON_EXE  = fullfile(getenv('USERPROFILE'),'.conda','envs','oasis','python.exe');
FALLBACK_FPS= 30;
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
warning('off','signal:findpeaks:largeMinPeakHeight');
outDir = fullfile(ROOT_DIR, OUT_SUB);
if ~isfolder(outDir), mkdir(outDir); end

d = dir(fullfile(ROOT_DIR, '**', '*_dFF.mat'));
% keep one dFF per folder
[~, ia] = unique({d.folder}, 'stable'); d = d(ia);
assert(~isempty(d), 'No *_dFF.mat under %s', ROOT_DIR);
fprintf('ca_spike_autodetect: %d session(s)\n\n', numel(d));

ALL = struct('session',{},'group',{},'roi',{},'regime',{},'g',{},'tau',{}, ...
             'sn',{},'th',{},'fdr',{},'n',{},'rate',{},'onsets',{},'fps',{},'T',{});

for i = 1:numel(d)
    fo  = d(i).folder;
    rel = fo(numel(ROOT_DIR)+2:end);
    grp = strtok(rel, filesep);
    if strcmp(grp,'Vglut2_test'), grp = 'Vglut2'; end

    L = load(fullfile(fo, d(i).name), 'dFF');
    if ~isfield(L,'dFF'), continue; end
    dFF = L.dFF; [T,N] = size(dFF);
    fps = detect_session_fps(fo, FALLBACK_FPS);
    md  = max(1, round(MIN_DIST_S*fps));

    % ---- per-ROI parameters, all label-free ----
    sn = zeros(1,N); gg = zeros(1,N); tau = nan(1,N);
    rhy = zeros(1,N); fpk = nan(1,N); asym = zeros(1,N); regime = cell(1,N);
    for k = 1:N
        x = dFF(:,k);
        sn(k) = median(abs(diff(x)))*1.4826/sqrt(2);
        if ~isfinite(sn(k)) || sn(k) <= 0, sn(k) = max(std(x),eps); end
        tau(k) = ac_tau(x, fps);
        if ~isfinite(tau(k)), tau(k) = 0.5; end
        gg(k)  = exp(-(1/fps)/tau(k));
        [fpk(k), rhy(k)] = rhythm_score(x, fps, [0.05 8]);
        z = (x-median(x))/sn(k);
        asym(k) = (numel(pk_idx(z,3,md))+1) / (numel(pk_idx(-z,3,md))+1);
        rt = numel(pk_idx(z,3,md))/(T/fps)*60;
        if asym(k) < ASYM_MIN
            regime{k} = 'skip';
        elseif rhy(k) > RHYTHM_MIN || (rt > RATE_MIN && rhy(k) > RHYTHM_MIN/2)
            regime{k} = 'rhythmic';
        else
            regime{k} = 'sparse';
        end
    end

    % ---- OASIS for the sparse ROIs, batched by quantised g ----
    isSp = strcmp(regime,'sparse');
    Spos = zeros(T,N); Sneg = zeros(T,N);
    if any(isSp)
        gq = round(gg*100)/100;
        for gv = unique(gq(isSp))
            sel = isSp & gq == gv;
            try
                Spos(:,sel) = run_oasis( dFF(:,sel), gv, PYTHON_EXE, repoRoot, sprintf('a%03dp%03d',i,round(gv*100)));
                Sneg(:,sel) = run_oasis(-dFF(:,sel), gv, PYTHON_EXE, repoRoot, sprintf('a%03dn%03d',i,round(gv*100)));
            catch ME
                fprintf(2,'  OASIS g=%.2f failed: %s\n', gv, ME.message);
                regime(sel) = {'skip'};
            end
        end
    end

    % ---- detect ----
    nDet = 0;
    for k = 1:N
        ons = []; th = NaN; fdr = NaN;
        switch regime{k}
            case 'sparse'
                zp = Spos(:,k)/sn(k); zn = Sneg(:,k)/sn(k);
                [th, fdr] = pick_threshold(zp, zn, md, TH_GRID, TARGET_FDR, TH_FLOOR);
                ons = pk_idx(zp, th, md);          % OASIS fires at the onset
            case 'rhythmic'
                [ons, th, fdr] = detect_cycles(dFF(:,k), fps, fpk(k), sn(k), ...
                                               TH_GRID, TARGET_FDR, TH_FLOOR, FOOT_FRAC);
        end
        ons = ons(ons>=1 & ons<=T);
        nDet = nDet + numel(ons);
        ALL(end+1) = struct('session',rel,'group',grp,'roi',k, ...
            'regime',regime{k}, 'g',gg(k), 'tau',tau(k), 'sn',sn(k), ...
            'th',th, 'fdr',fdr, 'n',numel(ons), 'rate',numel(ons)/(T/fps)*60, ...
            'onsets',ons, 'fps',fps, 'T',T); %#ok<AGROW>
    end

    % per-session save, alongside the data
    auto = ALL(strcmp({ALL.session}, rel));
    save(fullfile(fo,'ca_spike_auto.mat'), 'auto', 'TARGET_FDR', 'RHYTHM_MIN', '-v7.3');
    fprintf('%-50s %3d ROI (%2d rhy, %3d sp, %2d skip)  %5d events\n', ...
        rel(1:min(50,end)), N, sum(strcmp(regime,'rhythmic')), ...
        sum(strcmp(regime,'sparse')), sum(strcmp(regime,'skip')), nDet);
end

%% ------------------------------ summary ------------------------------------
reg = {ALL.regime};
fprintf('\n===================== SUMMARY =====================\n');
fprintf('%d ROIs: %d sparse, %d rhythmic, %d skipped (no calcium asymmetry)\n', ...
    numel(ALL), sum(strcmp(reg,'sparse')), sum(strcmp(reg,'rhythmic')), ...
    sum(strcmp(reg,'skip')));
act = ~strcmp(reg,'skip');
fprintf('events: %d total, %.1f /ROI-min median\n', sum([ALL.n]), median([ALL(act).rate]));
fprintf('chosen threshold (sn): 25/50/75 pct  %s\n', ...
    mat2str(round(prctile([ALL(act).th],[25 50 75]),2)));
fprintf('achieved sign-flipped FDR: 50/90/max  %s\n', ...
    mat2str(round(prctile([ALL(act).fdr],[50 90 100]),3)));
save(fullfile(outDir,'autodetect_all.mat'), 'ALL', '-v7.3');

%% ------------------------------ gallery ------------------------------------
sel = find(act);
rs0 = RandStream('mt19937ar','Seed',RNG_SEED);
pick = sel(randperm(rs0, numel(sel), min(GAL_N, numel(sel))));
draw_gallery(ALL, pick, ROOT_DIR, outDir, GAL_PER_PAGE, GAL_COLS, GAL_WIN_SEC);
fprintf('\nSaved: %s\n', outDir);
end

%% =============================== DETECTION ===============================
function [th, fdr] = pick_threshold(zp, zn, md, grid, targetFDR, floorTh)
% Most permissive threshold whose sign-flipped count stays under targetFDR.
% Scans high -> low and stops at the last threshold that still passes, so the
% result is the sensitivity ceiling consistent with the FDR budget.
    th = grid(1); fdr = 0;
    for t = grid
        np = numel(pk_idx(zp, t, md));
        nn = numel(pk_idx(zn, t, md));
        if np < 1, continue; end
        f = nn / np;
        if f > targetFDR, break; end
        th = t; fdr = f;
        if t <= floorTh, break; end
    end
    th = max(th, floorTh);
end

function [ons, th, fdr] = detect_cycles(x, fps, fPk, sn, grid, targetFDR, floorTh, footFrac)
% Rhythmic regime: one event per cycle, marked at the inspiration-side onset.
% Band-pass around the ROI's own rhythm to define cycles, then walk back on the
% RAW trace to the foot -- filtering rounds the sharp onset, which is the very
% feature being located.
    ons = []; th = NaN; fdr = NaN;
    if ~isfinite(fPk) || fPk <= 0, return; end
    bp = bandpass_simple(x, fps, [0.5 1.8]*fPk);
    md = max(1, round(0.6/fPk*fps));            % refractory = 60% of a cycle
    zb = bp / max(median(abs(diff(bp)))*1.4826/sqrt(2), eps);
    [th, fdr] = pick_threshold(zb, -zb, md, grid, targetFDR, floorTh);
    pks = pk_idx(zb, th, md);
    if isempty(pks), return; end

    ons = zeros(size(pks));
    for j = 1:numel(pks)
        p = pks(j);
        lo = 1; if j > 1, lo = pks(j-1); end
        seg = x(lo:p);
        if numel(seg) < 3, ons(j) = p; continue; end
        base = min(seg); amp = x(p) - base;
        if amp <= 0, ons(j) = p; continue; end
        % steepest rise, then walk back to the foot
        [~, is] = max(diff(seg));
        q = is;
        while q > 1 && seg(q) > base + footFrac*amp
            q = q - 1;
        end
        ons(j) = lo + q - 1;
    end
    ons = unique(ons(:));
end

function y = bandpass_simple(x, fps, band)
% Zero-phase band-pass without the Signal Processing toolbox design functions:
% difference of two moving averages, which is a clean enough band for cycle
% segmentation and cannot go unstable.
    band = max(band, 1e-3);
    wLo = max(3, round(fps/band(2)));   % keeps periods shorter than 1/band(2) out
    wHi = max(wLo+2, round(fps/band(1)));
    y = movmean(x, wLo) - movmean(x, wHi);
end

%% =============================== HELPERS ===============================
function loc = pk_idx(z, thr, minDistFr)
    [~, loc] = findpeaks(z, 'MinPeakHeight', thr, 'MinPeakDistance', minDistFr);
    loc = loc(:);
end

function tau = ac_tau(x, fps)
% AR(1) tau from the autocovariance at lags >= 1 (lag 0 carries the white noise).
    x = x - movmedian(x, round(20*fps));
    x = x - mean(x);
    L = min(max(3, round(0.6*fps)), floor(numel(x)/4));
    c = xcorr(x, L, 'biased'); c = c(L+1:end);
    c1 = c(2:end); ok = c1 > 0;
    if sum(ok) < 3, tau = NaN; return; end
    lags = (1:numel(c1))';
    p = polyfit(lags(ok), log(c1(ok)), 1);
    g = exp(p(1));
    if ~isfinite(g) || g <= 0 || g >= 1, tau = NaN; return; end
    tau = -1/(fps*log(g));
    if tau <= 0 || tau > 5, tau = NaN; end
end

function [fPk, score] = rhythm_score(x, fps, band)
    x = x - mean(x); n = numel(x);
    P = abs(fft(x.*hann(n))).^2; nf = floor(n/2)+1;
    P = P(1:nf); fr = (0:nf-1)'*(fps/n);
    m = fr >= band(1) & fr <= band(2);
    if ~any(m), fPk = NaN; score = NaN; return; end
    Pb = P; Pb(~m) = 0; [~,j] = max(Pb); fPk = fr(j);
    nb = fr >= 0.75*fPk & fr <= 1.25*fPk & m;
    score = sum(P(nb)) / max(sum(P(m)), eps);
end

function sp = run_oasis(Y, g, pyExe, repoRoot, tag)
    script = fullfile(repoRoot, 'ca_oasis_run.py');
    fin  = fullfile(tempdir, sprintf('auto_in_%s.mat',  tag));
    fout = fullfile(tempdir, sprintf('auto_out_%s.mat', tag));
    save(fin, 'Y', 'g', '-v7');
    [st, msg] = system(sprintf('"%s" "%s" --in "%s" --out "%s" --g %.6f --lam 0', ...
                               pyExe, script, fin, fout, g));
    if st ~= 0 || ~isfile(fout)
        rm(fin); rm(fout); error('ca_oasis_run failed (%d): %s', st, strtrim(msg));
    end
    O = load(fout,'S'); sp = O.S; rm(fin); rm(fout);
end

function rm(f)
    if isfile(f)
        try
            delete(f);
        catch
        end
    end
end

%% =============================== GALLERY ===============================
function draw_gallery(ALL, pick, ROOT_DIR, outDir, perPage, nCols, winSec)
    pdfPath = fullfile(outDir,'autodetect_gallery.pdf');
    if isfile(pdfPath), delete(pdfPath); end
    cache = containers.Map();
    nPage = ceil(numel(pick)/perPage);
    for p = 1:nPage
        idx = (p-1)*perPage + (1:perPage);
        idx = idx(idx <= numel(pick));
        nRow = ceil(numel(idx)/nCols);
        f = figure('Color','w','Position',[20 20 1820 1180],'Visible','off');
        tl = tiledlayout(f, nRow, nCols, 'TileSpacing','compact','Padding','compact');
        for j = 1:numel(idx)
            a = ALL(pick(idx(j)));
            key = a.session;
            if ~isKey(cache, key)
                fo = fullfile(ROOT_DIR, a.session);
                dh = dir(fullfile(fo,'*_dFF.mat'));
                Lz = load(fullfile(dh(1).folder,dh(1).name),'dFF');
                cache(key) = Lz.dFF;
            end
            D = cache(key);
            ax = nexttile(tl);
            draw_one(ax, D(:,a.roi), a, winSec);
            if j <= (nRow-1)*nCols, set(ax,'XTickLabel',[]); xlabel(ax,''); end
        end
        title(tl, sprintf(['Unsupervised detection (no hand labels) -- ' ...
            'threshold set per ROI by the sign-flipped null   |   page %d/%d'], ...
            p, nPage), 'FontWeight','bold');
        exportgraphics(f, fullfile(outDir, sprintf('autodetect_p%02d.png',p)), 'Resolution',150);
        exportgraphics(f, pdfPath, 'ContentType','vector', 'Append', p>1);
        close(f);
        fprintf('gallery page %d/%d\n', p, nPage);
    end
end

function draw_one(ax, x, a, winSec)
    T = numel(x); t = (0:T-1)'/a.fps; ons = a.onsets;
    if winSec > 0 && T/a.fps > winSec
        w = round(winSec*a.fps);
        if ~isempty(ons), c = ons(1+floor(numel(ons)/2)); else, c = round(T/2); end
        lo = max(1, min(T-w, c-round(w/2))); hi = lo + w;
    else
        lo = 1; hi = T;
    end
    o = ons(ons>=lo & ons<=hi);
    plot(ax, t(lo:hi), x(lo:hi), 'Color',[0.30 0.30 0.34], 'LineWidth',0.6); hold(ax,'on');
    yl = [min(x(lo:hi)) max(x(lo:hi))];
    if diff(yl)<=0, yl = yl+[-1 1]; end
    pad = 0.12*diff(yl); yl = yl+[-pad pad];
    if strcmp(a.regime,'rhythmic'), col = [0.85 0.35 0.05]; else, col = [0.10 0.55 0.85]; end
    if ~isempty(o)
        plot(ax, t(o), x(o), 'v', 'MarkerSize',5, 'MarkerFaceColor',col, 'MarkerEdgeColor','none');
    end
    ylim(ax,yl); xlim(ax,[t(lo) t(hi)]); grid(ax,'on'); box(ax,'off');
    ylabel(ax,'dF/F'); xlabel(ax,'time (s)');
    [~, short] = fileparts(fullfile(a.session));
    title(ax, sprintf('%s | %s ROI%d  [%s]  tau=%.2fs  th=%.1f sn  FDR=%.2f  n=%d', ...
        a.group, short(1:min(26,end)), a.roi, a.regime, a.tau, a.th, a.fdr, a.n), ...
        'FontSize',8,'Interpreter','none','FontWeight','normal');
end
