% whisk_imbalance_breath_xcorr_diag_persession_RZ.m
%
% Companion to whisk_imbalance_vs_envelope_breath_xcorr_RZ.m.
% Runs the SAME pipeline (fast-band L-R asymmetry & global whisking vs breath,
% masked pooled xcorr, full-trace shuffle), then for EVERY session:
%   - prints the asymmetry peak lag / peak |r| (+ significance, + global cols)
%   - draws ONE FIGURE PER SESSION:
%       left : asymmetry-vs-breath xcorr (purple) + shuffle 95% CI, peak marked;
%              global-whisking xcorr (teal) overlaid for context
%       right: longest whisking epoch -- z-scored asymmetry envelope vs z breath

clear; clc;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW    = 400;
LIK_THR = 0.6;

WHISK_BP  = [5 60];
BREATH_BP = [2 20];
INVERT_BREATH = false;
EXPOSURE_OFFSET_S = 0.36324e-3;

USE_EPOCHS = true;
ENV_THR    = 5;
ENV_SMOOTH = 0.05;
MIN_DUR    = 0.5;
MERGE_GAP  = 0.00;

MAXLAG_S    = 0.2;
NSHUF       = 200;
ALPHA       = 0.05;

NEAR_ZERO_S = 0.05;     % |lag| <= this flagged "near 0" in the table
SNIP_S      = 8;        % s, max length of the time snippet shown
EXCLUDE  = "5840027";   % animal id(s) to drop (whole animal)
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};  % {animal,run} drops

SAVE_FIG = false;       % true = also save each per-session figure as PNG
SAVE_DIR = "";          % folder for PNGs ("" = repo root)
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
if SAVE_FIG && strlength(SAVE_DIR) == 0, SAVE_DIR = string(repoRoot); end

ord = 4;
[b_wb, a_wb] = butter(3,   WHISK_BP/(fpsW/2),  'bandpass');
[b_br, a_br] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

maxLag  = round(MAXLAG_S * fpsW);
lagSec  = (-maxLag:maxLag)' / fpsW;
zsc = @(v) (v - mean(v,'omitnan')) ./ std(v,'omitnan');

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

R = struct('animal',{},'run',{},'imbLag',{},'imbR',{},'imbSgn',{},'imbSig',{}, ...
    'whLag',{},'whR',{},'nWhiskSec',{}, ...
    'r_imb',{},'ci',{},'r_whisk',{},'imbEnv',{},'brw',{},'segLong',{});

for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end
    if is_excluded_sess(animal, kRun, EXCLUDE_SESS), continue; end

    try
        M = dlc_gate_interp(pick_csv(whiskDir, ...
            sprintf('%s_whisk_n%d', animal, kRun)), LIK_THR);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)' / fpsW; N = numel(tW);

    La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0 - mean(La0,'omitnan');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0 - mean(Ra0,'omitnan');

    Lbp = filtfilt(b_wb,a_wb, La0);  Rbp = filtfilt(b_wb,a_wb, Ra0);
    whisk = 0.5 * (abs(hilbert(Lbp)) + abs(hilbert(Rbp)));
    imb   = abs(hilbert(Lbp - Rbp));

    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, ...
            sprintf('%s_nose_n%d', animal, kRun)), dataRoot);
        if ~isfile(Pn.breath), warning('no breath %s n%d', animal, kRun); continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb + EXPOSURE_OFFSET_S;
        brw = interp1(tB, br, tW, 'linear', NaN);
    catch ME
        warning('breath load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    keep = isfinite(brw);
    if nnz(keep) <= 2*maxLag + 1, continue; end
    imb = imb(keep); whisk = whisk(keep); brw = brw(keep); N = nnz(keep);
    if INVERT_BREATH, brw = -brw; end
    brw = filtfilt(b_br,a_br, brw);

    if USE_EPOCHS
        env = movmean(whisk, max(1, round(ENV_SMOOTH*fpsW)));
        ep  = bool_to_epochs(env > ENV_THR, (0:N-1)'/fpsW, MIN_DUR, MERGE_GAP);
        m = false(N,1); segLong = [1 1];
        for q = 1:size(ep,1)
            s_i = max(1, round(ep(q,1)*fpsW)+1); e_i = min(N, round(ep(q,2)*fpsW)+1);
            m(s_i:e_i) = true;
            if (e_i - s_i) > (segLong(2) - segLong(1)), segLong = [s_i e_i]; end
        end
    else
        m = true(N,1); segLong = [1 N];
    end
    if nnz(m) < 2*maxLag, continue; end

    [r_imb,   surImb] = masked_xcorr(imb,   brw, m, maxLag, NSHUF);
    [r_whisk, ~]      = masked_xcorr(whisk, brw, m, maxLag, NSHUF);
    [pi_R, pi_lag, pi_sgn] = peakAbs(r_imb,   lagSec);
    [pw_R, pw_lag, ~     ] = peakAbs(r_whisk, lagSec);
    nullPk = max(abs(surImb), [], 1);
    pi_sig = pi_R > prctile(nullPk, 100*(1-ALPHA));

    R(end+1) = struct('animal',animal,'run',kRun,'imbLag',pi_lag,'imbR',pi_R, ...
        'imbSgn',pi_sgn,'imbSig',pi_sig,'whLag',pw_lag,'whR',pw_R, ...
        'nWhiskSec',nnz(m)/fpsW,'r_imb',r_imb,'ci',prctile(surImb,[2.5 97.5],2), ...
        'r_whisk',r_whisk,'imbEnv',imb,'brw',brw,'segLong',segLong); %#ok<SAGROW>
end

assert(~isempty(R), 'no usable sessions');

% ============================== TABLE ==============================
fprintf('\n%-13s %9s %8s %5s   %9s %8s   %8s\n', ...
    'session','asymLag','asym|r|','sig','globLag','glob|r|','whisk(s)');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:numel(R)
    r = R(i);
    flag = tern(abs(r.imbLag) <= NEAR_ZERO_S, ' ', '*');   % * = off zero
    fprintf('%-8s n%-3d %+8.3f %8.2f %4s%s  %+8.3f %8.2f   %7.1f\n', ...
        r.animal, r.run, r.imbLag, r.imbR, tern(r.imbSig,'sig',' ns'), flag, ...
        r.whLag, r.whR, r.nWhiskSec);
end
fprintf('%s\n', repmat('-', 1, 70));
fprintf('(* = |asym lag| > %g s ; sig = asym peak clears %g shuffle CI)\n', ...
    NEAR_ZERO_S, 1-ALPHA);
offzero = sum(abs([R.imbLag]) > NEAR_ZERO_S);
fprintf('%d/%d sessions off-zero ; %d/%d asym peaks significant\n\n', ...
    offzero, numel(R), sum([R.imbSig]), numel(R));

% ===================== ONE FIGURE PER SESSION =====================
colImb = [0.55 0.10 0.65]; colWhisk = [0.00 0.55 0.45]; colBr = [0.20 0.40 0.80];

for i = 1:numel(R)
    r = R(i);
    f = figure('Color','w','Position',[60 60 1120 360], ...
        'Name', sprintf('%s n%d  asym-breath xcorr', r.animal, r.run));
    tlo = tiledlayout(f, 1, 2, 'TileSpacing','compact','Padding','compact');
    title(tlo, sprintf('%s  n%d    asym lag=%+.3fs  |r|=%.2f  (%s)    global |r|=%.2f @ %+.3fs', ...
        r.animal, r.run, r.imbLag, r.imbR, tern(r.imbSig,'SIG','ns'), r.whR, r.whLag), ...
        'FontWeight','bold','Interpreter','none');

    % --- left: xcorr curves ---
    axL = nexttile(tlo); hold(axL,'on'); grid(axL,'on');
    fill(axL, [lagSec; flipud(lagSec)], [r.ci(:,1); flipud(r.ci(:,2))], ...
        colImb, 'FaceAlpha',0.15, 'EdgeColor','none');
    plot(axL, lagSec, r.r_whisk, 'Color',colWhisk, 'LineWidth',1.2);
    plot(axL, lagSec, r.r_imb,   'Color',colImb,   'LineWidth',2);
    xline(axL,0,'k:'); yline(axL,0,'k:');
    plot(axL, r.imbLag, r.imbSgn, 'o', 'MarkerFaceColor',colImb, ...
        'MarkerEdgeColor','k','MarkerSize',8);
    xline(axL, r.imbLag, '--', 'Color',colImb);
    xlabel(axL,'lag (s)  (+ = whisk follows breath)'); ylabel(axL,'xcorr');
    legend(axL, {'shuffle CI','global','asymmetry'}, 'Location','northwest','Box','off');
    xlim(axL,[-MAXLAG_S MAXLAG_S]); box(axL,'off');
    title(axL, 'pooled masked xcorr');

    % --- right: time snippet over longest whisking epoch ---
    axR = nexttile(tlo); hold(axR,'on'); grid(axR,'on');
    w0 = r.segLong(1); w1 = min(r.segLong(2), w0 + round(SNIP_S*fpsW));
    tt = (w0:w1)'/fpsW;
    plot(axR, tt, zsc(r.imbEnv(w0:w1)), 'Color',colImb, 'LineWidth',1.2);
    plot(axR, tt, zsc(r.brw(w0:w1)),   'Color',colBr,  'LineWidth',1.0);
    xlabel(axR,'time (s)'); ylabel(axR,'z-score');
    legend(axR, {'asymmetry env','breath'}, 'Location','northeast','Box','off');
    title(axR, sprintf('longest whisk epoch  (%.1f s)', (w1-w0)/fpsW));
    xlim(axR,[tt(1) tt(end)]); box(axR,'off');

    if SAVE_FIG
        out = fullfile(char(SAVE_DIR), ...
            sprintf('%s_n%d_asym_breath_xcorr_diag.png', r.animal, r.run));
        exportgraphics(f, out, 'Resolution',150);
        fprintf('saved %s\n', out);
    end
end

% ============================= HELPERS =============================
function [r, sur] = masked_xcorr(x, b, m, maxLag, NSHUF)
% Exact within-whisking PEARSON cross-correlation: BOTH endpoints gated to the
% whisking mask m, so the breath endpoint is ALSO always whisking-breath (breath
% is never taken from quiet periods -> no margin leak). Pooled over all such
% pairs, FFT-fast. Matched null rotates breath WITHIN each whisking epoch, so the
% null breath stays breath-during-whisking (controls thermal-type artifacts).
%   r(tau) = (Sxb - Sx.*Sb./cnt) / sqrt((Sxx - Sx.^2./cnt).*(Sbb - Sb.^2./cnt))
% with cnt(tau)=#pairs both whisking; all S* whole-trace sliding sums (xcorr).
    x = x(:); b = b(:); m = double(m(:)); nLag = 2*maxLag + 1;
    xm = x .* m;  xxm = (x.^2) .* m;
    cnt = xcorr(m,   m, maxLag);     % # pairs with BOTH endpoints whisking
    Sx  = xcorr(xm,  m, maxLag);
    Sxx = xcorr(xxm, m, maxLag);
    r = pear_both(b .* m, m, xm, cnt, Sx, Sxx, maxLag, nLag);
    eidx = mask_runs(logical(m));
    sur = zeros(nLag, NSHUF);
    for s = 1:NSHUF
        bs = shift_within_epochs(b, eidx);            % within-epoch rotation
        sur(:,s) = pear_both(bs .* m, m, xm, cnt, Sx, Sxx, maxLag, nLag);
    end
end

function rr = pear_both(bm, m, xm, cnt, Sx, Sxx, maxLag, nLag)
    Sb  = xcorr(m,  bm,    maxLag);   % bm = b.*m (masked breath); bm.^2 = b.^2.*m
    Sbb = xcorr(m,  bm.^2, maxLag);
    Sxb = xcorr(xm, bm,    maxLag);
    num = Sxb - (Sx .* Sb) ./ cnt;
    dX  = Sxx - (Sx.^2)    ./ cnt;
    dB  = Sbb - (Sb.^2)    ./ cnt;
    rr  = zeros(nLag,1);
    ok  = cnt > 1 & dX > 0 & dB > 0;
    rr(ok) = num(ok) ./ sqrt(dX(ok) .* dB(ok));
end

function bs = shift_within_epochs(b, eidx)
    bs = b;
    for k = 1:size(eidx,1)
        i0 = eidx(k,1); i1 = eidx(k,2); L = i1 - i0 + 1;
        if L < 3, continue; end
        bs(i0:i1) = circshift(b(i0:i1), randi(L-1));   % random rotation within epoch
    end
end

function eidx = mask_runs(m)
    m = logical(m(:)); d = diff([false; m; false]);
    eidx = [find(d == 1), find(d == -1) - 1];
end

function tf = is_excluded_sess(animal, kRun, excl)
    tf = false;
    for q = 1:numel(excl)
        if strcmp(char(excl{q}{1}), char(animal)) && double(excl{q}{2}) == kRun
            tf = true; return;
        end
    end
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:)); d = diff([false; a; false]);
    s = find(d == 1); e = find(d == -1) - 1;
    ep = [t(s) t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap, ep2(end,2) = ep(i,2);
        else, ep2(end+1,:) = ep(i,:); end %#ok<AGROW>
    end
    ep = ep2; ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function [pk, lagAt, sgn] = peakAbs(r, lagSec)
    [pk, idx] = max(abs(r)); lagAt = lagSec(idx); sgn = r(idx);
end

function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));
    S = {}; key = {};
    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1} = k; %#ok<AGROW>
        S{end+1} = {tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]); S = S(o);
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
