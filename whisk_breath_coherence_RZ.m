% whisk_breath_coherence_RZ.m
%
% Continuous multitaper COHERENCE (Chronux) between RAW band-passed whisker and
% breathing traces (NOT Hilbert envelopes). Split out from the xcorr analysis
% (whisk_imbalance_vs_envelope_breath_xcorr_RZ.m) -- shares the same loading,
% cropping and whisk-epoch mask, but only does coherence here.
%
% Two scopes:
%   1) whole cropped overlap trace
%   2) whisk-epoch-only samples
% Each input is tiled into fixed COH_WIN_SEC windows = TRIALS (coherencyc
% trialave=1), so DOF = K * nWindows. Frequencies below the multitaper half-
% bandwidth W = TW/COH_WIN_SEC are unreliable and are omitted.
%
% Assumes the MODIFIED coherencyc:
%   output 1 = complex coherency Z12 ; output 2 = |C| ; output 7 = f ; output 8 = confC
%
% Pairs:  L vs R | L vs breath | R vs breath | W=0.5(L+R) vs breath
% For every |C| spectrum there is a paired dPhase spectrum.

clear; clc;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW    = 400;
LIK_THR = 0.6;

% --- filtering ---
WHISK_BP  = [5 50];       % Hz, whisker band
BREATH_BP = [1 30];       % Hz, breathing band
INVERT_BREATH = false;    % loaded breath already inhale-up

EXPOSURE_OFFSET_S = 0.36324e-3;   % thermal centroid later than Basler

% --- whisking epoch gating (for whisk-only scope; mask uses GLOBAL whisking) ---
USE_EPOCHS = true;
ENV_THR    = 5;           % deg, threshold on smoothed global whisking envelope
ENV_SMOOTH = 0.1;        % s
MIN_DUR    = 2;           % s, min whisk-epoch duration
MERGE_GAP  = 0.1;         % s

% --- coherence settings ---
DO_WHOLE     = true;      % whole-trace coherence
DO_WHISK     = true;      % whisk-epoch-only coherence
COH_P        = 0.01;      % p-value for coherence confidence threshold
COH_FPASS    = [0 40];
COH_FGRID    = (0:0.25:80)';
% whisk-only scope (short data): smaller TW -> finer floor, fewer tapers
COH_TW       = 2;  COH_K       = 3;  COH_WIN_SEC       = 2.0;  % floor = TW/win = 1.00 Hz
% whole-trace scope (~80s): bigger window + TW=3/K=5 tapers for DOF
COH_TW_WHOLE = 3;  COH_K_WHOLE = 5;  COH_WIN_SEC_WHOLE = 3.0;  % floor = 3/4  = 0.75 Hz
MIN_COH_S    = 3.0;       % s, minimum whisking data (>= 2 windows of COH_WIN_SEC)

% --- exclusions ---
EXCLUDE = "5840027";
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};

SAVE_FIG = false;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

fprintf('\n=== coherencyc path ===\n'); which coherencyc -all
fprintf(['p=%.3g | whisk: TW=%g K=%d win=%.1fs (floor %.2f Hz) | ' ...
    'whole: TW=%g K=%d win=%.1fs (floor %.2f Hz)\n\n'], COH_P, ...
    COH_TW, COH_K, COH_WIN_SEC, COH_TW/COH_WIN_SEC, ...
    COH_TW_WHOLE, COH_K_WHOLE, COH_WIN_SEC_WHOLE, COH_TW_WHOLE/COH_WIN_SEC_WHOLE);

ord = 4;
[b_wb, a_wb] = butter(3,   WHISK_BP/(fpsW/2),  'bandpass');
[b_br, a_br] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

cparams = struct('Fs',fpsW,'tapers',[COH_TW COH_K],'pad',0, ...
                 'fpass',COH_FPASS,'err',[2 COH_P],'trialave',1);
cparamsWhole = cparams; cparamsWhole.tapers = [COH_TW_WHOLE COH_K_WHOLE];
cohWinSamp      = round(COH_WIN_SEC*fpsW);        cohFmin      = COH_TW       / COH_WIN_SEC;
cohWinSampWhole = round(COH_WIN_SEC_WHOLE*fpsW);  cohFminWhole = COH_TW_WHOLE / COH_WIN_SEC_WHOLE;

fGrid = COH_FGRID(:); nF = numel(fGrid);

% --- accumulators ---
C_whole_LR=zeros(nF,0); Z_whole_LR=zeros(nF,0); C_whole_LB=zeros(nF,0); Z_whole_LB=zeros(nF,0);
C_whole_RB=zeros(nF,0); Z_whole_RB=zeros(nF,0); C_whole_WB=zeros(nF,0); Z_whole_WB=zeros(nF,0);
conf_whole=[];
C_whisk_LR=zeros(nF,0); Z_whisk_LR=zeros(nF,0); C_whisk_LB=zeros(nF,0); Z_whisk_LB=zeros(nF,0);
C_whisk_RB=zeros(nF,0); Z_whisk_RB=zeros(nF,0); C_whisk_WB=zeros(nF,0); Z_whisk_WB=zeros(nF,0);
conf_whisk=[];
nSess = 0; includedSessions = {};

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

fprintf('=== SESSION LOOP ===\n');
% ============================ SESSION LOOP ============================
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};

    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE))
        fprintf('EXCLUDING animal %s n%d\n', animal, kRun); continue;
    end
    if is_excluded_sess(animal, kRun, EXCLUDE_SESS)
        fprintf('EXCLUDING session %s n%d\n', animal, kRun); continue;
    end

    % -------------------- whisker DLC + raw band-passed --------------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), LIK_THR);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)' / fpsW;
    La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0 - mean(La0,'omitnan');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0 - mean(Ra0,'omitnan');

    LbpFull = filtfilt(b_wb,a_wb, La0);   RbpFull = filtfilt(b_wb,a_wb, Ra0);
    whiskAmpFull = 0.5*(abs(hilbert(LbpFull)) + abs(hilbert(RbpFull)));  % for the MASK only

    % -------------------- breath --------------------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
        if ~isfile(Pn.breath), warning('no breath %s n%d',animal,kRun); continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb + EXPOSURE_OFFSET_S;
        brwFull = interp1(tB, br, tW, 'linear', NaN);
    catch ME
        warning('breath load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    % -------------------- crop to shared overlap --------------------
    keep = isfinite(brwFull);
    if nnz(keep) <= cohWinSamp, warning('overlap too short %s n%d',animal,kRun); continue; end
    Lbp = LbpFull(keep);  Rbp = RbpFull(keep);  whiskAmp = whiskAmpFull(keep);
    brw = brwFull(keep);  if INVERT_BREATH, brw = -brw; end
    Bbp = filtfilt(b_br,a_br, brw);                 % raw band-passed breath
    N = numel(Bbp);

    % -------------------- whisking mask (smoothed GLOBAL whisking) --------------------
    if USE_EPOCHS
        env = movmean(whiskAmp, max(1, round(ENV_SMOOTH*fpsW)));
        ep  = bool_to_epochs(env > ENV_THR, (0:N-1)'/fpsW, MIN_DUR, MERGE_GAP);
        m = false(N,1);
        for q = 1:size(ep,1)
            s_i = max(1, round(ep(q,1)*fpsW)+1); e_i = min(N, round(ep(q,2)*fpsW)+1);
            m(s_i:e_i) = true;
        end
    else
        m = true(N,1);
    end
    nWhiskSec = nnz(m)/fpsW;

    % ============================ COHERENCE ============================
    nWhole = 0; nWhiskC = 0;
    if DO_WHOLE
        try
            OUT = run_coherence_set(Lbp, Rbp, Bbp, cparamsWhole, fGrid, cohWinSampWhole);
            C_whole_LR(:,end+1)=OUT.CmagLR; Z_whole_LR(:,end+1)=OUT.ZLR; %#ok<SAGROW>
            C_whole_LB(:,end+1)=OUT.CmagLB; Z_whole_LB(:,end+1)=OUT.ZLB; %#ok<SAGROW>
            C_whole_RB(:,end+1)=OUT.CmagRB; Z_whole_RB(:,end+1)=OUT.ZRB; %#ok<SAGROW>
            C_whole_WB(:,end+1)=OUT.CmagWB; Z_whole_WB(:,end+1)=OUT.ZWB; %#ok<SAGROW>
            conf_whole(end+1)=OUT.confC; %#ok<SAGROW>
            nWhole = N/fpsW;
        catch ME
            warning('WHOLE coherence failed %s n%d: %s', animal, kRun, ME.message);
        end
    end
    if DO_WHISK
        try
            % fixed-length windows tiled WITHIN each whisking epoch -- no window
            % crosses an epoch boundary (avoids artificial jumps from concatenation)
            [XLe, XRe, XBe] = make_epoch_windows(Lbp, Rbp, Bbp, m, cohWinSamp);
            if size(XLe,2) >= 2
                OUT = coh_pairs_from_windows(XLe, XRe, XBe, cparams, fGrid);
                C_whisk_LR(:,end+1)=OUT.CmagLR; Z_whisk_LR(:,end+1)=OUT.ZLR; %#ok<SAGROW>
                C_whisk_LB(:,end+1)=OUT.CmagLB; Z_whisk_LB(:,end+1)=OUT.ZLB; %#ok<SAGROW>
                C_whisk_RB(:,end+1)=OUT.CmagRB; Z_whisk_RB(:,end+1)=OUT.ZRB; %#ok<SAGROW>
                C_whisk_WB(:,end+1)=OUT.CmagWB; Z_whisk_WB(:,end+1)=OUT.ZWB; %#ok<SAGROW>
                conf_whisk(end+1)=OUT.confC; %#ok<SAGROW>
                nWhiskC = size(XLe,2)*COH_WIN_SEC;
            else
                warning('skip whisk-only coherence %s n%d: <2 epoch windows', animal, kRun);
            end
        catch ME
            warning('WHISK coherence failed %s n%d: %s', animal, kRun, ME.message);
        end
    end

    nSess = nSess + 1;
    includedSessions{end+1,1} = sprintf('%s n%d', animal, kRun); %#ok<SAGROW>
    fprintf('%s n%d: whisk=%.1fs | whole-coh %.1fs | whisk-coh %.1fs\n', ...
        animal, kRun, nWhiskSec, nWhole, nWhiskC);
end

assert(nSess > 0, 'no usable sessions');
fprintf('\n=== %d sessions ===\n', nSess);
for i = 1:numel(includedSessions), fprintf('  %2d. %s\n', i, includedSessions{i}); end

% ============================ FIGURES ============================
common = struct('fGrid',fGrid,'fMin',cohFmin,'WHISK_BP',WHISK_BP,'BREATH_BP',BREATH_BP, ...
    'SAVE_FIG',SAVE_FIG,'repoRoot',repoRoot,'COH_TW',COH_TW,'COH_K',COH_K,'COH_P',COH_P);

if DO_WHOLE && ~isempty(C_whole_LR)
    D = common; D.scopeName = 'whole cropped trace'; D.fMin = cohFminWhole;  % longer window -> lower floor
    D.COH_TW = COH_TW_WHOLE; D.COH_K = COH_K_WHOLE;
    D.confBySession = conf_whole;
    D.confLine = median(conf_whole,'omitnan'); D.nSess = size(C_whole_LR,2);
    D.CmagLR=C_whole_LR; D.ZLR=Z_whole_LR; D.CmagLB=C_whole_LB; D.ZLB=Z_whole_LB;
    D.CmagRB=C_whole_RB; D.ZRB=Z_whole_RB; D.CmagWB=C_whole_WB; D.ZWB=Z_whole_WB;
    fprintf('\nWHOLE TRACE coherence: %d sessions (confC~%.2f)\n', D.nSess, D.confLine);
    coherence_figures(D);
end
if DO_WHISK && ~isempty(C_whisk_LR)
    D = common; D.scopeName = 'whisk epochs, windowed trials';
    D.confBySession = conf_whisk;
    D.confLine = median(conf_whisk,'omitnan'); D.nSess = size(C_whisk_LR,2);
    D.CmagLR=C_whisk_LR; D.ZLR=Z_whisk_LR; D.CmagLB=C_whisk_LB; D.ZLB=Z_whisk_LB;
    D.CmagRB=C_whisk_RB; D.ZRB=Z_whisk_RB; D.CmagWB=C_whisk_WB; D.ZWB=Z_whisk_WB;
    fprintf('WHISK-ONLY coherence: %d sessions (confC~%.2f)\n', D.nSess, D.confLine);
    coherence_figures(D);
end

% ============================= HELPERS =============================
function OUT = run_coherence_set(Lbp, Rbp, Bbp, params, fGrid, winSamp)
    Wbp = 0.5 * (Lbp + Rbp);
    [zLR,cLR,fC,confC] = coherencyc_mod_pair(Lbp, Rbp, params, winSamp);
    [zLB,cLB,~, ~    ] = coherencyc_mod_pair(Lbp, Bbp, params, winSamp);
    [zRB,cRB,~, ~    ] = coherencyc_mod_pair(Rbp, Bbp, params, winSamp);
    [zWB,cWB,~, ~    ] = coherencyc_mod_pair(Wbp, Bbp, params, winSamp);
    OUT.CmagLR = interp1(fC, cLR, fGrid, 'linear', NaN);
    OUT.CmagLB = interp1(fC, cLB, fGrid, 'linear', NaN);
    OUT.CmagRB = interp1(fC, cRB, fGrid, 'linear', NaN);
    OUT.CmagWB = interp1(fC, cWB, fGrid, 'linear', NaN);
    OUT.ZLR = interp1(fC, zLR, fGrid, 'linear', NaN);
    OUT.ZLB = interp1(fC, zLB, fGrid, 'linear', NaN);
    OUT.ZRB = interp1(fC, zRB, fGrid, 'linear', NaN);
    OUT.ZWB = interp1(fC, zWB, fGrid, 'linear', NaN);
    OUT.confC = confC;
end

function [Z, Cmag, f, confC] = coherencyc_mod_pair(x, y, params, winSamp)
% Tiles continuous inputs into winSamp windows = TRIALS (coherencyc trialave=1).
% MODIFIED coherencyc output order: 1=Z(complex) 2=|C| 7=f 8=confC.
    x = x(:); y = y(:);
    ok = isfinite(x) & isfinite(y); x = x(ok); y = y(ok);
    nWin = floor(numel(x) / winSamp);
    if nWin < 2
        error('coherency: need >=2 windows (have %d samples, winSamp=%d)', numel(x), winSamp);
    end
    X = reshape(x(1:nWin*winSamp), winSamp, nWin);
    Y = reshape(y(1:nWin*winSamp), winSamp, nWin);
    X = X - mean(X,1);  Y = Y - mean(Y,1);
    [Z, Cmag, ~, ~, ~, ~, f, confC] = coherencyc(X, Y, params);
    Z = Z(:); Cmag = Cmag(:); f = f(:);
    if numel(Z) ~= numel(Cmag) || numel(Z) ~= numel(f)
        error('coherencyc output sizes mismatch. Check modified output order.');
    end
    if ~isreal(Cmag), warning('coherencyc output 2 not real; abs().'); Cmag = abs(Cmag); end
end

function coherence_figures(D)
% ONE figure per scope, 2x2:  left = |C|, right = dphi ;
%   row 1 = L vs R whisk ;  row 2 = W=0.5(L+R) vs breath.
    fGrid = D.fGrid; fMin = D.fMin;
    confLine = median(D.confBySession,'omitnan');   % dashed visual threshold
    below = fGrid < fMin;
    minSigSess = max(3, ceil(0.5 * D.nSess));

    rows = {
        'L vs R whisk',          D.CmagLR, D.ZLR, [0.10 0.30 0.85], [0 30]
        'W=0.5(L+R) vs breath',  D.CmagWB, D.ZWB, [0.55 0.10 0.65], [0 20]
    };

    f = figure('Color','w','Position',[60 60 1120 820], 'Name',['coherence: ' D.scopeName]);
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf(['%s    |C|: bold=mean, shaded=across-session bootstrap 95%% CI    ' ...
        'd\\phi: bold = >=%d sig sessions, dotted = exploratory    (p=%.3g, n=%d, <%.1f Hz omitted)'], ...
        D.scopeName, minSigSess, D.COH_P, D.nSess, fMin), 'FontWeight','bold','Interpreter','tex');

    for r = 1:2
        name = rows{r,1}; Cmat = rows{r,2}; Zmat = rows{r,3}; col = rows{r,4}; xl = rows{r,5};

        % ---- left: |C| (bootstrap CI band) ----
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
        Cm = Cmat; Cm(below,:) = NaN;
        plot(ax, fGrid, Cm, 'Color',[0.78 0.78 0.78], 'LineWidth',0.5);     % individual sessions
        [mC, loC, hiC] = coh_mean_boot(Cm, 1000);
        vv = isfinite(loC) & isfinite(hiC);
        if any(vv)
            fill(ax, [fGrid(vv); flipud(fGrid(vv))], [loC(vv); flipud(hiC(vv))], col, ...
                'FaceAlpha',0.18, 'EdgeColor','none');
        end
        plot(ax, fGrid, mC, 'Color',col, 'LineWidth',2.5);
        yline(ax, confLine, 'k--', sprintf('p=%.2g',D.COH_P), 'LabelHorizontalAlignment','left');
        xline(ax, fMin, ':', 'Color',[0.4 0.4 0.4]);
        xlim(ax,xl); ylim(ax,[0 1]);
        xlabel(ax,'frequency (Hz)'); ylabel(ax,'coherence |C|');
        title(ax, [name '    |C|'], 'Interpreter','none'); box(ax,'off');

        % ---- right: delta phase (bold only where enough sessions significant) ----
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
        ph = rad2deg(angle(mean(Zmat,2,'omitnan')));        % phase of complex mean
        nSig = sum(Cmat > D.confBySession(:)', 2, 'omitnan');
        showFaint = ~below & isfinite(ph);
        showBold  = ~below & nSig >= minSigSess & isfinite(ph);
        phFaint = ph; phFaint(~showFaint) = NaN;
        plot(ax, fGrid, phFaint, ':', 'Color', col + (1-col)*0.65, 'LineWidth', 1.0);
        phBold = ph; phBold(~showBold) = NaN;
        plot(ax, fGrid, phBold, '-', 'Color', col, 'LineWidth', 2.5);
        yline(ax,0,'k:'); xline(ax, fMin, ':', 'Color',[0.4 0.4 0.4]);
        xlim(ax,xl); ylim(ax,[-90 90]); yticks(ax,-90:45:90);   % +/- pi/2
        xlabel(ax,'frequency (Hz)'); ylabel(ax,'\Delta phase (deg)');
        title(ax, [name '    \Delta\phi'], 'Interpreter','tex'); box(ax,'off');
    end

    if D.SAVE_FIG
        out = fullfile(D.repoRoot, ['coherence_' regexprep(D.scopeName,'[^\w]+','_') '.png']);
        exportgraphics(f,out,'Resolution',150); fprintf('saved %s\n', out);
    end
end

function OUT = coh_pairs_from_windows(XL, XR, XB, params, fGrid)
% XL,XR,XB are [winSamp x nWin] TRIAL matrices (columns = independent windows).
% Computes the 4 coherence pairs (coherencyc trialave=1) and interps onto fGrid.
% MODIFIED coherencyc output order: 1=Z(complex) 2=|C| 7=f 8=confC.
    XW = 0.5 * (XL + XR);
    XL = XL - mean(XL,1);  XR = XR - mean(XR,1);
    XB = XB - mean(XB,1);  XW = XW - mean(XW,1);
    [zLR,cLR,~,~,~,~,fC,confC] = coherencyc(XL, XR, params);
    [zLB,cLB,~,~,~,~,~ ]       = coherencyc(XL, XB, params);
    [zRB,cRB,~,~,~,~,~ ]       = coherencyc(XR, XB, params);
    [zWB,cWB,~,~,~,~,~ ]       = coherencyc(XW, XB, params);
    fC = fC(:);
    OUT.CmagLR = interp1(fC, cLR(:), fGrid, 'linear', NaN);
    OUT.CmagLB = interp1(fC, cLB(:), fGrid, 'linear', NaN);
    OUT.CmagRB = interp1(fC, cRB(:), fGrid, 'linear', NaN);
    OUT.CmagWB = interp1(fC, cWB(:), fGrid, 'linear', NaN);
    OUT.ZLR = interp1(fC, zLR(:), fGrid, 'linear', NaN);
    OUT.ZLB = interp1(fC, zLB(:), fGrid, 'linear', NaN);
    OUT.ZRB = interp1(fC, zRB(:), fGrid, 'linear', NaN);
    OUT.ZWB = interp1(fC, zWB(:), fGrid, 'linear', NaN);
    OUT.confC = confC;
end

function [XL, XR, XB] = make_epoch_windows(Lbp, Rbp, Bbp, m, winSamp)
% Fixed-length winSamp windows tiled WITHIN each continuous whisking epoch.
% No window crosses an epoch boundary; sub-window remainders are dropped.
    runs = mask_runs(logical(m));
    XL = []; XR = []; XB = [];
    for q = 1:size(runs,1)
        i0 = runs(q,1); i1 = runs(q,2);
        Lseg = Lbp(i0:i1); Rseg = Rbp(i0:i1); Bseg = Bbp(i0:i1);
        nWin = floor(numel(Lseg) / winSamp);
        for w = 1:nWin
            a = (w-1)*winSamp + 1; b = w*winSamp;
            xL = Lseg(a:b); xR = Rseg(a:b); xB = Bseg(a:b);
            ok = isfinite(xL) & isfinite(xR) & isfinite(xB);
            if all(ok)
                XL(:,end+1) = xL(:); %#ok<AGROW>
                XR(:,end+1) = xR(:); %#ok<AGROW>
                XB(:,end+1) = xB(:); %#ok<AGROW>
            end
        end
    end
end

function eidx = mask_runs(m)
    m = logical(m(:));
    d = diff([false; m; false]);
    eidx = [find(d == 1), find(d == -1) - 1];
end

function [mC, loC, hiC] = coh_mean_boot(Cmat, nBoot)
% Across-session bootstrap CI for mean coherence magnitude. Cmat is [nFreq x nSess].
%   mC  = mean across sessions ; loC/hiC = 2.5/97.5 pctile bootstrap CI (clipped to [0,1]).
% NOTE: this is ACROSS-SESSION uncertainty of the mean, NOT Chronux within-session error.
    if nargin < 2, nBoot = 1000; end
    mC = mean(Cmat,2,'omitnan');
    nSess = size(Cmat,2);
    bootM = nan(size(Cmat,1), nBoot);
    for b = 1:nBoot
        idx = randi(nSess, [1 nSess]);          % resample sessions with replacement
        bootM(:,b) = mean(Cmat(:,idx),2,'omitnan');
    end
    loC = prctile(bootM, 2.5,  2);
    hiC = prctile(bootM, 97.5, 2);
    loC(loC < 0) = 0;   % clip to [0,1] but PRESERVE NaN (below-floor rows stay NaN,
    hiC(hiC > 1) = 1;   % so the band is not drawn there; max/min would turn NaN->0/1)
end

function plot_freq_background(ax, breathBand, whiskBand, ylims)
    patch(ax, [breathBand(1) breathBand(2) breathBand(2) breathBand(1)], ...
        [ylims(1) ylims(1) ylims(2) ylims(2)], [0.85 0.92 1.00], 'EdgeColor','none','FaceAlpha',0.35);
    patch(ax, [whiskBand(1) whiskBand(2) whiskBand(2) whiskBand(1)], ...
        [ylims(1) ylims(1) ylims(2) ylims(2)], [0.93 0.93 0.93], 'EdgeColor','none','FaceAlpha',0.30);
end

function shade_subres(ax, fMin, ylims)
    patch(ax, [0 fMin fMin 0], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
        [0.75 0.55 0.55], 'EdgeColor','none', 'FaceAlpha',0.30);
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:)); d = diff([false; a; false]);
    s = find(d == 1); e = find(d == -1) - 1; ep = [t(s) t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap, ep2(end,2) = ep(i,2);
        else, ep2(end+1,:) = ep(i,:); end %#ok<AGROW>
    end
    ep = ep2; ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function tf = is_excluded_sess(animal, kRun, excl)
    tf = false;
    for q = 1:numel(excl)
        if strcmp(char(excl{q}{1}), char(animal)) && double(excl{q}{2}) == kRun
            tf = true; return;
        end
    end
end

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
