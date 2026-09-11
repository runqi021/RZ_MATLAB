% breath_class_dff_group_260825.m
% -----------------------------------------------------------------------
%  GROUP COMPARISON: is the Ca response bigger on a GASP than on a normal breath?
%
%  Cycle classes come from your hand labelling (breath_cycle_class_pc1.mat,
%  written by breath_failure_gui_260825.m). Ca response is measured exactly as in
%  breath_amp_vs_dff_260825.m, so the graded and grouped analyses are the same
%  measurement scored two ways:
%
%      Y_i = max(dFF over [onset_i , onset_i + respWin])
%            - mean(dFF over [onset_i - basePre , onset_i])
%
%  FIXED absolute window, not the cycle. A gasp cycle is longer than a normal
%  one, so a cycle-length window would give the gasp's peak more time to
%  accumulate and manufacture the difference this script exists to test.
%
%  INDEX BASES -- the one thing that will silently ruin this
%  cyc_*_idx are RAW pc1 frames; the Ca traces have had nDrop frames tossed off
%  the front. So the Ca-frame onset is cyc_foot_idx - nDrop. Getting this wrong
%  shifts every window by a second and quietly destroys the result rather than
%  erroring, so nDrop is read from the label file's own params, not assumed.
%
%  THREE STATISTICS, because they fail differently
%    1. Mann-Whitney U over all cycles, plus Cliff's delta as the effect size.
%       Rank based, so it does not care that dF/F is skewed.
%    2. MATCHED comparison: each gasp against the median of the normal cycles
%       within +/- matchK cycles of it, then a signed-rank test on those
%       differences. Gasps are not spread evenly through a recording and dF/F
%       drifts (bleaching, arousal), so an unmatched test partly measures WHEN
%       the gasps happened. The matched version is local in time and immune to
%       that.
%    3. A per-cell circular-shift null on Cliff's delta: the dF/F trace is
%       shifted against the cycle labels, so the null keeps both the Ca
%       autocorrelation and the clustering of gasps in time.
%
%  POOLING ACROSS A CELL'S RECORDINGS: Y is z-scored WITHIN each recording
%  before pooling, so no absolute dF/F ever crosses a recording boundary. Most
%  of these cells are single-recording anyway.
%
%  WHAT THIS CANNOT TELL YOU. A gasp moves the animal ~2.2x more than a normal
%  breath (measured, in pixels of rigid displacement), so a bigger dF/F on gasps
%  is exactly what a breath-locked optical artifact would produce. The FOV
%  common-mode column in the CSV is the archive-only check: if the effect
%  survives removing the field mean it is at least private to the cell. The
%  decisive test is still a background annulus on the registered movie.
%
%  OUTPUT  <outDir>\
%     breath_class_dff_cells.csv     per cell: n, medians, delta, p, common-mode
%     class_triggered_avg.png/.pdf   per cell, normal vs gasp triggered average
%     class_response_dist.png/.pdf   per cell, per-cycle response distributions
%
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260806'));

%% ===================== USER-EDITABLE =====================
regFile  = 'D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat';
% The hand labels live with the acquisition copy, not the archive.
srcRoot  = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
outDir   = fullfile(srcRoot, 'breath_class_dff_260825');

cellList = [283 284 285 286 287 289 290 292 296 297];

respWin_sec = 0.50;    % fixed response window after inspiration onset
basePre_sec = 0.10;    % pre-onset baseline
fallback_fps = 30;

matchK      = 5;       % matched test: normals within +/- this many cycles
minPerClass = 5;       % need this many cycles of each class in a recording
nShuffle    = 1000;
shiftMinCyc = 3;
alphaSig    = 0.05;

doSave = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
if doSave && ~isfolder(outDir), mkdir(outDir); end
rng(260825);

NORM = int8(0); GASP = int8(1);

%% ===================== REGISTRY -> observations ======================
D = load(regFile,'CELL','OBS','REC');
obsOfCell = pooled_obs_260814(D.CELL, D.OBS);
fprintf('\n=========== breath_class_dff_group_260825 ===========\n');
fprintf('cells: %s\n', mat2str(cellList));

R = struct('cell',{},'nRec',{},'nNorm',{},'nGasp',{}, ...
           'medNorm',{},'medGasp',{},'delta',{},'p_mw',{},'p_shift',{}, ...
           'p_match',{},'medMatch',{},'delta_pop',{}, ...
           'Yn',{},'Yg',{},'avgN',{},'avgG',{},'seN',{},'seG',{},'tAx',{});

for c = cellList
    if c > numel(obsOfCell) || isempty(obsOfCell{c})
        fprintf(2,'  cell %d: not in the registry -- skipped\n', c); continue;
    end
    Yn = []; Yg = []; Pn = []; Pg = [];      % responses, and FOV common mode
    An = []; Ag = [];                        % triggered-average traces
    dMatch = []; nRec = 0; tAx = []; why = {};
    for o = obsOfCell{c}(:)'
        p = regexp(D.OBS(o).label,'/','split');
        recName = strjoin(p(3:end-1),'/');  roi = str2double(p{end});
        fp = fullfile(srcRoot, recName);
        Q = prep(fp, roi, respWin_sec, basePre_sec, fallback_fps);
        if isempty(Q)
            % Distinguish the two very different reasons a recording drops out.
            % An UNLABELLED recording is a gap in the labelling, not a property
            % of the cell; reporting both the same way sent one cell's result
            % the wrong way once already.
            if ~isfile(fullfile(fp,'breath_cycle_class_pc1.mat'))
                why{end+1} = sprintf('%s: NOT LABELLED', recName); %#ok<SAGROW>
            else
                why{end+1} = sprintf('%s: no dFF, or ROI %d out of range', recName, roi); %#ok<SAGROW>
            end
            continue;
        end

        gn = Q.cls == NORM;  gg = Q.cls == GASP;
        if nnz(gn) < minPerClass || nnz(gg) < minPerClass
            why{end+1} = sprintf('%s: only %d normal / %d gasp cycles (need %d each)', ...
                                 recName, nnz(gn), nnz(gg), minPerClass); %#ok<SAGROW>
            continue;
        end
        nRec = nRec + 1;

        % z within recording, so recordings pool without absolute dF/F crossing
        z  = @(v, ref) (v - mean(ref,'omitnan')) / max(std(ref,'omitnan'), eps);
        Yn = [Yn; z(Q.Y(gn), Q.Y)];   Yg = [Yg; z(Q.Y(gg), Q.Y)];   %#ok<AGROW>
        Pn = [Pn; z(Q.P(gn), Q.P)];   Pg = [Pg; z(Q.P(gg), Q.P)];   %#ok<AGROW>
        An = [An; Q.A(gn,:)];         Ag = [Ag; Q.A(gg,:)];         %#ok<AGROW>
        tAx = Q.tAx;

        % ---- matched: each gasp vs nearby normals, in cycle index --------
        ig = find(gg); inr = find(gn);
        for k = ig(:)'
            nb = inr(abs(inr - k) <= matchK);
            if isempty(nb), continue; end
            dMatch(end+1,1) = Q.Y(k) - median(Q.Y(nb)); %#ok<SAGROW>
        end
    end
    if nRec == 0
        fprintf(2,'  cell %d SKIPPED:\n', c);
        fprintf(2,'      %s\n', why{:});
        continue;
    end

    dlt  = cliffs(Yg, Yn);
    p_mw = ranksum_safe(Yg, Yn);
    dpop = cliffs(Pg, Pn);                    % same test on the FOV common mode
    p_mt = signrank_safe(dMatch);

    % ---- circular-shift null on the effect size ----------------------
    nullD = nan(nShuffle,1);
    for s = 1:nShuffle
        Ys = []; Yns = [];
        for o = obsOfCell{c}(:)'
            p = regexp(D.OBS(o).label,'/','split');
            fp = fullfile(srcRoot, strjoin(p(3:end-1),'/'));
            Q = prep(fp, str2double(p{end}), respWin_sec, basePre_sec, fallback_fps);
            if isempty(Q), continue; end
            sh = shiftDraw(Q.T, shiftMinCyc*Q.medIBIfr);
            Ysh = shiftedY(Q, sh);
            gn = Q.cls == NORM; gg = Q.cls == GASP;
            if nnz(gn) < minPerClass || nnz(gg) < minPerClass, continue; end
            Ys  = [Ys;  Ysh(gg)];  Yns = [Yns; Ysh(gn)]; %#ok<AGROW>
        end
        if ~isempty(Ys) && ~isempty(Yns), nullD(s) = cliffs(Ys, Yns); end
    end
    nv = nullD(~isnan(nullD));
    if isempty(nv), p_sh = NaN;
    else, p_sh = (1 + nnz(abs(nv) >= abs(dlt))) / (1 + numel(nv)); end

    R(end+1) = struct('cell',c,'nRec',nRec,'nNorm',numel(Yn),'nGasp',numel(Yg), ...
        'medNorm',median(Yn),'medGasp',median(Yg),'delta',dlt,'p_mw',p_mw, ...
        'p_shift',p_sh,'p_match',p_mt,'medMatch',median(dMatch), ...
        'delta_pop',dpop,'Yn',Yn,'Yg',Yg, ...
        'avgN',mean(An,1,'omitnan'),'avgG',mean(Ag,1,'omitnan'), ...
        'seN',std(An,0,1,'omitnan')/sqrt(size(An,1)), ...
        'seG',std(Ag,0,1,'omitnan')/sqrt(size(Ag,1)),'tAx',tAx); %#ok<SAGROW>

    fprintf(['  cell %3d  %d rec  n=%3d/%3d (norm/gasp)  delta %+.3f  ' ...
             'p_MW %.3g  p_shift %.3g  p_match %.3g\n'], ...
            c, nRec, numel(Yn), numel(Yg), dlt, p_mw, p_sh, p_mt);
end
assert(~isempty(R), 'no cell had enough labelled cycles of both classes');

%% ========================= SUMMARY ===================================
dlt = [R.delta]'; pS = [R.p_shift]'; dpop = [R.delta_pop]';
fprintf('\n--------------------------------------------------\n');
fprintf('%d cells: median Cliff''s delta %+.3f  (positive = bigger on GASP)\n', ...
        numel(R), median(dlt));
fprintf('significant by circular shift (p<%.2g): %d  (%d positive, %d negative)\n', ...
        alphaSig, nnz(pS<alphaSig), nnz(pS<alphaSig & dlt>0), nnz(pS<alphaSig & dlt<0));
fprintf('matched test significant: %d of %d\n', nnz([R.p_match]' < alphaSig), numel(R));
fprintf('sign test on delta vs 0: p = %.3g\n', signrank_safe(dlt));
fprintf('\nSAME test on the FOV common mode: median delta %+.3f\n', median(dpop));
if median(abs(dpop)) >= 0.6*median(abs(dlt))
    fprintf(2,['WARNING the field common mode shows a comparable effect. That is the\n' ...
               'signature of a breath-locked optical modulation, not private cellular\n' ...
               'tuning -- gasps displace the animal ~2.2x more than normal breaths.\n']);
end

%% ========================= FIGURES ===================================
fig_avg(R, outDir, doSave);
fig_dist(R, alphaSig, outDir, doSave);

if doSave
    fid = fopen(fullfile(outDir,'breath_class_dff_cells.csv'),'w');
    fprintf(fid,['cell,n_rec,n_normal,n_gasp,med_normal_z,med_gasp_z,cliffs_delta,' ...
                 'p_mannwhitney,p_circshift,p_matched,med_matched_dff,cliffs_delta_fovmean\n']);
    for i = 1:numel(R)
        fprintf(fid,'%d,%d,%d,%d,%.4f,%.4f,%.4f,%.6g,%.6g,%.6g,%.5f,%.4f\n', ...
            R(i).cell, R(i).nRec, R(i).nNorm, R(i).nGasp, R(i).medNorm, R(i).medGasp, ...
            R(i).delta, R(i).p_mw, R(i).p_shift, R(i).p_match, R(i).medMatch, R(i).delta_pop);
    end
    fclose(fid);
    fprintf('\nSaved to %s\n', outDir);
end
fprintf('Done.\n');

% =======================================================================
% =========================== LOCAL FUNCTIONS ===========================
% =======================================================================
function Q = prep(fp, roi, respWin_sec, basePre_sec, fallback_fps)
%PREP  Per-cycle Ca response and class labels for one (recording, ROI).
persistent K
if isempty(K), K = containers.Map('KeyType','char','ValueType','any'); end
key = sprintf('%s|%d', fp, roi);
if isKey(K, key), Q = K(key); return; end     % 1000 shuffles reuse this
Q = [];
try
    cf = fullfile(fp,'breath_cycle_class_pc1.mat');
    df = dir(fullfile(fp,'*_ch1_dFF.mat'));
    if ~isfile(cf) || isempty(df), return; end
    C  = load(cf);
    Dd = load(fullfile(df(1).folder, df(1).name),'dFF');
    dff = double(Dd.dFF);
    if roi < 1 || roi > size(dff,2), return; end

    fps = detect_session_fps(fp, fallback_fps);
    nDrop = 30;
    if isfield(C,'class_params') && isfield(C.class_params,'nDrop')
        nDrop = C.class_params.nDrop;
    end
    T = size(dff,1);
    w = max(2, round(respWin_sec*fps));
    b = max(1, round(basePre_sec*fps));

    on = C.cyc_foot_idx(:) - nDrop;          % RAW pc1 frames -> Ca frame base
    ok = on - b >= 1 & on + w - 1 <= T;
    on = on(ok);  cls = int8(C.class_final(ok));
    if numel(on) < 5, return; end

    x  = dff(:,roi);
    Y  = nan(numel(on),1);  A = nan(numel(on), w);
    for i = 1:numel(on)
        seg = x(on(i) : on(i)+w-1);
        base = mean(x(on(i)-b : on(i)-1), 'omitnan');
        Y(i) = max(seg) - base;
        A(i,:) = seg - base;
    end

    % FOV common mode: same response for EVERY ROI, z-scored then averaged. An
    % optical modulation is shared by the whole field, so this carries it.
    P = zeros(numel(on),1);  nR = size(dff,2);  Z = nan(numel(on), nR);
    for j = 1:nR
        xj = dff(:,j);
        for i = 1:numel(on)
            Z(i,j) = max(xj(on(i):on(i)+w-1)) - mean(xj(on(i)-b:on(i)-1),'omitnan');
        end
    end
    Z = (Z - mean(Z,1,'omitnan')) ./ max(std(Z,0,1,'omitnan'),eps);
    P = mean(Z,2,'omitnan');

    Q = struct('Y',Y,'A',A,'P',P,'cls',cls,'on',on,'T',T,'w',w,'b',b, ...
               'x',x,'fps',fps,'tAx',(0:w-1)/fps, ...
               'medIBIfr', median(diff(on)));
    K(key) = Q; %#ok<NASGU>
catch
    Q = [];
end
end

function Y = shiftedY(Q, sh)
%SHIFTEDY  Response at every cycle with the trace circularly shifted by sh.
idx = mod(Q.on - 1 + sh, Q.T) + 1;
Y = nan(numel(idx),1);
for i = 1:numel(idx)
    s = mod((idx(i):idx(i)+Q.w-1) - 1, Q.T) + 1;
    bslice = mod((idx(i)-Q.b:idx(i)-1) - 1, Q.T) + 1;
    Y(i) = max(Q.x(s)) - mean(Q.x(bslice),'omitnan');
end
end

function s = shiftDraw(T, minFr)
if 2*minFr >= T, s = randi(T)-1; return; end
s = round(minFr + (T - 2*minFr)*rand);
end

function d = cliffs(a, b)
%CLIFFS  Cliff's delta: P(a>b) - P(a<b). +1 = a always larger.
a = a(~isnan(a)); b = b(~isnan(b));
if isempty(a) || isempty(b), d = NaN; return; end
gt = 0; lt = 0;
for i = 1:numel(a)
    gt = gt + nnz(a(i) > b);
    lt = lt + nnz(a(i) < b);
end
d = (gt - lt) / (numel(a)*numel(b));
end

function p = ranksum_safe(a, b)
a = a(~isnan(a)); b = b(~isnan(b));
if numel(a) < 3 || numel(b) < 3, p = NaN; return; end
try, p = ranksum(a, b); catch, p = NaN; end
end

function p = signrank_safe(v)
v = v(~isnan(v));
if numel(v) < 6, p = NaN; return; end
try, p = signrank(v); catch, p = NaN; end
end

function fig_avg(R, outDir, doSave)
%FIG_AVG  Inspiration-triggered dF/F, normal vs gasp, mean +/- SEM.
n = numel(R); nc = 5; nr = ceil(n/nc);
f = figure('Color','w','Units','centimeters','Position',[2 2 4.8*nc, 4.2*nr+1]);
set(f,'DefaultAxesFontSize',7.5);
for i = 1:n
    ax = subplot(nr,nc,i); hold(ax,'on');
    t = R(i).tAx;
    band(ax, t, R(i).avgN, R(i).seN, [0.45 0.45 0.45]);
    band(ax, t, R(i).avgG, R(i).seG, [0.10 0.35 0.85]);
    plot(ax, t, R(i).avgN, '-','Color',[0.45 0.45 0.45],'LineWidth',1.2);
    plot(ax, t, R(i).avgG, '-','Color',[0.10 0.35 0.85],'LineWidth',1.2);
    xline(ax, 0, 'k:');
    title(ax, sprintf('cell %d   \\delta %+.2f', R(i).cell, R(i).delta), ...
          'FontWeight','normal');
    if mod(i-1,nc)==0, ylabel(ax,'\DeltaF/F - baseline'); end
    if i > n-nc, xlabel(ax,'s from insp onset'); end
    box(ax,'off'); xlim(ax,[t(1) t(end)]);
end
sgtitle('inspiration-triggered \DeltaF/F:  grey = normal breath,  blue = GASP', 'FontSize',9);
if doSave
    exportgraphics(f, fullfile(outDir,'class_triggered_avg.png'),'Resolution',300,'BackgroundColor','white');
    exportgraphics(f, fullfile(outDir,'class_triggered_avg.pdf'),'ContentType','vector','BackgroundColor','white');
end
close(f);
end

function band(ax, t, m, se, col)
if isempty(m) || all(isnan(m)), return; end
m = m(:)'; se = se(:)';
patch(ax, [t fliplr(t)], [m-se fliplr(m+se)], col, 'EdgeColor','none','FaceAlpha',0.20);
end

function fig_dist(R, alphaSig, outDir, doSave)
%FIG_DIST  Per-cycle response distributions, normal vs gasp, per cell.
n = numel(R); nc = 5; nr = ceil(n/nc);
f = figure('Color','w','Units','centimeters','Position',[2 2 4.2*nc, 4.2*nr+1]);
set(f,'DefaultAxesFontSize',7.5);
for i = 1:n
    ax = subplot(nr,nc,i); hold(ax,'on');
    strip(ax, 1, R(i).Yn, [0.45 0.45 0.45]);
    strip(ax, 2, R(i).Yg, [0.10 0.35 0.85]);
    set(ax,'XTick',[1 2],'XTickLabel',{'norm','gasp'});
    xlim(ax,[0.4 2.6]);
    mark = ''; if R(i).p_shift < alphaSig, mark = ' *'; end
    title(ax, sprintf('cell %d%s\n\\delta %+.2f  p %.3g', R(i).cell, mark, ...
          R(i).delta, R(i).p_shift), 'FontWeight','normal');
    if mod(i-1,nc)==0, ylabel(ax,'response (z within rec)'); end
    box(ax,'off');
end
sgtitle('per-cycle \DeltaF/F response by breath class  (p = circular-shift)', 'FontSize',9);
if doSave
    exportgraphics(f, fullfile(outDir,'class_response_dist.png'),'Resolution',300,'BackgroundColor','white');
    exportgraphics(f, fullfile(outDir,'class_response_dist.pdf'),'ContentType','vector','BackgroundColor','white');
end
close(f);
end

function strip(ax, x, v, col)
v = v(~isnan(v));
if isempty(v), return; end
scatter(ax, x + 0.10*randn(numel(v),1), v, 5, col, 'filled', 'MarkerFaceAlpha',0.35);
q = prctile(v,[25 50 75]);
plot(ax, x + [-0.28 0.28], [q(2) q(2)], '-','Color',col,'LineWidth',2);
plot(ax, [x x], q([1 3]), '-','Color',col,'LineWidth',1);
end
