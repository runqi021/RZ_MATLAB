% cycle_heatmap_events_260811.m
% -----------------------------------------------------------------------
%  EVENT-DENSITY HEATMAP over the breath cycle, one row per cell, for the
%  SIGNIFICANT cells (permutation p < alpha) and, separately, for the cells that
%  were tested and did NOT reach it.
%
%      x  = position in the breath cycle, 0..1, tiled TWICE (0..2) so a peak
%           sitting near the wrap point is readable instead of split in half
%      y  = one cell per row
%      c  = epc, events per cycle in that bin (NOT normalised per row -- a
%           bright row is a cell that fires a lot AND fires consistently)
%
%      LEFT  panel  ONSET-triggered
%      RIGHT panel  PEAK-triggered
%
%  SORTED BY TIME TO THE NEAREST DETECTED INSPIRATORY PEAK -- and BOTH panels
%  use that SAME order, deliberately. The two triggers are separated by the
%  onset-to-peak interval (244 ms median, but jittering breath to breath by
%  ~33 ms), so a shared sort makes the pair directly comparable: the diagonal
%  should appear in both, displaced by that interval, and any row where it does
%  not is a cell whose locking survives one landmark but not the other.
%
%  ANGLE IS A LINEAR TIME NORMALISATION -- latency divided by that cell's own
%  mean IBI. It is NOT the piecewise breath phase, which stretched inspiration
%  (~19% of the cycle) across half the axis and turned a fixed GCaMP lag into a
%  rate-dependent angle. A cell at 0.25 here fired a quarter of a cycle after
%  the trigger, in whatever milliseconds that is for that cell.
%
%  BIN COUNT IS CAPPED BY THE FRAME RATE, not by taste. The fastest breathers
%  here run 0.53 s cycles at 30 fps = 16 frames per cycle, so 16 positions is
%  all that exists for those cells and anything finer is empty resolution. 20
%  bins is a compromise set by the fastest session present; the slow sessions
%  (2.9 s cycles, 87 frames) could support far more.
%
%  Input : event_latency_260811\event_latency_data.mat
%  Output: cycle_heatmap_events_<sig|nonsig>.png / .pdf  + a combined figure
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath,'event_latency_260811','event_latency_data.mat');
outDir   = fullfile(rootPath,'event_latency_260811');
nBins    = 20;             % see header: capped by the fastest breathers
prcLim   = [0.5 99.5];     % colour clamp, over the drawn cells only
rowNorm  = false;          % true = each row scaled to its own max (timing only)
doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
D = load(inFile,'S','CELL','groups','prm');
S = D.S; CELL = D.CELL; groups = D.groups; prm = D.prm;
gi = [CELL.gi];
TRIG = {'onset','peak'};

tested = vertcat(S.tested);
sigM   = vertcat(S.sig);
anySig = any(sigM,2);
anyTst = any(tested,2);

SETS = { 'SIGNIFICANT',      find(anyTst &  anySig)
         'TESTED, NOT SIG',  find(anyTst & ~anySig) };

fprintf('%d cells | %d tested | %d significant (p<%.4g) | %d tested-not-sig\n', ...
        numel(S), nnz(anyTst), nnz(anySig), prm.alphaPerm, nnz(anyTst & ~anySig));

%% ---- build the two heatmaps per set ----
edges = linspace(0, 1, nBins+1);
ctrs  = edges(1:end-1) + diff(edges)/2;
x2    = [ctrs, ctrs + 1];               % the cycle, tiled twice

H = struct('name',{},'idx',{},'M',{},'ord',{},'sortVal',{});
for si = 1:size(SETS,1)
    idx = SETS{si,2};
    if isempty(idx), continue; end
    M = cell(1,2);
    for q = 1:2
        A = nan(numel(idx), nBins);
        for i = 1:numel(idx)
            k = idx(i);
            L = S(k).L{q};
            if isempty(L) || ~isfinite(S(k).ibi(q)), continue; end
            f = mod(L / (1000*S(k).ibi(q)), 1);         % LINEAR: latency / IBI
            A(i,:) = histcounts(f, edges) / max(S(k).nCyc(q),1);   % epc
        end
        M{q} = A;
    end

    % SORT: time to the nearest detected inspiratory PEAK, as a fraction of the
    % cycle. Cells never described (no mode) fall back to their circular median,
    % and then to the end, so nothing is silently dropped from the picture.
    sv = nan(numel(idx),1);
    for i = 1:numel(idx)
        k = idx(i);
        if     ~isnan(S(k).mode(2)), sv(i) = mod(S(k).mode(2)/(1000*S(k).ibi(2)),1);
        elseif ~isnan(S(k).med(2)),  sv(i) = mod(S(k).med(2) /(1000*S(k).ibi(2)),1);
        end
    end
    [~, ord] = sort(sv);                      % NaNs sort last in MATLAB
    H(end+1) = struct('name',SETS{si,1},'idx',idx,'M',{M},'ord',ord,'sortVal',sv); %#ok<AGROW>
    fprintf('  %-16s n=%3d | sort key spans %.3f to %.3f of a cycle\n', ...
            SETS{si,1}, numel(idx), min(sv), max(sv));
end

%% ---- colour scale: ONE scale across every panel, so brightness compares ----
allV = [];
for i = 1:numel(H), for q = 1:2, allV = [allV; H(i).M{q}(:)]; end, end %#ok<AGROW>
allV = allV(isfinite(allV));
cl = prctile(allV, prcLim);
if rowNorm, cl = [0 1]; end
fprintf('colour scale (epc): %.3f to %.3f\n', cl(1), cl(2));

%% ---- figure ----
fh = figure('Color','w','Position',[30 30 1400 950]);
tl = tiledlayout(fh, numel(H), 2, 'TileSpacing','compact','Padding','compact');
for i = 1:numel(H)
    for q = 1:2
        ax = nexttile(tl, (i-1)*2 + q);
        A = H(i).M{q}(H(i).ord, :);
        if rowNorm, A = A ./ max(A,[],2,'omitnan'); end
        imagesc(ax, x2, 1:size(A,1), [A A], 'AlphaData', ~isnan([A A]));
        set(ax,'CLim',cl,'YDir','reverse','TickDir','out','Layer','top','Color',[1 1 1]);
        colormap(ax, parula(256)); hold(ax,'on');
        xline(ax, 1, 'w:', 'LineWidth',1);                 % the wrap point
        xlim(ax,[0 2]);
        set(ax,'XTick',0:0.25:2);
        xlabel(ax, sprintf('cycle position  (latency / IBI, 0 = %s)', upper(TRIG{q})));
        if q == 1, ylabel(ax, sprintf('%s  (n=%d), sorted by time to nearest PEAK', ...
                                      H(i).name, numel(H(i).idx)));
        end
        title(ax, sprintf('%s  |  %s-triggered', H(i).name, upper(TRIG{q})), ...
              'FontWeight','bold','Interpreter','none');
        cb = colorbar(ax);
        if rowNorm, cb.Label.String = 'row-normalised'; else, cb.Label.String = 'epc'; end
    end
end

title(tl, {sprintf(['Event density over the breath cycle, one row per cell   |   %d bins   |   ' ...
                    'cycle tiled twice   |   both panels share ONE sort: time to the nearest detected PEAK'], nBins), ...
           sprintf(['SIG = permutation p < %.4g, RAW (no family-wise correction), %d shifts.   ' ...
                    'x is a LINEAR time normalisation (latency / that cell''s own IBI), NOT piecewise breath phase.   ' ...
                    'NO GCaMP lag compensation.'], prm.alphaPerm, prm.nPerm)}, ...
      'FontWeight','bold','FontSize',10,'Interpreter','none');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'cycle_heatmap_events');
    if rowNorm, base = [base '_rownorm']; end
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    T = table();
    for i = 1:numel(H)
        k = H(i).idx(H(i).ord);
        % repmat on a CHAR name gives an n x len array, and the two sets have
        % different name lengths -- vertcat then fails. Wrap in a cell.
        T = [T; table(repmat({H(i).name},numel(k),1), {CELL(k).key}', {CELL(k).label}', ...
             groups(gi(k))', H(i).sortVal(H(i).ord), (1:numel(k))', ...
             'VariableNames',{'set','cell_key','label','group','frac_to_peak','row'})]; %#ok<AGROW>
    end
    writetable(T,[base '.csv']);
    fprintf('saved %s.{png,pdf,csv}\n', base);
end

%% ---- what the rows show, in numbers ----
for i = 1:numel(H)
    sv = H(i).sortVal(isfinite(H(i).sortVal));
    if isempty(sv), continue; end
    fprintf('\n%s (n=%d): time to nearest peak, median %.3f of a cycle\n', ...
            H(i).name, numel(H(i).idx), median(sv));
    fprintf('   in 0.00-0.25 %d | 0.25-0.50 %d | 0.50-0.75 %d | 0.75-1.00 %d\n', ...
        nnz(sv<0.25), nnz(sv>=0.25&sv<0.5), nnz(sv>=0.5&sv<0.75), nnz(sv>=0.75));
end
