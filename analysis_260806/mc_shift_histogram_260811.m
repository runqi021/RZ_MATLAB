% mc_shift_histogram_260811.m
% -----------------------------------------------------------------------
%  How much did the brain actually MOVE, in MICRONS, across the whole ventral
%  surface archive?
%
%  This exists to bound the motion-artifact question. The breath-triggered
%  latency result (event_latency_260811) has an uncomfortable property: a
%  breathing-coupled motion artifact would predict exactly the observed answer,
%  events locked just after inspiration. Timing alone cannot separate them.
%  MAGNITUDE can bound it -- but only if the magnitude is measured correctly.
%
%  THERE ARE TWO MOTION-CORRECTION PASSES AND YOU MUST ADD THEM.
%      pass 1   raw -> *_ch1_preproc_MC.tif        options.max_shift = 25 px
%               shifts in `*_ch1_preproc_MC_output.mat`   (mcOut.shifts)
%      pass 2   *_preproc_MC -> *_preproc_MC_MC    options.max_shift =  6 px
%               shifts in `*_ch1_preproc_MC_MC_output.mat`
%  Both are RIGID (options.grid_size spans the whole frame), so the total
%  displacement of a frame is simply pass1 + pass2.
%
%  THE TRAP: `*_ch1_MCinfo.mat` carries `sx`,`sy`, which look like THE shifts and
%  are the obvious thing to grab -- but its own `mcOut.mc_path` points at
%  `*_preproc_MC_MC.tif`, so MCINFO IS PASS 2 ONLY, the RESIDUAL after the first
%  correction. Using it alone understates the motion by roughly 4x: measured over
%  this archive, the median per-recording MAX displacement is 0.89 um from
%  MCinfo alone versus 3.25 um for both passes summed, and the worst recording
%  goes from 3.75 um to 13.25 um -- a whole soma, not a fraction of a pixel.
%  This script reads BOTH and reports pass 1, pass 2 and the total separately so
%  the difference stays visible and this cannot be repeated silently.
%
%  Microns come from `pixelSize_um` in `*_meta.mat` (= 1.7778 / zoomFactor,
%  ScanImage), which is per-recording: zoom runs 1.2x to 12x here, so px->um
%  varies by 10x across the archive and a histogram in PIXELS would be
%  meaningless. That is the whole point of doing this in microns.
%
%  TWO DIFFERENT QUANTITIES, and they answer different questions:
%    DISPLACEMENT  |s - median(s)|, distance from the recording's own resting
%                  position. Slow drift over a session inflates this and is
%                  NOT what corrupts an event -- the movie is registered, and
%                  slow drift is tracked perfectly.
%    FRAME-TO-FRAME |diff(s)|, how far the brain moved between consecutive
%                  frames. THIS is the artifact-relevant number: a transient
%                  that lasts a few frames can only be faked by motion that
%                  happens on that timescale.
%  Both are plotted. The second is the one to quote.
%
%  WHAT THE SCALE BARS MEAN. Reference lines at 1 um and 10 um (a typical soma
%  is ~10-20 um across here). A frame-to-frame excursion well under a micron
%  cannot carry a cell off its mask; one approaching 10 um certainly can.
%
%  CAVEAT, STATED. These are the shifts that were DETECTED AND REMOVED, so they
%  are a lower bound on true motion: z-motion and non-rigid warping are invisible
%  to a rigid registration and do not appear here at all. This bounds the
%  in-plane rigid component only -- which is also the component already measured
%  NOT to predict events (motion-speed AUC 0.50 across 25 sessions).
%
%  Output: D:\Ventral_surface_summary\mc_shift_260811\
%
%  Runqi Zhang / 2026-08-11

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
outDir    = fullfile(rootPath,'mc_shift_260811');
scan_dirs = {'ChAT','Vglut2','Vgat','Sst','Sert'};
groups    = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};
gCol      = [0.35 0.35 0.35;    % IO     grey
             0.85 0.10 0.10;    % ChAT   red
             0.10 0.65 0.20;    % Vglut2 green
             0.10 0.30 0.85;    % Vgat   blue
             0.55 0.20 0.75;    % Sst    purple
             0.90 0.45 0.10];   % Sert   orange
applyExclusions = true;   % same coh_cfg data criteria as every other figure
refLines_um = [1 10];     % 1 um, and a typical soma diameter
nBins       = 80;
doSave      = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
excl = {};
if applyExclusions, excl = coh_cfg_260727().excludeRecordings; end

R = struct('rec',{},'grp',{},'gi',{},'umpx',{},'zoom',{},'nFr',{},'nPass',{}, ...
           'disp',{},'step',{},'d1',{},'d2',{}, ...
           'medD',{},'p95D',{},'maxD',{},'medS',{},'p95S',{},'maxS',{},'max1',{},'max2',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), continue; end
    hits = dir(fullfile(gdir,'**','*_ch1_preproc_MC_output.mat'));
    fprintf('=== %-7s %d pass-1 files\n', sname, numel(hits));
    for h = 1:numel(hits)
        fp  = hits(h).folder;
        rec = folder_basename(fp);
        if any(strcmp(rec, excl)), fprintf('  skip (excluded): %s\n', rec); continue; end
        gname = sname;
        rel   = regexprep(strrep(fp,gdir,''),'^[\\/]+','');
        parts = regexp(rel,'[\\/]','split');
        if numel(parts) >= 2 && ~isempty(regexpi(parts{2},'IO','once')), gname = 'IO'; end
        gi = find(strcmp(groups,gname),1);

        try
            S1 = shifts_from(fullfile(fp, hits(h).name));
            if isempty(S1), continue; end
            p2 = dir(fullfile(fp,'*_ch1_preproc_MC_MC_output.mat'));
            S2 = [];  nPass = 1;
            if ~isempty(p2)
                S2 = shifts_from(fullfile(p2(1).folder, p2(1).name));
                if ~isempty(S2), nPass = 2; end
            end
            if isempty(S2), S2 = zeros(size(S1)); end
            n  = min(size(S1,1), size(S2,1));
            S1 = S1(1:n,:);  S2 = S2(1:n,:);
            ST = S1 + S2;                    % both passes are RIGID, so they add

            % microns per pixel, per recording -- zoom varies 1.2x to 12x here
            [~, meta] = detect_session_fps(fp, 30);
            if ~isfield(meta,'pixelSize_um') || ~isfinite(meta.pixelSize_um) || meta.pixelSize_um <= 0
                fprintf('  skip (no pixelSize_um): %s\n', rec); continue;
            end
            umpx = meta.pixelSize_um;  zm = meta.zoomFactor;

            % displacement from the recording's OWN resting position, and the
            % frame-to-frame step. Median, not mean: a slow drift should not
            % define "centre".
            dmag = @(S) hypot(S(:,1)-median(S(:,1)), S(:,2)-median(S(:,2))) * umpx;
            d  = dmag(ST);
            s  = hypot(diff(ST(:,1)), diff(ST(:,2))) * umpx;

            R(end+1) = struct('rec',rec,'grp',gname,'gi',gi,'umpx',umpx,'zoom',zm, ...
                'nFr',n,'nPass',nPass,'disp',d,'step',s,'d1',dmag(S1),'d2',dmag(S2), ...
                'medD',median(d),'p95D',prctile(d,95),'maxD',max(d), ...
                'medS',median(s),'p95S',prctile(s,95),'maxS',max(s), ...
                'max1',max(dmag(S1)),'max2',max(dmag(S2))); %#ok<AGROW>
        catch ME
            fprintf(2,'  ERROR %s: %s\n', rec, ME.message);
        end
    end
end
assert(~isempty(R), 'no *_ch1_preproc_MC_output.mat found under %s', rootPath);

allD = vertcat(R.disp);  allS = vertcat(R.step);
fprintf('\n%d recordings (%d with both passes, %d pass-1 only) | %d frames\n', ...
    numel(R), nnz([R.nPass]==2), nnz([R.nPass]==1), numel(allD));
fprintf('px size %.3f-%.3f um (zoom %.1f-%.1fx)\n', ...
    min([R.umpx]), max([R.umpx]), min([R.zoom]), max([R.zoom]));
fprintf('\nper-recording MAX displacement (medians over recordings):\n');
fprintf('   pass1 %.2f um | pass2 (MCinfo alone) %.2f um | TOTAL %.2f um\n', ...
    median([R.max1]), median([R.max2]), median([R.maxD]));
fprintf('   worst recording TOTAL %.2f um | %d of %d recordings exceed 10 um\n', ...
    max([R.maxD]), nnz([R.maxD]>10), numel(R));
fprintf('\npooled over all frames:\n');
fprintf('  DISPLACEMENT from rest : median %.2f um | 95th %.2f | max %.2f\n', ...
    median(allD), prctile(allD,95), max(allD));
fprintf('  FRAME-TO-FRAME step    : median %.3f um | 95th %.2f | max %.2f\n', ...
    median(allS), prctile(allS,95), max(allS));
fprintf('  frames with a step > 1 um: %.2f%% | > 10 um: %.3f%%\n', ...
    100*mean(allS>1), 100*mean(allS>10));
fprintf('\n%-7s %5s %10s %10s | %10s %10s\n','group','nRec','medDisp','p95Disp','medStep','p95Step');
for g = 1:numel(groups)
    m = [R.gi]==g; if ~any(m), continue; end
    fprintf('%-7s %5d %10.2f %10.2f | %10.3f %10.2f\n', groups{g}, nnz(m), ...
        median(vertcat(R(m).disp)), prctile(vertcat(R(m).disp),95), ...
        median(vertcat(R(m).step)), prctile(vertcat(R(m).step),95));
end

%% ---- figure ----
fh = figure('Color','w','Position',[40 40 1500 900]);
tl = tiledlayout(fh,2,2,'TileSpacing','compact','Padding','compact');

% (1) TOTAL displacement, per genotype
ax = nexttile(tl,1); hold(ax,'on');
edD = linspace(0, prctile(allD,99.5), nBins);
hLg = gobjects(0); lLg = {};
for g = 1:numel(groups)
    m = [R.gi]==g; if ~any(m), continue; end
    hh = histogram(ax, vertcat(R(m).disp), edD, 'Normalization','probability', ...
                   'FaceColor',gCol(g,:), 'FaceAlpha',0.45, 'EdgeColor','none');
    hLg(end+1) = hh;  lLg{end+1} = sprintf('%s (%d rec)', groups{g}, nnz(m)); %#ok<AGROW>
end
ref_lines(ax, refLines_um);
xlabel(ax,'TOTAL displacement from rest (\mum, both passes)'); ylabel(ax,'fraction of frames');
title(ax, sprintf('displacement  |  median %.2f, 95th %.2f \\mum', median(allD), prctile(allD,95)));
set(ax,'TickDir','out'); box(ax,'off');
legend(ax, hLg, lLg, 'Box','off', 'FontSize',8);

% (2) TOTAL frame-to-frame step -- the artifact-relevant one
ax = nexttile(tl,2); hold(ax,'on');
edS = linspace(0, prctile(allS,99.5), nBins);
for g = 1:numel(groups)
    m = [R.gi]==g; if ~any(m), continue; end
    histogram(ax, vertcat(R(m).step), edS, 'Normalization','probability', ...
              'FaceColor',gCol(g,:), 'FaceAlpha',0.45, 'EdgeColor','none');
end
ref_lines(ax, refLines_um);
xlabel(ax,'TOTAL frame-to-frame step (\mum, both passes)'); ylabel(ax,'fraction of frames');
title(ax, sprintf('FRAME-TO-FRAME  |  median %.3f, 95th %.2f \\mum  (quote this one)', ...
      median(allS), prctile(allS,95)), 'FontWeight','bold');
set(ax,'TickDir','out'); box(ax,'off');

% (3) cumulative, log x -- reads off "what fraction of frames exceed X um"
ax = nexttile(tl,3); hold(ax,'on');
for g = 1:numel(groups)
    m = [R.gi]==g; if ~any(m), continue; end
    v = sort(vertcat(R(m).step)); v = v(v>0);
    if isempty(v), continue; end
    plot(ax, v, (1:numel(v))/numel(v), '-', 'Color',gCol(g,:), 'LineWidth',1.4);
end
set(ax,'XScale','log','TickDir','out'); box(ax,'off');
ref_lines(ax, refLines_um);
yline(ax,0.95,':','Color',[0.5 0.5 0.5]);
xlabel(ax,'frame-to-frame step (\mum, log)'); ylabel(ax,'cumulative fraction of frames');
title(ax,'cumulative step size');

% (4) per-recording MAX displacement, and WHY BOTH PASSES MATTER.
% Open circle = pass 2 alone (what *_MCinfo.mat gives you); filled = total.
% The line between them is the motion an MCinfo-only version would miss.
ax = nexttile(tl,4); hold(ax,'on');
for g = 1:numel(groups)
    m = find([R.gi]==g); if isempty(m), continue; end
    x = g + 0.30*(rand(numel(m),1)-0.5);
    for i = 1:numel(m)
        plot(ax, [x(i) x(i)], [R(m(i)).max2 R(m(i)).maxD], '-', ...
             'Color',[gCol(g,:) 0.35], 'LineWidth',0.8);
    end
    scatter(ax, x, [R(m).max2], 22, gCol(g,:), 'MarkerEdgeColor',gCol(g,:), 'LineWidth',0.8);
    scatter(ax, x, [R(m).maxD], 36, gCol(g,:), 'filled', 'MarkerFaceAlpha',0.8, ...
            'MarkerEdgeColor','w');
    plot(ax, g+[-0.24 0.24], median([R(m).maxD])*[1 1], '-', 'Color',gCol(g,:), 'LineWidth',2.4);
end
ref_lines(ax, refLines_um, true);
set(ax,'XTick',1:numel(groups),'XTickLabel',groups,'TickDir','out','YScale','log');
ylabel(ax,'MAX displacement (\mum)'); box(ax,'off');
title(ax,'per recording: open = pass 2 only (MCinfo), filled = BOTH passes');

title(tl, {sprintf(['Motion correction actually applied, in MICRONS   |   %d recordings, %d frames   |   ' ...
                    'px size %.2f-%.2f \\mum (zoom %.1f-%.1fx)'], numel(R), numel(allD), ...
                    min([R.umpx]), max([R.umpx]), min([R.zoom]), max([R.zoom])), ...
           ['BOTH rigid passes summed (pass1 *_preproc_MC_output + pass2 *_preproc_MC_MC_output). ' ...
            '*_MCinfo.mat is PASS 2 ONLY and understates this ~4x. Still a LOWER BOUND: z-motion is invisible to rigid registration.']}, ...
      'FontWeight','bold','FontSize',10,'Interpreter','tex');

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'mc_shift_um');
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    writetable(table({R.rec}',{R.grp}',[R.zoom]',[R.umpx]',[R.nFr]',[R.nPass]', ...
        [R.max1]',[R.max2]',[R.maxD]', ...
        [R.medD]',[R.p95D]',[R.medS]',[R.p95S]',[R.maxS]', ...
        'VariableNames',{'recording','group','zoom','um_per_px','n_frames','n_passes', ...
        'max_disp_pass1_um','max_disp_pass2_um','max_disp_total_um', ...
        'median_disp_um','p95_disp_um','median_step_um','p95_step_um','max_step_um'}), ...
        [base '.csv']);
    fprintf('\nsaved %s.{png,pdf,csv}\n', base);
end

%% ---- local ----
function S = shifts_from(matPath)
% NoRMCorre stores one struct per frame with .shifts sized [1 1 1 2] for a rigid
% run (row shift, col shift) or [gy gx 1 2] for a non-rigid grid. Flatten to
% [T x 2]; for a grid, average the patches to get the global translation, which
% is what a displacement in microns means.
S = [];
A = load(matPath,'mcOut');
if ~isfield(A,'mcOut') || ~isfield(A.mcOut,'shifts'), return; end
sh = A.mcOut.shifts;
if isempty(sh), return; end
S = zeros(numel(sh), 2);
for i = 1:numel(sh)
    v = sh(i).shifts;
    if ndims(v) == 4 && size(v,4) >= 2
        S(i,1) = mean(reshape(v(:,:,:,1),[],1));
        S(i,2) = mean(reshape(v(:,:,:,2),[],1));
    else
        v = reshape(v,1,[]);
        S(i,1) = v(1);  S(i,2) = v(min(2,numel(v)));
    end
end
end

function ref_lines(ax, um, horiz)
if nargin < 3, horiz = false; end
lbl = {'1 \mum','soma ~10 \mum'};
for i = 1:numel(um)
    if horiz, yline(ax, um(i), '--', lbl{min(i,end)}, 'Color',[0.4 0.4 0.4], ...
                    'LineWidth',1, 'LabelHorizontalAlignment','left', 'FontSize',8);
    else,     xline(ax, um(i), '--', lbl{min(i,end)}, 'Color',[0.4 0.4 0.4], ...
                    'LineWidth',1, 'LabelOrientation','horizontal', 'FontSize',8);
    end
end
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; rebuild the full last segment
% for folders with dots in the name (e.g. "1.7x", "15.5lp").
p = char(p);
while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
[~,n,e] = fileparts(p);  name = [n e];
end
