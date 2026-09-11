% stack_dff_logz_260812.m
% -----------------------------------------------------------------------
%  Sheared "diagonal wave" stack of dF/F traces for the PHASE-LOCKED cells,
%  coloured by genotype, N chunks per cell -- one figure per significance
%  threshold.
%
%  Descendant of Breath_summary_timeNphase_260331/stackDFF_chunkALL_260413.m
%  (same parallelogram geometry), but the cell list now comes from the CURRENT
%  Rayleigh analysis instead of the old coherence-based fov_map.mat.
%
%  THRESHOLDS.  logZ here is the NATURAL log of the Rayleigh Z, not log10 --
%  the upstream heatmap script pins it: logZ 1.097 <-> p .05 and 1.933 <-> p
%  .001, and exp(-exp(1.097)) = 0.05. So the analytic Rayleigh p is
%
%       p = exp(-Z) = exp(-exp(logZ))        <->    logZ = log(-log(p))
%
%  which is what converts P_LIST into cut-offs. NOTE this is the ANALYTIC
%  Rayleigh p; the phases it is computed on are occupancy-corrected and
%  weighted upstream, so treat these p values as labels on a ranking rather
%  than exact false-positive rates.
%
%  Cell selection:  rayZ >= cut  from breath_trig_heatmap_GENOTYPES_<trig>.mat.
%  roiLbl is 'MMDD/FOV/roi' -- the SITE level is not in the label, so each ROI
%  is resolved by globbing <genotypeRoot>/MMDD/*/FOV/. That recovered site is
%  what identifies the IO recordings: IO is a site inside the genotype folders
%  (Vglut2/1124/IO, ChAT/0521/IO), not a genotype of its own.
%
%  Traces are loaded ONCE at the most permissive threshold; every stricter
%  figure is a subset of those same traces, so the alignment and the chunk
%  choice are identical across thresholds and the figures are comparable.
%
%  Colours are the SAME group_colors used by the polar figures, read straight
%  out of polar_coh_vs_rayleigh_data.mat, except IO which is drawn grey
%  because it is a site rather than a genotype.
%
%  Output: stack_dff_<tag>_<trig>.png / .pdf  next to the input .mat
%
%  Runqi Zhang / 2026-08-12

clear; clc; close all;

%% ===================== USER-EDITABLE =====================
archRoot   = 'D:\Ventral_surface_summary';
dataFile   = fullfile(archRoot,'breath_trig_heatmap_260806','breath_trig_heatmap_GENOTYPES_peak.mat');
rayFile    = fullfile(archRoot,'polar_coh_vs_rayleigh_260808','polar_coh_vs_rayleigh_data.mat');

Z_LIST     = [1 2 3];            % Rayleigh logZ cut-offs
P_LIST     = [0.05 0.01 0.001];  % analytic-p cut-offs -> converted to logZ

nChunk     = 3;          % chunks per cell
chunk_sec  = 15;         % chunk duration (s)
CHUNK_PICK = 'events';   % 'events' | 'first' | 'even'

BaselineWinSec = 20;
fps_fallback   = 30;

align_xcorr    = true;
max_shift_sec  = 5;
n_align_passes = 3;

dFF_scale           = 0.3;
spacingFrac         = 0.15;
shear_from_vert_deg = 12;
ax_w_cm             = 36;
ax_h_cm_per_row     = 0.16;

% ---- cells from sessions that are NOT in dataFile -------------------------
% dataFile is the 260806 heatmap output and predates these two vagotomised
% sessions, so their cells cannot be selected from it. They are listed here
% explicitly and appended to the selection AFTER it is made. Nothing else
% changes: same colours, same geometry, same ordering rule.
%
% logZ is each cell's own occupancy-corrected Rayleigh logZ -- for Vglut2/0824
% from polar_coh_rayleigh_260824\polar_scores_percell.csv, for ChAT/0826 the
% figure_pooled_logZ written by per_cell_summary_260812. They are on the same
% footing as S.rayZ, so these cells face the SAME thresholds as every other
% cell rather than being force-included.
%
%   { genotype, MMDD, folder, roi, logZ, pooled_id }
EXTRA_CELLS = { ...
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell1','roi1_z-10_10x_3000f_16lp_00001'), 4, 2.8064, 283
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell1','roi1_z-13_8x_3000f_19lp_00001'), 5, 3.0010, 284
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell1','roi1_z-13_8x_3000f_19lp_00001'), 8, 3.5993, 285
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell1','roi1_z0_8x_3000f_12lp_00001'),   6, 3.2437, 287
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell1','roi1_z3_6x_3000f_11lp_00001'),   4, 2.1505, 289
 'Vglut2','0824', fullfile(archRoot,'Vglut2','0824','cell4','roi2_z-25_6x_3000f_14lp_00001'), 1, 2.8882, 297
 'ChAT',  '0826', fullfile(archRoot,'ChAT','0826','roi1_5x_6000f_00001'),                     3, 2.9239, 298
 };

doSave     = true;
% =========================================================

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontSize',7);

%% ---- thresholds ----
THR = struct('tag',{},'cut',{},'p',{},'label',{});
for z = Z_LIST
    THR(end+1) = struct('tag', sprintf('logZ%g',z), 'cut', z, ...
        'p', exp(-exp(z)), 'label', sprintf('logZ >= %g', z));           %#ok<SAGROW>
end
for pv = P_LIST
    z = log(-log(pv));
    THR(end+1) = struct('tag', sprintf('p%g',pv), 'cut', z, ...
        'p', pv, 'label', sprintf('p < %g  (logZ >= %.3f)', pv, z));     %#ok<SAGROW>
end
[~, o] = sort([THR.cut]); THR = THR(o);

fprintf('thresholds:\n');
for i = 1:numel(THR)
    fprintf('  %-8s cut logZ %.3f   p %.3g\n', THR(i).tag, THR(i).cut, THR(i).p);
end
zAll = min([THR.cut]);

%% ---- colours ----
RC = load(rayFile, 'groups', 'group_colors');
GENO_ORDER = RC.groups(:).';
GENO_COLOR = RC.group_colors;
iIO = find(strcmpi(GENO_ORDER,'IO'), 1);
if ~isempty(iIO); GENO_COLOR(iIO,:) = [0.55 0.55 0.55]; end   % IO = site, not genotype

%% ---- select cells at the MOST PERMISSIVE cut ----
S = load(dataFile,'roiLbl','rayZ','grpOf','GROUPS','TRIGGER');
sel = find(S.rayZ >= zAll);
fprintf('\n%d cells at the most permissive cut (logZ >= %.3f), trigger %s\n', ...
        numel(sel), zAll, S.TRIGGER);

cells = struct('grp',{},'site',{},'mmdd',{},'fov',{},'roi',{},'folder',{},'logZ',{});
for i = 1:numel(sel)
    lb   = strsplit(S.roiLbl{sel(i)}, '/');
    root = S.GROUPS{S.grpOf(sel(i)),2};
    h = dir(fullfile(root, lb{1}, '*', lb{2}, '*_cpSAM_output.mat'));
    if isempty(h); h = dir(fullfile(root, lb{1}, lb{2}, '*_cpSAM_output.mat')); end
    if isempty(h)
        warning('unresolved: %s/%s', S.GROUPS{S.grpOf(sel(i)),1}, S.roiLbl{sel(i)});
        continue;
    end
    parts = strsplit(h(1).folder, filesep);
    site  = parts{end-1};
    gname = S.GROUPS{S.grpOf(sel(i)),1};
    if strcmpi(site,'IO'); gname = 'IO'; end
    cells(end+1) = struct('grp',gname,'site',site,'mmdd',lb{1},'fov',lb{2}, ...
        'roi',str2double(lb{3}),'folder',h(1).folder,'logZ',S.rayZ(sel(i))); %#ok<SAGROW>
end
fprintf('resolved %d / %d\n', numel(cells), numel(sel));

% ---- append the sessions dataFile does not cover --------------------------
nAdd = 0;
for e = 1:size(EXTRA_CELLS,1)
    gname = EXTRA_CELLS{e,1};  mmdd = EXTRA_CELLS{e,2};
    fold  = EXTRA_CELLS{e,3};  r    = EXTRA_CELLS{e,4};
    zz    = EXTRA_CELLS{e,5};
    if ~isfolder(fold) || isempty(dir(fullfile(fold,'*_cpSAM_output.mat')))
        warning('extra cell %d: no cpSAM output under %s -- skipped', EXTRA_CELLS{e,6}, fold);
        continue;
    end
    pp = strsplit(regexprep(fold,'[\/]+$',''), filesep);
    cells(end+1) = struct('grp',gname,'site',pp{end-1},'mmdd',mmdd, ...
        'fov',pp{end},'roi',r,'folder',fold,'logZ',zz); %#ok<SAGROW>
    nAdd = nAdd + 1;
end
if nAdd > 0
    fprintf('added %d cell(s) from sessions outside dataFile (pooled ids %s)\n', ...
            nAdd, mat2str([EXTRA_CELLS{:,6}]));
end

gIdx = cellfun(@(g) find(strcmp(GENO_ORDER,g),1), {cells.grp});
[~, ord] = sortrows([gIdx(:), -[cells.logZ].'], [1 2]);
cells = cells(ord);  gIdx = gIdx(ord);

%% ---- load chunks ONCE ----
L_target  = round(chunk_sec * fps_fallback);
stack_dff = zeros(L_target, 0);
stack_lbl = {};  stack_g = [];  stack_z = [];

[uFov, ~, fovOf] = unique({cells.folder});
fprintf('\nloading %d FOVs...\n', numel(uFov));

for f = 1:numel(uFov)
    folder = uFov{f};
    mine   = find(fovOf == f);
    try
        [fps, ~] = detect_session_fps(folder);
    catch
        fps = fps_fallback;
    end
    if ~isfinite(fps) || fps <= 0; fps = fps_fallback; end

    h  = dir(fullfile(folder,'*_cpSAM_output.mat'));
    sd = load(fullfile(h(1).folder, h(1).name), 'F');
    dout = helper.dFF_RZ(double(sd.F), 'FPS', fps, 'BaselineWinSec', BaselineWinSec);
    D    = dout.dFF;

    spk = [];
    sp  = fullfile(folder,'ca_spike_data.mat');
    if strcmpi(CHUNK_PICK,'events') && isfile(sp)
        try
            q = load(sp); fn = fieldnames(q);
            for a = 1:numel(fn)
                if isstruct(q.(fn{a})) && numel(q.(fn{a})) >= 1
                    spk = q.(fn{a}); break;
                end
            end
        catch
            spk = [];
        end
    end

    L = round(chunk_sec * fps);
    for c = mine(:).'
        r = cells(c).roi;
        if r < 1 || r > size(D,2); continue; end
        y = D(:, r);
        nFull = floor(numel(y) / L);
        if nFull < 1; continue; end

        score = zeros(nFull,1);
        for b = 1:nFull
            seg = y((b-1)*L+1 : b*L);
            score(b) = max(seg) - median(seg);
        end
        if ~isempty(spk) && numel(spk) >= r
            try
                si = spk(r).spike_idx;
                for b = 1:nFull
                    score(b) = nnz(si > (b-1)*L & si <= b*L);
                end
            catch
                % keep the amplitude-based score
            end
        end

        switch lower(CHUNK_PICK)
            case 'first', pick = 1:min(nChunk,nFull);
            case 'even',  pick = unique(round(linspace(1,nFull,min(nChunk,nFull))));
            otherwise
                [~, sIdx] = sort(score,'descend');
                pick = sort(sIdx(1:min(nChunk,nFull)).');
        end

        for b = pick
            seg  = y((b-1)*L+1 : b*L);
            segR = interp1(linspace(0,1,numel(seg)).', seg, ...
                           linspace(0,1,L_target).', 'linear');
            stack_dff(:, end+1) = segR;                                    %#ok<SAGROW>
            stack_lbl{end+1}    = sprintf('%s %s/%s/%d c%d', cells(c).grp, ...
                cells(c).mmdd, cells(c).fov, cells(c).roi, b);             %#ok<SAGROW>
            stack_g(end+1)      = gIdx(c);                                 %#ok<SAGROW>
            stack_z(end+1)      = cells(c).logZ;                           %#ok<SAGROW>
        end
    end
    if mod(f,10)==0 || f==numel(uFov)
        fprintf('  %d/%d FOVs\n', f, numel(uFov));
    end
end

n_all = size(stack_dff,2);
fprintf('\n%d traces from %d cells\n', n_all, numel(cells));
assert(n_all > 0, 'No traces collected.');

%% ---- align ONCE, so every threshold figure shares the alignment ----
if align_xcorr && n_all > 1
    max_shift = round(max_shift_sec * fps_fallback);
    Z = stack_dff - mean(stack_dff,1,'omitnan');  Z(isnan(Z)) = 0;
    template = mean(Z,2);
    for pass = 1:n_align_passes
        shifts = zeros(1,n_all);
        for k = 1:n_all
            [cc, lags] = xcorr(Z(:,k), template, max_shift, 'coeff');
            [~, mi] = max(cc);  shifts(k) = lags(mi);
        end
        for k = 1:n_all
            Z(:,k)         = circshift(Z(:,k),         -shifts(k));
            stack_dff(:,k) = circshift(stack_dff(:,k), -shifts(k));
        end
        template = mean(Z,2);
        fprintf('align pass %d/%d: |shift| mean %.2f frames\n', ...
                pass, n_align_passes, mean(abs(shifts)));
    end
end

%% ---- one figure per threshold ----
fprintf('\n');
summary = cell(numel(THR),1);
for i = 1:numel(THR)
    keep  = find(stack_z >= THR(i).cut);
    nCell = nnz([cells.logZ] >= THR(i).cut);   % one entry per cell already

    if isempty(keep)
        fprintf('%-8s  no cells above cut -- skipped\n', THR(i).tag);
        summary{i} = sprintf('%-8s  0 cells', THR(i).tag);
        continue;
    end

    [~, ordK] = sortrows([stack_g(keep).', (1:numel(keep)).']);
    keep = keep(ordK);

    base = draw_stack(stack_dff(:,keep), stack_g(keep), gIdx, cells, ...
        GENO_ORDER, GENO_COLOR, THR(i), nCell, S.TRIGGER, ...
        L_target, fps_fallback, chunk_sec, nChunk, CHUNK_PICK, ...
        dFF_scale, spacingFrac, shear_from_vert_deg, ax_w_cm, ax_h_cm_per_row, ...
        fileparts(dataFile), doSave);

    summary{i} = sprintf('%-8s  %3d cells / %3d traces   %s', ...
        THR(i).tag, nCell, numel(keep), base);
    fprintf('%s\n', summary{i});
end

if doSave
    save(fullfile(fileparts(dataFile), 'stack_dff_thresholds_cells.mat'), ...
         'cells','stack_lbl','stack_g','stack_z','THR','nChunk','chunk_sec','GENO_ORDER');
end

fprintf('\n---- summary ----\n');
for i = 1:numel(summary); fprintf('%s\n', summary{i}); end


%% ========================================================================
function base = draw_stack(dff, gsub, gIdxAll, cells, GENO_ORDER, GENO_COLOR, ...
    thr, nCell, TRIGGER, L_target, fps, chunk_sec, nChunk, CHUNK_PICK, ...
    dFF_scale, spacingFrac, shear_deg, ax_w_cm, ax_h_cm_per_row, outDir, doSave)

    n_kept = size(dff,2);

    amp_g = max(dff(:)) - min(dff(:));
    if amp_g == 0; amp_g = 1; end
    gap_g   = spacingFrac * amp_g;
    offsets = gap_g * (0:n_kept-1);
    y_lo    = min(dff(:)) - gap_g;
    y_hi    = offsets(end) + max(dff(:)) + 4*gap_g;
    y_range = y_hi - y_lo;

    ax_h_cm  = max(10, ax_h_cm_per_row * n_kept);
    LEFT_CM  = 5.5;  BOT_CM = 2.0;
    fig_w_cm = ax_w_cm + LEFT_CM + 2.5;
    fig_h_cm = ax_h_cm + BOT_CM + 2.5;

    % A figure canvas is capped at the screen no matter what Position reports,
    % and an axes taller than the canvas is silently cut -- print() then
    % stretches the cut canvas onto the paper. Scale to fit; DPI restores detail.
    old = get(0,'Units'); set(0,'Units','centimeters');
    scr = get(0,'ScreenSize'); set(0,'Units',old);
    sFit = min([1, (scr(3)-2)/fig_w_cm, (scr(4)-3)/fig_h_cm]);
    if sFit < 1
        ax_w_cm = ax_w_cm*sFit;  ax_h_cm = ax_h_cm*sFit;
        LEFT_CM = LEFT_CM*sFit;  BOT_CM  = BOT_CM*sFit;
        fig_w_cm = fig_w_cm*sFit; fig_h_cm = fig_h_cm*sFit;
    end

    gap_cm = (gap_g / y_range) * ax_h_cm;
    xs_cm  = gap_cm * tand(shear_deg);
    denom  = ax_w_cm - xs_cm * max(n_kept-1,0);
    if denom <= 0
        xs_cm = 0.8 * ax_w_cm / max(n_kept-1,1);
        denom = ax_w_cm - xs_cm * max(n_kept-1,0);
    end
    xs_sec  = xs_cm * chunk_sec / denom;
    tot_sec = xs_sec * max(n_kept-1,0);

    fh = figure('Color','w','Visible','off','Units','centimeters', ...
                'Position',[1 1 fig_w_cm fig_h_cm]);
    ax = axes(fh,'Units','centimeters','Position',[LEFT_CM BOT_CM ax_w_cm ax_h_cm]);
    hold(ax,'on');

    t_base = (0:L_target-1)'/fps;
    for k = 1:n_kept
        plot(ax, t_base + (k-1)*xs_sec, dff(:,k) + offsets(k), ...
             'Color', GENO_COLOR(gsub(k),:), 'LineWidth', 0.5);
    end

    x_bar = chunk_sec + tot_sec;
    y0_b  = offsets(end) + 0.5*amp_g - dFF_scale/2;
    plot(ax, [x_bar x_bar], [y0_b y0_b+dFF_scale], 'k', 'LineWidth', 2);
    text(ax, x_bar, y0_b+dFF_scale, sprintf('  %.2g \\DeltaF/F', dFF_scale), ...
         'HorizontalAlignment','left','VerticalAlignment','top','FontSize',8);

    tickY = []; tickL = {};
    for k = 1:numel(GENO_ORDER)
        ii = find(gsub == k);
        if isempty(ii); continue; end
        nc = nnz(gIdxAll==k & [cells.logZ] >= thr.cut);
        tickY(end+1) = offsets(round(median(ii)));                       %#ok<AGROW>
        tickL{end+1} = sprintf('%s  %d cells / %d traces', ...
                               GENO_ORDER{k}, nc, numel(ii));            %#ok<AGROW>
    end
    [tickY, iT] = sort(tickY);
    set(ax,'YTick',tickY,'YTickLabel',tickL(iT));

    hL = gobjects(0); lbl = {};
    for k = 1:numel(GENO_ORDER)
        if any(gsub == k)
            hL(end+1) = plot(ax, nan, nan, '-', 'Color', GENO_COLOR(k,:), ...
                             'LineWidth', 2);                            %#ok<AGROW>
            lbl{end+1} = GENO_ORDER{k};                                  %#ok<AGROW>
        end
    end
    lg = legend(ax, hL, lbl, 'Box','off', 'FontSize', 9);
    lg.Units = 'centimeters';
    lg.Position(1) = LEFT_CM + ax_w_cm - lg.Position(3) - 0.6;
    lg.Position(2) = BOT_CM + 0.6;

    xlim(ax, [0, x_bar + 1]);  ylim(ax, [y_lo, y_hi]);
    xlabel(ax, 'Time (s)');
    title(ax, sprintf(['%s  (p = %.3g)  |  %d cells, %d traces  |  ' ...
                       '%d x %g s chunks (%s)  |  trigger %s'], ...
          thr.label, thr.p, nCell, n_kept, nChunk, chunk_sec, CHUNK_PICK, TRIGGER), ...
          'Interpreter','none');
    box(ax,'off'); hold(ax,'off');

    base = fullfile(outDir, sprintf('stack_dff_%s_%s', thr.tag, lower(TRIGGER)));
    if doSave
        set(fh,'PaperUnits','centimeters','PaperPosition',[0 0 fig_w_cm fig_h_cm], ...
               'PaperSize',[fig_w_cm fig_h_cm],'PaperPositionMode','manual');
        print(fh, [base '.png'], '-dpng', '-r400');
        set(fh,'Color','none','InvertHardcopy','off');
        print(fh, [base '.pdf'], '-dpdf', '-painters');
    end
    close(fh);
end
