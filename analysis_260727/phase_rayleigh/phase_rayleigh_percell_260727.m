function phase_rayleigh_percell_260727()
%% phase_rayleigh_percell_260727  One QC panel per cell, showing the whole chain
%  from raw event phases to the significance call, so no number in the summary
%  has to be taken on trust.
% -----------------------------------------------------------------------
% Five panels per cell:
%   1  raw event-phase histogram with the TIME-OCCUPANCY curve overlaid. If the
%      two have the same shape, the cell is not modulated -- it is just firing
%      where the time is. This is the panel that makes the occupancy problem
%      visible for that individual cell.
%   2  occupancy-normalised rate profile (events per second in each phase bin)
%      with the preferred-phase arrow. This is the biologically meaningful tuning
%      curve, and where th_rate and modulation depth come from.
%   3  ECDF-corrected phase histogram, which is what the Rayleigh statistic
%      actually sees. Flat = unmodulated under the corrected null.
%   4  the cycle-shuffle null distribution of logZ with the observed value marked.
%      At the event counts here (median 70) this is the only honest significance
%      statement, so it is plotted rather than summarised.
%   5  a text block: counts, both preferred phases, depth, both reliabilities,
%      both p-values and the FDR q.
%
% Input : phase_rayleigh_data.mat  (phase_rayleigh_260727.m, which must have been
%         run AFTER null storage was added -- it saves null_logZ_breath/uniform)
% Output: <phys>\analysis_260727\phase_rayleigh\percell\cell_###.png
%
% Runqi Zhang / 2026-07-27
close all;

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir), 'coh_ca_breath'));

cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
which_cells = 'tested';   % 'sig' | 'tested' | 'all'   (all = every cell with >=1 event)
NULL        = 'breath';   % which null to display / rank by: 'breath' | 'uniform'
showFirst   = 6;          % open this many on screen; the rest saved silently
% =====================================================================

dataFile = fullfile(cfg.phaseDir,'phase_rayleigh_data.mat');
assert(isfile(dataFile), ['phase_rayleigh_data.mat not found:\n  %s\nRun phase_rayleigh_260727.m first.'], dataFile);
S = load(dataFile,'Cell','E','binEdges','binCtrs','sigB','sigU','qB','qU','params');
Cell = S.Cell; binEdges = S.binEdges; binCtrs = S.binCtrs; params = S.params;
switch lower(NULL)
    case 'uniform', sig = S.sigU; q = S.qU; nullName = 'uniform shift';
    otherwise,      sig = S.sigB; q = S.qB; nullName = 'whole-breath shift';
end

assert(isfield(Cell,'null_logZ_breath'), ...
    ['This phase_rayleigh_data.mat predates null storage. Re-run ' ...
     'phase_rayleigh_260727.m so the null distributions are saved.']);

switch lower(which_cells)
    case 'sig',    sel = find(sig);
    case 'all',    sel = find([Cell.n_events] > 0);
    otherwise,     sel = find([Cell.tested]);
end
assert(~isempty(sel), 'No cells selected (which_cells = %s).', which_cells);
[~,o] = sort([Cell(sel).logZ],'descend','MissingPlacement','last');  sel = sel(o);

outSub = fullfile(cfg.phaseDir,'percell');
if ~isfolder(outSub), mkdir(outSub);
else, old = dir(fullfile(outSub,'cell_*.png')); for k=1:numel(old), delete(fullfile(old(k).folder,old(k).name)); end
end

col = cfg.genotype_color;
hOff = figure('Color','w','Visible','off');
fprintf('rendering %d cells (%s, ranked by logZ)...\n', numel(sel), which_cells);

for si = 1:numel(sel)
    c = sel(si);  C = Cell(c);
    onScreen = si <= showFirst;
    if onScreen, hf = figure('Color','w','Visible','on'); else, hf = hOff; clf(hf); end
    set(hf,'Units','centimeters','Position',[2 2 34 9]);

    occFrac  = C.occ_prof / max(sum(C.occ_prof), eps);
    evCount  = histcounts(C.phi_raw, binEdges);
    evFrac   = evCount / max(sum(evCount), eps);

    % --- 1. raw event phases vs time occupancy ---
    a1 = subplot(1,5,1); hold(a1,'on'); box(a1,'on');
    bar(a1, binCtrs, evFrac, 1, 'FaceColor',col, 'EdgeColor','none', 'FaceAlpha',0.75);
    plot(a1, binCtrs, occFrac, 'k-', 'LineWidth',1.6);
    yline(a1, 1/numel(binCtrs), 'k:');
    xlim(a1,[0 2*pi]); set(a1,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
    xlabel(a1,'raw phase'); ylabel(a1,'fraction');
    title(a1,{'events (bars) vs TIME occupancy (line)','same shape = not modulated'},'FontSize',7);

    % --- 2. occupancy-normalised rate profile ---
    a2 = subplot(1,5,2); hold(a2,'on'); box(a2,'on');
    bar(a2, binCtrs, C.rate_prof, 1, 'FaceColor',[0 .45 .74], 'EdgeColor','none');
    yl = ylim(a2);
    plot(a2, [C.th_rate C.th_rate], yl, '-', 'Color',[.85 .10 .10], 'LineWidth',2);
    xlim(a2,[0 2*pi]); set(a2,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
    xlabel(a2,'raw phase'); ylabel(a2,'events / s in bin');
    % cosine fit overlaid: rate ~ a + b*cos(theta - mu), depth = b/a
    if isfinite(C.mod_depth_cos) && isfinite(C.th_cos)
        aFit = mean(C.rate_prof(isfinite(C.rate_prof)));
        yFit = aFit * (1 + C.mod_depth_cos*cos(binCtrs - C.th_cos));
        plot(a2, binCtrs, yFit, '-', 'Color',[.85 .10 .10], 'LineWidth',1.2);
    end
    title(a2,sprintf('rate profile   b/a=%.2f (R^2=%.2f)\nred = preferred %.0f deg', ...
          C.mod_depth_cos, C.r2_cos, mod(rad2deg(C.th_rate),360)),'FontSize',7);

    % --- 3. ECDF-corrected phases = what Rayleigh sees ---
    a3 = subplot(1,5,3); hold(a3,'on'); box(a3,'on');
    histogram(a3, C.psi, binEdges, 'Normalization','probability', ...
              'FaceColor',[.20 .60 .30], 'EdgeColor','none');
    yline(a3, 1/numel(binCtrs), 'k:','LineWidth',1);
    xlim(a3,[0 2*pi]); set(a3,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
    xlabel(a3,'corrected phase \psi'); ylabel(a3,'fraction');
    title(a3,sprintf('occupancy-corrected\nrbar=%.3f  logZ=%.2f', C.rbar, C.logZ),'FontSize',7);

    % --- 4. the shuffle null ---
    a4 = subplot(1,5,4); hold(a4,'on'); box(a4,'on');
    if strcmpi(NULL,'uniform'), nl = double(C.null_logZ_uniform); else, nl = double(C.null_logZ_breath); end
    if ~isempty(nl)
        histogram(a4, nl, 30, 'FaceColor',[.72 .72 .72], 'EdgeColor','none');
        yl = ylim(a4);
        plot(a4, [C.logZ C.logZ], yl, '-', 'Color',col, 'LineWidth',2);
        plot(a4, [prctile(nl,95) prctile(nl,95)], yl, 'k--', 'LineWidth',1);
    end
    xlabel(a4,'logZ under the null'); ylabel(a4,'count');
    title(a4,sprintf('%s null (%d iter)\ncolored = observed, dashed = 95th pct', ...
          nullName, params.nShuffle),'FontSize',7);

    % --- 5. the numbers ---
    a5 = subplot(1,5,5); axis(a5,'off');
    pv = C.p_breath; if strcmpi(NULL,'uniform'), pv = C.p_uniform; end
    txt = {
        sprintf('\\bfcell %d\\rm', C.cell_id)
        sprintf('%d recording(s)', C.n_obs)
        ''
        sprintf('events        %d', C.n_events)
        sprintf('duration      %.0f s', C.dur_s)
        sprintf('rate          %.2f /min', C.rate_perMin)
        sprintf('breaths       %d', C.n_breaths)
        sprintf('  with event  %d', C.n_breaths_with_event)
        ''
        sprintf('logZ          %.2f', C.logZ)
        sprintf('rbar          %.3f', C.rbar)
        sprintf('depth b/a     %.2f  (R2 %.2f)', C.mod_depth_cos, C.r2_cos)
        ''
        sprintf('pref (rate)   %.0f deg', mod(rad2deg(C.th_rate),360))
        sprintf('pref (raw)    %.0f deg', mod(rad2deg(C.th_raw),360))
        ''
        sprintf('recruitment   %.3f', C.recruitment)
        sprintf('precision     %.3f', C.precision)
        ''
        sprintf('p (breath)    %.4f', C.p_breath)
        sprintf('p (uniform)   %.4f', C.p_uniform)
        sprintf('q (BH)        %.4f', q(c))
        sprintf('\\bf%s\\rm', ternary(sig(c), 'SIGNIFICANT', 'not significant'))
        };
    text(a5, 0, 1, txt, 'Units','normalized','VerticalAlignment','top', ...
         'FontName','FixedWidth','FontSize',7.5, 'Interpreter','tex');

    sgtitle(sprintf('cell %d  |  %s', C.cell_id, strjoin(cellstr(C.rec_names),'  ')), ...
            'Interpreter','none','FontSize',8);

    exportgraphics(hf, fullfile(outSub, sprintf('cell_%03d_logZ%05.2f.png', C.cell_id, max(C.logZ,0))), ...
                   'Resolution',150);
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d per-cell panels to %s\n', numel(sel), outSub);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
