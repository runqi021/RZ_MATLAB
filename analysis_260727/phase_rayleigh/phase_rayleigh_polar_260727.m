function phase_rayleigh_polar_260727()
%% phase_rayleigh_polar_260727  The polar plot, on Rayleigh modulation instead of coherence.
% -----------------------------------------------------------------------
%   angle  = preferred phase (occupancy-normalised rate profile -- see below)
%   radius = logZ = log(R^2/n), the Sirota modulation score
%   fill   = significant after the cycle shuffle + Benjamini-Hochberg FDR
%
% RADIUS CAVEAT. logZ is unbounded and goes NEGATIVE whenever n*rbar^2 < 1, which
% happens for genuinely unmodulated cells. A polar radius cannot be negative, so
% those are drawn clamped at 0 and marked with an inward tick; their true value is
% in the CSV. The radial axis therefore starts at logZ = 0 (Z = 1), which is the
% expected value of Z under the null -- so the origin is "no modulation" and
% distance from the origin is evidence, which reads correctly.
%
% ANGLE. Preferred phase is taken from the occupancy-normalised rate profile, not
% from the raw resultant. With 15% of each cycle spent in [0,pi) on this data, the
% raw resultant points into expiration for every cell regardless of tuning. Panel 2
% plots raw against rate-normalised so the size of that bias is visible rather
% than assumed.
%
% Convention: 0 = inspiration onset (foot), pi = inspiratory peak / start of
% expiration, 2pi = next inspiration onset.
%
% Input : phase_rayleigh_data.mat  (phase_rayleigh_260727.m)
% Output: <phys>\analysis_260727\phase_rayleigh\phase_rayleigh_polar.png/.pdf
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir), 'coh_ca_breath'));

cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
NULL       = 'breath';    % which null decides significance: 'breath' | 'uniform'
sizeByN    = true;        % marker size scales with event count
doSave     = true;
% =====================================================================

dataFile = fullfile(cfg.phaseDir,'phase_rayleigh_data.mat');
assert(isfile(dataFile), ['phase_rayleigh_data.mat not found:\n  %s\nRun phase_rayleigh_260727.m first.'], dataFile);
S = load(dataFile, 'Cell','sigB','sigU','qB','qU','params','binCtrs');
Cell = S.Cell;  params = S.params;  binCtrs = S.binCtrs;
switch lower(NULL)
    case 'uniform', sig = S.sigU;  q = S.qU;  nullName = 'uniform-shift';
    otherwise,      sig = S.sigB;  q = S.qB;  nullName = 'whole-breath-shift';
end

tested = [Cell.tested];
idx    = find(tested);
assert(~isempty(idx), 'No tested cells in phase_rayleigh_data.mat.');

logZ  = [Cell(idx).logZ]';
thR   = [Cell(idx).th_rate]';
thRaw = [Cell(idx).th_raw]';
nEv   = [Cell(idx).n_events]';
nObs  = [Cell(idx).n_obs]';
depth = [Cell(idx).mod_depth_cos]';   % cosine-fit b/a. NOT mod_depth: that one
                                      % saturates at 1 whenever a bin is empty,
                                      % which is nearly every cell at these counts.
recr  = [Cell(idx).recruitment]';
prec  = [Cell(idx).precision]';
p95   = [Cell(idx).logZ_null_p95_breath]';
if strcmpi(NULL,'uniform'), p95 = [Cell(idx).logZ_null_p95_uniform]'; end
sg    = sig(idx)';

rPlot   = max(logZ, 0);                 % clamp negatives to the origin
clamped = logZ < 0;
rMax    = max([rPlot; p95; 1]) * 1.08;
col     = cfg.genotype_color;

fig = figure('Color','w','Name','Rayleigh phase modulation (polar)', ...
             'Units','centimeters','Position',[2 2 32 13]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);

%% ---- panel 1: the polar ----
ax = polaraxes(fig,'Position',[0.03 0.09 0.30 0.78]); hold(ax,'on');
thc = linspace(0,2*pi,360);
polarplot(ax, thc, repmat(median(p95,'omitnan'),1,360), 'k--','LineWidth',1);
for k = 1:numel(idx)
    ms = 5; if sizeByN, ms = 4 + 4*min(log10(max(nEv(k),1))/log10(500), 1); end
    if sg(k)
        polarplot(ax, thR(k), rPlot(k), 'o', 'MarkerFaceColor',col, ...
                  'MarkerEdgeColor','k', 'MarkerSize',ms, 'LineWidth',0.4);
    else
        polarplot(ax, thR(k), rPlot(k), 'o', 'MarkerFaceColor','none', ...
                  'MarkerEdgeColor',col+(1-col)*0.45, 'MarkerSize',ms, 'LineWidth',0.5);
    end
    if clamped(k)   % logZ < 0: drawn at the origin, flagged
        polarplot(ax, thR(k), 0, 'x', 'Color',[.55 .55 .55], 'MarkerSize',4, 'LineWidth',0.5);
    end
end
ax.RLim = [0 rMax];
ax.ThetaZeroLocation = 'right'; ax.ThetaDir = 'counterclockwise'; ax.RAxisLocation = 180;
ax.ThetaTick = [0 90 180 270];
ax.ThetaTickLabel = {'0 insp onset','\pi/2','\pi peak','3\pi/2'};
title(ax, {sprintf('%s: breathing modulation of Ca events', cfg.genotype), ...
           'radius = logZ   filled = sig   dashed = median null 95th pct'}, 'FontSize',8);

%% ---- panel 2: how much occupancy biases the raw preferred phase ----
ax2 = axes(fig,'Position',[0.395 0.60 0.24 0.31]); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, [0 360],[0 360],'k:');
scatter(ax2, mod(rad2deg(thR),360), mod(rad2deg(thRaw),360), 22, col, 'filled', 'MarkerFaceAlpha',0.6);
xlim(ax2,[0 360]); ylim(ax2,[0 360]);
set(ax2,'XTick',0:90:360,'YTick',0:90:360);
xlabel(ax2,'preferred phase, rate-normalised (deg)'); ylabel(ax2,'raw (deg)');
title(ax2,'raw preferred phase is dragged toward expiration','FontSize',7);

%% ---- panel 3: logZ vs event count, with the null band ----
ax3 = axes(fig,'Position',[0.395 0.11 0.24 0.31]); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
scatter(ax3, nEv(~sg), logZ(~sg), 20, [.72 .72 .72], 'filled');
scatter(ax3, nEv(sg),  logZ(sg),  26, col, 'filled');
scatter(ax3, nEv, p95, 12, 'k', '+');
set(ax3,'XScale','log'); yline(ax3, 0, 'k:');
xlabel(ax3,'# events pooled into the cell'); ylabel(ax3,'logZ');
title(ax3,'+ = that cell''s own null 95th pct (it moves with n)','FontSize',7);

%% ---- panel 4: modulation depth vs logZ ----
ax4 = axes(fig,'Position',[0.71 0.60 0.26 0.31]); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
scatter(ax4, depth(~sg), logZ(~sg), 20, [.72 .72 .72], 'filled');
scatter(ax4, depth(sg),  logZ(sg),  26, col, 'filled');
set(ax4,'XScale','log');
xlabel(ax4,'cosine-fit modulation depth  b/a'); ylabel(ax4,'logZ');
title(ax4,'effect size vs evidence','FontSize',7);

%% ---- panel 5: cycle reliability, both denominators ----
ax5 = axes(fig,'Position',[0.71 0.11 0.26 0.31]); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');
scatter(ax5, recr(~sg), prec(~sg), 20, [.72 .72 .72], 'filled');
scatter(ax5, recr(sg),  prec(sg),  26, col, 'filled');
xlabel(ax5,sprintf('recruitment (breaths near pref / ALL breaths)'));
ylabel(ax5,'precision (near pref / event-bearing breaths)');
ylim(ax5,[0 1]);
title(ax5,sprintf('window \\pm%.0f deg; low recruit + high precision = rare but reliable', ...
                  rad2deg(params.relWin_rad)),'FontSize',7);

sgtitle(sprintf(['%s  |  %d cells tested (rate \\geq %g/min)  |  %d significant, %s null, BH q<%.2f' ...
                 '  |  %d shuffles  |  %d with logZ<0 drawn at origin'], ...
        cfg.genotype, numel(idx), params.minRate_perMin, nnz(sg), nullName, params.fdr_q, ...
        params.nShuffle, nnz(clamped)), 'FontSize',9);

%% ---- save ----
if doSave
    if ~isfolder(cfg.phaseDir), mkdir(cfg.phaseDir); end
    exportgraphics(fig, fullfile(cfg.phaseDir,'phase_rayleigh_polar.png'), ...
                   'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(cfg.phaseDir,'phase_rayleigh_polar.pdf'), ...
                   'ContentType','vector','BackgroundColor','white');
    fprintf('Saved phase_rayleigh_polar.png/.pdf to %s\n', cfg.phaseDir);
end

fprintf('\n%d/%d cells significant (%s null, BH q<%.2f)\n', nnz(sg), numel(idx), nullName, params.fdr_q);
if any(sg)
    [~,o] = sort(logZ(sg),'descend');  si = idx(sg);  si = si(o);
    fprintf('top modulated cells:\n');
    for k = si(1:min(10,numel(si)))
        fprintf('  cell %3d  logZ=%5.2f  pref=%4.0f deg  n=%4d ev  b/a=%.2f (R2 %.2f)  prec=%.2f  q=%.4f\n', ...
            Cell(k).cell_id, Cell(k).logZ, mod(rad2deg(Cell(k).th_rate),360), ...
            Cell(k).n_events, Cell(k).mod_depth_cos, Cell(k).r2_cos, Cell(k).precision, q(k));
    end
end
end
