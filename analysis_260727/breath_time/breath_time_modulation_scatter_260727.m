function breath_time_modulation_scatter_260727(triggerIn)
%% breath_time_modulation_scatter_260727  Population summary: one point per cell.
% -----------------------------------------------------------------------
% TWO scatters, because they answer different questions and only one of them is
% guaranteed to carry two independent dimensions.
%
% A. MODULATION vs PRECISION      x = M_exc      y = precision_z
%       upper right : strongly recruited AND precisely timed
%       upper left  : quiet or weakly recruited, but precisely locked when active
%       lower right : broad modulation with poor temporal precision
%       lower left  : no detectable respiratory relationship
%    CAVEAT, and it is checked and printed rather than assumed: both axes are
%    shuffle-standardised statistics whose power scales with event count, so they
%    can collapse onto a single axis. The correlation is computed and, if it is
%    high, the figure says so on its face and defers to scatter B.
%
% B. RECRUITMENT vs PRECISION     x = recruitment   y = 1 / latency MAD
%    These are raw biological quantities with no shuffle standardisation:
%    "what fraction of breaths recruit this cell" against "when it does fire, how
%    reproducible is the timing". They dissociate by construction -- a cell can be
%    recruited on 5% of breaths and still fire at the same latency every time --
%    so this plot keeps its two dimensions whatever the z-scores do.
%
% Encodings on both: marker size = event count, colour = preferred latency,
% edge = which tests survived FDR.
%
% Input : breath_time_peth_data.mat
% Output: <dataset>\analysis_260727\breath_time\breath_time_modulation_scatter.png/.pdf
%
% Runqi Zhang / 2026-07-27
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();

doSave = true;

TRIGGER  = 'onset';   % 'onset' | 'peak' -- must match the run you want to plot
if nargin >= 1 && ~isempty(triggerIn), TRIGGER = char(triggerIn); end
outDir   = fullfile(cfg.outRoot,'breath_time',lower(TRIGGER));
dataFile = fullfile(outDir,'breath_time_peth_data.mat');
assert(isfile(dataFile), ['breath_time_peth_data.mat not found:\n  %s\n' ...
    'Run breath_time_peth_260727.m first.'], dataFile);
S = load(dataFile,'R','sig_exc','sig_sup','sig_pre','params');
R = S.R; params = S.params;
ti = find([R.tested]);
assert(~isempty(ti), 'No tested cells.');

M    = [R(ti).mod_exc_z]';
Pz   = [R(ti).precision_z]';
recr = [R(ti).recruitment]';
mad_ = [R(ti).latency_mad_s]';
lat  = 1000*[R(ti).preferred_latency_s]';
nEv  = [R(ti).n_events]';
sE   = S.sig_exc(ti)';   sP = S.sig_pre(ti)';   sS = S.sig_sup(ti)';

ok = isfinite(M) & isfinite(Pz);
rr = NaN;
if nnz(ok) > 3, rr = corr(M(ok), Pz(ok)); end
oneAxis = abs(rr) > 0.8;
fprintf('\ncorr(M_exc, precision_z) = %.2f over %d cells\n', rr, nnz(ok));
if oneAxis
    fprintf('  -> these are effectively ONE axis. Scatter B (recruitment vs precision)\n');
    fprintf('     is the one that carries independent information; lead with it.\n');
else
    fprintf('  -> genuinely two axes; scatter A is informative.\n');
end

sz = 20 + 60*min(log10(max(nEv,1))/log10(500),1);

fig = figure('Color','w','Name','breath time modulation vs precision', ...
             'Units','centimeters','Position',[2 2 32 13]);
set(fig,'DefaultAxesFontSize',8);

% ---- A ----
a1 = axes(fig,'Position',[0.07 0.13 0.36 0.72]); hold(a1,'on'); box(a1,'on'); grid(a1,'on');
scatter(a1, M, Pz, sz, lat, 'filled', 'MarkerFaceAlpha',0.85);
gi = sE & sP;
if any(gi),  scatter(a1, M(gi), Pz(gi), sz(gi), 'k', 'LineWidth',1.1); end
yline(a1,0,'k:'); xline(a1,0,'k:');
colorbar(a1);   % NB: do NOT touch cb.Label / ylabel(cb) -- in R2021b that getter fires
                % a listener that errors on a scatter with CData. Label goes in the title.
colormap(a1, parula);
xlabel(a1,'M_{exc}  (shuffle-standardised PETH modulation)');
ylabel(a1,'precision_z  (shuffle-standardised latency concentration)');
ttl = sprintf('A. modulation vs precision   r = %.2f', rr);
if oneAxis, ttl = [ttl '   -- nearly ONE axis, prefer B']; end
title(a1, {ttl,'colour = preferred latency (ms); size = # events; circled = sig on both'}, 'FontSize',8);

% ---- B ----
a2 = axes(fig,'Position',[0.57 0.13 0.36 0.72]); hold(a2,'on'); box(a2,'on'); grid(a2,'on');
invmad = 1./max(mad_, params.binWidth_s);      % floor at the bin width: latency cannot
                                               % be resolved finer than one bin
scatter(a2, recr, invmad, sz, lat, 'filled', 'MarkerFaceAlpha',0.85);
if any(sE), scatter(a2, recr(sE), invmad(sE), sz(sE), 'k', 'LineWidth',1.1); end
colorbar(a2);
colormap(a2, parula);
xlabel(a2,'recruitment  (accepted breaths with an event near preferred latency)');
ylabel(a2,'1 / latency MAD  (s^{-1})');
title(a2,{'B. recruitment vs timing precision   (colour = preferred latency, ms)', ...
          'upper left = fires rarely but always at the same moment'},'FontSize',8);

sgtitle(sprintf(['%s  |  %d cells tested  |  exc %d, sup %d, precision %d significant (BH q<%.2f)  ' ...
                 '|  %d shuffles, %.0f ms bins, \\pm%.0f ms precision window'], ...
        cfg.genotype, numel(ti), nnz(sE), nnz(sS), nnz(sP), params.fdr_q, ...
        params.nShuffle, 1000*params.binWidth_s, 1000*params.relWin_s), 'FontSize',9);

if doSave
    exportgraphics(fig, fullfile(outDir,'breath_time_modulation_scatter.png'), ...
                   'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'breath_time_modulation_scatter.pdf'), ...
                   'ContentType','vector','BackgroundColor','white');
    fprintf('Saved breath_time_modulation_scatter.png/.pdf to %s\n', outDir);
end

% ---- a short ranked listing ----
if any(sE)
    [~,o] = sort(M,'descend'); o = o(sE(o));
    fprintf('\ntop excited cells:\n');
    for k = o(1:min(10,numel(o)))'
        c = ti(k);
        fprintf('  cell %3d  M_exc=%5.2f  lat=%4.0f ms  n=%4d ev  recruit=%.3f  prec=%.3f  MAD=%3.0f ms\n', ...
            R(c).cell_id, R(c).mod_exc_z, 1000*R(c).preferred_latency_s, R(c).n_events, ...
            R(c).recruitment, R(c).precision_fraction, 1000*R(c).latency_mad_s);
    end
end
end
