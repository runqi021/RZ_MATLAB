function breath_time_trigger_compare_260727()
%% breath_time_trigger_compare_260727  Onset-locked or peak-locked?
% -----------------------------------------------------------------------
% Runs after breath_time_peth_260727.m has been run BOTH ways (TRIGGER='onset'
% and TRIGGER='peak'). Compares the two alignments cell by cell.
%
% WHY THIS IS A RESULT, NOT A ROBUSTNESS CHECK
% Inspiration onset and the inspiratory peak are two different physical events,
% separated here by ~267 ms with CV ~0.15 (SD ~40 ms). That separation is not
% constant, so a cell locked to the PEAK is smeared by the onset->peak jitter when
% aligned on the ONSET, and vice versa. Whichever alignment gives the TIGHTER
% latency distribution is the landmark the cell is actually locked to.
%
% The primary discriminator is latency MAD, not M_exc: MAD measures dispersion
% directly, whereas M_exc mixes dispersion with firing rate and with how much
% evidence the cell provided, so it barely moves between alignments.
%
% CONSISTENCY CHECK, printed: median onset latency minus the median onset->peak
% interval should equal the median peak latency. If it does not, the two runs are
% not describing the same event and something is wrong upstream.
%
% Output: <dataset>\analysis_260727\breath_time\trigger_compare.csv / .png / .pdf
%
% Runqi Zhang / 2026-07-27
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();
base = fullfile(cfg.outRoot,'breath_time');

fA = fullfile(base,'onset','breath_time_peth_data.mat');
fB = fullfile(base,'peak' ,'breath_time_peth_data.mat');
for f = {fA,fB}
    assert(isfile(f{1}), ['Missing %s\nRun breath_time_peth_260727.m with ' ...
        'TRIGGER = ''onset'' AND with TRIGGER = ''peak''.'], f{1});
end
A = load(fA,'R','sig_exc','params');  B = load(fB,'R','sig_exc','params');

n = min(numel(A.R), numel(B.R));
ids = []; madA = []; madB = []; mA = []; mB = []; latA = []; latB = []; nEv = []; sig = [];
for c = 1:n
    if ~A.R(c).tested || ~B.R(c).tested, continue; end
    ids(end+1,1)  = A.R(c).cell_id; %#ok<AGROW>
    madA(end+1,1) = A.R(c).latency_mad_s; madB(end+1,1) = B.R(c).latency_mad_s; %#ok<AGROW>
    mA(end+1,1)   = A.R(c).mod_exc_z;     mB(end+1,1)   = B.R(c).mod_exc_z; %#ok<AGROW>
    latA(end+1,1) = A.R(c).preferred_latency_s; latB(end+1,1) = B.R(c).preferred_latency_s; %#ok<AGROW>
    nEv(end+1,1)  = A.R(c).n_events; %#ok<AGROW>
    sig(end+1,1)  = A.sig_exc(c) | B.sig_exc(c); %#ok<AGROW>
end
assert(~isempty(ids), 'No cell was tested in both runs.');
sig = logical(sig);

tol         = 1e-9;
tied        = abs(madA - madB) <= tol;   % latency MAD is quantised, so exact ties
lockedOnset = madA < madB - tol;         % are common and must not be silently
lockedPeak  = madB < madA - tol;         % broken in favour of whichever is tested first
fprintf('\n======= onset vs peak alignment =======\n');
fprintf('cells tested in both : %d (%d significant in at least one)\n', numel(ids), nnz(sig));
fprintf('latency MAD  : onset median %.0f ms | peak median %.0f ms\n', ...
        1000*median(madA), 1000*median(madB));
fprintf('M_exc        : onset median %.2f    | peak median %.2f\n', median(mA), median(mB));
fprintf('tighter on ONSET : %d   tighter on PEAK : %d   TIED : %d   (of %d)\n', ...
        nnz(lockedOnset), nnz(lockedPeak), nnz(tied), numel(ids));
if any(sig)
    fprintf('among significant cells only: ONSET %d, PEAK %d, TIED %d\n', ...
            nnz(lockedOnset & sig), nnz(lockedPeak & sig), nnz(tied & sig));
end
d = median(latA) - median(latB);
fprintf(['\nconsistency: median onset latency %.0f ms - median peak latency %.0f ms = %.0f ms\n' ...
         '             this should equal the median onset->peak interval; if it does,\n' ...
         '             both runs are describing the same event in two reference frames.\n'], ...
        1000*median(latA), 1000*median(latB), 1000*d);

%% ---- figure ----
fig = figure('Color','w','Name','onset vs peak alignment','Units','centimeters','Position',[2 2 30 11]);
set(fig,'DefaultAxesFontSize',8);
sz = 20 + 60*min(log10(max(nEv,1))/log10(500),1);

a1 = axes(fig,'Position',[0.07 0.14 0.25 0.70]); hold(a1,'on'); box(a1,'on'); grid(a1,'on');
mx = 1000*max([madA;madB])*1.05;
plot(a1,[0 mx],[0 mx],'k:');
scatter(a1, 1000*madA(~sig), 1000*madB(~sig), sz(~sig), [.72 .72 .72],'filled');
scatter(a1, 1000*madA(sig),  1000*madB(sig),  sz(sig),  cfg.genotype_color,'filled');
xlim(a1,[0 mx]); ylim(a1,[0 mx]);
xlabel(a1,'latency MAD, ONSET-aligned (ms)'); ylabel(a1,'latency MAD, PEAK-aligned (ms)');
title(a1,{'below the line = tighter when peak-aligned','above = tighter when onset-aligned'},'FontSize',7);

a2 = axes(fig,'Position',[0.40 0.14 0.25 0.70]); hold(a2,'on'); box(a2,'on'); grid(a2,'on');
mm = max([mA;mB])*1.05; mn = min([0;mA;mB]);
plot(a2,[mn mm],[mn mm],'k:');
scatter(a2, mA(~sig), mB(~sig), sz(~sig), [.72 .72 .72],'filled');
scatter(a2, mA(sig),  mB(sig),  sz(sig),  cfg.genotype_color,'filled');
xlabel(a2,'M_{exc}, ONSET-aligned'); ylabel(a2,'M_{exc}, PEAK-aligned');
title(a2,{'M_{exc} barely discriminates:','it mixes dispersion with rate and evidence'},'FontSize',7);

a3 = axes(fig,'Position',[0.73 0.14 0.24 0.70]); hold(a3,'on'); box(a3,'on'); grid(a3,'on');
dd = 1000*(madB - madA);                       % >0 means onset is tighter
histogram(a3, dd, 15, 'FaceColor',cfg.genotype_color,'EdgeColor','none');
xline(a3,0,'k--','LineWidth',1.5);
xlabel(a3,'MAD_{peak} - MAD_{onset}  (ms)'); ylabel(a3,'# cells');
title(a3,{'right of 0 = ONSET-locked','left of 0 = PEAK-locked'},'FontSize',7);

sgtitle(sprintf(['%s  |  %d cells in both  |  tighter on onset %d, on peak %d  |  ' ...
                 'onset lat %.0f ms, peak lat %.0f ms (difference %.0f ms)'], ...
        cfg.genotype, numel(ids), nnz(lockedOnset), nnz(lockedPeak), ...
        1000*median(latA), 1000*median(latB), 1000*d),'FontSize',9);

exportgraphics(fig, fullfile(base,'trigger_compare.png'),'Resolution',200,'BackgroundColor','white');
exportgraphics(fig, fullfile(base,'trigger_compare.pdf'),'ContentType','vector','BackgroundColor','white');

T = table(ids, nEv, sig, 1000*madA, 1000*madB, mA, mB, 1000*latA, 1000*latB, ...
          repmat("",numel(ids),1), ...
    'VariableNames',{'cell_id','n_events','sig_either','mad_onset_ms','mad_peak_ms', ...
                     'mod_exc_z_onset','mod_exc_z_peak','latency_onset_ms','latency_peak_ms','locked_to'});
T.locked_to(lockedOnset) = "onset";
T.locked_to(lockedPeak)  = "peak";
T.locked_to(tied)        = "tied";
writetable(T, fullfile(base,'trigger_compare.csv'));
fprintf('\nSaved trigger_compare.csv/.png/.pdf to %s\n', base);
end
