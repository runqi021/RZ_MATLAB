% pop_features_260813.m
% -----------------------------------------------------------------------
%  Build the per-CELL feature matrix for the population analysis.
%
%  Step 1 of the two-space plan:
%     FUNCTIONAL SPACE  cluster on dynamics alone (this file feeds it)
%     then overlay genotype / xy / breath stats as SUPPLEMENTARY variables.
%
%  THE ACTIVE FEATURE IS A TIME-DOMAIN WAVEFORM: the peak-triggered average dF/F
%  over a FIXED +/-1 s window, every cycle of every recording of that cell
%  averaged together. Time, not phase, because GCaMP rise and decay are properties
%  of the indicator, not of the breath: warping the axis by each cell's IBI would
%  stretch the same kinetics differently for a 0.4 Hz and a 2 Hz breather and put
%  identical cells in different places. The phase profiles are still exported, as
%  supplementary descriptors.
%
%  WHAT IS AN ACTIVE FEATURE AND WHAT IS SUPPLEMENTARY
%    ACTIVE (defines the space)   the dF/F phase profile and the event-rate phase
%                                 profile, each normalised to a SHAPE per cell
%    SUPPLEMENTARY (projected in) logZ, Rbar, preferred phase, event rate, IBI,
%                                 xy, genotype, session
%  logZ is a summary of the same phase tuning that the profiles already carry, so
%  using it as an active feature would double-count the signal and manufacture
%  structure. It is kept as a validation variable. p-values are not exported as
%  coordinates at all: ours is a permutation p with a floor at 1/(nShuffle+1), so
%  -log10(p) is censored for every strongly modulated cell.
%
%  NESTING. 276 cells come from 19 sessions, and IO_1124 alone is 121 of them
%  (44%). Session is exported so it can be checked as a confound -- any clustering
%  here must be tested against "did I just recover the big session?".
%
%  Output: <rootPath>\pop_analysis_260813\pop_features.mat  (+ a flat CSV of the
%  scalar columns; the profiles stay in the .mat).
%
%  Runqi Zhang / 2026-08-13
% -----------------------------------------------------------------------

clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
inFile   = fullfile(rootPath, 'event_latency_260811', 'event_latency_data.mat');
rayCsv   = fullfile(rootPath, 'breath_time_summary_260808', 'breath_time_summary_cells.csv');
outDir   = fullfile(rootPath, 'pop_analysis_260813');

nPhaseBins = 36;            % 10 deg bins for the supplementary phase profiles
featWin_sec = 1.0;          % ACTIVE feature: peak-triggered average over +/- this

% ACTIVE = pooled firing rate >= 2 events per minute.
%   Replaces the archive's old ">5 events" count, which was duration-dependent:
%   recordings run 49 s to 995 s here, so a flat 6-event floor meant anywhere from
%   7.3 to 0.36 events/min depending only on how long the recording happened to be.
%   Rate is pooled over the cell's recordings (nSpikes / durTot), like every other
%   pooled quantity. MUST match per_cell_summary_260812.m, or T and the rendered
%   figure folder describe different populations and the IBI split mismatches.
activeMinRateHz = 2/60;
P = struct();
P.doCoh=false; P.nDrop=30; P.fallback_fps=30; P.TW_spec=6; P.alpha_sig=0.01;
P.minSpikes=2; P.ca_lag_sec=0; P.f_breath_search=[0.2 4]; P.fwhm_factor=0.6;
P.min_bw=0.05; P.fmin=0.05; P.fmax=15; P.trigWin_sec=[]; P.trigWinIBI=2;
P.ylim_dff=[]; P.ylim_epc=[]; P.histBinFrames=2; P.shiftMinCyc=3;
% 1200 shuffles: the permutation p must MATCH the figures, so the count must too.
P.nShuffle=1200;
P.pad_um=20; P.clip_pct=[0.5 99.9]; P.scalebar_um=50;
P.rayPhaseBins=nPhaseBins; P.featWin_sec=featWin_sec; P.featFs=30;
% gamma is unused here (statsOnly draws nothing) but is kept at 1 to match the
% figures, so nobody later copies 0.6 out of this file believing it is current.
P.gamma_val=1; P.PixelSizeBase=1.7778;
P.outlineLW=0.8; P.sortMode='none';
P.dffColor=[0.2 0.7 0.2]; P.onsetCol=[0.9 0.1 0.1]; P.peakCol=[0.35 0.75 1];
P.statsOnly = true;
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end

Din  = load(inFile, 'CELL', 'OBS', 'REC');
CELL = Din.CELL; OBS = Din.OBS; REC = Din.REC;
rayArchive = nan(numel(CELL),1);
if isfile(rayCsv)
    Bt = readtable(rayCsv);
    [tf, loc] = ismember(string({CELL.key}), string(Bt.cell_key));
    rayArchive(tf) = Bt.ray_logZ(loc(tf));
end

%% ---- cell identity corrections (must match per_cell_summary_260812.m) ----
% Same merge list, applied the same way, so T.cell indexes the same cells as the
% rendered figure filenames -- the IBI split copies figures by matching T.cell
% against '_cell%03d.png', and any divergence here silently mismatches them.
% A merged cell takes the LOWEST of its original numbers; archiveLogZ goes NaN
% because the archive CSV holds one value per OLD cell.
[obsMerged, mergeInfo, origCells] = ...
        apply_cell_merges_260814(CELL, OBS, cell_merge_overrides_260814());
for q = 1:numel(mergeInfo)
    if ~strcmp(mergeInfo(q).status,'merged')
        fprintf(2,'  merge group %d: %s\n', mergeInfo(q).group, mergeInfo(q).status);
    end
end
CELL2 = struct('key',{},'obs',{});
cellNum = zeros(numel(obsMerged),1);
ray2    = nan(numel(obsMerged),1);
for i = 1:numel(obsMerged)
    cellNum(i)   = min(origCells{i});
    CELL2(i).obs = obsMerged{i};
    CELL2(i).key = CELL(cellNum(i)).key;
    if numel(origCells{i}) == 1, ray2(i) = rayArchive(origCells{i}); end
end
fprintf('cell merges: %d groups applied, %d cells -> %d\n', ...
        nnz(strcmp({mergeInfo.status},'merged')), numel(CELL), numel(CELL2));
CELL = CELL2;  rayArchive = ray2;

%% ===================== EXTRACT =====================
n = numel(CELL);
dffPhase  = nan(n, nPhaseBins);
ratePhase = nan(n, nPhaseBins);
occPhase  = nan(n, nPhaseBins);
featPeak  = [];  featOnset = [];  featTau = [];   % sized from the first cell
S = struct('cell',{},'key',{},'group',{},'session',{},'nRec',{},'nEvents',{}, ...
           'nSpikes',{},'rateHz',{},'durTot',{},'IBI',{},'fps',{},'logZ',{}, ...
           'Rbar',{},'nEff',{},'prefPhase',{},'x',{},'y',{},'archiveLogZ',{}, ...
           'pOnset',{},'pPeak',{});
tic;
for c = 1:n
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = CELL(c).obs(:)'
        p = regexp(OBS(o).label,'/','split');
        if numel(p) < 3, continue; end
        fp = REC(OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp,'roi',str2double(p{end}), ...
            'recName',strjoin(p(3:end-1),'/'),'group',p{1},'recDate',p{2}); %#ok<SAGROW>
    end
    if isempty(O), continue; end
    try
        [~,~,st] = temporal_phase_cell_fig_260812(O, P);
    catch ME
        fprintf(2,'  cell %d failed: %s\n', c, ME.message); continue;
    end
    dffPhase(c,:)  = st.dffPhase;
    ratePhase(c,:) = st.ratePhase;
    occPhase(c,:)  = st.occPhase;
    if isempty(featPeak)
        featTau   = st.featTau;
        featPeak  = nan(n, numel(featTau));
        featOnset = nan(n, numel(featTau));
    end
    featOnset(c,:) = st.featOnsetDff; %#ok<SAGROW>
    % Lengths are now identical for every cell (fixed P.featFs). Assert rather
    % than skip: silently dropping a mismatched row is what lost 24 cells.
    assert(numel(st.featPeakDff) == size(featPeak,2), ...
        'cell %d returned %d feature samples, expected %d', ...
        c, numel(st.featPeakDff), size(featPeak,2));
    featPeak(c,:) = st.featPeakDff; %#ok<SAGROW>
    S(end+1) = struct('cell',cellNum(c), 'key',CELL(c).key, 'group',O(1).group, ...
        'session',[O(1).group '_' O(1).recDate], 'nRec',numel(O), ...
        'nEvents',st.nEvents, 'nSpikes',st.nSpikes, 'rateHz',st.rateHz, ...
        'durTot',st.durTot, 'IBI',st.IBI, 'fps',st.fps, 'logZ',st.logZ, ...
        'Rbar',st.Rbar, 'nEff',st.nEff, 'prefPhase',st.prefPhase, ...
        'x',st.xy(1), 'y',st.xy(2), 'archiveLogZ',rayArchive(c), ...
        'pOnset',st.pOnset, 'pPeak',st.pPeak); %#ok<SAGROW>
    if mod(c,25)==0, fprintf('  %d/%d  (%.1f s)\n', c, n, toc); end
end
fprintf('extracted %d of %d cells in %.1f min\n', numel(S), n, toc/60);

keepRow   = ~all(isnan(dffPhase),2);
dffPhase  = dffPhase(keepRow,:);
ratePhase = ratePhase(keepRow,:);
occPhase  = occPhase(keepRow,:);
featPeak  = featPeak(keepRow,:);
featOnset = featOnset(keepRow,:);
phaseCtrs = (linspace(0,2*pi,nPhaseBins+1));
phaseCtrs = (phaseCtrs(1:end-1)+phaseCtrs(2:end))/2;

%% ===================== SHAPE NORMALISATION =====================
% Cluster on SHAPE, not amplitude: a bright cell and a dim cell with the same
% respiratory pattern belong together, and amplitude is kept separately as a
% scalar (modDepth / rateHz) so nothing is lost.
zscore_rows = @(M) (M - mean(M,2,'omitnan')) ./ max(std(M,0,2,'omitnan'), eps);
featShape  = zscore_rows(featPeak);       % ACTIVE block: peak-triggered
featShapeOn= zscore_rows(featOnset);     % same, onset-triggered
dffShape  = zscore_rows(dffPhase);
rateShape = zscore_rows(ratePhase);
featAmp   = max(featPeak,[],2) - min(featPeak,[],2);   % kept separately
modDepthDff  = (max(dffPhase,[],2)  - min(dffPhase,[],2));
modDepthRate = (max(ratePhase,[],2) - min(ratePhase,[],2)) ./ ...
                max(mean(ratePhase,2,'omitnan'), eps);

T = struct2table(S);
T.modDepthDff  = modDepthDff;
T.modDepthRate = modDepthRate;
T.prefPhaseDeg = mod(rad2deg(T.prefPhase),360);
T.featAmp      = featAmp;

%% ===================== ACTIVE GATE =====================
% Applied to T and to every feature matrix in one place, because they are row-
% aligned and a filter applied to some but not all of them silently mispairs every
% cell after the first drop.
%
% NOTE this can only ever REMOVE cells. The old ">5 events" rule was already
% applied upstream when event_latency_data.mat was built (326 of 1955 observations
% became cell members), so ROIs that fired e.g. 4 events in a 49 s recording --
% 4.9 events/min, active by this definition -- were discarded before we see them.
% What comes out here is the INTERSECTION of the two definitions, not a clean
% rate-defined population. Rebuilding the registry from OBS would be needed for
% that, and is deliberately not done here.
isActive = T.rateHz >= activeMinRateHz;
fprintf('\nACTIVE gate (>= %.3g events/min): %d of %d cells pass, %d dropped\n', ...
        activeMinRateHz*60, nnz(isActive), numel(isActive), nnz(~isActive));
if any(~isActive)
    for i = find(~isActive)'
        fprintf('   dropped cell %3d  %-10s %5.2f ev/min  (%d ev in %.0f s)\n', ...
                T.cell(i), T.group{i}, T.rateHz(i)*60, T.nSpikes(i), T.durTot(i));
    end
end
T         = T(isActive,:);
dffPhase  = dffPhase(isActive,:);   ratePhase = ratePhase(isActive,:);
occPhase  = occPhase(isActive,:);   featPeak  = featPeak(isActive,:);
featOnset = featOnset(isActive,:);  featShape = featShape(isActive,:);
featShapeOn = featShapeOn(isActive,:);
dffShape  = dffShape(isActive,:);   rateShape = rateShape(isActive,:);
assert(height(T) == size(featPeak,1) && height(T) == size(dffPhase,1), ...
       'row alignment broken by the active gate');

save(fullfile(outDir,'pop_features.mat'), 'T','featPeak','featShape','featTau', ...
     'featOnset','featShapeOn', ...
     'dffPhase','ratePhase','occPhase','dffShape','rateShape','phaseCtrs', ...
     'nPhaseBins','featWin_sec','P','-v7.3');
writetable(T, fullfile(outDir,'pop_features_scalars.csv'));

fprintf('\n---- feature matrix ----\n');
fprintf('  cells            : %d\n', height(T));
fprintf('  ACTIVE feature   : peak-triggered dF/F, +/-%.1f s, %d samples\n', ...
        featWin_sec, numel(featTau));
fprintf('  phase bins       : %d  (%.0f deg)  [supplementary]\n', nPhaseBins, 360/nPhaseBins);
fprintf('  sessions         : %d\n', numel(unique(T.session)));
fprintf('  groups           : %s\n', strjoin(unique(T.group)', ', '));
fprintf('  IBI range        : %.2f - %.2f s\n', min(T.IBI), max(T.IBI));
fprintf('  events per cell  : %d - %d\n', min(T.nSpikes), max(T.nSpikes));
fprintf('  saved            : %s\n', fullfile(outDir,'pop_features.mat'));
