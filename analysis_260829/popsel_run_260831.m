% popsel_run_260831.m
% -----------------------------------------------------------------------
%  ONE ENTRY POINT for the population include/exclude curation.
%
%  Edit the four lines in the USER block, press run. It builds the cache if it
%  is missing and then opens the GUI. Nothing else to remember, and no shared
%  file has to be edited to point at a different group.
%
%  WHY THIS EXISTS RATHER THAN EDITING popsel_gui_260816.m EACH TIME
%    popsel_precompute_260816.m picks cells with strcmp(p{1}, GROUP) -- the
%    GENOTYPE alone. That is wrong wherever one genotype spans an experimental
%    manipulation: 'Sert' pools Sert/0721 (intact) with 0828+0829 (vagotomised),
%    and 'Vglut2' pools 0224/0728/0810/1124 with 0824 (vagotomised). Here a group
%    is a NAME plus a genotype plus a list of session dates, so the population is
%    whatever you say it is and the cache and decisions CSV are named after it.
%
%  ONE CACHE AND ONE CSV PER GROUP NAME
%    popsel_cache_<GROUP>.mat and popsel_decisions_<GROUP>.csv. Pick a new NAME
%    and you get a fresh curation; reuse a name and you resume the old one. The
%    GUI resumes from the CSV automatically, so stopping half way is safe.
%
%  Runqi Zhang / 2026-08-31
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));
addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(fullfile(repoRoot,'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
% The seven populations. {name, genotype, {session dates}} -- {} = every session
% of that genotype. Vagotomy is the independent variable, so the vagotomised
% sessions are their own group and are NOT pooled with the intact ones.
GROUPS = { ...
    'Vglut2',              'Vglut2', {'0224','0728','0810','1124'}   % intact
    'Vglut2_vagotomized',  'Vglut2', {'0824'}
    'Sert',                'Sert',   {'0721'}                        % intact
    'Sert_vagotomized',    'Sert',   {'0828','0829'}
    'Sst',                 'Sst',    {}
    'Vgat',                'Vgat',   {}           % 0323 + 0730, see below
    'IO',                  'IO',     {} };
% Vgat/0730 is INCLUDED, by RZ decision 2026-08-31, having been briefly excluded
% and put back. Know what it is: C:\260730_shiver_vgat_g8m -- a SHIVERER animal
% expressing GCaMP8m, not the GCaMP8s wildtype every other session uses, and the
% only session in the archive imaged deliberately deep. So the Vgat group mixes
% two genotypes AND two indicators, and 8 of its 12 cells come from that one
% animal. It is also where the entire "Vgat I" class came from. Without it Vgat
% is 4 cells from 0323, which is too thin to plot. Say which it is in any caption.

OPEN     = 'Sert';   % which group to open in the GUI. '' = build caches only.
TRIG     = 'onset';     % 'peak' | 'onset' -- which trigger the population panels use
REBUILD  = false;      % true = rebuild every cache even if it already exists
sumRoot  = 'D:\Ventral_surface_summary';
%% =========================================================

outDir  = fullfile(sumRoot,'popsel_260816');
regFile = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
figRoot = fullfile(sumRoot,'per-cell-summary_active_260812');
if ~isfolder(outDir), mkdir(outDir); end

for gi = 1:size(GROUPS,1)
    nm = GROUPS{gi,1};
    cf = fullfile(outDir, sprintf('popsel_cache_%s.mat', nm));
    if isfile(cf) && ~REBUILD
        fprintf('%-20s cache exists, reusing\n', nm);
    else
        build_cache(regFile, figRoot, cf, nm, GROUPS{gi,2}, GROUPS{gi,3});
    end
end

fprintf('\n---- decisions so far ----\n');
for gi = 1:size(GROUPS,1)
    nm = GROUPS{gi,1};
    df = fullfile(outDir, sprintf('popsel_decisions_%s.csv', nm));
    if isfile(df)
        Td = readtable(df,'TextType','string');
        fprintf('%-20s %3d include / %3d exclude / %3d undecided\n', nm, ...
            nnz(Td.decision=="include"), nnz(Td.decision=="exclude"), nnz(Td.decision=="undecided"));
    else
        K0 = load(fullfile(outDir, sprintf('popsel_cache_%s.mat', nm)), 'C');
        fprintf('%-20s not started (%d cells)\n', nm, numel(K0.C));
    end
end

if isempty(OPEN)
    fprintf('\nOPEN is empty -- caches built, GUI not launched.\n');
    return;
end
NAME   = OPEN;
cacheF = fullfile(outDir, sprintf('popsel_cache_%s.mat', NAME));
assert(isfile(cacheF), 'no cache for "%s" -- is it in GROUPS?', NAME);
fprintf('\nopening %s\n', NAME);

K = load(cacheF);
altFigDirs = { fullfile(figRoot,'_below_active_gate_2evmin')
               fullfile(figRoot,'_superseded_by_merge_260814') };
popsel_launch_260816(K.C, K.tauN, K.ctrN, NAME, TRIG, outDir, altFigDirs);


% =========================================================================
function build_cache(regFile, figRoot, cacheOut, NAME, GENOTYPE, DATES)
% Body copied from popsel_precompute_260816.m so the cached quantities are
% identical; only the cell SELECTION differs (genotype + dates, not genotype).
nTauN = 201; nHistN = 41; nShuffle = 0;

D = load(regFile,'CELL','OBS','REC');
obsOf = pooled_obs_260814(D.CELL, D.OBS);

P = struct('doCoh',false,'nDrop',30,'fallback_fps',30,'TW_spec',6,'alpha_sig',0.01, ...
    'minSpikes',2,'ca_lag_sec',0,'f_breath_search',[0.2 4],'fwhm_factor',0.6, ...
    'min_bw',0.05,'fmin',0.05,'fmax',15,'trigWin_sec',[],'trigWinIBI',2, ...
    'ylim_dff',[],'ylim_epc',[],'histBinFrames',2,'nShuffle',nShuffle,'shiftMinCyc',3, ...
    'pad_um',20,'clip_pct',[0.5 99.9],'scalebar_um',50,'rayPhaseBins',36,'featFs',30, ...
    'gamma_val',1,'PixelSizeBase',1.7778,'outlineLW',0.8,'sortMode','none', ...
    'dffColor',[0.2 0.7 0.2],'onsetCol',[0.9 0.1 0.1],'peakCol',[0.35 0.75 1], ...
    'statsOnly',true);

tauN  = linspace(-1, 1, nTauN);
edgeN = linspace(-1, 1, nHistN+1);
ctrN  = (edgeN(1:end-1) + edgeN(2:end))/2;

% ---- ABSOLUTE-TIME axis, in seconds ------------------------------------------
% The normalised tau/IBI axis above is what the GUI averages on. The population
% figure is asked for in real seconds, so a SECOND pass runs the same function
% with trigWin_sec set, and the curves land on this grid.
%
% trigWin_sec is a TOTAL width (temporal_phase_cell_fig_260812:63-69 halves it),
% so +/-WIN_S needs 2*WIN_S.
WIN_S = 1.5;
tauS  = linspace(-WIN_S, WIN_S, 301);
edgeS = linspace(-WIN_S, WIN_S, 61);
ctrS  = (edgeS(1:end-1) + edgeS(2:end))/2;
Ps = P;  Ps.trigWin_sec = 2*WIN_S;

% ---- which cells: genotype AND date ----
ids = [];
for c = 1:numel(obsOf)
    if isempty(obsOf{c}), continue; end
    p = regexp(D.OBS(obsOf{c}(1)).label,'/','split');
    if ~strcmp(p{1}, GENOTYPE), continue; end
    if ~isempty(DATES) && ~any(strcmp(p{2}, DATES)), continue; end
    ids(end+1) = c; %#ok<AGROW>
end
if isempty(DATES), dStr = 'all sessions'; else, dStr = strjoin(DATES, ', '); end
fprintf('%s: genotype %s, %s -> %d cells\n', NAME, GENOTYPE, dStr, numel(ids));
assert(~isempty(ids), '%s selected no cells -- check GENOTYPE and DATES', NAME);

% per-session breakdown, so a missing session is visible before curating
sess = strings(numel(ids),1);
for k = 1:numel(ids)
    p = regexp(D.OBS(obsOf{ids(k)}(1)).label,'/','split');  sess(k) = p{2};
end
for u = unique(sess)'
    fprintf('    %s/%s : %d cells\n', GENOTYPE, u, nnz(sess==u));
end

C = struct('cell',{},'group',{},'date',{},'label',{},'png',{},'nRec',{}, ...
           'IBI',{},'fps',{},'nSpikes',{},'rateHz',{},'logZ',{}, ...
           'pOnset',{},'pPeak',{},'dffPeakN',{},'dffOnsetN',{}, ...
           'dffPeakZN',{},'dffOnsetZN',{}, ...
           'histPeakN',{},'histOnsetN',{},'evSum',{},'occ',{},'phaseCtrs',{},'dOnPk',{}, ...
           'dffPeakS',{},'dffOnsetS',{},'dffPeakZS',{},'dffOnsetZS',{}, ...
           'histPeakS',{},'histOnsetS',{},'histWinHalfS',{});

tic;
for k = 1:numel(ids)
    c = ids(k);
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = obsOf{c}(:)'
        p = regexp(D.OBS(o).label,'/','split');
        fp = D.REC(D.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp,'roi',str2double(p{end}), ...
            'recName',strjoin(p(3:end-1),'/'),'group',p{1},'recDate',p{2}); %#ok<AGROW>
    end
    if isempty(O), fprintf(2,'  cell %d: no usable recording\n', c); continue; end
    try
        [~,~,st] = temporal_phase_cell_fig_260812(O, P);
    catch ME
        fprintf(2,'  cell %d failed: %s\n', c, ME.message);  continue;
    end

    tI = st.tau / st.IBI;
    dffPeakN  = interp1(tI, st.muPeak,  tauN, 'linear', NaN);
    dffOnsetN = interp1(tI, st.muOnset, tauN, 'linear', NaN);
    % Z-scored versions, from the SAME function that draws the per-cell figure,
    % rather than z-scored here from the raw curve. Raw dF/F is dominated by the
    % brightest cells; the z-scored trace gives every cell equal weight and shows
    % the SHAPE of the modulation. Both panels get plotted.
    dffPeakZN  = interp1(tI, st.muPeakZ,  tauN, 'linear', NaN);
    dffOnsetZN = interp1(tI, st.muOnsetZ, tauN, 'linear', NaN);
    hI = st.histCtrs / st.IBI;
    histPeakN  = interp1(hI, st.histPeak,  ctrN, 'linear', NaN);
    histOnsetN = interp1(hI, st.histOnset, ctrN, 'linear', NaN);

    % ---- absolute-time pass -------------------------------------------------
    % Same function, same cell, window fixed in SECONDS instead of in cycles.
    % The dF/F curves fill +/-WIN_S exactly. The HISTOGRAM does not: line 358 of
    % temporal_phase_cell_fig_260812 hardcodes its half-window to the cell's own
    % IBI and ignores trigWin_sec, so a cell with IBI 0.4 s carries histogram
    % only out to +/-0.4 s and is NaN beyond. histWinHalfS records that reach so
    % the figure can say how many cells actually support each part of the axis.
    try
        [~,~,ss] = temporal_phase_cell_fig_260812(O, Ps);
        dffOnsetS  = interp1(ss.tau, ss.muOnset,  tauS, 'linear', NaN);
        dffPeakS   = interp1(ss.tau, ss.muPeak,   tauS, 'linear', NaN);
        dffOnsetZS = interp1(ss.tau, ss.muOnsetZ, tauS, 'linear', NaN);
        dffPeakZS  = interp1(ss.tau, ss.muPeakZ,  tauS, 'linear', NaN);
        histOnsetS = interp1(ss.histCtrs, ss.histOnset, ctrS, 'linear', NaN);
        histPeakS  = interp1(ss.histCtrs, ss.histPeak,  ctrS, 'linear', NaN);
        histWinHalfS = ss.histWinHalf;
    catch ME
        fprintf(2,'  cell %d: absolute-time pass failed (%s)\n', c, ME.message);
        dffOnsetS = nan(size(tauS)); dffPeakS = nan(size(tauS));
        dffOnsetZS = nan(size(tauS)); dffPeakZS = nan(size(tauS));
        histOnsetS = nan(size(ctrS)); histPeakS = nan(size(ctrS));
        histWinHalfS = NaN;
    end

    occ   = st.occPhase(:).';
    evSum = st.ratePhase(:).' .* occ / st.fps;
    evSum(~isfinite(evSum)) = 0;
    occ(~isfinite(occ) | occ <= 0) = NaN;

    pngf = fullfile(figRoot, sprintf('%s_%s_cell%03d.png', O(1).group, O(1).recDate, c));
    C(end+1) = struct('cell',c,'group',O(1).group,'date',O(1).recDate, ...
        'label',D.OBS(obsOf{c}(1)).label,'png',pngf,'nRec',numel(O), ...
        'IBI',st.IBI,'fps',st.fps,'nSpikes',st.nSpikes,'rateHz',st.rateHz, ...
        'logZ',st.logZ,'pOnset',st.pOnset,'pPeak',st.pPeak, ...
        'dffPeakN',dffPeakN,'dffOnsetN',dffOnsetN, ...
        'dffPeakZN',dffPeakZN,'dffOnsetZN',dffOnsetZN, ...
        'histPeakN',histPeakN,'histOnsetN',histOnsetN, ...
        'evSum',evSum,'occ',occ,'phaseCtrs',st.phaseCtrs(:).', ...
        'dOnPk',st.dOnPk, ...
        'dffPeakS',dffPeakS,'dffOnsetS',dffOnsetS, ...
        'dffPeakZS',dffPeakZS,'dffOnsetZS',dffOnsetZS, ...
        'histPeakS',histPeakS,'histOnsetS',histOnsetS, ...
        'histWinHalfS',histWinHalfS); %#ok<AGROW>
    if mod(k,10)==0, fprintf('  %d/%d  (%.0f s)\n', k, numel(ids), toc); end
end
fprintf('cached %d of %d cells in %.1f min\n', numel(C), numel(ids), toc/60);

nMissPng = sum(~arrayfun(@(x) isfile(x.png), C));
if nMissPng > 0
    fprintf(2,['%d cell(s) have no rendered figure -- the GUI shows a placeholder.\n' ...
               '   Re-run per_cell_summary_260812.m for them.\n'], nMissPng);
end

GROUP = NAME; %#ok<NASGU>  kept for compatibility with popsel_cache_* readers
save(cacheOut, 'C','tauN','ctrN','edgeN','GROUP','P', ...
     'tauS','ctrS','edgeS','WIN_S','-v7.3');
fprintf('saved -> %s\n', cacheOut);
end
