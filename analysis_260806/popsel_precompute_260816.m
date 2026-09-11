% popsel_precompute_260816.m
% -----------------------------------------------------------------------
%  Cache everything the include/exclude GUI needs, so that toggling a cell
%  re-averages instantly instead of recomputing.
%
%  NORMALISED CYCLE AXIS. Every triggered curve is resampled onto tau/IBI, i.e.
%  -1 .. +1 in units of the cell's OWN inter-breath interval. This is not a
%  cosmetic choice: IBI spans 0.43-2.83 s across the archive (1.74-2.63 s within
%  Sert alone), so averaging in seconds lets slow breathers stretch the window
%  and smears the transient. It is the same "windows in breath cycles" rule the
%  260806 heatmaps use.
%
%  POOLED RAYLEIGH, PRECOMPUTED PER CELL. The GUI needs a population statistic
%  that updates on every click, so the per-cell pieces are cached in a form that
%  POOLS BY SUMMATION:
%      evSum(b)  events in phase bin b        (from ratePhase .* occPhase / fps)
%      occ(b)    frames the recording spent in bin b
%  Each event carries weight 1/occ(b) from ITS OWN recording, exactly as the
%  ventral analysis defines it, so the pooled resultant is just a sum over the
%  selected cells and Kish's n_eff follows from the same sums. Nothing has to be
%  recomputed when the selection changes.
%
%  Runqi Zhang / 2026-08-16
% -----------------------------------------------------------------------

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));
addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(fullfile(repoRoot,'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
sumRoot  = 'D:\Ventral_surface_summary';
regFile  = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
figRoot  = fullfile(sumRoot,'per-cell-summary_active_260812');
outDir   = fullfile(sumRoot,'popsel_260816');

% One cache per group, written to popsel_cache_<GROUP>.mat. Listing several here
% caches them all in one pass; the GUI opens whichever one it is pointed at.
GROUPS   = {'IO'};        % e.g. {'Sert','IO','Vglut2','Vgat','Sst','ChAT'}
nTauN    = 201;           % points on the normalised -1..+1 IBI grid
nHistN   = 41;            % bins for the normalised event histogram
nShuffle = 0;             % per-cell permutation p is already in pop_features
overwrite = false;        % false = skip a group whose cache already exists
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end

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

tauN  = linspace(-1, 1, nTauN);                       % units of IBI
edgeN = linspace(-1, 1, nHistN+1);
ctrN  = (edgeN(1:end-1) + edgeN(2:end))/2;

for gI = 1:numel(GROUPS)
GROUP = GROUPS{gI};
cacheOut = fullfile(outDir, sprintf('popsel_cache_%s.mat', GROUP));
if isfile(cacheOut) && ~overwrite
    fprintf('%s: cache exists, skipping (set overwrite = true to rebuild)\n', GROUP);
    continue;
end

%% ---- which cells ----
ids = [];
for c = 1:numel(obsOf)
    if isempty(obsOf{c}), continue; end
    p = regexp(D.OBS(obsOf{c}(1)).label,'/','split');
    if strcmp(p{1}, GROUP), ids(end+1) = c; end %#ok<SAGROW>
end
fprintf('\n%s: %d cells in the registry\n', GROUP, numel(ids));

C = struct('cell',{},'group',{},'date',{},'label',{},'png',{},'nRec',{}, ...
           'IBI',{},'fps',{},'nSpikes',{},'rateHz',{},'logZ',{}, ...
           'pOnset',{},'pPeak',{},'dffPeakN',{},'dffOnsetN',{}, ...
           'histPeakN',{},'histOnsetN',{},'evSum',{},'occ',{},'phaseCtrs',{});

tic;
for k = 1:numel(ids)
    c = ids(k);
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = obsOf{c}(:)'
        p = regexp(D.OBS(o).label,'/','split');
        fp = D.REC(D.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp,'roi',str2double(p{end}), ...
            'recName',strjoin(p(3:end-1),'/'),'group',p{1},'recDate',p{2}); %#ok<SAGROW>
    end
    if isempty(O), fprintf(2,'  cell %d: no usable recording\n', c); continue; end
    try
        [~,~,st] = temporal_phase_cell_fig_260812(O, P);
    catch ME
        fprintf(2,'  cell %d failed: %s\n', c, ME.message);  continue;
    end

    % --- resample onto the normalised cycle axis ---
    tI = st.tau / st.IBI;                       % seconds -> IBI units
    dffPeakN  = interp1(tI, st.muPeak,  tauN, 'linear', NaN);
    dffOnsetN = interp1(tI, st.muOnset, tauN, 'linear', NaN);
    hI = st.histCtrs / st.IBI;
    histPeakN  = interp1(hI, st.histPeak,  ctrN, 'linear', NaN);
    histOnsetN = interp1(hI, st.histOnset, ctrN, 'linear', NaN);

    % --- pieces for the pooled occupancy-weighted Rayleigh ---
    occ   = st.occPhase(:).';
    evSum = st.ratePhase(:).' .* occ / st.fps;   % events per phase bin
    evSum(~isfinite(evSum)) = 0;
    occ(~isfinite(occ) | occ <= 0) = NaN;        % empty bins carry no weight

    pngf = fullfile(figRoot, sprintf('%s_%s_cell%03d.png', O(1).group, O(1).recDate, c));
    C(end+1) = struct('cell',c,'group',O(1).group,'date',O(1).recDate, ...
        'label',D.OBS(obsOf{c}(1)).label,'png',pngf,'nRec',numel(O), ...
        'IBI',st.IBI,'fps',st.fps,'nSpikes',st.nSpikes,'rateHz',st.rateHz, ...
        'logZ',st.logZ,'pOnset',st.pOnset,'pPeak',st.pPeak, ...
        'dffPeakN',dffPeakN,'dffOnsetN',dffOnsetN, ...
        'histPeakN',histPeakN,'histOnsetN',histOnsetN, ...
        'evSum',evSum,'occ',occ,'phaseCtrs',st.phaseCtrs(:).'); %#ok<SAGROW>
    if mod(k,10)==0, fprintf('  %d/%d  (%.0f s)\n', k, numel(ids), toc); end
end
fprintf('cached %d of %d cells in %.1f min\n', numel(C), numel(ids), toc/60);

nMissPng = sum(~arrayfun(@(x) isfile(x.png), C));
if nMissPng > 0
    fprintf(2,'%d cell(s) have no rendered figure -- the GUI will show a placeholder\n', nMissPng);
end

save(cacheOut, 'C','tauN','ctrN','edgeN','GROUP','P','-v7.3');
fprintf('saved -> %s\n', cacheOut);
end   % GROUPS loop
