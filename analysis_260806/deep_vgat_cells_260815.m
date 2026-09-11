% deep_vgat_cells_260815.m
% -----------------------------------------------------------------------
%  The DEEP recordings, analysed SEPARATELY from the surface population.
%
%  WHY SEPARATE, AND WHY NOT POOLED
%  event_latency_data.mat contains only recordings at z <= 100 um; the archive
%  also holds deep ones (z 250-320 um in Vgat/0730\deep and \cell2) that were
%  never registered. They are analysed here with the SAME statistics as the
%  surface cells, but they are not merged into that population, for two reasons:
%
%   1. The breath-locked OPTICAL ARTIFACT was characterised in exactly this deep
%      Vgat data: the top SVD mode of the movie tracks chest breathing at r=0.96,
%      survives motion correction, and is ~4x stronger in background than in
%      somata. It produces breath phase-locking that is not neural. Any
%      phase-locking reported below has to be defended against that first --
%      pooling it with surface cells would launder it into the main result.
%   2. Depth changes the measurement, not just the sample: penetration falls off
%      with an exponential length constant (WT ~25 um), so at 250-320 um the SNR
%      and the point spread are not comparable to a 0-100 um recording.
%
%  CELL IDENTITY: no cross-FOV matcher covers these recordings, so one cell per
%  (recording, ROI). Where the same neuron appears in two of them this OVERCOUNTS
%  -- stated rather than silently assumed.
%
%  Runqi Zhang / 2026-08-15
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
sumRoot   = 'D:\Ventral_surface_summary';
outDir    = fullfile(sumRoot,'deep_cells_260815');
depthCut  = 120;              % um. Above this the registry has nothing.
activeMinRateHz = 2/60;       % same gate as the surface population
nShuffle  = 1200;             % same as the surface figures, so p is comparable

% Deep recordings, given explicitly rather than globbed: the same recording is
% duplicated under \deep\ and \cell2\, and globbing would analyse it twice.
DEEP = { 'Vgat','0730', fullfile(sumRoot,'Vgat','0730','deep','roi1_2.4x_x1300y900_z250_3000f_30lp_00001'),  250
         'Vgat','0730', fullfile(sumRoot,'Vgat','0730','deep','roi1_3x_x1350y850_z265_3000f_31lp_00001'),    265
         'Vgat','0730', fullfile(sumRoot,'Vgat','0730','deep','roi1_2.4x_x1000y1000_z320_3000f_32lp_00001'), 320 };
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

P = struct('doCoh',false,'nDrop',30,'fallback_fps',30,'TW_spec',6,'alpha_sig',0.01, ...
    'minSpikes',2,'ca_lag_sec',0,'f_breath_search',[0.2 4],'fwhm_factor',0.6, ...
    'min_bw',0.05,'fmin',0.05,'fmax',15,'trigWin_sec',[],'trigWinIBI',2, ...
    'ylim_dff',[],'ylim_epc',[],'histBinFrames',2,'nShuffle',nShuffle,'shiftMinCyc',3, ...
    'pad_um',20,'clip_pct',[0.5 99.9],'scalebar_um',50,'rayPhaseBins',36,'featFs',30, ...
    'gamma_val',1,'PixelSizeBase',1.7778,'outlineLW',0.8,'sortMode','none', ...
    'dffColor',[0.2 0.7 0.2],'onsetCol',[0.9 0.1 0.1],'peakCol',[0.35 0.75 1], ...
    'statsOnly',true);

%% ===================== MEASURE =====================
Row = struct('group',{},'date',{},'rec',{},'z_um',{},'roi',{}, ...
             'nSpikes',{},'durTot',{},'rateHz',{},'IBI',{},'logZ',{}, ...
             'pOnset',{},'pPeak',{},'nEvents',{});
tAll = tic;
for r = 1:size(DEEP,1)
    fp = DEEP{r,3};
    if ~isfolder(fp), fprintf(2,'missing: %s\n', fp); continue; end
    sam = dir(fullfile(fp,'*_cpSAM_output.mat'));
    if isempty(sam), fprintf(2,'no cpSAM: %s\n', fp); continue; end
    S = load(fullfile(sam(1).folder,sam(1).name),'maskL');
    rois = setdiff(unique(S.maskL(:)),0).';
    fprintf('%s/%s  z=%d  %s : %d ROIs\n', DEEP{r,1}, DEEP{r,2}, DEEP{r,4}, ...
            'deep', numel(rois));
    for q = rois
        O = struct('folder',fp,'roi',double(q),'recName',DEEP{r,3}, ...
                   'group',DEEP{r,1},'recDate',DEEP{r,2});
        try
            [~,~,st] = temporal_phase_cell_fig_260812(O, P);
        catch ME
            fprintf(2,'   ROI %d failed: %s\n', q, ME.message);  continue;
        end
        Row(end+1) = struct('group',DEEP{r,1},'date',DEEP{r,2}, ...
            'rec',string(DEEP{r,3}),'z_um',DEEP{r,4},'roi',double(q), ...
            'nSpikes',st.nSpikes,'durTot',st.durTot,'rateHz',st.rateHz, ...
            'IBI',st.IBI,'logZ',st.logZ,'pOnset',st.pOnset,'pPeak',st.pPeak, ...
            'nEvents',st.nEvents); %#ok<SAGROW>
    end
end
fprintf('\nmeasured %d ROIs in %.1f min\n\n', numel(Row), toc(tAll)/60);
assert(~isempty(Row), 'nothing measured');

T = struct2table(Row);
T.active = T.rateHz >= activeMinRateHz;
T.sig05  = (T.pOnset < 0.05)  | (T.pPeak < 0.05);
T.sig01  = (T.pOnset < 0.01)  | (T.pPeak < 0.01);
T.sigRay = T.logZ >= 1.93;                       % alpha = 0.001
writetable(T, fullfile(outDir,'deep_cells_stats.csv'));

%% ===================== REPORT =====================
uz = unique(T.z_um);
fprintf('%-8s %-6s %6s %8s %8s %10s %10s %10s\n', ...
        'group','z(um)','nROI','active','%active','sig p<.05','sig p<.01','logZ>=1.93');
fprintf('%s\n', repmat('-',1,74));
for i = 1:numel(uz)
    m = T.z_um == uz(i);
    a = m & T.active;
    fprintf('%-8s %-6d %6d %8d %7.1f%% %10d %10d %10d\n', T.group{find(m,1)}, uz(i), ...
        nnz(m), nnz(a), 100*nnz(a)/max(nnz(m),1), ...
        nnz(a & T.sig05), nnz(a & T.sig01), nnz(a & T.sigRay));
end
a = T.active;
fprintf('%s\n', repmat('-',1,74));
fprintf('%-8s %-6s %6d %8d %7.1f%% %10d %10d %10d\n','ALL DEEP','', height(T), nnz(a), ...
    100*nnz(a)/height(T), nnz(a & T.sig05), nnz(a & T.sig01), nnz(a & T.sigRay));

fprintf(['\nSIGNIFICANCE HERE IS NOT EVIDENCE OF NEURAL PHASE-LOCKING.\n' ...
    'These are the depths where the breath-locked optical artifact was\n' ...
    'characterised (top SVD mode r=0.96 with chest breathing, survives motion\n' ...
    'correction, ~4x stronger in background than in somata). Any locking above\n' ...
    'has to clear a spatial artifact test -- core/ring contrast, or background\n' ...
    'ROIs showing the same phase -- before it means anything.\n']);

save(fullfile(outDir,'deep_cells_260815.mat'),'T','DEEP','activeMinRateHz','nShuffle');
fprintf('\nsaved -> %s\n', outDir);
