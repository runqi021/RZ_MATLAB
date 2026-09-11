% popsel_group_figures_260831.m
% -----------------------------------------------------------------------
%  ONE FIGURE PER POPSEL GROUP, in the same 2x3 style as
%  popsel_population_figure_260831.m, with that group's own trigger-averaged
%  BREATHING WAVEFORM overlaid on every panel.
%
%      row 1  INSPIRATION-ONSET triggered      row 2  BREATH-PEAK triggered
%      col 1  raw dF/F    col 2  z-scored dF/F    col 3  event histogram
%      grey   breathing waveform, min-max normalised to 0-1, right-hand axis
%
%  X axis is ABSOLUTE TIME, +/-1.5 s, matching the combined figure.
%
%  WHY THE BREATH IS NORMALISED 0-1. The breath trace is SVD PC1 (or the
%  fixed-metric displacement), whose units are arbitrary and whose SIGN is
%  arbitrary per video -- see project_breath_sign_is_analysis_immune. Nothing
%  here depends on its amplitude, only on its timing relative to the trigger, so
%  it is min-max scaled over the plotted window and drawn on its own axis. Do not
%  read amplitude off the grey curve.
%
%  HOW THE BREATH AVERAGE IS BUILT. Per RECORDING, not per cell: a recording
%  contributes one waveform however many included cells it holds, so a FOV with
%  eight cells does not outvote one with a single cell. The trace, the event
%  indices and the truncation are loaded exactly as
%  temporal_phase_cell_fig_260812.m's load_obs_local does it, because a breath
%  curve aligned by a different rule than the calcium would be a lie about the
%  relative timing:
%     * bw = detrend(BP.breath), first P.nDrop samples dropped, then de-meaned
%     * PEAK   = insp_onset_idx  from breath_peak_pc1.mat      (yes, that name)
%     * ONSET  = insp_start_idx  from breath_insp_start_pc1.mat (the feet)
%     * T = min(dFF, breath, events), applied before any epoching
%     * Vglut2/1124 used a rising-edge trigger, so its breath and events are one
%       frame late and are shifted back. That gate is on the GENOTYPE FOLDER, not
%       the group label: the IO sites of that session live under Vglut2\1124\IO\
%       and carry group 'IO', so testing the label alone would skip exactly those.
%
%  Runqi Zhang / 2026-08-31
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ===================== USER-EDITABLE =====================
sumRoot = 'D:\Ventral_surface_summary';
inDir   = fullfile(sumRoot,'popsel_260816');
outDir  = fullfile(sumRoot,'popsel_population_260831');
regFile = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');

GRP = { ...
  'IO',                  'IO',                 [0    0    0   ]
  'Vglut2',              'Vglut2',             [0.10 0.65 0.20]
  'Vglut2_vagotomized',  'Vglut2 vagotomised', [0.10 0.65 0.20]
  'Sert',                'Sert',               [0.90 0.45 0.10]
  'Sert_vagotomized',    'Sert vagotomised',   [0.90 0.45 0.10]
  'Sst',                 'Sst',                [0.55 0.20 0.75]
  'Vgat',                'Vgat',               [0.10 0.30 0.85] };

XLIM   = [-1.5 1.5];
nDrop  = 30;            % must match P.nDrop in the cache build
fbFps  = 30;
LW     = 2.0;
BCOL   = [0.45 0.45 0.45];   % breath grey
%% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
D = load(regFile,'CELL','OBS','REC');
obsOf = pooled_obs_260814(D.CELL, D.OBS);

summary = {};
for gi = 1:size(GRP,1)
    nm = GRP{gi,1};
    cf = fullfile(inDir, sprintf('popsel_cache_%s.mat', nm));
    df = fullfile(inDir, sprintf('popsel_decisions_%s.csv', nm));
    if ~isfile(cf) || ~isfile(df)
        fprintf(2,'%-20s missing cache or decisions, skipped\n', nm); continue;
    end
    K = load(cf);  C = K.C;
    T = readtable(df,'TextType','string');
    assert(isfield(C,'dffOnsetS'), ...
        '%s cache has no absolute-time curves -- rebuild with popsel_run_260831.m', nm);

    inc = T.cell(T.decision == "include");
    [tf, loc] = ismember(inc, [C.cell].');
    sel = loc(tf);
    if isempty(sel), fprintf(2,'%-20s no included cells, skipped\n', nm); continue; end
    cells = [C(sel).cell];

    %% ---- breath average over the recordings those cells came from ----------
    folders = strings(0,1); labs = strings(0,1);
    for c = cells
        for o = obsOf{c}(:)'
            folders(end+1,1) = string(D.REC(D.OBS(o).rec).folder); %#ok<SAGROW>
            labs(end+1,1)    = string(D.OBS(o).label);             %#ok<SAGROW>
        end
    end
    [folders, ia] = unique(folders); labs = labs(ia);

    BOn = []; BPk = []; nRecUsed = 0;
    for r = 1:numel(folders)
        [bOn, bPk] = breath_trig_avg(char(folders(r)), char(labs(r)), K.tauS, nDrop, fbFps);
        if isempty(bOn) && isempty(bPk), continue; end
        if ~isempty(bOn), BOn = [BOn; bOn]; end %#ok<AGROW>
        if ~isempty(bPk), BPk = [BPk; bPk]; end %#ok<AGROW>
        nRecUsed = nRecUsed + 1;
    end
    bOnMu = norm01(mean(BOn,1,'omitnan'));
    bPkMu = norm01(mean(BPk,1,'omitnan'));

    cov = min([C(sel).histWinHalfS]);
    fprintf('%-20s %3d cells, %3d recordings   hist full to +/-%.2f s\n', ...
            nm, numel(sel), nRecUsed, cov);

    %% ---- figure ------------------------------------------------------------
    col = GRP{gi,3};
    fig = figure('Color','w','Units','pixels','Position',[60 60 1450 820],'Visible','off');
    tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

    PAN = { 'dffOnsetS','mean dF/F','raw dF/F','tauS'
            'dffOnsetZS','mean dF/F (z-scored)','z-scored dF/F','tauS'
            'histOnsetS','events per cycle','event histogram','ctrS'
            'dffPeakS','mean dF/F','raw dF/F','tauS'
            'dffPeakZS','mean dF/F (z-scored)','z-scored dF/F','tauS'
            'histPeakS','events per cycle','event histogram','ctrS' };
    ROWLAB = {'inspiration onset','breath peak'};

    for p = 1:6
        ax = nexttile(tl,p); row = ceil(p/3);
        M = cell2mat({C(sel).(PAN{p,1})}.');
        x = K.(PAN{p,4});
        mu = mean(M,1,'omitnan'); sd = std(M,0,1,'omitnan');
        nn = sum(~isnan(M),1);    se = sd ./ max(sqrt(nn),1);

        yyaxis(ax,'left'); hold(ax,'on');
        good = ~isnan(mu) & nn > 1;
        if any(good)
            fill(ax,[x(good) fliplr(x(good))],[mu(good)+se(good) fliplr(mu(good)-se(good))], ...
                 col,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
        end
        plot(ax,x,mu,'-','Color',col,'LineWidth',LW,'DisplayName',sprintf('%s (n=%d)',GRP{gi,2},numel(sel)));
        ax.YColor = [0.15 0.15 0.15];
        ylabel(ax,PAN{p,2});
        yl = ylim(ax);
        tc = [0.90 0.10 0.10]; if row == 2, tc = [0.20 0.50 0.95]; end
        plot(ax,[0 0],yl,'-','Color',tc,'LineWidth',1,'HandleVisibility','off');
        ylim(ax,yl);

        % ---- breath, normalised 0-1, own axis ----
        yyaxis(ax,'right');
        if row == 1, b = bOnMu; else, b = bPkMu; end
        if ~isempty(b) && any(isfinite(b))
            plot(ax,K.tauS,b,'-','Color',BCOL,'LineWidth',1.3,'DisplayName','breath (0-1)');
        end
        ylim(ax,[-0.05 1.35]); ax.YColor = BCOL;
        if mod(p,3) == 0, ylabel(ax,'breath (norm.)'); else, set(ax,'YTickLabel',[]); end

        yyaxis(ax,'left');
        % shade where the histogram no longer has every cell behind it
        if strcmp(PAN{p,4},'ctrS') && isfinite(cov) && cov < XLIM(2)
            for s = [-1 1]
                xa = sort([s*cov, s*XLIM(2)]);
                fill(ax,[xa(1) xa(2) xa(2) xa(1)],[yl(1) yl(1) yl(2) yl(2)], ...
                     [0.5 0.5 0.5],'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
            end
            text(ax,0,yl(2),sprintf('all cells only within \\pm%.2f s',cov), ...
                 'HorizontalAlignment','center','VerticalAlignment','top', ...
                 'FontSize',8,'Color',[0.35 0.35 0.35]);
            ylim(ax,yl);
        end
        xlim(ax,XLIM);
        xlabel(ax,sprintf('time from %s (s)',ROWLAB{row}));
        title(ax,sprintf('%s  -  %s triggered',PAN{p,3},ROWLAB{row}),'FontWeight','normal');
        grid(ax,'on'); ax.GridAlpha = 0.12; box(ax,'off');
        if p == 1, legend(ax,'Location','northwest','Box','off','FontSize',8); end
    end

    title(tl, sprintf(['%s   |   %d cells, %d recordings   |   grey = breathing, ' ...
        'min-max normalised (amplitude not meaningful)'], GRP{gi,2}, numel(sel), nRecUsed), ...
        'FontWeight','normal');

    stem = fullfile(outDir, sprintf('popsel_group_%s', nm));
    exportgraphics(fig,[stem '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','none');
    dd = dir([stem '.png']);
    if dd.bytes < 20000
        fprintf(2,'  %s PNG only %d bytes -- possible OpenGL blank stub\n', nm, dd.bytes);
    end
    close(fig);
    summary(end+1,:) = {nm, numel(sel), nRecUsed, cov, dd.bytes/1024}; %#ok<SAGROW>
end

Ts = cell2table(summary,'VariableNames', ...
    {'group','n_cells','n_recordings','hist_full_coverage_s','png_KB'});
writetable(Ts, fullfile(outDir,'popsel_group_figures.csv'));
fprintf('\nsaved %d figures -> %s\n', height(Ts), outDir);
disp(Ts);


% =========================================================================
function y = norm01(x)
%NORM01  Min-max to 0-1 over the finite samples. Sign and scale of the breath
% trace are arbitrary, so only the shape survives this and that is intended.
y = x;
lo = min(x(isfinite(x))); hi = max(x(isfinite(x)));
if isempty(lo) || hi <= lo, return; end
y = (x - lo) / (hi - lo);
end

% =========================================================================
function [bOn, bPk] = breath_trig_avg(folder, label, tauS, nDrop, fbFps)
%BREATH_TRIG_AVG  Onset- and peak-triggered mean breath waveform for one
% recording, on the tauS grid. Mirrors load_obs_local in
% temporal_phase_cell_fig_260812.m -- see the header of this file.
bOn = []; bPk = [];
bp = dir(fullfile(folder,'breath_peak_pc1.mat'));
ip = dir(fullfile(folder,'breath_insp_start_pc1.mat'));
df = dir(fullfile(folder,'*_ch1_dFF.mat'));
if isempty(bp) || isempty(df), return; end

fps = detect_session_fps(folder, fbFps);
BP  = load(fullfile(bp(1).folder, bp(1).name));
Dd  = load(fullfile(df(1).folder, df(1).name),'dFF');

bw = detrend(double(BP.breath(:)));
bw(1:min(nDrop,numel(bw))) = [];
bw = bw - mean(bw);
nB = numel(BP.breath);

ev = zeros(nB,1);
oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
ev(1:min(nDrop,numel(ev))) = [];

ev_foot = [];
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    ev_foot = zeros(nB,1); fi = round(IP.insp_start_idx(:));
    ev_foot(fi(fi>=1 & fi<=nB)) = 1;
    ev_foot(1:min(nDrop,numel(ev_foot))) = [];
end

% Vglut2/1124 rising-edge trigger: breath and events are one frame late. Gate on
% the genotype FOLDER, because that session's IO sites carry the label 'IO'.
p = regexp(label,'/','split');
geno = p{1}; dateStr = p{2};
if strcmpi(geno,'IO') && contains(folder,[filesep 'Vglut2' filesep],'IgnoreCase',true)
    geno = 'Vglut2';
end
if strcmpi(geno,'Vglut2') && strcmp(dateStr,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end

T = min([size(Dd.dFF,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
if T < 10, return; end
bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); else, ev_foot = zeros(T,1); end

% z-score, as the per-cell breath strip does (bwz in that file)
bwz = (bw - mean(bw)) / max(std(bw), eps);

w  = max(1, ceil(max(abs(tauS))*fps) + 1);
tk = (-w:w)/fps;
bPk = epoch_mean(bwz, find(ev>0),      w, tk, tauS);
bOn = epoch_mean(bwz, find(ev_foot>0), w, tk, tauS);
end

% =========================================================================
function mu = epoch_mean(x, idx, w, tk, tauS)
%EPOCH_MEAN  Mean of x cut around idx, regridded onto tauS. Cut wider than the
% target and interpolate inwards so interp1 never extrapolates at the edges.
mu = [];
idx = idx(idx-w >= 1 & idx+w <= numel(x));
if isempty(idx), return; end
E = zeros(numel(idx), 2*w+1);
for k = 1:numel(idx)
    E(k,:) = x(idx(k)-w : idx(k)+w);
end
mu = interp1(tk, mean(E,1,'omitnan'), tauS, 'linear', NaN);
end
