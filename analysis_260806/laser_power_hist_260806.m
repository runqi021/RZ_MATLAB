% laser_power_hist_260806.m
% -----------------------------------------------------------------------
%  Distribution of imaging laser power over every recording in this analysis,
%  in MILLIWATTS at the sample -- not in ScanImage percent.
%
%  WHY mW AND NOT %.  ScanImage % is a Pockels setpoint; the %->mW curve is a
%  property of (rig, date), and the rig was re-measured on 2026-07-23 when the
%  curve MOVED.  The two tables differ by ~4x at 1 %.  Two-photon signal goes as
%  P^2, so quoting % across sessions that straddle that date is not a scale error,
%  it is a different curve.  Every conversion here goes through
%  laser_power_calibration(pct, acqDate), which picks the table BY DATE and uses
%  pchip (the curve is strongly superlinear at the low end, so linear interpolation
%  is wrong there).
%
%  The acquisition date comes from the session folder recorded in the meta's
%  source_tif ('...\260721_Sert_soma_G8s\...' -> 2026-07-21), NOT from file
%  timestamps, which move whenever data is copied between drives.
%
%  Output: laser_power_hist.png / .pdf  (+ a per-recording csv)
%
%  Runqi Zhang / 2026-08-06

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot, fullfile(repoRoot,'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
GROUPS = { ...
    'Vglut2',     'D:\Ventral_surface_summary\Vglut2'     ; ...
    'Vgat',       'D:\Ventral_surface_summary\Vgat'       ; ...
    'Sst',        'D:\Ventral_surface_summary\Sst'        ; ...
    'ChAT',       'D:\Ventral_surface_summary\ChAT'       ; ...
    'Sert',       'D:\Ventral_surface_summary\Sert'       };
% The 260806/260807 Sst sessions are NOT listed separately any more: since
% archive_sst_260806_260807 copied them to Sst\0806 and Sst\0807 they arrive
% through the Sst root, and listing the C: paths as well counted them twice.
outDir  = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806';
binW_mW = 5;               % histogram bin width, mW
doSave  = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

pct = []; mW = []; grp = {}; rec = {}; dt = datetime.empty(0,1); tbl = {};
for g = 1:size(GROUPS,1)
    hits = dir(fullfile(GROUPS{g,2}, '**', 'breath_peak_pc1.mat'));
    for h = 1:numel(hits)
        fp = hits(h).folder;
        mh = dir(fullfile(fp,'*_ch1_meta.mat'));
        if isempty(mh) || isempty(dir(fullfile(fp,'*_ch1_dFF.mat'))), continue; end
        M = load(fullfile(mh(1).folder, mh(1).name));
        if ~isfield(M,'laserPower_pct') || isempty(M.laserPower_pct), continue; end

        % acquisition date from the session folder named in source_tif
        d = NaT;
        if isfield(M,'source_tif')
            tok = regexp(char(M.source_tif), '[\\/](\d{6})_', 'tokens', 'once');
            if ~isempty(tok), d = datetime(tok{1}, 'InputFormat','yyMMdd'); end
        end
        if isnat(d)     % fall back to the archive's own MMDD folder
            pp = strsplit(fp, filesep);
            k  = find(~cellfun(@isempty, regexp(pp,'^\d{4}$','once')), 1, 'last');
            if ~isempty(k), d = datetime(['26' pp{k}], 'InputFormat','yyMMdd'); end
        end
        if isnat(d), warning('no date for %s, skipped', fp); continue; end

        [w, cal] = laser_power_calibration(double(M.laserPower_pct), d);
        pct(end+1) = double(M.laserPower_pct); %#ok<SAGROW>
        mW(end+1)  = w;                        %#ok<SAGROW>
        grp{end+1} = GROUPS{g,1};              %#ok<SAGROW>
        [~,rn]     = fileparts(fp); rec{end+1} = rn; %#ok<SAGROW>
        dt(end+1)  = d;                        %#ok<SAGROW>
        tbl{end+1} = cal.name;                 %#ok<SAGROW>
    end
end
pct = pct(:); mW = mW(:); dt = dt(:);
isPost = strcmp(tbl(:), 'post_260723');
fprintf('%d recordings | %d on post_260723 table, %d on pre_260723\n', ...
        numel(mW), nnz(isPost), nnz(~isPost));
fprintf('percent : %.1f - %.1f (median %.1f)\n', min(pct), max(pct), median(pct));
fprintf('mW      : %.1f - %.1f (median %.1f)\n', min(mW),  max(mW),  median(mW));

%% ---- figure ----
fh = figure('Color','w','Position',[60 80 1250 480]);
tl = tiledlayout(fh, 1, 3, 'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1);
histogram(ax, pct, 'BinWidth',1, 'FaceColor',[0.4 0.4 0.4]);
xlabel(ax,'ScanImage setpoint (%)'); ylabel(ax,'recordings'); box(ax,'off');
set(ax,'TickDir','out'); title(ax,'raw setpoint (NOT comparable across dates)');

ax = nexttile(tl,2); hold(ax,'on');
edges = 0:binW_mW:ceil(max(mW)/binW_mW)*binW_mW;
histogram(ax, mW(~isPost), 'BinEdges',edges, 'FaceColor',[0.85 0.3 0.1], 'FaceAlpha',0.7);
histogram(ax, mW( isPost), 'BinEdges',edges, 'FaceColor',[0.1 0.4 0.85], 'FaceAlpha',0.7);
xlabel(ax,'power at sample (mW)'); ylabel(ax,'recordings'); box(ax,'off');
set(ax,'TickDir','out');
legend(ax, {sprintf('pre 260723 (%d)',nnz(~isPost)), sprintf('post 260723 (%d)',nnz(isPost))}, ...
       'Box','off','Location','northeast');
title(ax,'converted to mW, per-date calibration');

ax = nexttile(tl,3); hold(ax,'on');
uG = unique(grp,'stable');
for k = 1:numel(uG)
    m = strcmp(grp, uG{k});
    scatter(ax, dt(m), mW(m), 26, 'filled', 'MarkerFaceAlpha',0.75);
end
xline(ax, datetime(2026,7,23), 'k--', 'LineWidth',1);
ylabel(ax,'mW'); box(ax,'off'); set(ax,'TickDir','out');
legend(ax, uG, 'Box','off','Location','northwest','Interpreter','none','FontSize',7);
title(ax,'by session (dashed = calibration change)');

title(tl, sprintf('imaging laser power, %d recordings   |   median %.0f mW (range %.0f-%.0f)', ...
      numel(mW), median(mW), min(mW), max(mW)), 'FontWeight','bold');

%% ---- the histogram on its own ----
fh2 = figure('Color','w','Position',[80 120 720 520]);
ax  = axes(fh2); hold(ax,'on');
histogram(ax, mW, 'BinEdges',edges, 'FaceColor',[0.30 0.45 0.70], 'FaceAlpha',0.9, ...
          'EdgeColor','w');
xline(ax, median(mW), 'r-', 'LineWidth', 1.5);
text(ax, median(mW), 0, sprintf('  median %.0f mW', median(mW)), ...
     'Color','r', 'VerticalAlignment','bottom', 'FontWeight','bold');
xlabel(ax, 'power at the sample (mW)'); ylabel(ax, 'recordings');
set(ax, 'TickDir','out'); box(ax,'off');
title(ax, sprintf(['imaging power, %d recordings   (%%->mW per acquisition date: ' ...
                   '%d post-260723, %d pre)'], numel(mW), nnz(isPost), nnz(~isPost)));

%% ---- save ----
if doSave
    exportgraphics(fh2, fullfile(outDir,'laser_power_hist_mW.png'), 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh2, fullfile(outDir,'laser_power_hist_mW.pdf'), 'ContentType','vector', 'BackgroundColor','white');
end
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir, 'laser_power_hist');
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    writetable(table(grp(:), rec(:), dt, pct, mW, tbl(:), ...
               'VariableNames',{'group','recording','date','pct','mW','cal_table'}), ...
               [base '.csv']);
    fprintf('saved %s.{png,pdf,csv}\n', base);
end
