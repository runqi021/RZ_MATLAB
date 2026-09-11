% ventral_pie_summary_260811.m
% -----------------------------------------------------------------------
%  Nine pies on one page. EVERY pie is a breakdown BY GENOTYPE -- same wedges,
%  same colours, one question per pie -- so the colour key is read once:
%
%      [1] animals imaged        [2] recorded cells      [3] active cells
%      [4] active, logZ >= 3     [5] active, logZ 2-3    [6] active, logZ 1-2
%      [7] active, logZ < 1      [8] phase-locked (>=1)  [9] legend
%
%  Every wedge is labelled with its COUNT -- N for animals, n for cells. No
%  percentages are drawn anywhere.
%
%  DEFINITIONS, inherited from the pipeline rather than invented here:
%    animal   = one session date (Genotype/MMDD).  One date = one animal.
%    cell     = cell_link identity where a matcher was run, else one cell per
%               ROI.  Masks marked `tossed` in curation are NOT cells.
%    recorded = every cell in a recording that has ca_spike_data, after the
%               coh_cfg exclusions (stage-zero FOV + the Vgat 120 um depth bar).
%    active   = pooled nnz(spike_train>0) > 5, the archive's own criterion.
%    logZ     = occupancy-weighted Rayleigh from the per-cell summary, which
%               computes it for exactly the active cells it draws.
%
%  IO IS NOT A SEPARATE ANIMAL SET. It is a SITE inside ChAT and Vglut2
%  sessions, so its animals are also counted under those genotypes and pie 1
%  sums to more than the number of mice. Stated on the figure.
%
%  Output: ventral_pie_summary.png / .pdf / .csv
%
%  Runqi Zhang / 2026-08-11

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot, fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
sumCsv    = fullfile(rootPath,'breath_time_summary_260808','breath_time_summary_cells.csv');
outDir    = fullfile(rootPath,'breath_trig_heatmap_260806');
scan_dirs = {'ChAT','Vglut2','Vgat','Sst','Sert'};
groups    = {'ChAT','Sst','Vglut2','Vgat','Sert','IO'};
gCol      = [0.85 0.10 0.10;    % ChAT   red
             0.55 0.20 0.75;    % Sst    purple
             0.10 0.65 0.20;    % Vglut2 green
             0.10 0.30 0.85;    % Vgat   blue
             0.90 0.45 0.10;    % Sert   orange
             0.35 0.35 0.35];   % IO     grey
zCuts     = [3 2 1];
% (no percentage labels anywhere: every wedge is labelled with its COUNT)
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
excl = coh_cfg_260727().excludeRecordings;

%% ---- scan the archive: recorded cells + animals, per genotype ----
cellKey = strings(0,1); cellGrp = strings(0,1); animKey = strings(0,1);
linkCache = containers.Map('KeyType','char','ValueType','any');
for s = 1:numel(scan_dirs)
    sname = scan_dirs{s};
    hits  = dir(fullfile(rootPath, sname, '**', 'ca_spike_data.mat'));
    for h = 1:numel(hits)
        fp  = hits(h).folder;
        rel = strrep(fp, fullfile(rootPath,sname), ''); rel = regexprep(rel,'^[\\/]+','');
        parts = regexp(rel,'[\\/]','split');
        if numel(parts) < 3, continue; end
        dateStr = parts{1};  site = parts{2};  rec = parts{end};
        if any(strcmp(rec, excl)), continue; end
        gname = sname;
        if ~isempty(regexpi(site,'IO','once')), gname = 'IO'; end

        C = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
        nROI = numel(C.roi_spikes);

        key = sprintf('%s/%s', sname, dateStr);
        if ~isKey(linkCache, key)
            lk = fullfile(rootPath, sname, dateStr, 'cell_pooled', 'cell_link.mat');
            if isfile(lk), Lk = load(lk,'link'); linkCache(key) = Lk.link.obsT;
            else,          linkCache(key) = [];
            end
        end
        obsT = linkCache(key);

        for r = 1:nROI
            ck = sprintf('roi:%s/%s/%s/%d', gname, dateStr, rec, r);
            if ~isempty(obsT)
                j = find(strcmp(string(obsT.rec_name), rec) & obsT.roi_index == r, 1);
                if ~isempty(j)
                    if strcmp(string(obsT.status(j)),'tossed'), continue; end
                    ck = sprintf('%s/%s#c%d', sname, dateStr, obsT.cell_id(j));
                end
            end
            cellKey(end+1,1) = string(ck);    %#ok<SAGROW>
            cellGrp(end+1,1) = string(gname); %#ok<SAGROW>
        end
        animKey(end+1,1) = string(sprintf('%s/%s', gname, dateStr)); %#ok<SAGROW>
    end
end
[uk, ia] = unique(cellKey);  ukGrp = cellGrp(ia);

%% ---- active cells + logZ ----
T = readtable(sumCsv);  aGrp = string(T.group);  aZ = T.ray_logZ;

nG = numel(groups);
S = struct('grp',{},'animals',{},'recorded',{},'active',{},'bands',{});
fprintf('%-7s %8s %10s %8s | %6s %6s %6s %6s\n','group','animals','recorded','active','>=3','2-3','1-2','<1');
for g = 1:nG
    m  = ukGrp == groups{g};
    am = startsWith(animKey, groups{g} + "/");
    a  = aGrp == groups{g};   z = aZ(a);
    bands = [nnz(z>=zCuts(1)), nnz(z>=zCuts(2) & z<zCuts(1)), ...
             nnz(z>=zCuts(3) & z<zCuts(2)), nnz(z<zCuts(3) | isnan(z))];
    S(g) = struct('grp',groups{g},'animals',numel(unique(animKey(am))), ...
                  'recorded',nnz(m),'active',nnz(a),'bands',bands);
    fprintf('%-7s %8d %10d %8d | %6d %6d %6d %6d\n', groups{g}, S(g).animals, ...
            S(g).recorded, S(g).active, bands);
end
B = vertcat(S.bands);
fprintf('%-7s %8d %10d %8d | %6d %6d %6d %6d\n','ALL', sum([S.animals]), ...
        sum([S.recorded]), sum([S.active]), sum(B,1));

%% ---- figure: 9 pies, every one split BY GENOTYPE ----
PIES = { 'animals imaged',            [S.animals],  'N'
         'recorded cells',            [S.recorded], '%'
         'active cells',              [S.active],   '%'
         'active, logZ \geq 3',       B(:,1).',     '%'
         'active, logZ 2-3',          B(:,2).',     '%'
         'active, logZ 1-2',          B(:,3).',     '%'
         'active, logZ < 1',          B(:,4).',     '%'
         'phase-locked (logZ \geq 1)',(B(:,1)+B(:,2)+B(:,3)).', '%' };

fh = figure('Color','w','Position',[30 30 1350 1250]);
tl = tiledlayout(fh, 3, 3, 'TileSpacing','compact','Padding','compact');
for k = 1:size(PIES,1)
    ax = nexttile(tl,k);  v = PIES{k,2};
    if strcmp(PIES{k,3},'N')
        lbl = arrayfun(@(g) sprintf('%s (N=%d)', S(g).grp, v(g)), 1:nG, 'uni',0);
    else
        lbl = arrayfun(@(g) sprintf('%s (n=%d)', S(g).grp, v(g)), 1:nG, 'uni',0);
    end
    pie_lbl(ax, v, gCol, lbl);   %% counts only -- NO percentages anywhere
    title(ax, sprintf('%s  (n=%d)', PIES{k,1}, sum(v)), 'FontWeight','bold');
end

ax = nexttile(tl,9); axis(ax,'off'); hold(ax,'on');
for g = 1:nG
    patch(ax, [0 .10 .10 0], 0.92-0.115*g + [0 0 .07 .07], gCol(g,:), 'EdgeColor','none');
    text(ax, 0.14, 0.955-0.115*g, ...
         sprintf('%s   %d animals, %d cells, %d active', S(g).grp, S(g).animals, S(g).recorded, S(g).active), ...
         'FontSize',9, 'VerticalAlignment','middle');
end
xlim(ax,[0 1]); ylim(ax,[0 1]);
text(ax, 0, 0.24, {'every pie is split BY GENOTYPE; the colours are the same throughout', ...
                   'animal = one session date', ...
                   'cell = cell\_link identity, curation-tossed masks excluded', ...
                   'active = pooled >5 events;  logZ = occupancy-weighted Rayleigh', ...
                   'IO is a SITE inside ChAT/Vglut2 sessions, so its animals count there too'}, ...
     'FontSize',8.5, 'Color',[0.3 0.3 0.3], 'VerticalAlignment','top');
title(ax,'legend');

title(tl, 'Ventral surface archive by genotype: animals, cells, and breath phase locking', ...
      'FontWeight','bold','FontSize',13);

%% ---- save ----
if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    base = fullfile(outDir,'ventral_pie_summary');
    exportgraphics(fh,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    writetable(table({S.grp}', [S.animals]', [S.recorded]', [S.active]', ...
                     B(:,1), B(:,2), B(:,3), B(:,4), ...
        'VariableNames',{'group','animals','recorded_cells','active_cells', ...
                         'logZ_ge3','logZ_2to3','logZ_1to2','logZ_lt1'}), [base '.csv']);
    fprintf('saved %s.{png,pdf,csv}\n', base);
end

%% ---- local ----
function pie_pct(ax, v, col, minPct)
v = v(:).'; keep = v > 0;
if ~any(keep), axis(ax,'off'); text(ax,0.5,0.5,'none','HorizontalAlignment','center'); return; end
p  = pie(ax, v(keep));                      % no labels -> MATLAB writes percentages
ip = arrayfun(@(x) isgraphics(x,'patch'), p);
hp = p(ip);  cc = col(keep,:);
for k = 1:numel(hp), set(hp(k),'FaceColor',cc(k,:),'EdgeColor','w','LineWidth',0.8); end
ht = p(~ip);
for k = 1:numel(ht)
    set(ht(k),'FontSize',9,'FontWeight','bold');
    if str2double(erase(ht(k).String,'%')) < minPct, ht(k).String = ''; end
end
end

function pie_lbl(ax, v, col, lbls)
v = v(:).'; keep = v > 0;
if ~any(keep), axis(ax,'off'); return; end
p  = pie(ax, v(keep), lbls(keep));
ip = arrayfun(@(x) isgraphics(x,'patch'), p);
hp = p(ip);  cc = col(keep,:);
for k = 1:numel(hp), set(hp(k),'FaceColor',cc(k,:),'EdgeColor','w','LineWidth',0.8); end
set(p(~ip),'FontSize',9);
end
