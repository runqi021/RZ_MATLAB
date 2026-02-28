plot_random5_onepanel_per_z

%%
function plot_random5_onepanel_per_z()
% plot_random5_onepanel_per_z
% ------------------------------------------------------------
% Single-panel summary:
% - finds z-folders (z0, z20, ...) under outRoot dynamically
% - loads each slice's *_dFF.mat (variable name: dFF)
% - randomly selects 5 ROI columns per slice
% - plots ALL slices in ONE axes, stacked vertically
% - draws ONE yline per slice (center of that slice block)
% - y-ticks at those lines, labeled by slice z (um)
% - one global dFF scale bar
%
% Outputs in outRoot:
%   random5_onepanel_perZ.pdf
%   random5_onepanel_perZ.png
%   random5_onepanel_perZ_picks.mat
% ------------------------------------------------------------

%% ================= USER =================
outRoot   = "C:\Users\Admin\Desktop\260210_chat_soma_g8s\roi1_-1050_400_z-25-125_IO_00001_ZSTACK_PIPELINE";
nPick     = 10;
seed      = 6;      % [] => new random each run; number => reproducible

fps       = 30;

% same "stackDFF vibe"
gain      = 2.5;
minOffset = 1.0;
offsetMult= 12;

% single global dFF scale bar
barAmp    = 0.20;   % 20% dFF
gapMult   = 1.2;    % vertical gap between slices in units of offset
%% =======================================

assert(isfolder(outRoot), "Folder not found: %s", outRoot);
if ~isempty(seed), rng(seed); end

%% ---- find z folders dynamically ----
dd = dir(fullfile(outRoot, "z*"));
dd = dd([dd.isdir]);

Z = [];
Zdir = strings(0,1);
for i = 1:numel(dd)
    name = string(dd(i).name);
    tok = regexp(name, "^z(-?\d+)$", "tokens", "once");
    if isempty(tok), continue; end
    Z(end+1,1) = str2double(tok{1}); %#ok<AGROW>
    Zdir(end+1,1) = string(fullfile(dd(i).folder, dd(i).name)); %#ok<AGROW>
end
assert(~isempty(Z), "No z-folders found under: %s", outRoot);

% z0 should be TOP, z_end BOTTOM
[Z, ord] = sort(Z, "ascend");
Zdir = Zdir(ord);

%% ---- locate dFF mats per z ----
DFFMAT = strings(numel(Z),1);
ok = false(numel(Z),1);

for i = 1:numel(Z)
    subdir = Zdir(i);

    d1 = dir(fullfile(subdir, "*cpSAM_output_dFF.mat"));
    if ~isempty(d1)
        [~,ix] = max([d1.datenum]);
        DFFMAT(i) = string(fullfile(d1(ix).folder, d1(ix).name));
        ok(i) = true;
        continue
    end

    d2 = dir(fullfile(subdir, "*_dFF.mat"));
    if ~isempty(d2)
        [~,ix] = max([d2.datenum]);
        DFFMAT(i) = string(fullfile(d2(ix).folder, d2(ix).name));
        ok(i) = true;
    end
end

Z      = Z(ok);
Zdir   = Zdir(ok);
DFFMAT = DFFMAT(ok);
assert(~isempty(Z), "No *_dFF.mat found inside z folders.");

nZ = numel(Z);

%% ---- pick ROIs + load 5 traces per slice (baseline removed + gain applied) ----
picked = cell(nZ,1);
Dcell  = cell(nZ,1);
Tlen   = zeros(nZ,1);

for i = 1:nZ
    pth = char(DFFMAT(i));

    vinfo = whos('-file', pth, 'dFF');
    assert(~isempty(vinfo), "No variable 'dFF' in: %s", pth);
    sz = vinfo.size;  % [T N]
    T = sz(1); N = sz(2);

    kPick = min(nPick, N);
    sel = randperm(N, kPick);
    picked{i} = sel;

    % partial load (your dFF mats are -v7.3, so this should work)
    try
        mf = matfile(pth);
        d = mf.dFF(:, sel);
    catch
        S = load(pth, 'dFF');
        d = S.dFF(:, sel);
    end

    % baseline remove per trace + gain
    for j = 1:size(d,2)
        d(:,j) = d(:,j) - median(d(:,j), "omitnan");
    end
    d = gain * d;

    Dcell{i} = d;
    Tlen(i)  = size(d,1);
end

% global spacing scale (one scale for whole figure)
allVals = vertcat(Dcell{:});
s = mad(allVals(:), 1);
if ~isfinite(s) || s==0
    s = std(allVals(:), "omitnan");
end
if ~isfinite(s) || s==0
    s = 0.05 * gain;
end
offset = max(minOffset, offsetMult * s);
gap    = gapMult * offset;

% time axis (use max length)
Tmax = max(Tlen);
tmax = (Tmax-1) / fps;

% each slice block height
nPickEff = cellfun(@(d)size(d,2), Dcell);   % some slices may have <5 ROI
blockSpan = (max(nPickEff)-1) * offset;     % reserve vertical space uniformly

%% ---- plot: ONE axes ----
rowPx = 120;                 % tweak if you want taller
figH  = max(800, rowPx*nZ);
figW  = 900;

fig = figure("Color","w","Units","pixels","Position",[30 30 figW figH]);
ax  = axes(fig); hold(ax,"on");

yLinePos = zeros(nZ,1);   % y positions for yline + ytick labels

for i = 1:nZ
    d = Dcell{i};
    Ti = size(d,1);
    t  = (0:Ti-1)/fps;

    yBase = (i-1) * (blockSpan + gap);

    % stack the picked traces within slice block
    for j = 1:size(d,2)
        y0 = yBase + (j-1)*offset;
        plot(ax, t, -d(:,j) + y0, "k", "LineWidth", 0.8);
    end

    % one separator line "in the middle" of this slice block
    yMid = yBase + blockSpan/2;
    yLinePos(i) = yMid;
    %yl = yline(ax, yMid, "-", "LineWidth", 0.5);
    %yl.Color = [0.6 0.6 0.6];
end

% formatting
xlim(ax, [0 tmax]);
ylim(ax, [-0.2*gap, (nZ-1)*(blockSpan+gap) + blockSpan + 0.2*gap]);

set(ax, "YDir","reverse");   % z0 at top, z_end at bottom
xlabel(ax, "Time (s)");
ylabel(ax, "Depth (um)");

% y ticks/labels correspond to the ylines
yticks(ax, yLinePos);
yticklabels(ax, "z" + string(Z));

box(ax,"off");

%% ---- single global dFF scale bar (one) ----
barH = barAmp * gain;

% put the bar in the first gap region (usually empty)
if nZ >= 2
    yBar0 = (blockSpan + gap) - 0.65*gap;   % inside gap after slice 1
else
    yBar0 = 0.15*gap;
end
xBar  = 0.92 * tmax;

plot(ax, [xBar xBar], [yBar0 yBar0+barH], "k", "LineWidth", 2);
text(ax, xBar, yBar0+barH/2, sprintf("  %d%% dFF", round(barAmp*100)), ...
    "Color","k","FontSize",10,"FontName","Arial", ...
    "HorizontalAlignment","left","VerticalAlignment","middle");

%% ---- save ----
outPdf = fullfile(outRoot, "random5_onepanel_perZ.pdf");
outPng = fullfile(outRoot, "random5_onepanel_perZ.png");
outMat = fullfile(outRoot, "random5_onepanel_perZ_picks.mat");

exportgraphics(fig, outPdf, "ContentType","vector");
exportgraphics(fig, outPng, "Resolution", 220);
close(fig);

save(outMat, "Z", "Zdir", "DFFMAT", "picked", "nPick", "seed", "fps", ...
    "gain","minOffset","offsetMult","barAmp","offset","gap","blockSpan");

fprintf("Saved:\n  %s\n  %s\n  %s\n", outPdf, outPng, outMat);
end
