function vessel_hist_raw_260729(vesselDir, roiTag)
% VESSEL_HIST_RAW_260729  Raw intensity distributions per z and per laser power.
%
% Debug companion to vessel_mode_tail_vs_z_260728.m.  That script reduces each
% slice to a single "mode" number and the result looks wrong, so this one throws
% the reduction away and shows the thing the reduction was computed from: the
% full pixel-value distribution of every slice, on ONE fixed intensity axis
% shared by every slice and every power round.
%
% WHY A FIXED AXIS IS THE WHOLE POINT
%   vessel_mode_tail_vs_z_260728 builds its histogram between prctile(I,0.5)
%   and prctile(I,99.5) with a fixed 512 bins.  That makes the bin WIDTH a
%   function of the slice: a dim deep slice gets sub-ADU bins, a bright
%   saturating surface slice gets ~60 ADU bins, and if >99.5% of the pixels are
%   bright the histogram window never even contains the background peak -- the
%   "first peak from the left" is then a peak of the signal, not of the noise
%   floor.  Nothing about that is visible in a mode-vs-z line plot.  Here every
%   histogram uses the same edges, so slices and rounds are literally comparable
%   and the failure (if that is the failure) is visible as a ridge that jumps.
%
% WHAT IS PLOTTED, per power round
%   1. log10 count heatmap, z (rows) x intensity (columns), FULL ADU range
%   2. the same heatmap zoomed onto the low-intensity end where the background
%      / noise-floor peak lives, with several candidate "mode" estimators drawn
%      on top so they can be compared against the actual ridge
%   3. the same, but pooled over the RAW frames instead of the frame average
%      (averaging framesPerSlice frames narrows the noise peak by ~sqrt(n), so
%      the two are NOT the same distribution and the mode differs)
%   4. individual histogram curves at a handful of depths, log y, with every
%      candidate estimator marked
%   and finally one summary figure with all rounds' estimators on shared axes.
%
% CANDIDATE ESTIMATORS (all on the fixed-edge histogram unless noted)
%   asScript  - exact reproduction of vessel_mode_tail_vs_z_260728's mode,
%               percentile-windowed 512 bins, leftmost findpeaks peak
%   leftPeak  - leftmost findpeaks peak of the FIXED-edge histogram
%   argmax    - the plain, unconditional most-common bin
%   p1/p5/p50 - 1st, 5th, 50th percentile of the slice
%
% CACHING
%   Reading is ~7.6 GB, so all histograms are cached to
%   vessel_hist_raw_<roiTag>.mat.  If that file exists the script only replots.
%   To recompute: delete the .mat.
%
% USAGE
%   vessel_hist_raw_260729
%   vessel_hist_raw_260729('D:\260728_vglut2_soma-g8s\vessel', 'roi1')

%% --------------------------- USER PARAMETERS -------------------------------
VESSEL_DIR = 'D:\260728_vglut2_soma-g8s\vessel';
ROI_TAG    = 'roi1';

Z_TOP      = 10;        % keep slices with SI zs <= this
BIN_W      = 8;         % fixed histogram bin width, ADU
ADU_MIN    = -1024;     % fixed histogram lower edge (ScanImage int16 can be < 0)
LOW_MAX    = 4000;      % upper limit of the "zoomed" low-intensity panels, ADU
SUBSAMPLE  = 2;         % pixel stride, matches vessel_mode_tail_vs_z_260728
NBINS_SCR  = 512;       % script-replica histogram bins
SMOOTH_W   = 9;
MINPROM_FR = 0.01;
NSHOW      = 6;         % how many individual-depth histogram curves per round
SAT_LEVEL  = [];        % [] = auto from pixel class
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(vesselDir), VESSEL_DIR = vesselDir; end
if nargin >= 2 && ~isempty(roiTag),    ROI_TAG    = roiTag;    end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

cacheFile = fullfile(VESSEL_DIR, sprintf('vessel_hist_raw_%s.mat', ROI_TAG));

if exist(cacheFile, 'file')
    fprintf('loading cached histograms: %s\n', cacheFile);
    fprintf('  (delete that file to force a re-read of the TIFFs)\n');
    S = load(cacheFile);
    D = S.D;  edges = S.edges;  centers = S.centers;
else
    [D, edges, centers] = compute_hists(VESSEL_DIR, ROI_TAG, Z_TOP, BIN_W, ...
        ADU_MIN, SUBSAMPLE, NBINS_SCR, SMOOTH_W, MINPROM_FR, SAT_LEVEL);
    save(cacheFile, 'D', 'edges', 'centers', 'BIN_W', 'ADU_MIN', 'Z_TOP', ...
        'VESSEL_DIR', 'ROI_TAG', '-v7.3');
    fprintf('cached -> %s\n', cacheFile);
end

%% --- per-round figures -------------------------------------------------------
for i = 1:numel(D)
    plot_round(D(i), centers, LOW_MAX, NSHOW, VESSEL_DIR, ROI_TAG, SAVE_FIG);
end

%% --- summary across rounds ---------------------------------------------------
plot_summary(D, VESSEL_DIR, ROI_TAG, SAVE_FIG);

%% --- text report -------------------------------------------------------------
fprintf(['\nSANITY: where does the script-replica mode disagree with the plain\n' ...
         'argmax of the fixed-edge histogram?  Large disagreement = the leftmost-\n' ...
         'peak rule picked something that is not the dominant background peak.\n\n']);
fprintf('  %-10s %6s %8s %8s %9s %9s %8s\n', ...
    'round', 'z', 'asScript', 'argmax', 'ratio', 'binW_scr', 'sat%');
fprintf('  %s\n', repmat('-', 1, 70));
for i = 1:numel(D)
    d = D(i);
    bad = abs(log2(max(d.asScript, 1) ./ max(d.argmax, 1))) > 1;   % >2x apart
    idx = find(bad);
    if isempty(idx)
        fprintf('  %-10s  (all slices agree within 2x)\n', d.label);
    else
        show = idx(round(linspace(1, numel(idx), min(6, numel(idx)))));
        for k = show(:).'
            fprintf('  %-10s %6.1f %8.1f %8.1f %9.2f %9.2f %7.2f\n', ...
                d.label, d.z(k), d.asScript(k), d.argmax(k), ...
                d.asScript(k) / max(d.argmax(k), 1), d.binW_scr(k), ...
                100 * d.satFrac(k));
        end
        fprintf('  %-10s  ^ %d of %d slices disagree by >2x\n', ...
            d.label, numel(idx), numel(d.z));
    end
end

end

% ============================== compute =====================================

function [D, edges, centers] = compute_hists(VESSEL_DIR, ROI_TAG, Z_TOP, ...
    BIN_W, ADU_MIN, SUBSAMPLE, NBINS_SCR, SMOOTH_W, MINPROM_FR, SAT_LEVEL)

files = dir(fullfile(VESSEL_DIR, [ROI_TAG '*.tif']));
assert(~isempty(files), 'No %s*.tif under %s', ROI_TAG, VESSEL_DIR);

D = struct('name', {}, 'label', {}, 'power', {}, 'z', {});
meta_all = cell(numel(files), 1);
pw0 = nan(numel(files), 1);
for i = 1:numel(files)
    info = imfinfo(fullfile(files(i).folder, files(i).name));
    meta_all{i} = si_meta(info);
    pw0(i) = si_num(meta_all{i}, "SI.hBeams.stackStartPower");
    meta_all{i} = {meta_all{i}, info};
end
[~, ord] = sort(pw0);
files = files(ord);  meta_all = meta_all(ord);

edges   = ADU_MIN : BIN_W : 32768;
centers = (edges(1:end-1) + edges(2:end)) / 2;

for i = 1:numel(files)
    meta = meta_all{i}{1};  info = meta_all{i}{2};
    tp   = fullfile(files(i).folder, files(i).name);

    zs   = si_vec(meta, "SI.hStackManager.zs");
    nSl  = si_num(meta, "SI.hStackManager.numSlices");
    fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
    P0   = si_num(meta, "SI.hBeams.stackStartPower");
    nCh  = max(numel(si_vec(meta, "SI.hChannels.channelSave")), 1);
    assert(nCh == 1, '%s: %d channels saved; this script assumes 1', ...
        files(i).name, nCh);
    assert(numel(zs) == nSl, '%s: zs/numSlices mismatch', files(i).name);
    assert(numel(info) == nSl * fpsl, '%s: page count mismatch', files(i).name);

    H = info(1).Height;  W = info(1).Width;
    idxZ = find(zs <= Z_TOP);
    nK   = numel(idxZ);
    nB   = numel(centers);

    fprintf('\n[%d/%d] %s  %g%% lp  %d slices (of %d) x %d frames  %dx%d\n', ...
        i, numel(files), files(i).name, P0, nK, nSl, fpsl, H, W);

    Havg = zeros(nK, nB, 'single');
    Hraw = zeros(nK, nB, 'single');
    asScript = nan(nK, 1);  binW_scr = nan(nK, 1);
    argmaxV  = nan(nK, 1);  leftPeak = nan(nK, 1);
    p1 = nan(nK, 1);  p5 = nan(nK, 1);  p50 = nan(nK, 1);
    satFrac  = nan(nK, 1);
    meanV    = nan(nK, 1);

    t  = Tiff(tp, 'r');
    cl = onCleanup(@() t.close());
    if isempty(SAT_LEVEL)
        setDirectory(t, 1);
        satLev = double(intmax(class(t.read())));
    else
        satLev = SAT_LEVEL;
    end

    tStack = tic;
    for k = 1:nK
        acc = zeros(H, W);  nSat = 0;
        for fr = 1:fpsl
            setDirectory(t, (idxZ(k) - 1) * fpsl + fr);
            F = double(t.read());
            nSat = nSat + nnz(F >= satLev);
            acc  = acc + F;
            Hraw(k, :) = Hraw(k, :) + single(histcounts(F(1:SUBSAMPLE:end, 1:SUBSAMPLE:end), edges));
        end
        I = acc / fpsl;
        satFrac(k) = nSat / (fpsl * H * W);

        Is = I(1:SUBSAMPLE:end, 1:SUBSAMPLE:end);
        v  = Is(:);
        Havg(k, :) = single(histcounts(v, edges));

        % --- fixed-edge estimators ---
        [~, im]     = max(Havg(k, :));
        argmaxV(k)  = centers(im);
        cs          = smoothdata(double(Havg(k, :)), 'gaussian', SMOOTH_W);
        [pk, lc]    = findpeaks(cs, centers, 'MinPeakProminence', max(cs) * MINPROM_FR);
        if isempty(pk), leftPeak(k) = argmaxV(k); else, leftPeak(k) = min(lc); end
        q = prctile(v, [1 5 50]);
        p1(k) = q(1);  p5(k) = q(2);  p50(k) = q(3);
        meanV(k) = mean(v);

        % --- exact replica of vessel_mode_tail_vs_z_260728's mode ---
        lo = prctile(v, 0.5);  hi = prctile(v, 99.5);
        if hi <= lo, hi = lo + 1; end
        binW_scr(k) = (hi - lo) / NBINS_SCR;
        e2 = linspace(lo, hi, NBINS_SCR + 1);
        c2 = (e2(1:end-1) + e2(2:end)) / 2;
        s2 = smoothdata(histcounts(v, e2), 'gaussian', SMOOTH_W);
        [pk2, lc2] = findpeaks(s2, c2, 'MinPeakProminence', max(s2) * MINPROM_FR);
        if isempty(pk2)
            [~, i2] = max(s2);  asScript(k) = c2(i2);
        else
            asScript(k) = min(lc2);
        end

        if mod(k, max(1, floor(nK / 6))) == 0 || k == nK
            fprintf('    z %6.1f (%3d/%3d)  asScript %8.1f  argmax %8.1f  p50 %8.1f  binW_scr %6.2f  sat %.3f%%\n', ...
                zs(idxZ(k)), k, nK, asScript(k), argmaxV(k), p50(k), ...
                binW_scr(k), 100 * satFrac(k));
        end
    end
    clear cl
    fprintf('    %.1f s\n', toc(tStack));

    n = numel(D) + 1;
    D(n).name  = files(i).name;
    D(n).label = sprintf('%g%% lp', P0);
    D(n).power = P0;
    D(n).z     = zs(idxZ);
    D(n).Havg  = Havg;
    D(n).Hraw  = Hraw;
    D(n).asScript = asScript;
    D(n).leftPeak = leftPeak;
    D(n).argmax   = argmaxV;
    D(n).p1  = p1;   D(n).p5 = p5;   D(n).p50 = p50;
    D(n).mean = meanV;
    D(n).satFrac  = satFrac;
    D(n).binW_scr = binW_scr;
    D(n).satLevel = satLev;
    D(n).fpsl = fpsl;
end
end

% ============================== plotting ====================================

function plot_round(d, centers, LOW_MAX, NSHOW, VESSEL_DIR, ROI_TAG, SAVE_FIG)

f = figure('Color', 'w', 'Position', [40 40 1500 950], ...
    'Name', sprintf('raw hist %s', d.label));
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

zEdge = [d.z(1), d.z(end)];
lowIdx = centers <= LOW_MAX;

% --- (1,1) full-range heatmap, frame-averaged image -------------------------
ax = nexttile;
imagesc(ax, centers, zEdge, log10(double(d.Havg) + 1));
set(ax, 'YDir', 'normal');  colormap(ax, turbo);
cb = colorbar(ax);  cb.Label.String = 'log_{10}(count + 1)';
xlabel(ax, 'pixel value (ADU), frame-averaged');
ylabel(ax, 'SI z (\mum)');
title(ax, 'FULL range -- every slice, same bins', 'FontWeight', 'normal');
xline(ax, d.satLevel, 'w--', 'ADC ceiling', 'LabelVerticalAlignment', 'bottom');

% --- (1,2) zoomed heatmap + estimators --------------------------------------
ax = nexttile;
imagesc(ax, centers(lowIdx), zEdge, log10(double(d.Havg(:, lowIdx)) + 1));
set(ax, 'YDir', 'normal');  colormap(ax, turbo);  hold(ax, 'on');
cb = colorbar(ax);  cb.Label.String = 'log_{10}(count + 1)';
plot(ax, d.argmax,   d.z, 'w-',  'LineWidth', 2.0, 'DisplayName', 'argmax bin');
plot(ax, d.asScript, d.z, 'r-',  'LineWidth', 1.6, 'DisplayName', 'asScript mode');
plot(ax, d.leftPeak, d.z, 'm--', 'LineWidth', 1.2, 'DisplayName', 'leftmost peak (fixed bins)');
plot(ax, d.p50,      d.z, 'k-',  'LineWidth', 1.0, 'DisplayName', 'median');
plot(ax, d.p5,       d.z, 'k:',  'LineWidth', 1.0, 'DisplayName', '5th pct');
legend(ax, 'Location', 'southeast', 'TextColor', 'w', 'Color', 'k', ...
    'EdgeColor', 'none');
xlabel(ax, 'pixel value (ADU), frame-averaged');
ylabel(ax, 'SI z (\mum)');
title(ax, sprintf('ZOOM 0..%g ADU -- does the red line follow the ridge?', LOW_MAX), ...
    'FontWeight', 'normal');

% --- (2,1) zoomed heatmap, RAW frames ---------------------------------------
ax = nexttile;
imagesc(ax, centers(lowIdx), zEdge, log10(double(d.Hraw(:, lowIdx)) + 1));
set(ax, 'YDir', 'normal');  colormap(ax, turbo);
cb = colorbar(ax);  cb.Label.String = 'log_{10}(count + 1)';
xlabel(ax, 'pixel value (ADU), RAW single frames pooled');
ylabel(ax, 'SI z (\mum)');
title(ax, sprintf('same zoom, un-averaged (%d frames/slice) -- wider noise peak', ...
    d.fpsl), 'FontWeight', 'normal');

% --- (2,2) individual histograms at selected depths --------------------------
ax = nexttile; hold(ax, 'on');
pick = round(linspace(1, numel(d.z), NSHOW));
cc = turbo(NSHOW);
for j = 1:NSHOW
    k = pick(j);
    plot(ax, centers(lowIdx), double(d.Havg(k, lowIdx)) + 0.1, '-', ...
        'Color', cc(j,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('z = %.0f  (sat %.2f%%)', d.z(k), 100*d.satFrac(k)));
    plot(ax, d.asScript(k), 1, 'v', 'Color', cc(j,:), 'MarkerFaceColor', cc(j,:), ...
        'HandleVisibility', 'off');
end
set(ax, 'YScale', 'log');  grid(ax, 'on');  box(ax, 'on');
xlim(ax, [centers(1) LOW_MAX]);
xlabel(ax, 'pixel value (ADU)');  ylabel(ax, 'count');
title(ax, 'individual slices (\nabla = asScript mode)', 'FontWeight', 'normal');
legend(ax, 'Location', 'northeast', 'Box', 'off');

title(tl, sprintf('%s  --  %s  --  raw pixel-value distributions', d.name, d.label), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

if SAVE_FIG
    base = fullfile(VESSEL_DIR, sprintf('vessel_hist_raw_%s_%s', ROI_TAG, ...
        matlab.lang.makeValidName(d.label)));
    exportgraphics(f, [base '.png'], 'Resolution', 180, 'BackgroundColor', 'white');
    fprintf('saved %s.png\n', base);
end
end

function plot_summary(D, VESSEL_DIR, ROI_TAG, SAVE_FIG)
cols = lines(max(numel(D), 4));
f = figure('Color', 'w', 'Position', [60 60 1400 800], 'Name', 'estimators vs z');
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

fields = {'asScript', 'argmax', 'p50', 'satFrac'};
names  = {'asScript mode (current metric)', 'argmax bin (plain mode)', ...
          'median', 'saturated pixel fraction'};
for q = 1:4
    ax = nexttile; hold(ax, 'on');
    for i = 1:numel(D)
        y = D(i).(fields{q});
        if strcmp(fields{q}, 'satFrac'), y = 100 * y; end
        plot(ax, D(i).z, y, '-', 'Color', cols(i,:), 'LineWidth', 1.4, ...
            'DisplayName', D(i).label);
    end
    set(ax, 'XDir', 'reverse');  grid(ax, 'on');  box(ax, 'on');
    if q < 4, set(ax, 'YScale', 'log'); end
    xlabel(ax, 'SI z (\mum)');
    ylabel(ax, tern(q == 4, '% pixels at ADC ceiling', 'ADU'));
    title(ax, names{q}, 'FontWeight', 'normal');
    if q == 1, legend(ax, 'Location', 'best', 'Box', 'off'); end
end
title(tl, sprintf('%s  --  candidate background estimators vs depth', ROI_TAG), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

if SAVE_FIG
    base = fullfile(VESSEL_DIR, sprintf('vessel_hist_raw_%s_summary', ROI_TAG));
    exportgraphics(f, [base '.png'], 'Resolution', 180, 'BackgroundColor', 'white');
    fprintf('saved %s.png\n', base);
end
end

% ========================= ScanImage metadata helpers ========================

function meta = si_meta(info)
meta = "";
if isfield(info(1), 'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1), 'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
end
assert(strlength(meta) > 0, 'No ScanImage metadata in Software/ImageDescription');
meta = replace(meta, char(8629), newline);
end

function v = si_num(meta, key)
ln = si_line(meta, key);
if strlength(ln) == 0, v = NaN; return; end
tok = regexp(ln, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = si_vec(meta, key)
vec = [];
ln = si_line(meta, key);
if strlength(ln) == 0, return; end
m = regexp(ln, '\[(.*)\]', 'tokens', 'once');
if isempty(m), src = ln; else, src = m{1}; end
nums = regexp(src, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end

function ln = si_line(meta, key)
ln = "";
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
ln = extractBefore(s, newline);
if isempty(ln) || ismissing(ln), ln = s; end
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
