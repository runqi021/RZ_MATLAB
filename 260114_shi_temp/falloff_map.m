
optical_falloff_hist_tail_mode_preavg_nometa()
%%
function optical_falloff_hist_tail_mode_preavg_nometa()
% optical_falloff_hist_tail_mode_preavg_nometa
% ------------------------------------------------------------
% INPUT TIFF: already Z-averaged stack (each page = one Z slice avg)
% OUTPUTS:
%   (1) Overlaid UNFILLED histograms for selected Z slices
%   (2) 2x2 panels:
%        TL: Tail (raw) + power schedule (yyaxis right)
%        TR: Tail / P^2 (power-corrected) + exp fit
%        BL: Mode (first peak from left) + power schedule
%        BR: Mode / P^2 (power-corrected) + exp fit
%
% No TIFF metadata used. You manually specify:
%   - z axis: zStart, zEnd, zStep   (for plotting AND power schedule)
%   - power schedule: P0, P1, LC_um
% ------------------------------------------------------------

%% ---------------- USER ----------------
tiffPath   = "E:\live_251104_definitely_not_sst_tdtomato_fitc\raw\groupZ\matlab_stitch\stitched_stack_imagejstyle_nometa.tif";

% Histogram/metrics
subsample   = 1;        % <-- 1 = use ALL pixels per Z (you asked for this)
nBins       = 512;
smoothW     = 9;
minPromFrac = 0.01;     % findpeaks MinPeakProminence = max(counts)*minPromFrac
Nbright     = 1000;     % brightest N pixels for tail metric

% Z axis (KEEP your convention)
zStart = -20;  zEnd = 200;  zStep = 10;

% Power schedule (MANUAL, no metadata)
P0_pct = 7;            % start power (%)
P1_pct = 26;            % end power (%)
LC_um  = 167.7;           % length constant (um) in the exponential power ramp

% Option: force exact end match (so power hits exactly P1 at last slice)
force_end_match = false;

% Plot controls
sliceList = [3 13 23];
skipN     = 2;          % skip first N slices for exp fit (use skipN+1:end)

%% ---------------- BASIC TIFF INFO ----------------
assert(isfile(tiffPath), "File not found: %s", tiffPath);

info   = imfinfo(tiffPath);
nPages = numel(info);
H = info(1).Height;
W = info(1).Width;

% For pre-averaged stack: each page is one Z
nZ = nPages;

% ----- YOUR z labeling convention (do NOT change) -----
z_expected = (zStart:zStep:zEnd);
if numel(z_expected) == nZ
    z_um = z_expected(:);
else
    warning("Expected %d z from [%g:%g:%g], got %d pages. Using linspace.", ...
        numel(z_expected), zStart, zStep, zEnd, nZ);
    z_um = linspace(zStart, zEnd, nZ).';
end

%% ---------------- POWER SCHEDULE (MANUAL) ----------------
% Model: P(z) = P0 * exp((z - z0)/LC)
z0 = z_um(1);
p_pct = P0_pct .* exp((z_um - z0) ./ LC_um);

if force_end_match
    % enforce exact endpoint by log-linear interpolation between endpoints
    p_pct = exp(interp1([z_um(1) z_um(end)], log([P0_pct P1_pct]), z_um, 'linear', 'extrap'));
end

p_frac = p_pct / 100;

fprintf('[input] nZ=%d pages (assumed 1 page per Z)\n', nZ);
fprintf('[input] z range: %.3g -> %.3g (dz=%.3g)\n', z_um(1), z_um(end), mean(diff(z_um)));
fprintf('[input] power: P0=%.2f%%, P(end)=%.2f%% (target P1=%.2f%%), LC=%.4g um\n', ...
    P0_pct, p_pct(end), P1_pct, LC_um);

%% ---------------- PER-Z METRICS ----------------
modeVal = zeros(nZ,1);
minVal  = zeros(nZ,1);
tailVal = zeros(nZ,1);

t = Tiff(tiffPath,'r');
cleanup = onCleanup(@() t.close()); %#ok<NASGU>

for iz = 1:nZ
    I = read_page_as_double(t, iz);

    % collect pixels for histogram (ALL pixels if subsample=1)
    Is   = I(1:subsample:end, 1:subsample:end);
    vals = Is(:);
    vals = vals(isfinite(vals));

    % Min = low-edge proxy
    mn = prctile(vals, 0.5);
    minVal(iz) = mn;

    % Histogram in robust range
    lo = prctile(vals, 0.5);
    hi = prctile(vals, 99.5);
    if hi <= lo, hi = lo + 1; end

    edges    = linspace(lo, hi, nBins+1);
    counts   = histcounts(vals, edges);
    centers  = (edges(1:end-1) + edges(2:end))/2;
    counts_s = smoothdata(counts, 'gaussian', smoothW);

    % Mode = FIRST PEAK FROM LEFT
    [pks, locs] = findpeaks(counts_s, centers, ...
        'MinPeakProminence', max(counts_s) * minPromFrac);

    if isempty(pks)
        [~, imax] = max(counts_s);
        modeVal(iz) = centers(imax);
    else
        modeVal(iz) = min(locs);
    end

    % Tail = mean of brightest N pixels of (I - Min)
    vv = I(:) - mn;
    vv = vv(isfinite(vv));
    Nb = min(Nbright, numel(vv));
    tailVal(iz) = mean(maxk(vv, Nb));

    if mod(iz, max(1,floor(nZ/10)))==0 || iz==nZ
        fprintf("z %d/%d: Min=%.2f  Mode(firstPeak)=%.2f  Tail=%.2f  P=%.2f%%\n", ...
            iz, nZ, minVal(iz), modeVal(iz), tailVal(iz), p_pct(iz));
    end
end

%% ---------------- POWER CORRECTION (2p ~ P^2) ----------------
mode_corr = modeVal ./ max(p_frac.^2, eps);
tail_corr = tailVal ./ max(p_frac.^2, eps);

%% ---------------- FIT EXP ON (skipN+1:end) ----------------
idxFit = (skipN+1):nZ;
z_plot = z_um(idxFit);

[L_tail_c, A_tail_c] = local_fit_exp(z_plot, tail_corr(idxFit));
[L_mode_c, A_mode_c] = local_fit_exp(z_plot, mode_corr(idxFit));

z0fit = z_plot(1);
tail_fit = A_tail_c * exp(-(z_plot - z0fit)/L_tail_c);
mode_fit = A_mode_c * exp(-(z_plot - z0fit)/L_mode_c);

%% ============================================================
% (1) FIGURE: OVERLAID HISTOGRAMS (UNFILLED)
% ============================================================
loAll = inf; hiAll = -inf;

tH = Tiff(tiffPath,'r');
cH = onCleanup(@() tH.close()); %#ok<NASGU>

for ii = 1:numel(sliceList)
    iz = sliceList(ii);
    if iz < 1 || iz > nZ
        warning('slice %d out of range (nZ=%d), skipping.', iz, nZ);
        continue;
    end
    I = read_page_as_double(tH, iz);
    Is = I(1:subsample:end, 1:subsample:end);
    vals = Is(:); vals = vals(isfinite(vals));
    loAll = min(loAll, prctile(vals, 0.5));
    hiAll = max(hiAll, prctile(vals, 99.5));
end

if ~isfinite(loAll) || ~isfinite(hiAll) || hiAll <= loAll
    error('Could not determine histogram range for selected slices.');
end

edges_common   = linspace(loAll, hiAll, nBins+1);
centers_common = (edges_common(1:end-1) + edges_common(2:end))/2;

figure('Color','w','Position',[120 120 900 500]); hold on;
tH = Tiff(tiffPath,'r');
cH = onCleanup(@() tH.close()); %#ok<NASGU>

leg = {};
for ii = 1:numel(sliceList)
    iz = sliceList(ii);
    if iz < 1 || iz > nZ, continue; end
    I = read_page_as_double(tH, iz);
    Is = I(1:subsample:end, 1:subsample:end);
    vals = Is(:); vals = vals(isfinite(vals));
    counts = histcounts(vals, edges_common);
    plot(centers_common, counts, 'LineWidth', 1);   % unfilled overlay
    leg{end+1} = sprintf('%.0f \\mum', z_um(iz)); %#ok<AGROW>
end

set(gca,'YScale','log');
xlabel('Intensity (ADU)');
xlim([0, 32700]);   % keep your convention
ylabel('Counts (log)');
title('F Histogram (pre-avg Z stack, no metadata)');
legend(leg, 'Location','northeast');
axis square;

%% ============================================================
% (2) FIGURE: 2x2 PANELS (KEEP YOUR STYLE)
% ============================================================
figure('Color','w','Position',[80 80 1200 800]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% --- (TL) Tail raw + power curve (yyaxis) ---
ax = nexttile; hold(ax,'on');
yyaxis(ax,'left');
plot(ax, z_plot, tailVal(idxFit), 'k-', 'LineWidth', 1);
ylabel(ax, 'Average Intensity in Tail');
ylim([0 32767]);

yyaxis(ax,'right'); ylim([0 100]);
plot(ax, z_plot, p_pct(idxFit), 'r--', 'LineWidth', 1);
ylabel(ax, 'Laser power (%)');

xlabel(ax, 'Depth (\mum)');
axis square;
xlim([0 700]);
title(ax, 'Tail');
grid(ax,'on');

% --- (TR) Tail corrected + exp fit ---
ax = nexttile; hold(ax,'on');
semilogy(ax, z_plot, tail_corr(idxFit), 'k-', 'LineWidth', 1);
semilogy(ax, z_plot, tail_fit, 'r--', 'LineWidth', 1);
xlabel(ax, 'Depth (\mum)');
ylabel(ax, 'Tail / P^2');
xlim([0 700]);
grid(ax,'on');
set(ax, 'YScale', 'log');
axis square;

% --- (BL) Mode raw + power curve (yyaxis) ---
ax = nexttile; hold(ax,'on');
yyaxis(ax,'left');
plot(ax, z_plot, modeVal(idxFit), 'k-', 'LineWidth', 1);
ylim([0 2200]);
ylabel(ax, 'Intensity at Mode');

yyaxis(ax,'right'); ylim([0 100]);
plot(ax, z_plot, p_pct(idxFit), 'r--', 'LineWidth', 1);
ylabel(ax, 'Laser power (%)');

xlabel(ax, 'Depth (\mum)');
xlim([0 700]);
axis square;
title(ax, 'Mode');
grid(ax,'on');

% --- (BR) Mode corrected + exp fit ---
ax = nexttile; hold(ax,'on');
semilogy(ax, z_plot, mode_corr(idxFit), 'k-', 'LineWidth', 1);
semilogy(ax, z_plot, mode_fit, 'r--', 'LineWidth', 1);
xlabel(ax, 'Depth (\mum)');
ylabel(ax, 'Mode / P^2');
grid(ax,'on');
set(ax, 'YScale', 'log');
xlim([0 700]);
axis square;

sgtitle(tl, sprintf('Optical falloff (pre-avg Z, skip first %d). Power: %.1f%%→%.1f%% (LC=%.1f \\mum)', ...
    skipN, p_pct(1), p_pct(end), LC_um));

end

%% ================= HELPERS =================
function I = read_page_as_double(tiffObj, pageIdx)
setDirectory(tiffObj, pageIdx);
I = double(tiffObj.read());
end

function [L, A] = local_fit_exp(z, y)
% Fit y ~ A * exp(-(z - z0)/L) via polyfit on log(y)
m = isfinite(y) & (y > 0) & isfinite(z);
z = z(m); y = y(m);
if numel(z) < 5
    L = NaN; A = NaN; return;
end
z0 = z(1);
p = polyfit(z - z0, log(y), 1);    % log(y)=log(A) - (z-z0)/L
A = exp(p(2));
L = -1/p(1);
end
