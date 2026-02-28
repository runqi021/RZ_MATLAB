optical_falloff_hist_tail_mode

%%
function optical_falloff_hist_tail_mode()
% optical_falloff_hist_tail_mode
% ------------------------------------------------------------
% Two outputs:
%   (1) Figure: overlaid UNFILLED histograms for selected Z slices
%   (2) Figure: 2x2 panels:
%        TL: Tail (raw) + laser power schedule (yyaxis right)
%        TR: Tail / P^2 (power-corrected) + exp fit
%        BL: Mode (raw; FIRST peak from left) + laser power schedule (yyaxis right)
%        BR: Mode / P^2 (power-corrected) + exp fit
%
% This version:
%   - Keeps YOUR plotting convention (z_um from zStart:zStep:zEnd, xlim=[0 700], stairs+plot, etc.)
%   - Uses ScanImage TIFF metadata ONLY for the power schedule (P0/P1/LC and zs)
% ------------------------------------------------------------

%% ---------------- USER ----------------
tiffPath   = "C:\Users\Admin\Downloads\shi_250115\FN_1200_1600_z50-700_15-55lp_00002.tif";

framesPerZ = 50;          % must divide nPages
subsample  = 2;           % histogram speed; 2–4 is fine
nBins      = 512;
smoothW    = 9;
minPromFrac = 0.01;       % findpeaks MinPeakProminence = max(counts)*minPromFrac

Nbright    = 1000;        % brightest N pixels for tail metric

% Z axis (for labeling only)  <-- KEEP YOUR CONVENTION
zStart = -40; zEnd = 710; zStep = 10;

% (kept for title text only; power used in computation comes from TIFF meta)
pMin = 15; pMax = 55;

sliceList = [15 35 55];
skipN = 4;    % skip 1..5, use 6..end

%% ---------------- META (basic TIFF) ----------------
assert(isfile(tiffPath), "File not found: %s", tiffPath);

info   = imfinfo(tiffPath);
nPages = numel(info);
H = info(1).Height;
W = info(1).Width;

if mod(nPages, framesPerZ) ~= 0
    error("nPages=%d not divisible by framesPerZ=%d", nPages, framesPerZ);
end
nZ = nPages / framesPerZ;

% ----- YOUR z labeling convention (do NOT change) -----
z_expected = (zStart:zStep:zEnd);
if numel(z_expected) == nZ
    z_um = z_expected(:);
else
    warning("Expected %d z from [%g:%g:%g], got %d. Using linspace.", ...
        numel(z_expected), zStart, zStep, zEnd, nZ);
    z_um = linspace(zStart, zEnd, nZ).';
end

%% ---------------- READ ScanImage META (POWER ONLY) ----------------
% Robust metadata extraction (Software -> ImageDescription -> low-level tag)
meta = "";
if isfield(info(1),'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
else
    t0 = Tiff(tiffPath,'r');
    c0 = onCleanup(@() t0.close());
    try
        meta = string(t0.getTag('ImageDescription'));
    catch
        try
            meta = string(t0.getTag('Software'));
        catch
            error("No ScanImage metadata found in Software/ImageDescription tags.");
        end
    end
end
meta = replace(meta, "↵", newline);

getNum = @(key) local_get_num(meta, key);
getVec = @(key) local_get_vec(meta, key);

LC_SI = getNum("SI.hBeams.lengthConstants");
P0_SI = getNum("SI.hBeams.stackStartPower");   % percent
P1_SI = getNum("SI.hBeams.stackEndPower");     % percent

zs_SI = getVec("SI.hStackManager.zs");         % motor Z list (preferred)
z0_SI = getNum("SI.hStackManager.stackZStartPos");
z1_SI = getNum("SI.hStackManager.stackZEndPos");
dz_SI = getNum("SI.hStackManager.stackZStepSize");

% Build motor-Z list in slice order (length must be nZ)
if ~isempty(zs_SI)
    z_motor = zs_SI(:);
else
    if ~all(isfinite([z0_SI z1_SI dz_SI])) || dz_SI==0
        error("Could not read SI.hStackManager.zs or stackZStart/End/Step.");
    end
    z_motor = (z0_SI:dz_SI:z1_SI).';
end

if numel(z_motor) ~= nZ
    warning("ScanImage zs count (%d) != nZ (%d). Using linspace for motor Z.", numel(z_motor), nZ);
    z_motor = linspace(z_motor(1), z_motor(end), nZ).';
end

% Power schedule from ScanImage (aligned by SLICE INDEX; independent of your z_um labeling)
% Model: P(z_motor) = P0 * exp((z_motor - z0)/LC)
p_pct  = P0_SI .* exp((z_motor - z_motor(1)) ./ LC_SI);

% Optional: exact endpoint match (if you want strict P0->P1 regardless of LC rounding)
% p_pct = exp(interp1([z_motor(1), z_motor(end)], log([P0_SI,P1_SI]), z_motor, 'linear','extrap'));

p_frac = p_pct / 100;

% sanity print
depth_meta = z_motor(1) - z_motor;  % depth-like axis from SI (not used for your plot)
if any(diff(depth_meta) <= 0), depth_meta = -depth_meta; end

fprintf('[meta] nZ=%d, framesPerZ=%d\n', nZ, framesPerZ);
fprintf('[meta] SI motorZ: %g -> %g (dz~%g)\n', z_motor(1), z_motor(end), mean(diff(z_motor)));
fprintf('[meta] SI depth : %g -> %g\n', depth_meta(1), depth_meta(end));
fprintf('[meta] SI power : %.2f%% -> %.2f%%, LC=%.4g um (implied LC=%.4g)\n', ...
    p_pct(1), p_pct(end), LC_SI, (z_motor(end)-z_motor(1))/log(P1_SI/P0_SI));

%% ---------------- PER-Z METRICS ----------------
modeVal = zeros(nZ,1);
minVal  = zeros(nZ,1);
tailVal = zeros(nZ,1);

t = Tiff(tiffPath,'r');
cleanup = onCleanup(@() t.close()); %#ok<NASGU>

for iz = 1:nZ
    I = avg_image_for_z(t, iz, framesPerZ, H, W);

    % subsample pixels for histogram
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

[L_tail_c, A_tail_c] = local_fit_exp(z_um(idxFit), tail_corr(idxFit));
[L_mode_c, A_mode_c] = local_fit_exp(z_um(idxFit), mode_corr(idxFit));

% Fit curves (anchored at first fit depth for nicer overlay)
z_plot = z_um(idxFit);
z0 = z_plot(1);

tail_fit = A_tail_c * exp(-(z_plot - z0)/L_tail_c);
mode_fit = A_mode_c * exp(-(z_plot - z0)/L_mode_c);

%% ============================================================
% (1) FIGURE: OVERLAID HISTOGRAMS (UNFILLED), FULL RANGE
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

    I = avg_image_for_z(tH, iz, framesPerZ, H, W);
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

    I = avg_image_for_z(tH, iz, framesPerZ, H, W);
    Is = I(1:subsample:end, 1:subsample:end);
    vals = Is(:); vals = vals(isfinite(vals));

    counts = histcounts(vals, edges_common);
    plot(centers_common, counts, 'LineWidth', 1);   % unfilled overlay
    leg{end+1} = sprintf('%.0f \\mum', z_um(iz)); %#ok<AGROW>
end

set(gca,'YScale','log');
xlabel('Intensity (ADU)');
xlim([0,32700]);
ylabel('Counts (log)');
title('F Histogram');
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
%stairs(z_plot, tailVal(idxFit), 'k-', 'LineWidth', 1);
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
%stairs(z_plot, modeVal(idxFit), 'k-', 'LineWidth', 1);
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

% keep your title formatting, but now show REAL SI power endpoints
sgtitle(tl, sprintf('Optical falloff (skip first %d slices). Power ramp: %.1f%%→%.1f%% exp', ...
    skipN, p_pct(1), p_pct(end)));

end

%% ================= HELPERS =================
function I = avg_image_for_z(tiffObj, iz, framesPerZ, H, W)
acc = zeros(H,W,'double');
for f = 1:framesPerZ
    page = (iz-1)*framesPerZ + f;
    setDirectory(tiffObj, page);
    acc = acc + double(tiffObj.read());
end
I = acc / framesPerZ;
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

function v = local_get_num(meta, key)
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), v = NaN; return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
tok = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = local_get_vec(meta, key)
vec = [];
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
m = regexp(line, '\[(.*)\]', 'tokens', 'once');
if isempty(m), return; end
nums = regexp(m{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end
