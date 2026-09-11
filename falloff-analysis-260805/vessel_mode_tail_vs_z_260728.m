function vessel_mode_tail_vs_z_260728(vesselDir, roiTag)
% VESSEL_MODE_TAIL_VS_Z_260728  Mode + tail vs depth, one colour per laser-power round.
%
% A depth series acquired as several OVERLAPPING stacks, each at a different
% fixed laser power, so no single stack has to span the whole depth range at
% one power.  This plots mode and tail on a single shared z axis with one
% colour per power round, so the rounds can be compared where they overlap.
%
% Metric definitions follow 260114_shi_temp/F_attenuation.m:
%   mode = FIRST histogram peak from the left of the frame-averaged image
%          (the background / noise-floor peak, not the global max)
%   tail = mean of the brightest NBRIGHT pixels of (I - prctile(I, 0.5))
%
% TWO THINGS THIS SCRIPT REFUSES TO TAKE FROM THE FILENAME
%
% 1. Z.  The filename z is DEPTH, positive downward, sign dropped:
%        roi1_z40-160_lp21  ->  SI.hStackManager.zs runs -40 to -160.
%    ScanImage z is motor position, negative downward.  Only `zs` is used
%    here; the filename z is parsed solely to print the mismatch.
%
% 2. POWER.  Read from SI.hBeams, not from the "lp" token.
%    On this session pzAdjust = false and stackStartPower == stackEndPower,
%    so power is CONSTANT within each stack.  SI.hBeams.lengthConstants is
%    still present (-374.237) but is stale and must NOT be used to build a
%    within-stack ramp -- F_attenuation.m applies P0*exp((z-z0)/LC)
%    unconditionally, which here would invent a ~35%% rise by the bottom of a
%    flat-power stack.  This script asserts the power really is flat and
%    falls back to the ramp only when ScanImage says the ramp is on.
%
% Two-photon signal goes as P^2, so the P^2-corrected columns are what make
% rounds at different powers comparable in their shared z range.
%
% USAGE
%   vessel_mode_tail_vs_z_260728                                  % defaults below
%   vessel_mode_tail_vs_z_260728('D:\260728_vglut2_soma-g8s\vessel', 'roi1')
%
% OUTPUT (written into vesselDir)
%   vessel_mode_tail_vs_z_<roiTag>.png / .pdf / .mat

%% --------------------------- USER PARAMETERS -------------------------------
VESSEL_DIR = 'D:\260728_vglut2_soma-g8s\vessel';
ROI_TAG    = 'roi1';       % which acquisition series; 'roi2' is unfinished (1 stack)

Z_TOP      = 10;           % keep slices with zs <= this, in SI motor convention
NBRIGHT    = 1000;         % brightest N pixels for the tail metric
SUBSAMPLE  = 2;            % pixel stride for the histogram (mode) only
NBINS      = 512;
SMOOTH_W   = 9;
MINPROM_FR = 0.01;         % findpeaks prominence, as a fraction of max count
FRAME_STRIDE = 1;          % 1 = average all framesPerSlice frames (as F_attenuation)

% Saturation.  The tail is the mean of the brightest pixels, so as soon as the
% PMT/digitiser clips, the tail stops measuring brightness and starts reporting
% the ADC ceiling -- it goes flat, which looks like a plateau rather than like
% a failure.  On this session the top slices of the 21/35/48%% rounds all clip.
% Saturation is counted on the RAW frames (before averaging, which would hide
% partial clipping) and flagged slices are drawn dotted and dropped from the
% overlap statistics.
SAT_LEVEL    = [];         % [] = auto from the pixel class (int16 -> 32767)
SAT_FRAC_MAX = 1e-3;       % slice is "saturated" above this fraction of clipped px

% ScanImage "%" is a Pockels setpoint, NOT power -- the %-to-mW curve is
% strongly superlinear at low % and rolls over above ~70%.  The plotted P^2
% panels use the raw % (no external assumption), but the overlap check also
% reports what the ratios become under the measured calibration below, copied
% from laser_power_vs_z_260530.m.  If that calibration is current for this rig
% it is the better denominator; if it is stale the two columns bracket the
% answer.  Verify before treating either as absolute.
CAL_PCT = [0.10 1 3 5 8 10 12 15 18 20 22 25 27 30 32 35 37 40 ...
           45 50 55 60 65 70 75 80 85 90];
CAL_MW  = [2.3  2.3 3.4 5.7 11.2 16.5 23 34.4 48.8 59.6 70.5 89.5 103 ...
           123 138 161 177 200 240 276 312 341 365 383 392 392 385 370];
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(vesselDir), VESSEL_DIR = vesselDir; end
if nargin >= 2 && ~isempty(roiTag),    ROI_TAG    = roiTag;    end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

files = dir(fullfile(VESSEL_DIR, [ROI_TAG '*.tif']));
assert(~isempty(files), 'No %s*.tif under %s', ROI_TAG, VESSEL_DIR);

fprintf('vessel_mode_tail_vs_z: %d stack(s) matching %s* in %s\n', ...
    numel(files), ROI_TAG, VESSEL_DIR);

%% --- pass 1: metadata only, so the z/power layout is visible before the slow part
R = struct('name', {}, 'path', {}, 'zs', {}, 'power', {}, 'nSlice', {}, ...
           'framesPerSlice', {}, 'H', {}, 'W', {}, 'nameZ', {}, 'flat', {});
for i = 1:numel(files)
    tp   = fullfile(files(i).folder, files(i).name);
    info = imfinfo(tp);
    meta = si_meta(info);

    zs   = si_vec(meta, "SI.hStackManager.zs");
    nSl  = si_num(meta, "SI.hStackManager.numSlices");
    fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
    P0   = si_num(meta, "SI.hBeams.stackStartPower");
    P1   = si_num(meta, "SI.hBeams.stackEndPower");
    LC   = si_num(meta, "SI.hBeams.lengthConstants");
    pzOn = si_flag(meta, "SI.hBeams.pzAdjust");
    nCh  = max(numel(si_vec(meta, "SI.hChannels.channelSave")), 1);

    if isempty(zs)
        z0 = si_num(meta, "SI.hStackManager.stackZStartPos");
        z1 = si_num(meta, "SI.hStackManager.stackZEndPos");
        dz = si_num(meta, "SI.hStackManager.stackZStepSize");
        assert(all(isfinite([z0 z1 dz])) && dz ~= 0, ...
            'No zs and no usable stackZStart/End/Step in %s', files(i).name);
        zs = (z0:dz:z1).';
    end
    assert(numel(zs) == nSl, 'zs has %d entries but numSlices = %g in %s', ...
        numel(zs), nSl, files(i).name);
    assert(numel(info) == nSl * fpsl * nCh, ...
        '%s: %d pages != numSlices*framesPerSlice*nCh = %g', ...
        files(i).name, numel(info), nSl * fpsl * nCh);

    % Constant power, or a genuine ScanImage depth ramp?  Decided by ScanImage,
    % never by the presence of lengthConstants.
    flat = ~pzOn || (isfinite(P0) && isfinite(P1) && abs(P1 - P0) < 1e-9);
    if flat
        pw = repmat(P0, numel(zs), 1);
    else
        assert(isfinite(LC) && LC ~= 0, ...
            '%s: pzAdjust is on but lengthConstants is unusable', files(i).name);
        pw = P0 .* exp((zs - zs(1)) ./ LC);
    end

    n = numel(R) + 1;
    R(n).name  = files(i).name;
    R(n).path  = tp;
    R(n).zs    = zs(:);
    R(n).power = pw(:);
    R(n).nSlice = nSl;
    R(n).framesPerSlice = fpsl;
    R(n).H = info(1).Height;
    R(n).W = info(1).Width;
    R(n).flat = flat;
    R(n).nameZ = parse_name_z(files(i).name);
end

% Order rounds by power so the colour ramp is meaningful.
[~, ord] = sort(cellfun(@(p) p(1), {R.power}));
R = R(ord);

fprintf('\n  %-34s %8s %8s %22s %22s\n', ...
    'file', 'power%', 'nSlice', 'zs range (SI motor)', 'filename z (depth)');
fprintf('  %s\n', repmat('-', 1, 100));
for i = 1:numel(R)
    fprintf('  %-34s %8.3g %8d %10.1f .. %-10.1f %10s .. %-10s%s\n', ...
        R(i).name, R(i).power(1), R(i).nSlice, R(i).zs(1), R(i).zs(end), ...
        num2str(R(i).nameZ(1)), num2str(R(i).nameZ(2)), ...
        tern(R(i).flat, '', '  [POWER RAMPED]'));
end
fprintf(['  note: filename z is DEPTH (positive down); SI zs is motor position\n' ...
         '        (negative down).  Everything below uses zs.\n']);

%% --- pass 2: per-slice mode and tail ----------------------------------------
tAll = tic;
for i = 1:numel(R)
    keep = R(i).zs <= Z_TOP;
    idxZ = find(keep);
    nK   = numel(idxZ);
    fprintf('\n[%d/%d] %s -- %d of %d slices at zs <= %g\n', ...
        i, numel(R), R(i).name, nK, R(i).nSlice, Z_TOP);

    modeVal = nan(nK, 1);
    tailVal = nan(nK, 1);
    minVal  = nan(nK, 1);
    satFrac = nan(nK, 1);

    t = Tiff(R(i).path, 'r');
    cl = onCleanup(@() t.close());
    if isempty(SAT_LEVEL)
        setDirectory(t, 1);
        satLev = double(intmax(class(t.read())));
    else
        satLev = SAT_LEVEL;
    end
    tStack = tic;
    for k = 1:nK
        [I, satFrac(k)] = avg_slice(t, idxZ(k), R(i).framesPerSlice, ...
                                    R(i).H, R(i).W, FRAME_STRIDE, satLev);

        % --- mode: leftmost peak of the smoothed histogram -------------------
        Is   = I(1:SUBSAMPLE:end, 1:SUBSAMPLE:end);
        vals = Is(:);  vals = vals(isfinite(vals));
        lo = prctile(vals, 0.5);
        hi = prctile(vals, 99.5);
        if hi <= lo, hi = lo + 1; end
        edges   = linspace(lo, hi, NBINS + 1);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        cs      = smoothdata(histcounts(vals, edges), 'gaussian', SMOOTH_W);
        [pks, locs] = findpeaks(cs, centers, 'MinPeakProminence', max(cs) * MINPROM_FR);
        if isempty(pks)
            [~, im] = max(cs);
            modeVal(k) = centers(im);
        else
            modeVal(k) = min(locs);     % FIRST peak from the left
        end

        % --- tail: brightest NBRIGHT pixels above the low-edge proxy ---------
        mn = lo;
        vv = I(:) - mn;
        vv = vv(isfinite(vv));
        tailVal(k) = mean(maxk(vv, min(NBRIGHT, numel(vv))));
        minVal(k)  = mn;

        if mod(k, max(1, floor(nK/6))) == 0 || k == nK
            fprintf('    z %6.1f  (%3d/%3d)  mode %8.1f  tail %9.1f  sat %.3f%%\n', ...
                R(i).zs(idxZ(k)), k, nK, modeVal(k), tailVal(k), 100*satFrac(k));
        end
    end
    clear cl
    fprintf('    %.1f s\n', toc(tStack));

    R(i).z    = R(i).zs(idxZ);
    R(i).p    = R(i).power(idxZ);
    R(i).mode = modeVal;
    R(i).tail = tailVal;
    R(i).min  = minVal;
    R(i).satFrac = satFrac;
    R(i).sat  = satFrac > SAT_FRAC_MAX;
    R(i).satLevel = satLev;
    if any(R(i).sat)
        zsat = R(i).z(R(i).sat);
        fprintf('    SATURATED: %d of %d slices (z %.0f .. %.0f) -- tail invalid there\n', ...
            sum(R(i).sat), nK, max(zsat), min(zsat));
    end
    % Two-photon fluorescence goes as the square of the excitation power, so
    % this is what puts rounds acquired at different powers on one scale.
    pf = R(i).p / 100;
    R(i).mode_c = modeVal ./ max(pf.^2, eps);
    R(i).tail_c = tailVal ./ max(pf.^2, eps);
end
fprintf('\ntotal read+metric time: %.1f s\n', toc(tAll));

%% --- overlap agreement -------------------------------------------------------
% The point of overlapping z at different powers: after the P^2 correction the
% rounds should agree where they share depth.  A systematic offset there means
% the correction is not capturing everything (PMT nonlinearity, bleaching,
% or a power calibration that is not actually quadratic).
fprintf(['\noverlap check (median ratio of P^2-corrected values, later/earlier round)\n' ...
         '  1.00 = the two rounds agree once scaled by power^2.  Saturated slices are\n' ...
         '  excluded from the tail ratio: a clipped tail reads the ADC ceiling, not the\n' ...
         '  sample, and would drag the ratio toward an arbitrary number.\n']);
for i = 1:numel(R)-1
    a = R(i); b = R(i+1);
    zc = intersect(round(a.z, 6), round(b.z, 6));
    if isempty(zc)
        fprintf('  %s vs %s: no shared z\n', pw_lbl(a), pw_lbl(b)); continue;
    end
    [~, ia] = ismember(zc, round(a.z, 6));
    [~, ib] = ismember(zc, round(b.z, 6));
    okT = ~a.sat(ia) & ~b.sat(ib);
    rM  = median(b.mode_c(ib) ./ a.mode_c(ia), 'omitnan');
    % Same ratio re-expressed with mW instead of % in the denominator.  Only
    % the P^2 factor changes, so this is a pure rescale of the number above.
    mwA = interp1(CAL_PCT, CAL_MW, min(max(a.p(1), min(CAL_PCT)), max(CAL_PCT)), 'pchip');
    mwB = interp1(CAL_PCT, CAL_MW, min(max(b.p(1), min(CAL_PCT)), max(CAL_PCT)), 'pchip');
    k   = ((b.p(1)/a.p(1)) / (mwB/mwA))^2;
    if any(okT)
        rT = median(b.tail_c(ib(okT)) ./ a.tail_c(ia(okT)), 'omitnan');
        sT = sprintf('x%.2f (%d clean z)', rT, sum(okT));
        sM = sprintf('x%.2f', rT * k);
    else
        sT = 'ALL SATURATED';  sM = '--';
    end
    fprintf('  %-9s vs %-9s %3d shared z (%4.0f..%4.0f)  tail %-22s mode x%.2f   [mW-cal: %.0f->%.0f mW, tail %s]\n', ...
        pw_lbl(a), pw_lbl(b), numel(zc), max(zc), min(zc), sT, rM, mwA, mwB, sM);
end

%% --- figure -----------------------------------------------------------------
cols = [0.10 0.45 0.80;    % increasing power -> cool to warm
        0.15 0.65 0.45;
        0.90 0.60 0.10;
        0.80 0.20 0.20];
if numel(R) > size(cols, 1), cols = turbo(numel(R)); end

f = figure('Color', 'w', 'Position', [60 60 1300 860], 'Name', 'mode/tail vs z');
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

panels = { 'tail',   'Tail (mean of brightest %d px)', false
           'tail_c', 'Tail / P^2',                     true
           'mode',   'Mode (first histogram peak)',    false
           'mode_c', 'Mode / P^2',                     true };

for q = 1:size(panels, 1)
    ax = nexttile; hold(ax, 'on');
    isTail = contains(panels{q,1}, 'tail');
    for i = 1:numel(R)
        y = R(i).(panels{q,1});
        if isTail && any(R(i).sat)
            % Saturated stretch drawn dotted and kept out of the legend: it is
            % the ADC ceiling, not a measurement, and it is exactly the part
            % that looks most like a clean plateau.
            ys = y;  ys(~R(i).sat) = NaN;
            yc = y;  yc(R(i).sat)  = NaN;
            plot(ax, R(i).z, ys, ':', 'Color', cols(i,:), 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
            plot(ax, R(i).z, yc, '-', 'Color', cols(i,:), 'LineWidth', 1.4, ...
                'DisplayName', pw_lbl(R(i)));
        else
            plot(ax, R(i).z, y, '-', 'Color', cols(i,:), 'LineWidth', 1.4, ...
                'DisplayName', pw_lbl(R(i)));
        end
    end
    grid(ax, 'on'); box(ax, 'on');
    % zs is negative downward, so reversing X puts the surface on the left and
    % increasing depth to the right while the tick labels stay in SI convention.
    set(ax, 'XDir', 'reverse');
    xlabel(ax, 'ScanImage z (\mum, SI.hStackManager.zs)');
    if panels{q,3}
        set(ax, 'YScale', 'log');
    end
    if isTail
        ylabel(ax, sprintf('tail%s', tern(panels{q,3}, ' / P^2', ' (ADU)')));
        title(ax, sprintf([panels{q,2} '  --  dotted = SATURATED'], NBRIGHT), ...
            'FontWeight', 'normal');
    else
        ylabel(ax, sprintf('mode%s', tern(panels{q,3}, ' / P^2', ' (ADU)')));
        title(ax, panels{q,2}, 'FontWeight', 'normal');
    end
    if q == 1
        legend(ax, 'Location', 'northeast', 'Box', 'off');
    end
end

title(tl, sprintf('%s  --  %s  --  mode & tail vs depth, one colour per laser-power round', ...
    VESSEL_DIR, ROI_TAG), 'FontWeight', 'bold', 'Interpreter', 'none');

if SAVE_FIG
    base = fullfile(VESSEL_DIR, sprintf('vessel_mode_tail_vs_z_%s', ROI_TAG));
    try
        exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor', 'white');
        exportgraphics(f, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    catch
        print(f, [base '.png'], '-dpng', '-r200');
    end
    save([base '.mat'], 'R', 'Z_TOP', 'NBRIGHT', 'VESSEL_DIR', 'ROI_TAG');
    fprintf('\nsaved %s.png/.pdf/.mat\n', base);
end

end

% ============================== local functions ==============================

function [I, satFrac] = avg_slice(t, iz, framesPerSlice, H, W, stride, satLev)
% Mean image of one z slice, plus the fraction of RAW samples at the ADC
% ceiling.  Averaging before the histogram is deliberate and matches
% F_attenuation.m -- the mode of an averaged image is not the mode of the
% single frames, so the two conventions are not interchangeable.  Saturation
% must be counted before that average, though: a pixel clipped in 10 of 30
% frames averages well below the ceiling and would otherwise look clean.
acc = zeros(H, W);
n = 0; nSat = 0;
for fr = 1:stride:framesPerSlice
    setDirectory(t, (iz - 1) * framesPerSlice + fr);
    F = t.read();
    nSat = nSat + nnz(double(F) >= satLev);
    acc = acc + double(F);
    n = n + 1;
end
I = acc / n;
satFrac = nSat / (n * H * W);
end

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

function tf = si_flag(meta, key)
ln = lower(strtrim(si_line(meta, key)));
tf = startsWith(ln, "true");
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

function z = parse_name_z(nm)
% The z token in the filename, kept ONLY so the report can show how far it is
% from the metadata.  "roi1_z-10-100_lp11" and "roi1_z40-160_lp21" are both
% depth ranges with the sign dropped; nothing downstream uses this.
z = [NaN NaN];
tok = regexp(nm, '_z(-?\d+)-(\d+)_', 'tokens', 'once');
if isempty(tok), return; end
z = [str2double(tok{1}), str2double(tok{2})];
end

function s = pw_lbl(r)
s = sprintf('%g%% lp', r.power(1));
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
