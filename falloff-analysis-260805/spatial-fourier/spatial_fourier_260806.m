function spatial_fourier_260806(dataDir, roiTag, zSurface, powerLut, depthMode, distSpec, channel)
% SPATIAL_FOURIER_260806  Radial spatial power spectrum vs depth.
%
% Companion to falloff_260804.m.  That one asks how much signal survives with
% depth; this one asks what happens to the STRUCTURE -- scattering blurs the
% focus, so fine detail should disappear faster than coarse detail, and the
% spectrum should lose its high-frequency end with depth.
%
% For every z slice: 2-D FFT of the frame-averaged image, power |F|^2, averaged
% in radial bins of spatial frequency k.  Result is P(k, depth).
%
% READS THE RAW TIFFS AT NATIVE PIXEL SIZE, not the 1 um isotropic volume the
% projection script caches.  Resampling 0.8889 -> 1 um is itself a low-pass
% step (Nyquist 0.5625 -> 0.5 cyc/um, and bicubic attenuates the top of the
% band), which would corrupt exactly the measurement being made here.
%
% CONVENTIONS shared with falloff_260804.m: overlap between stacks is split at
% its midpoint; depth comes from DEPTH_MODE; power from the date-aware LUT.
%
% NORMALISATION.  I ~ mW^2, so P = |FFT|^2 ~ mW^4.  Spectra are divided by
% mW^4 to put stacks acquired at different powers on one scale.  The
% shape-normalised panel divides each spectrum by its own low-k value instead,
% which removes brightness entirely and shows only the change in SHAPE.
%
% USAGE
%   spatial_fourier_260806('C:\fall-off\260728_vglut2_vessel','roi1',0,'pre_260723','filename')
%
% OUTPUT (into dataDir)
%   spfft_<roiTag>.png / .pdf / .mat      figure + spectra
%   spfft_<roiTag>_cache.mat              spectra cache; delete to re-read

%% --------------------------- USER PARAMETERS -------------------------------
DATA_DIR   = 'C:\fall-off\260728_vglut2_vessel';
ROI_TAG    = 'roi1';
POWER_LUT  = 'pre_260723';
DEPTH_MODE = 'filename';
Z_SURFACE  = 0;
P_EXPONENT = 2;             % signal ~ mW^P; power spectrum ~ mW^(2P)
CHANNEL    = [];            % SI channel carrying the vessel dye when more than
                            % one is saved (260806: channelSave [1;3], vessels on
                            % 3).  Pages interleave, so the wrong channel gives a
                            % spectrum of the other label with no error.

NK         = 120;           % radial frequency bins
% [first depth, step] for the single-plane spectra.  Same convention as the
% falloff figure: [25 50] wildtype, [50 100] shiverer.
DIST_SPEC   = [50 100];
DEPTH_CLIM = [0 600];       % FIXED jet colourbar, same as every other figure
PERIOD_LIM = [2 150];       % plotted band, spatial PERIOD in um.  The FFT can
                            % only measure periods CROP_UM/m for whole m, so the
                            % coarse end is intrinsically sparse: 400, 200, 133,
                            % 100, 80 ... Stopping at 100 keeps the band where
                            % several measurable periods actually exist.
CROP_UM    = 400;           % CENTRE CROP in MICRONS, not pixels.  400 um is the
                            % largest square that fits inside the smallest FOV in
                            % the project (455 um at zoom 2), so every session can
                            % supply it and all of them get the SAME spatial
                            % aperture.  Without that, sessions at different zoom
                            % have a different lowest measurable frequency and
                            % their spectra are not comparable.
                            % The aperture sets the longest measurable period:
                            % 1 cycle = CROP_UM (a gradient, not a periodicity),
                            % 2 cycles = CROP_UM/2.  Past ~CROP_UM/2 the curve is
                            % the Hann window and any slow illumination gradient,
                            % not structure.  That limit is drawn on the plots.
SCALES_UM  = [3 10 30 100];  % spatial periods for the printed per-scale fit
BAND_UM    = [10 30];        % panel 3: band to integrate, um period
FIT_D      = [10 200];       % depth range for the per-scale exponential fit
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(dataDir),  DATA_DIR   = dataDir;  end
if nargin >= 2 && ~isempty(roiTag),   ROI_TAG    = roiTag;   end
if nargin >= 3 && ~isempty(zSurface), Z_SURFACE  = zSurface; end
if nargin >= 4 && ~isempty(powerLut), POWER_LUT  = powerLut; end
if nargin >= 5 && ~isempty(depthMode),DEPTH_MODE = depthMode;end
if nargin >= 6 && ~isempty(distSpec), DIST_SPEC  = distSpec; end
if nargin >= 7 && ~isempty(channel),  CHANNEL    = channel;  end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fileparts(thisDir));          % falloff-analysis-260805, for the LUT
addpath(fileparts(fileparts(thisDir)));   % repo root

TAG = regexprep(ROI_TAG, '[\*\.\s/\\]+', '_');
TAG = regexprep(TAG, '^_|_$', '');

cacheFile = fullfile(DATA_DIR, sprintf('spfft_%s_cache.mat', TAG));
useCache = false;
if isfile(cacheFile)
    S = load(cacheFile);
    cCh = [];  if isfield(S, 'CHANNEL'), cCh = S.CHANNEL; end
    useCache = isequal(cCh, CHANNEL);   % another channel is another label
    if ~useCache
        fprintf('cache was channel %s, now asking for %s -- re-reading\n', ...
            mat2str(cCh), mat2str(CHANNEL));
    end
end
if useCache
    fprintf('using cached spectra: %s  (delete to re-read)\n', cacheFile);
    R = S.R;  kAx = S.kAx;  pxUm = S.pxUm;  cropPx = S.cropPx;  kCnt = S.kCnt;
else
    [R, kAx, pxUm, cropPx, kCnt] = read_spectra(DATA_DIR, ROI_TAG, Z_SURFACE, DEPTH_MODE, NK, CROP_UM, CHANNEL);
    save(cacheFile, 'R', 'kAx', 'pxUm', 'cropPx', 'kCnt', 'DATA_DIR', 'ROI_TAG', ...
        'NK', 'CROP_UM', 'CHANNEL');
    fprintf('cached -> %s\n', cacheFile);
end
fprintf('pixel size %.4f um  ->  Nyquist %.4f cyc/um\n', pxUm, 1/(2*pxUm));

%% --- power normalisation + overlap ownership ---------------------------------
for i = 1:numel(R)
    [R(i).mW, cal] = laser_power_calibration(R(i).pctZ, POWER_LUT);
    R(i).lutName = cal.name;
    % P = |FFT(I)|^2 and I ~ mW^P_EXPONENT, so P ~ mW^(2*P_EXPONENT)
    R(i).specN = R(i).spec ./ (R(i).mW(:) .^ (2*P_EXPONENT));
end

[~, ord] = sort(arrayfun(@(r) min(r.depth), R));
R = R(ord);
lo = -inf(numel(R),1);  hi = inf(numel(R),1);
for i = 1:numel(R)-1
    a = max(min(R(i).depth), min(R(i+1).depth));
    b = min(max(R(i).depth), max(R(i+1).depth));
    if b > a, m = (a+b)/2; else, m = (max(R(i).depth)+min(R(i+1).depth))/2; end
    hi(i) = m;  lo(i+1) = m;
end
fprintf('\n  %-12s %10s %18s %8s\n', 'stack', 'mW', 'owns depths', 'nSlice');
for i = 1:numel(R)
    R(i).own = R(i).depth > lo(i) & R(i).depth <= hi(i);
    fprintf('  %-12s %10.1f %7.0f -> %-8.0f %8d\n', R(i).label, R(i).mW(1), ...
        max(lo(i),min(R(i).depth)), min(hi(i),max(R(i).depth)), sum(R(i).own));
end

% one continuous depth axis from the owned slices
dAll = [];  Sall = [];
for i = 1:numel(R)
    dAll = [dAll; R(i).depth(R(i).own)];        %#ok<AGROW>
    Sall = [Sall; R(i).specN(R(i).own, :)];     %#ok<AGROW>
end
[dAll, o] = sort(dAll);  Sall = Sall(o, :);

%% --- resolution cutoff vs depth ----------------------------------------------
% Where each spectrum drops to a fixed fraction of its own low-k value.  That is
% a brightness-free number, so it isolates loss of DETAIL from loss of signal.
kLoIdx = kAx > 0.01 & kAx < 0.05;         % coarse-structure reference band
FRAC   = 0.01;
kCut = nan(numel(dAll),1);
for j = 1:numel(dAll)
    p = Sall(j,:);
    ref = mean(p(kLoIdx));
    if ~isfinite(ref) || ref <= 0, continue; end
    r = p / ref;
    idx = find(r < FRAC & kAx(:).' > 0.05, 1);
    if isempty(idx) || idx < 2, continue; end
    % linear interpolation in log power for a sub-bin cutoff
    k1 = kAx(idx-1);  k2 = kAx(idx);
    r1 = log(max(r(idx-1),eps));  r2 = log(max(r(idx),eps));
    kCut(j) = k1 + (log(FRAC)-r1)/(r2-r1) * (k2-k1);
end
fprintf(['\ncutoff = k where P falls to %g%% of its 0.01-0.05 cyc/um value\n' ...
         '  %.4f -> %.4f cyc/um  (period %.1f -> %.1f um) over depth %g -> %g um\n'], ...
    100*FRAC, kCut(find(isfinite(kCut),1)), kCut(find(isfinite(kCut),1,'last')), ...
    1/kCut(find(isfinite(kCut),1)), 1/kCut(find(isfinite(kCut),1,'last')), ...
    dAll(find(isfinite(kCut),1)), dAll(find(isfinite(kCut),1,'last')));

%% --- attenuation length per spatial scale (printed, not plotted) -------------
% Power goes as intensity^2, so these should be HALF the falloff length from
% falloff_260804.m.  That is the cross-check between the two analyses.
fprintf('\nattenuation length per spatial scale (fit %g-%g um depth)\n', FIT_D);
fprintf('  %10s %10s %10s\n', 'period um', 'l (um)', 'R^2');
for j = 1:numel(SCALES_UM)
    [~, ik] = min(abs((1./kAx) - SCALES_UM(j)));
    y = Sall(:, ik);
    m = dAll >= FIT_D(1) & dAll <= FIT_D(2) & y > 0;
    if nnz(m) < 6, continue; end
    pf  = polyfit(dAll(m), log(y(m)), 1);
    res = log(y(m)) - polyval(pf, dAll(m));
    r2  = 1 - sum(res.^2)/sum((log(y(m)) - mean(log(y(m)))).^2);
    fprintf('  %10g %10.1f %10.3f\n', SCALES_UM(j), -1/pf(1), r2);
end

%% --- figure -------------------------------------------------------------------
% Two panels: the power-corrected spectrum, and the same curves each rescaled
% to 0-1 on its own min/max.  x = spatial PERIOD, log scale.
nC = 256;  cmap = jet(nC);
ramp = @(v,lim) cmap(min(nC,max(1,round(1+(nC-1)*(v-lim(1))/diff(lim)))), :);
kNy   = 1/(2*pxUm);
fovUm = pxUm * cropPx;
pAx   = 1 ./ kAx;

ib = find(pAx >= PERIOD_LIM(1) & pAx <= PERIOD_LIM(2));
% Strictly inside the band, and never bin 1 (k starts at 0, so bin 1 has
% infinite period width).  No expansion past PERIOD_LIM(2): bins are uniform in
% k, so at the long-period end they are sparse -- on the coarser-pixel sessions
% the next bin out sits at 237 um, and including it drags a line clean off the
% axis and makes it look as if there were data out there.
ib = max(min(ib), 2) : max(ib);
[pB, iSrt] = sort(pAx(ib));
idxB = ib(iSrt);
pTick = [2 5 10 20 50 100 150];

fprintf(['crop %d x %d px = %.0f x %.0f um\n' ...
         'band plotted %g - %g um period; bins used %.1f - %.1f um (Nyquist %.2f)\n'], ...
    cropPx, cropPx, fovUm, fovUm, PERIOD_LIM, min(pB), max(pB), 1/kNy);

maxD  = max(dAll);
dSpec = DIST_SPEC;
dWant = dSpec(1) : dSpec(2) : maxD;

f = figure('Color','w','Position',[40 80 1700 500]);
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact','Padding','compact');
tl.Title.FontSize = 10;

% (1) laser-power corrected
ax = nexttile; hold(ax,'on');
for j = 1:numel(dWant)
    [~, q] = min(abs(dAll - dWant(j)));
    plot(ax, pB, Sall(q, idxB), '-', 'Color', ramp(dAll(q), DEPTH_CLIM), 'LineWidth', 1.3);
end
set(ax,'XScale','log','YScale','log'); grid(ax,'on'); box(ax,'on');
xlim(ax, PERIOD_LIM); set(ax,'XTick',pTick);
colormap(ax, cmap); caxis(ax, DEPTH_CLIM);
cb = colorbar(ax); cb.Label.String = 'depth below pia (\mum)'; cb.Direction = 'reverse';
xlabel(ax, 'spatial period (\mum)');
ylabel(ax, sprintf('power / mW^%g', 2*P_EXPONENT));
title(ax, sprintf('laser-power corrected, %g:%g:%g um depth', ...
    dSpec(1), dSpec(2), max(dWant)), 'FontWeight','normal');

% (2) each curve rescaled to 0-1 on its OWN min and max
ax = nexttile; hold(ax,'on');
for j = 1:numel(dWant)
    [~, q] = min(abs(dAll - dWant(j)));
    y = Sall(q, idxB);
    y = (y - min(y)) / max(max(y) - min(y), eps);
    plot(ax, pB, y, '-', 'Color', ramp(dAll(q), DEPTH_CLIM), 'LineWidth', 1.3);
end
set(ax,'XScale','log','YScale','linear'); grid(ax,'on'); box(ax,'on');
xlim(ax, PERIOD_LIM); set(ax,'XTick',pTick);  ylim(ax, [0 1]);
xlabel(ax, 'spatial period (\mum)');
ylabel(ax, 'normalised power (0-1)');
title(ax, 'normalised individually', 'FontWeight','normal');

% (3) area under the NORMALISED spectrum over BAND_UM, at every depth.
% Same 0-1 normalisation as panel 2 (each slice on its own min/max), then
% trapezoidal integration between BAND_UM(1) and BAND_UM(2).
inBnd = pB >= BAND_UM(1) & pB <= BAND_UM(2);
bandArea = nan(numel(dAll),1);
for j = 1:numel(dAll)
    y = Sall(j, idxB);
    y = (y - min(y)) / max(max(y) - min(y), eps);
    bandArea(j) = trapz(pB(inBnd), y(inBnd));
end

ax = nexttile; hold(ax,'on');
yPlot = bandArea;  yPlot(yPlot <= 0) = NaN;   % log axis cannot take the floor
plot(ax, dAll, yPlot, 'k-', 'LineWidth', 1.5);
set(ax,'YScale','log'); grid(ax,'on'); box(ax,'on');
xlim(ax, [-50 600]); set(ax,'XTick',0:100:600);
xlabel(ax, 'depth (\mum)');
ylabel(ax, sprintf('integrated band power, %g-%g \\mum', BAND_UM));
title(ax, sprintf('integrated %g-%g \\mum band', BAND_UM), 'FontWeight','normal');
fprintf('\narea under normalised spectrum %g-%g um: %.1f at depth 0, %.1f at depth %g\n', ...
    BAND_UM, bandArea(find(dAll>=0,1)), bandArea(end), dAll(end));

title(tl, {sprintf('%s  --  %s  --  spatial power spectrum vs depth', DATA_DIR, ROI_TAG), ...
    sprintf('native %.4f um/px (Nyquist period %.2f um)   |   %.0f um crop   |   Hann   |   LUT %s', ...
        pxUm, 1/kNy, fovUm, R(1).lutName), ''}, ...
    'FontWeight','bold','Interpreter','none');

if SAVE_FIG
    base = fullfile(DATA_DIR, sprintf('spfft_%s', TAG));
    exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'R', 'kAx', 'dAll', 'Sall', 'kCut', 'pxUm', 'POWER_LUT');
    fprintf('\nsaved %s.png/.pdf/.mat\n', base);
end
end

% ============================== read =========================================

function [R, kAx, pxUm, cropPx, kCnt] = read_spectra(DATA_DIR, ROI_TAG, Z_SURFACE, DEPTH_MODE, NK, CROP_UM, CHANNEL)
if contains(ROI_TAG,'*') || endsWith(lower(ROI_TAG),'.tif')
    pat = ROI_TAG;  if ~endsWith(lower(pat),'.tif'), pat = [pat '.tif']; end
else
    pat = [ROI_TAG '*.tif'];
end
files = dir(fullfile(DATA_DIR, pat));
assert(~isempty(files), 'No stacks matching ''%s'' under %s', ROI_TAG, DATA_DIR);

R = struct('name',{},'label',{},'depth',{},'pctZ',{},'spec',{});
P0all = nan(numel(files),1);  M = cell(numel(files),1);
for i = 1:numel(files)
    info = imfinfo(fullfile(files(i).folder, files(i).name));
    M{i} = {si_meta(info), info};
    P0all(i) = si_num(M{i}{1}, "SI.hBeams.stackStartPower");
end
[~, ord] = sort(P0all);  files = files(ord);  M = M(ord);

kAx = [];  pxUm = [];  cropPx = [];  kCnt = [];
for i = 1:numel(files)
    meta = M{i}{1};  info = M{i}{2};
    tp   = fullfile(files(i).folder, files(i).name);

    n    = si_num(meta, "SI.hStackManager.numSlices");
    fpsl = si_num(meta, "SI.hStackManager.framesPerSlice");
    dzS  = si_num(meta, "SI.hStackManager.stackZStepSize");
    z0S  = si_num(meta, "SI.hStackManager.stackZStartPos");
    P0   = si_num(meta, "SI.hBeams.stackStartPower");
    P1   = si_num(meta, "SI.hBeams.stackEndPower");
    LC   = si_num(meta, "SI.hBeams.lengthConstants");
    pzOn = startsWith(lower(strtrim(si_line(meta,"SI.hBeams.pzAdjust"))), "true");
    zoom = si_num(meta, "SI.hRoiManager.scanZoomFactor");
    ppl  = si_num(meta, "SI.hRoiManager.pixelsPerLine");
    zs   = si_vec(meta, "SI.hStackManager.zs");
    chSave = si_vec(meta, "SI.hChannels.channelSave");
    nCh    = max(numel(chSave), 1);
    chIdx  = pick_channel(chSave, CHANNEL, files(i).name);
    assert(numel(info) == n*fpsl*nCh, '%s: page count mismatch', files(i).name);

    H = info(1).Height;  W = info(1).Width;
    if ~isfinite(ppl) || ppl <= 0, ppl = W; end
    px = 1.7778 * (512/ppl) / zoom;
    if isempty(pxUm), pxUm = px; end
    assert(abs(px - pxUm) < 1e-9, 'stacks have different pixel sizes');

    % depth axis
    switch lower(DEPTH_MODE)
    case 'filename'
        t = regexp(files(i).name, '_z(-?\d+)-(\d+)_', 'tokens', 'once');
        assert(~isempty(t), 'No _z<start>-<end>_ token in %s', files(i).name);
        d = str2double(t{1}) + (0:n-1).' * abs(dzS);
    case 'motor'
        d = Z_SURFACE - (z0S + (0:n-1).' * dzS);
    otherwise
        error('DEPTH_MODE must be ''filename'' or ''motor''');
    end

    % per-slice power setpoint
    if numel(zs) ~= n, zs = (z0S + (0:n-1)*dzS).'; end
    if pzOn && ~(isfinite(P0) && isfinite(P1) && abs(P1-P0) < 1e-9)
        zRef = si_num(meta, "SI.hStackManager.zPowerReference");
        if ~isfinite(zRef), zRef = zs(1); end
        pctZ = P0 .* exp((zs(:) - zRef) ./ LC);
        lbl  = sprintf('%g-%g%% lp', P0, P1);
    else
        pctZ = repmat(P0, n, 1);
        lbl  = sprintf('%g%% lp', P0);
    end

    % CENTRE CROP to a fixed PHYSICAL size, so every session is analysed through
    % the same aperture regardless of zoom.
    cp = min([round(CROP_UM/pxUm), H, W]);
    r0 = floor((H - cp)/2) + 1;  c0 = floor((W - cp)/2) + 1;
    rIdx = r0:r0+cp-1;  cIdx = c0:c0+cp-1;
    if isempty(cropPx), cropPx = cp; end
    assert(cp == cropPx, 'stacks give different crop sizes');
    Hc = cp;  Wc = cp;

    % frequency grid + radial bins (same for every slice, built once)
    if isempty(kAx)
        fy = (-floor(Hc/2):ceil(Hc/2)-1).' / (Hc*pxUm);
        fx = (-floor(Wc/2):ceil(Wc/2)-1)   / (Wc*pxUm);
        KR = hypot(repmat(fx, Hc, 1), repmat(fy, 1, Wc));
        kMax = min(max(abs(fx)), max(abs(fy)));
        kEdge = linspace(0, kMax, NK+1);
        kAx = (kEdge(1:end-1) + kEdge(2:end)) / 2;
        bin = discretize(KR(:), kEdge);
        keep = ~isnan(bin);
        binK = bin(keep);
        cnt  = accumarray(binK, 1, [NK 1]);
        kCnt = cnt;   % coefficients per annulus -- annuli have very different areas
    end
    % Hann window: the FOV edges are a hard discontinuity and without a window
    % they inject a k^-2 cross that swamps the real high-frequency content.
    win = hann(Hc) * hann(Wc).';

    fprintf('\n[%d/%d] %s  %s  %d slices x %d frames  %.4f um/px\n', ...
        i, numel(files), files(i).name, lbl, n, fpsl, pxUm);
    spec = zeros(n, NK);
    t = Tiff(tp,'r');  cl = onCleanup(@() t.close());
    tS = tic;
    for k = 1:n
        acc = zeros(H, W);
        for fr = 1:fpsl
            % Channels interleave page by page within each frame.
            setDirectory(t, ((k-1)*fpsl + fr-1)*nCh + chIdx);
            acc = acc + double(t.read());
        end
        I = acc(rIdx, cIdx) / fpsl;          % centre crop, fixed physical size
        J = (I - mean(I(:))) .* win;
        P = abs(fftshift(fft2(J))).^2 / (Hc*Wc);
        s = accumarray(binK, P(keep), [NK 1]) ./ max(cnt, 1);
        spec(k,:) = s.';
        if mod(k, max(1,floor(n/4))) == 0 || k == n
            fprintf('    %3d/%3d  depth %6.1f um\n', k, n, d(k));
        end
    end
    clear cl
    fprintf('    %.1f s\n', toc(tS));

    q = numel(R)+1;
    R(q).name = files(i).name;  R(q).label = lbl;
    R(q).depth = d;  R(q).pctZ = pctZ;  R(q).spec = spec;
end
end

% ========================= ScanImage metadata helpers ========================

function chIdx = pick_channel(chSave, CHANNEL, name)
% Position within the interleaved pages of the requested SI channel.
nCh = max(numel(chSave), 1);
if nCh <= 1
    chIdx = 1;
    if ~isempty(CHANNEL) && ~isempty(chSave) && chSave(1) ~= CHANNEL
        warning('%s saves only channel %g; CHANNEL = %g ignored', name, chSave(1), CHANNEL);
    end
    return
end
assert(~isempty(CHANNEL), ['%s saves %d channels (%s) -- set CHANNEL to the ' ...
    'one carrying the vessel dye'], name, nCh, mat2str(chSave(:).'));
chIdx = find(chSave == CHANNEL, 1);
assert(~isempty(chIdx), '%s: channel %g not saved (saved: %s)', ...
    name, CHANNEL, mat2str(chSave(:).'));
end

function meta = si_meta(info)
meta = "";
if isfield(info(1),'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
end
assert(strlength(meta) > 0, 'No ScanImage metadata');
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
