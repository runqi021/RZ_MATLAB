function R = stitch_shift_null_260910(coordFile, varargin)
%STITCH_SHIFT_NULL_260910  Is the curated stitch better than a random placement?
%
%   R = stitch_shift_null_260910(coordFile)
%   R = stitch_shift_null_260910(coordFile, 'Ch', 3, 'NShift', 500, ...)
%
% Builds a null distribution for the manual alignment: take the curated coords,
% displace one tile of each seam by a uniform random (dx, dy), score the
% registration, and repeat.  The curated placement is then read against that
% null.  3 tiles = 2 seams, so NShift draws per seam.
%
% THE SHIFT IS DRAWN FROM AN ANNULUS, not a disc: |dx| and |dy| are each uniform
% over the integers in [ShiftMin, ShiftMax] with a random sign.  ShiftMin exists
% because a 1 um draw is not a misalignment -- it is the curated placement plus
% noise, and such draws would score almost as well as the curated one and pile
% up in the null's upper tail, making a good alignment look unremarkable.  The
% null has to represent placements a human would visibly reject.
%
% THE METRIC.  "Yellow = registered" in the seam view, summed over the overlap:
%
%       YF = sum( min(a,b) ) / sum( max(a,b) )
%
% a soft intersection-over-union of the two overlap images.  Structure that both
% tiles see lands in the numerator AND the denominator; structure only one of
% them sees lands in the denominator alone.  1 = identical, 0 = disjoint.
%
% BACKGROUND SUBTRACTION IS NOT OPTIONAL HERE, it is the whole reason the metric
% works.  A flat background sits in min() and max() equally, so on raw data YF
% climbs toward 1 for ANY shift and the alignment signal is a ripple on top of a
% large constant.  Each patch therefore has its own 25th percentile removed and
% is divided by its own (99th - 25th) before the ratio -- per tile, so a
% brightness difference between tiles cannot register as a difference in
% structure either.
%
% Pearson NCC over the same pixels is reported next to it.  YF is the quantity
% asked for and is the one the eye is using; NCC is the standard registration
% score and exists here to catch YF misbehaving, not because two metrics are
% better than one.  If they disagree, distrust YF.
%
% THE COMPARISON WINDOW IS FIXED, and that matters more than the metric choice.
% Recomputing the overlap at each shift would change both its SIZE and its
% CONTENT, so the score would partly be measuring how many pixels happened to be
% compared -- a shift that shrinks the overlap onto a bright vessel could score
% well for no good reason.  Instead the curated overlap is eroded by ShiftMax on
% every side, tile A's patch is held FIXED, and only tile B's sampling moves.
% Every draw then scores the same pixel count against the same reference.
% Erosion by exactly ShiftMax is what guarantees tile B always has real data to
% offer: the shifted window can reach the original overlap boundary but never
% cross it.
%
% Shifts are INTEGER microns, drawn independently in x and y.  Integer because
% the voxel is 1 um and the tiles are placed on the pixel grid; a sub-pixel shift
% would need interpolation, which would itself blur the patch and confound the
% score with a smoothing.  Note that forcing BOTH axes past ShiftMin means a pure
% x or pure y misalignment is never sampled -- set ShiftMin to 0 if that failure
% mode is the one of interest.
%
% Z IS HELD AT THE CURATED VALUE.  Only x and y are randomised.
%
% OUTPUT  -> <coord dir>\shift_null_ch<N>\
%   shift_null_ch<N>.png/.pdf    per seam: YF null + curated, YF vs |shift|, NCC null
%   shift_null_ch<N>.mat         every draw, both metrics, and the curated values
%   shift_null_ch<N>.csv         the summary table printed to the console
%
% Runqi Zhang / 2026-09-10

p = inputParser;
p.addParameter('Ch',        1);
p.addParameter('NShift',    500, @(x)isscalar(x) && x >= 1);
p.addParameter('ShiftMin',  3,   @(x)isscalar(x) && x >= 0);
p.addParameter('ShiftMax',  30,  @(x)isscalar(x) && x >= 1);
p.addParameter('Seed',      260910);
p.addParameter('ProjMode',  'max');    % matches the seam view's projection
p.addParameter('OutDir',    '');
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(opt.Seed);

C = load(coordFile, 'stitch'); St = C.stitch;
coordDir = fileparts(coordFile);
outDir = opt.OutDir;
if isempty(outDir), outDir = fullfile(coordDir, sprintf('shift_null_ch%d', opt.Ch)); end
if ~isfolder(outDir), mkdir(outDir); end

n  = numel(St.metaFiles);
SM = opt.ShiftMax;
assert(opt.ShiftMin <= SM, 'ShiftMin (%g) exceeds ShiftMax (%g)', opt.ShiftMin, SM);

fprintf('\n============== stitch_shift_null_260910 ==============\n');
fprintf('coords   : %s   (saved %s)\n', coordFile, St.savedOn);
fprintf('channel  : %d   projection %s\n', opt.Ch, opt.ProjMode);
fprintf('shifts   : %d per seam, integer |d| uniform in [%d %d] um with random sign, x and y\n', ...
        opt.NShift, opt.ShiftMin, SM);

%% ------------------------------ tile images --------------------------------
XY = cell(n,1);
for i = 1:n
    M = load(St.metaFiles{i}, 'meta');
    XY{i} = tile_xy(M.meta, opt.Ch, opt.ProjMode);
    fprintf('  loaded %-22s %s\n', St.names{i}, mat2str(size(XY{i})));
end

ox = round(St.x(:) + St.manX(:));  oy = round(St.y(:) + St.manY(:));
ox = ox - min(ox);                 oy = oy - min(oy);
h = St.tileH(:); w = St.tileW(:);

%% -------------------------------- seams ------------------------------------
R = struct('a',{},'b',{},'yf0',{},'ncc0',{},'yf',{},'ncc',{},'dx',{},'dy',{},'npx',{});
for i = 1:n
    for j = i+1:n
        c0 = max(ox(i), ox(j));  c1 = min(ox(i)+w(i), ox(j)+w(j)) - 1;
        r0 = max(oy(i), oy(j));  r1 = min(oy(i)+h(i), oy(j)+h(j)) - 1;
        if c1 < c0 || r1 < r0, continue; end
        % erode by the maximum shift, so a displaced sample never leaves tile j
        c0 = c0 + SM; c1 = c1 - SM;  r0 = r0 + SM; r1 = r1 - SM;
        if c1 < c0 || r1 < r0
            warning('seam %d-%d: overlap %d px is too small to erode by %d -- skipped', ...
                    i, j, min(c1-c0, r1-r0) + 2*SM, SM);
            continue
        end
        ri = (r0:r1) - oy(i) + 1;  ci = (c0:c1) - ox(i) + 1;
        rj = (r0:r1) - oy(j) + 1;  cj = (c0:c1) - ox(j) + 1;

        A = XY{i}(ri, ci);
        a = prep(A);

        [yf0, ncc0] = score(a, prep(XY{j}(rj, cj)));

        dx = signed_shift(opt.NShift, opt.ShiftMin, SM);
        dy = signed_shift(opt.NShift, opt.ShiftMin, SM);
        yf = zeros(opt.NShift,1); ncc = zeros(opt.NShift,1);
        for k = 1:opt.NShift
            B = XY{j}(rj + dy(k), cj + dx(k));
            [yf(k), ncc(k)] = score(a, prep(B));
        end

        R(end+1) = struct('a', i, 'b', j, 'yf0', yf0, 'ncc0', ncc0, ...
                          'yf', yf, 'ncc', ncc, 'dx', dx, 'dy', dy, ...
                          'npx', numel(ri)*numel(ci)); %#ok<AGROW>
        fprintf(['\nseam %d-%d  (%s vs %s)  window %d x %d px\n' ...
                 '  curated   YF %.4f   NCC %+.4f\n' ...
                 '  null      YF %.4f +/- %.4f   (min %.4f  max %.4f)\n' ...
                 '  curated is at percentile %.2f of the YF null, z = %+.2f\n' ...
                 '  draws scoring >= curated: %d / %d\n'], ...
                 i, j, St.names{i}, St.names{j}, numel(ri), numel(ci), ...
                 yf0, ncc0, mean(yf), std(yf), min(yf), max(yf), ...
                 100*mean(yf < yf0), (yf0 - mean(yf))/max(eps,std(yf)), ...
                 nnz(yf >= yf0), opt.NShift);
    end
end
assert(~isempty(R), 'no usable seam found');

%% -------------------------------- figure -----------------------------------
nS = numel(R);
f = figure('Color','w','Units','pixels','Position',[40 40 1350 380*nS]);
drawnow; got = get(f,'Position'); FS = min(1, got(3)/1350);
set(f,'DefaultAxesFontSize',max(7,10*FS));
tl = tiledlayout(f, nS, 3, 'TileSpacing','compact','Padding','compact');

for s = 1:nS
    r = R(s);
    nm = sprintf('seam %d-%d', r.a, r.b);

    ax = nexttile(tl, (s-1)*3 + 1);
    histogram(ax, r.yf, 30, 'FaceColor',[.6 .6 .6], 'EdgeColor','none'); hold(ax,'on');
    yl = ylim(ax);
    plot(ax, [r.yf0 r.yf0], yl, '-', 'Color',[.85 .2 .2], 'LineWidth',2);
    text(ax, r.yf0, yl(2), '  curated', 'Color',[.85 .2 .2], 'VerticalAlignment','top');
    xlabel(ax, 'YF = \Sigma min(a,b) / \Sigma max(a,b)'); ylabel(ax, 'draws');
    title(ax, sprintf('%s   YF null, n = %d', nm, numel(r.yf)), 'FontWeight','normal');
    box(ax,'off'); hold(ax,'off');

    ax = nexttile(tl, (s-1)*3 + 2);
    d = hypot(r.dx, r.dy);
    scatter(ax, d, r.yf, 9, [.45 .45 .45], 'filled', 'MarkerFaceAlpha',0.5); hold(ax,'on');
    plot(ax, 0, r.yf0, 'o', 'MarkerFaceColor',[.85 .2 .2], 'MarkerEdgeColor','none', 'MarkerSize',8);
    xlabel(ax, 'shift magnitude |(dx,dy)|  (\mum)'); ylabel(ax, 'YF');
    title(ax, sprintf('%s   YF vs displacement', nm), 'FontWeight','normal');
    box(ax,'off'); hold(ax,'off');

    ax = nexttile(tl, (s-1)*3 + 3);
    histogram(ax, r.ncc, 30, 'FaceColor',[.6 .7 .85], 'EdgeColor','none'); hold(ax,'on');
    yl = ylim(ax);
    plot(ax, [r.ncc0 r.ncc0], yl, '-', 'Color',[.85 .2 .2], 'LineWidth',2);
    xlabel(ax, 'Pearson NCC'); ylabel(ax, 'draws');
    title(ax, sprintf('%s   NCC null (cross-check)', nm), 'FontWeight','normal');
    box(ax,'off'); hold(ax,'off');
end

% Interpreter is 'none' here (the title carries a Windows path, where a
% backslash is a TeX control sequence and % starts a comment), so the unit has
% to be a literal -- a \mum would print as the four characters \mum.
title(tl, {sprintf(['stitch shift null  --  ch %d, %s projection  --  %d shifts per seam, ' ...
        '|dx| and |dy| uniform integer in [%d, %d] um with random sign'], ...
        opt.Ch, opt.ProjMode, opt.NShift, opt.ShiftMin, SM), ...
        sprintf('%s   (curated %s)   |   window eroded by %d px, tile A fixed, per-patch background removed', ...
        coordFile, St.savedOn, SM)}, ...
      'FontWeight','bold','Interpreter','none','FontSize',max(7,9*FS));

base = fullfile(outDir, sprintf('shift_null_ch%d', opt.Ch));
exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');

T = table(arrayfun(@(r) r.a, R).', arrayfun(@(r) r.b, R).', ...
          arrayfun(@(r) r.npx, R).', arrayfun(@(r) r.yf0, R).', ...
          arrayfun(@(r) mean(r.yf), R).', arrayfun(@(r) std(r.yf), R).', ...
          arrayfun(@(r) max(r.yf), R).', ...
          arrayfun(@(r) 100*mean(r.yf < r.yf0), R).', ...
          arrayfun(@(r) (r.yf0-mean(r.yf))/max(eps,std(r.yf)), R).', ...
          arrayfun(@(r) r.ncc0, R).', arrayfun(@(r) mean(r.ncc), R).', ...
          'VariableNames', {'tileA','tileB','nPixels','YF_curated','YF_null_mean', ...
                            'YF_null_sd','YF_null_max','YF_percentile','YF_z','NCC_curated','NCC_null_mean'});
writetable(T, [base '.csv']);
save([base '.mat'], 'R', 'St', 'opt');
fprintf('\nsaved %s.png / .pdf / .csv / .mat\n', base);
end

% =============================================================================
function d = signed_shift(n, lo, hi)
%SIGNED_SHIFT  n integers uniform over {-hi..-lo} union {lo..hi}.
if lo <= 0, d = randi([-hi hi], n, 1); return; end
d = randi([lo hi], n, 1) .* (2*randi([0 1], n, 1) - 1);
end

function y = prep(X)
%PREP  Background out, scale normalised, per patch.
% The 25th percentile is the background estimate: on a sparse label most pixels
% ARE background, so it sits safely below the signal, and unlike the mean it is
% not dragged up by a few bright somata.  Clipping at 0 afterwards is what makes
% background pixels drop out of both min() and max() instead of inflating both.
v = single(X(:));
p = prctile(v, [25 99]);
if p(2) <= p(1), p(2) = p(1) + 1; end
y = max(0, min(1, (single(X) - p(1)) / (p(2) - p(1))));
end

function [yf, ncc] = score(a, b)
mn = min(a, b);  mx = max(a, b);
yf = sum(mn(:)) / max(eps, sum(mx(:)));
av = a(:) - mean(a(:));  bv = b(:) - mean(b(:));
den = sqrt(sum(av.^2) * sum(bv.^2));
if den <= 0, ncc = 0; else, ncc = sum(av .* bv) / den; end
end

function XY = tile_xy(meta, ch, mode)
ci = find(meta.channels == ch, 1);
assert(~isempty(ci), 'channel %d is not in this tile', ch);
tif = meta.outTif{ci};
if ~isfile(tif)
    [~, nm, ex] = fileparts(tif);
    tif = fullfile(fileparts(meta.tifPath), [nm ex]);
end
t = Tiff(tif, 'r'); c = onCleanup(@() close(t));
A = t.read(); nz = numel(meta.depth);
if strcmpi(mode,'avg'), XY = zeros(size(A),'single'); else, XY = -inf(size(A),'single'); end
t.setDirectory(1);
for k = 1:nz
    P = single(t.read());
    if strcmpi(mode,'avg'), XY = XY + P; else, XY = max(XY, P); end
    if k < nz, t.nextDirectory(); end
end
if strcmpi(mode,'avg'), XY = XY / nz; end
end
