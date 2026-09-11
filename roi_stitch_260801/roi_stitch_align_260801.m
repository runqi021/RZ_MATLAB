function roi_stitch_align_260801()
%% roi_stitch_align_260801  Manually stitch overlapping FOVs; write fov_offsets.mat
% -----------------------------------------------------------------------
% Step 1 of the stitch-first pipeline. You align the FOVs by eye; the matcher
% then treats "same place" as the definition of "same neuron", so the quality of
% the whole result is set here.
%
% HOW IT WORKS. FOVs are placed one at a time. The first (largest, most connected)
% is the anchor and never moves. Every later FOV is shown against the already-
% placed neighbours it overlaps, seeded with an FFT phase-correlation guess, and
% you nudge it until the vasculature and somata line up.
%
% Phase correlation seeds the guess rather than deciding it: it whitens magnitude
% and keeps only phase, so it locks onto structure rather than brightness, but it
% assumes pure translation and one dominant shift. Where that assumption breaks
% (little overlap, different depths, tissue moved between recordings) it returns a
% confident wrong answer -- which is why you confirm every pair.
%
% CONTROLS
%   arrow keys      nudge STEP_UM (0.5 um)
%   shift + arrow   nudge 10x
%   r               reset this FOV to the raw stage position
%   p               re-seed from phase correlation
%   b               cycle which already-placed neighbour is shown underneath
%   t               toggle red/green <-> checkerboard
%   Enter           accept and move to the next FOV
%   s               skip (keep raw stage position)
%   q               save what is done so far and quit
%
% STANDALONE: writes only <outDir>\fov_offsets.mat. Touches no other pipeline.

cfg = stitch_cfg_260801();
tf  = fullfile(cfg.outDir,'fov_table.mat');
assert(isfile(tf), 'Run roi_stitch_scan_260801 first (%s missing).', tf);
S = load(tf); F = S.F; nF = numel(F);

STEP_UM = 0.5;

% ---- which FOVs can be aligned against which ------------------------------
% Two gates: footprint overlap (there must be shared tissue to align on) and dz
% (fields far apart in depth show different neurons -- nothing to match).
ov = zeros(nF); dz = zeros(nF);
for a = 1:nF
    for b = 1:nF
        if a == b, continue; end
        ax = [F(a).motor(1)-F(a).halfW_um, F(a).motor(1)+F(a).halfW_um];
        ay = [F(a).motor(2)-F(a).halfH_um, F(a).motor(2)+F(a).halfH_um];
        bx = [F(b).motor(1)-F(b).halfW_um, F(b).motor(1)+F(b).halfW_um];
        by = [F(b).motor(2)-F(b).halfH_um, F(b).motor(2)+F(b).halfH_um];
        w = max(0, min(ax(2),bx(2)) - max(ax(1),bx(1)));
        h = max(0, min(ay(2),by(2)) - max(ay(1),by(1)));
        aA = diff(ax)*diff(ay); bA = diff(bx)*diff(by);
        ov(a,b) = (w*h) / max(min(aA,bA), eps);
        dz(a,b) = abs(F(a).motor(3) - F(b).motor(3));
    end
end
enoughOv = ov >= cfg.minOverlapFrac;
linkCell = enoughOv & dz <= cfg.maxDz_align_um;                       % cell regime
linkVess = enoughOv & dz >  cfg.maxDz_align_um & dz <= cfg.maxDz_vessel_um;

% The vessel regime is a rescue, not a default: it is offered only for FOVs that
% have no cell-regime partner at all. Otherwise a 20 um pair would compete with a
% 2 um pair and you would be asked to align planes that share no somata.
orphan  = ~any(linkCell,2);
linkVess = linkVess & (orphan | orphan');
link     = linkCell | linkVess;

% Pair QUALITY. Overlap alone is the wrong score here: six FOVs in this dataset
% sit at overlap 1.00, so overlap cannot break the tie and the order fell through
% to FOV index -- which is how a dz=19.4 um pair got presented before a dz=2.2 um
% pair of identical overlap. Penalising dz makes the best-focused pair come first.
qual = ov .* max(0, 1 - dz/cfg.maxDz_vessel_um) .* link;

fprintf('=== roi_stitch_align_260801 ===\n');
fprintf('dataset: %s\n%d FOVs\n', cfg.datasetPath, nF);
fprintf('pairs, cell regime   (ov >= %.2f, |dz| <= %g um): %d\n', ...
    cfg.minOverlapFrac, cfg.maxDz_align_um, nnz(triu(linkCell,1)));
fprintf('pairs, vessel rescue (%g < |dz| <= %g um, orphans only): %d\n\n', ...
    cfg.maxDz_align_um, cfg.maxDz_vessel_um, nnz(triu(linkVess,1)));

nLink = sum(link,2);
if ~any(nLink)
    fprintf(2,'No FOV pair overlaps enough to stitch. Nothing to do here --\n');
    fprintf(2,'run roi_stitch_match_260801 on raw stage coordinates instead.\n');
    return;
end

% ---- resume ---------------------------------------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); placed = false(nF,1); skipped = false(nF,1);
if isfile(cfg.offsetFile)
    O = load(cfg.offsetFile);
    if isfield(O,'fovTable')
        nm = string({F.name});
        for i = 1:nF
            j = find(O.fovTable.fov == nm(i), 1);
            if ~isempty(j)
                tx(i)=O.fovTable.tx_um(j); ty(i)=O.fovTable.ty_um(j);
                placed(i)=O.fovTable.placed(j);
            end
        end
        fprintf('resumed: %d FOVs already placed\n\n', nnz(placed));
    end
end

% ---- order: anchor = highest total pair quality, not merely most links -------
% The anchor is recomputed on every run, including a resumed one, so that a
% resumed session still has a fixed reference. Previously the anchor was only
% chosen when nothing was placed, which left a resumed run with no anchor at all.
[~, anchor] = max(sum(qual,2));
if ~any(placed)
    placed(anchor) = true;
end
fprintf('anchor (never moves): %s   z %.1f\n\n', F(anchor).name, F(anchor).motor(3));

% Work list = every linked FOV except the anchor, whether or not it is already
% placed. Resuming pre-loads its offset and lets you confirm or adjust; it does
% NOT remove it from the queue. The old behaviour queued only unplaced FOVs, so
% once a dataset had been stitched the tool exited immediately with "All linked
% FOVs handled" and there was no way to re-align anything.
% Set SKIP_PLACED = true to get the old fill-in-the-gaps behaviour.
SKIP_PLACED = false;
visited = false(nF,1);
nAlready = nnz(placed) - 1;
if nAlready > 0 && ~SKIP_PLACED
    fprintf('%d FOV(s) already placed will be REVISITED (Enter keeps as-is, s skips)\n\n', nAlready);
end

fig = figure('Name','roi_stitch_align_260801','Color','w', ...
             'NumberTitle','off','WindowState','maximized');
mode = 1;  % 1 = red/green, 2 = checkerboard

while true
    elig = ~visited & ~skipped & (1:nF)' ~= anchor & any(link,2) & ...
           (~SKIP_PLACED | ~placed);
    % a FOV is only workable once something it links to is placed
    cand = find(elig & any(link(:,placed),2));
    if isempty(cand)
        stranded = find(elig);
        if isempty(stranded), fprintf('\nAll linked FOVs handled.\n');
        else
            fprintf('\n%d FOV(s) link only to unplaced FOVs -- nothing to align against:\n', numel(stranded));
            for i = stranded(:)', fprintf('  %s\n', F(i).name); end
        end
        break;
    end
    sc = max(qual(cand, placed), [], 2);        % best available pair, not sum:
    [~,k] = max(sc); mv = cand(k);              %   one good partner is what matters
    nb  = find(link(mv,:) & placed');
    [~,so] = sort(qual(mv,nb),'descend'); nb = nb(so);   % best-focused neighbour first
    ib = 1;

    Am = loadavg(F(mv));
    % An already-placed FOV keeps its existing offset as the starting point --
    % re-seeding from phase correlation would silently discard work you did by
    % hand. Only a never-placed FOV gets the automatic seed.
    if ~placed(mv)
        tx(mv) = 0; ty(mv) = 0;
        [tx(mv), ty(mv)] = phasecorr_seed(F, mv, nb(1), tx, ty, Am);
    end

    % The canvas and the FIXED image are built once per (moving, neighbour) pair,
    % never per keypress. Nudging only rewrites the pixels of the moving FOV, so
    % the axes are untouched and your zoom/pan survives every arrow key.
    V = build_view(fig, F, mv, nb(ib), tx, ty);
    update_view(V, F, mv, nb(ib), tx, ty, Am, mode, STEP_UM, numel(nb), ib);

    done = false; act = ''; rebuild = false;
    while ~done
        bi = nb(ib);
        w = waitforbuttonpress; if w ~= 1, continue; end
        c = get(fig,'CurrentCharacter'); m = get(fig,'CurrentModifier');
        big = any(strcmp(m,'shift')); st = STEP_UM * (1 + 9*big);
        if isempty(c), continue; end
        switch double(c)
            case 28, tx(mv) = tx(mv) - st;
            case 29, tx(mv) = tx(mv) + st;
            case 30, ty(mv) = ty(mv) - st;
            case 31, ty(mv) = ty(mv) + st;
            case 13, act = 'accept'; done = true;
            otherwise
                switch lower(c)
                    case 'r', tx(mv)=0; ty(mv)=0;
                    case 'p', [tx(mv),ty(mv)] = phasecorr_seed(F,mv,bi,tx,ty,Am);
                    case 'b', ib = mod(ib, numel(nb)) + 1; rebuild = true;
                    case 't', mode = 3 - mode;
                    case 's', act = 'skip'; done = true;
                    case 'q', act = 'quit'; done = true;
                    otherwise, continue;   % unmapped key: do not redraw at all
                end
        end
        if done, break; end
        if rebuild                                   % only a neighbour switch
            V = build_view(fig, F, mv, nb(ib), tx, ty); rebuild = false;
        end
        update_view(V, F, mv, nb(ib), tx, ty, Am, mode, STEP_UM, numel(nb), ib);
    end

    visited(mv) = true;
    switch act
        case 'accept'
            placed(mv) = true;
            fprintf('  placed  %-44s  dx %+7.1f  dy %+7.1f um\n', ...
                F(mv).name(1:min(44,end)), tx(mv), ty(mv));
        case 'skip'
            skipped(mv) = true; tx(mv)=0; ty(mv)=0; placed(mv) = false;
            fprintf('  skipped %-44s  (raw stage position kept)\n', F(mv).name(1:min(44,end)));
        case 'quit'
            fprintf('\nquit -- saving progress\n'); break;
    end
end
if isvalid(fig), close(fig); end

% ---- save ------------------------------------------------------------------
fovTable = table(string({F.name})', tx, ty, placed, ...
                 arrayfun(@(f) f.motor(3), F)', ...
                 'VariableNames', {'fov','tx_um','ty_um','placed','z_um'});
save(cfg.offsetFile, 'fovTable', 'cfg');
fprintf('\n%d/%d FOVs placed | max correction %.1f um\n', nnz(placed), nF, max(hypot(tx,ty)));
fprintf('Saved %s\n', cfg.offsetFile);
fprintf('\nNext: roi_stitch_match_260801\n');

un = find(~placed);
if ~isempty(un)
    fprintf(2,'\n%d FOV(s) keep their raw stage position (no overlap, or skipped):\n', numel(un));
    for i = un(:)', fprintf(2,'  %s\n', F(i).name); end
end
end

%% ---------------------------------------------------------------------------
function A = loadavg(f)
if isempty(f.avgPath) || ~isfile(f.avgPath), A = zeros(f.H,f.W); return; end
A = double(imread(f.avgPath));
A = A - prctile(A(:),1); A = A / max(prctile(A(:),99.5), eps);
A = min(max(A,0),1);
% Pixels must go through the same transform as the coordinates, or the picture is
% drawn in a different orientation from the ROI markers laid over it.
A = roi_stitch_orient_260801(A);
end

function [dx, dy] = phasecorr_seed(F, mv, bi, tx, ty, Am)
% FFT phase correlation of the two FOVs placed at their STAGE positions on a
% canvas covering the union of both.
%
% The previous version cropped both images from their top-left corner to a common
% size and correlated those, which throws away the stage offset entirely and
% compares whichever tissue happens to sit at the corner of each field. It
% produced a 205 um "seed" on Vglut2/0224 for a tile whose true correction is
% ~10 um -- and because a seed looks authoritative, it got accepted.
%
% Union-canvas placement is the same arrangement roi_stitch_check_axes_260801
% validates and roi_stitch_map_260801 uses. Still only a seed: confirm every pair.
dx = 0; dy = 0;
try
    Ab = loadavg(F(bi));
    px = max(F(mv).px_um, F(bi).px_um);
    A = imresize(Am, F(mv).px_um/px);  B = imresize(Ab, F(bi).px_um/px);
    ax = F(mv).motor(1)/px;               ay = F(mv).motor(2)/px;
    bx = (F(bi).motor(1)+tx(bi))/px;      by = (F(bi).motor(2)+ty(bi))/px;
    ha = [size(A,2) size(A,1)]/2;  hb = [size(B,2) size(B,1)]/2;
    x0 = min(ax-ha(1), bx-hb(1)); x1 = max(ax+ha(1), bx+hb(1));
    y0 = min(ay-ha(2), by-hb(2)); y1 = max(ay+ha(2), by+hb(2));
    W = round(x1-x0); H = round(y1-y0);
    if W<32 || H<32 || W>4000 || H>4000, return; end
    CA = pasteMax(zeros(H,W), A, round(ay-ha(2)-y0)+1, round(ax-ha(1)-x0)+1);
    CB = pasteMax(zeros(H,W), B, round(by-hb(2)-y0)+1, round(bx-hb(1)-x0)+1);
    if nnz(CA)==0 || nnz(CB)==0, return; end
    CA = CA-mean(CA(:)); CB = CB-mean(CB(:));
    w = hann(H)*hann(W)';
    R = fft2(CA.*w) .* conj(fft2(CB.*w)); R = R./max(abs(R),eps);
    c = fftshift(real(ifft2(R)));
    [~,ix] = max(c(:)); [r0,c0] = ind2sub(size(c), ix);
    dx = -(c0-floor(W/2)-1)*px;  dy = -(r0-floor(H/2)-1)*px;
    if hypot(dx,dy) > 100, dx = 0; dy = 0; end   % implausible -> no seed
catch
    dx = 0; dy = 0;
end
end

function C = pasteMax(C, A, r0, c0)
[h,w] = size(A);
dr = max(r0,1):min(r0+h-1,size(C,1)); dc = max(c0,1):min(c0+w-1,size(C,2));
if isempty(dr)||isempty(dc), return; end
C(dr,dc) = max(C(dr,dc), A(dr-r0+1, dc-c0+1));
end

function V = build_view(fig, F, mv, bi, tx, ty)
% Build the canvas, the FIXED background and every graphics handle ONCE.
% Called only when the pair changes -- never on a nudge.
figure(fig); clf;
cfg = stitch_cfg_260801();
Ab  = loadavg(F(bi));

% Margin so the moving FOV stays inside the canvas as you nudge; without it the
% canvas would have to be resized (and the axes rebuilt) on every keypress.
MARG = 150;
px = min(F(mv).px_um, F(bi).px_um);
x0 = min(F(mv).motor(1)-F(mv).halfW_um, F(bi).motor(1)+tx(bi)-F(bi).halfW_um) - MARG;
x1 = max(F(mv).motor(1)+F(mv).halfW_um, F(bi).motor(1)+tx(bi)+F(bi).halfW_um) + MARG;
y0 = min(F(mv).motor(2)-F(mv).halfH_um, F(bi).motor(2)+ty(bi)-F(bi).halfH_um) - MARG;
y1 = max(F(mv).motor(2)+F(mv).halfH_um, F(bi).motor(2)+ty(bi)+F(bi).halfH_um) + MARG;
px = max(px, max(x1-x0, y1-y0)/2500);      % keep the canvas a sane size
W  = round((x1-x0)/px); H = round((y1-y0)/px);

C1 = place(Ab, F(bi), tx(bi), ty(bi), x0, y0, px, H, W);

ax  = axes('Parent',fig);
him = imshow(zeros(H,W,3),'Parent',ax); axis(ax,'image'); hold(ax,'on');
set(ax,'YDir','normal');   % Y up, as everywhere else (imshow defaults to reverse)

% ROI overlay. Fixed-FOV ROIs never move, so they are drawn once.
fx = (F(bi).roi_x_um + tx(bi) - x0)/px + 1;
fy = (F(bi).roi_y_um + ty(bi) - y0)/px + 1;
plot(ax, fx, fy, 'o', 'MarkerEdgeColor',[0 1 0], 'MarkerSize',9, 'LineWidth',1.2);
hMov  = plot(ax, nan, nan, 'o', 'MarkerEdgeColor',[1 .2 .2], 'MarkerSize',9, 'LineWidth',1.2);
hLink = plot(ax, nan, nan, '-', 'Color',[1 1 0], 'LineWidth',2);
hHit  = plot(ax, nan, nan, 'o', 'MarkerEdgeColor',[1 1 0], 'MarkerSize',16, 'LineWidth',2);
ht    = title(ax,'','Interpreter','none','FontSize',10);

V = struct('ax',ax,'him',him,'hMov',hMov,'hLink',hLink,'hHit',hHit,'ht',ht, ...
           'C1',C1,'x0',x0,'y0',y0,'px',px,'H',H,'W',W, ...
           'fx',fx,'fy',fy,'maxD',cfg.maxMatchDist_um,'dzCell',cfg.maxDz_align_um);
end

function update_view(V, F, mv, bi, tx, ty, Am, mode, step, nnb, ib)
% Nudge path: rewrite CData and a few line handles. The axes object is never
% touched, so zoom, pan and figure layout all survive.
C2 = place(Am, F(mv), tx(mv), ty(mv), V.x0, V.y0, V.px, V.H, V.W);
if mode == 1
    RGB = cat(3, C2, V.C1, zeros(V.H,V.W));
else
    blk = 40; [gx,gy] = meshgrid(1:V.W, 1:V.H);
    m = mod(floor(gx/blk)+floor(gy/blk),2) == 0;
    v = V.C1; v(m) = C2(m); RGB = repmat(v,1,1,3);
end
set(V.him, 'CData', RGB);

% ---- ROI overlay + live match count ---------------------------------------
% The alignment criterion is not just "the vasculature looks right": it is how
% many ROIs actually pair up under the matcher's own distance gate. Nudge until
% this count peaks and the residual bottoms out, and you have aligned the FOVs
% by the same rule the matcher will later apply.
mx = (F(mv).roi_x_um + tx(mv) - V.x0)/V.px + 1;
my = (F(mv).roi_y_um + ty(mv) - V.y0)/V.px + 1;
set(V.hMov, 'XData', mx, 'YData', my);

thr = V.maxD / V.px;                      % match gate in canvas pixels
lx = []; ly = []; hx = []; hy = []; res = [];
for i = 1:numel(mx)
    d = hypot(V.fx - mx(i), V.fy - my(i));
    [dm, j] = min(d);
    if ~isempty(dm) && dm <= thr
        lx = [lx, mx(i), V.fx(j), NaN]; %#ok<AGROW>
        ly = [ly, my(i), V.fy(j), NaN]; %#ok<AGROW>
        hx = [hx, (mx(i)+V.fx(j))/2];   %#ok<AGROW>
        hy = [hy, (my(i)+V.fy(j))/2];   %#ok<AGROW>
        res = [res, dm*V.px];           %#ok<AGROW>
    end
end
set(V.hLink,'XData',lx,'YData',ly);
set(V.hHit, 'XData',hx,'YData',hy);

if isempty(res), ms = 'no ROI pairs within gate';
else, ms = sprintf('%d ROI pairs matched  |  median residual %.1f um', numel(res), median(res));
end

% Say WHICH evidence this pair can be aligned on. A vessel-regime pair showing no
% yellow is expected, not a malfunction, and the operator should be told that
% rather than left to wonder why the cells do not correspond.
adz = abs(F(mv).motor(3) - F(bi).motor(3));
if adz <= V.dzCell
    regime = sprintf('CELL regime (dz %.1f um): somata should correspond -- nudge to maximise yellow', adz);
else
    regime = sprintf(['VESSEL regime (dz %.1f um > %g): planes too far apart for somata to recur. ' ...
                      'Align on VASCULATURE; expect little or no yellow.'], adz, V.dzCell);
end

set(V.ht,'String',sprintf([ ...
    'MOVING (red circles): %s   z=%.0f\n' ...
    'FIXED (green circles): %s   z=%.0f      [neighbour %d/%d]\n' ...
    '%s\n' ...
    'YELLOW = ROI pairing within %.0f um gate  ->  %s\n' ...
    'dx %+.1f  dy %+.1f um   |   step %.1f um (shift = x10)   |   ' ...
    'arrows nudge, Enter accept, s skip, r reset, p re-seed, b next nb, t view, q quit'], ...
    F(mv).name, F(mv).motor(3), F(bi).name, F(bi).motor(3), ib, nnb, ...
    regime, V.maxD, ms, tx(mv), ty(mv), step));
drawnow limitrate;
end

function C = place(A, f, tx, ty, x0, y0, px, H, W)
% Paste one FOV into the shared canvas, clipping whatever falls outside.
% Written as explicit source/destination ranges rather than clever offset
% arithmetic -- an off-by-one here shifts the image under the operator's eye and
% silently corrupts every offset they then record.
C  = zeros(H,W);
As = imresize(A, f.px_um/px);
[h,w] = size(As);
r0 = round((f.motor(2)+ty-f.halfH_um - y0)/px) + 1;
c0 = round((f.motor(1)+tx-f.halfW_um - x0)/px) + 1;

dr = max(r0,1) : min(r0+h-1, H);        % destination rows
dc = max(c0,1) : min(c0+w-1, W);
if isempty(dr) || isempty(dc), return; end
sr = dr - r0 + 1;                        % matching source rows
sc = dc - c0 + 1;
C(dr,dc) = As(sr,sc);
end
