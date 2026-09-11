function Ventral_surface_manual_align_260731()
%% Ventral_surface_manual_align_260731  Nudge overlapping FOVs into register by hand.
% -----------------------------------------------------------------------
% Suggests every pair of FOVs that actually OVERLAP, shows them as a two-colour
% overlay, and lets you slide one with the arrow keys until the anatomy lines up.
%
% WHY PAIRS THAT OVERLAP, not pairs that are merely close. The site merge uses
% centre distance (mergeSiteUm), which says nothing about whether the two fields
% share any tissue -- and you can only align by eye what visibly overlaps. Here a
% pair is offered only if its footprints intersect by MIN_OVERLAP of the smaller
% one.
%
% ORDER. Within a session the FOV with the most ROIs is the anchor and never
% moves. After that the pair with the largest overlap involving one already-placed
% FOV is offered next, so every FOV is aligned exactly once against something
% already fixed -- a spanning walk, the same way manual mosaicking works. Aligning
% arbitrary pairs instead lets later pairs silently undo earlier ones.
%
% SEEDED FROM THE AUTOMATIC FIT. If ROI_on_cartoon_registered.mat exists its
% tie-point offsets are the starting position, so this is a correction pass over
% the automatic result rather than starting from stage coordinates. FOVs with no
% tie-point start where the stage put them.
%
% COLOUR. Fixed FOV in green, moving FOV in magenta. Misalignment reads as colour
% fringing; alignment goes grey/white. That is far easier to judge than two
% side-by-side panels.
%
% KEYS
%   arrows        nudge  STEP_UM        shift+arrows  nudge  10x
%   ENTER / y     accept this FOV's position and move on
%   n / space     skip (keep the seeded position)
%   r             reset this FOV to its seed
%   b             back one FOV
%   q             save and quit
%
% Output: <Desktop>\ROI_on_cartoon_manual_offsets.mat  -- same schema as the
% tie-point fovTable, so the map's USE_REGISTERED path can read either.

%% ----------------------------- USER PARAMETERS -----------------------------
DATA_MAT    = 'C:\Users\Admin\Desktop\ROI_on_cartoon_data.mat';       % stage-coord run
SEED_MAT    = 'C:\Users\Admin\Desktop\ROI_on_cartoon_registered.mat'; % '' = no seed
OUT_MAT     = 'C:\Users\Admin\Desktop\ROI_on_cartoon_manual_offsets.mat';
MIN_OVERLAP = 0.10;    % pair offered only if footprints share this fraction of the smaller
STEP_UM     = 0.5;     % arrow key step (um); shift = 10x
CLIP_PCT    = [2 99.5];
GAMMA       = 0.55;
%% ---------------------------------------------------------------------------

assert(isfile(DATA_MAT), 'Not found: %s', DATA_MAT);
S = load(DATA_MAT); F = S.F; R = S.R;   % R carries the per-ROI cell_key
nF = numel(F);
fovKey = string({F.sess}) + "/" + string({F.fov});

% ---- seed from the automatic tie-point fit ---------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); seeded = false(nF,1);
if ~isempty(SEED_MAT) && isfile(SEED_MAT)
    Q = load(SEED_MAT,'fovTable'); T = Q.fovTable;
    for i = 1:nF
        j = find(T.fov == fovKey(i), 1);
        if ~isempty(j) && T.linked(j)
            tx(i) = T.tx_um(j); ty(i) = T.ty_um(j); seeded(i) = true;
        end
    end
    fprintf('seeded %d FOVs from the tie-point fit\n', sum(seeded));
end
tx0 = tx; ty0 = ty;

% ---- load and pre-render each FOV image ------------------------------------
IM = cell(nF,1); EXT = nan(nF,4);   % [halfW halfH] in um, plus size
for i = 1:nF
    if isempty(F(i).avgPath) || ~isfile(F(i).avgPath), continue; end
    A = double(imread(F(i).avgPath));
    lo = prctile(A(:),CLIP_PCT(1)); hi = prctile(A(:),CLIP_PCT(2));
    aN = max(0,min(1,(A-lo)/max(hi-lo,eps))).^GAMMA;
    if F(i).rot90, aN = rot90(aN,-1); end
    IM{i} = aN;
    [hp,wp] = size(aN);
    EXT(i,:) = [(wp/2)*F(i).px_um_fov, (hp/2)*F(i).px_um_fov, hp, wp];
end
has = ~cellfun(@isempty, IM);
fprintf('%d of %d FOVs have an average projection to align\n', sum(has), nF);

% ---- candidate pairs: footprints must actually intersect --------------------
sess = string({F.sess});
ov = zeros(nF);
for a = 1:nF
    if ~has(a), continue; end
    for b = a+1:nF
        if ~has(b) || sess(a) ~= sess(b), continue; end
        ovx = min(F(a).x_um+EXT(a,1), F(b).x_um+EXT(b,1)) - max(F(a).x_um-EXT(a,1), F(b).x_um-EXT(b,1));
        ovy = min(F(a).y_um+EXT(a,2), F(b).y_um+EXT(b,2)) - max(F(a).y_um-EXT(a,2), F(b).y_um-EXT(b,2));
        if ovx <= 0 || ovy <= 0, continue; end
        aSm = min(4*EXT(a,1)*EXT(a,2), 4*EXT(b,1)*EXT(b,2));
        f = (ovx*ovy)/aSm;
        if f >= MIN_OVERLAP, ov(a,b) = f; ov(b,a) = f; end
    end
end
fprintf('overlapping pairs offered: %d\n', nnz(triu(ov)));

% ---- how many MATCHED CELLS each pair shares -------------------------------
% This is NOT a selection criterion -- measured on this dataset, every pair that
% shares a cell already passes the footprint test, so adding it would offer
% nothing new. It is evidence ABOUT a pair: a shared cell is proof the two fields
% contain the same tissue, whereas footprint overlap is only the coordinates
% claiming so. 204 of 228 pairs here have no such proof, and on those the manual
% alignment is judged purely on image content with nothing to check it against.
% So it is displayed, and used to order the walk: pairs with shared cells are
% aligned first, because they make the most trustworthy anchors for the rest.
nShared = zeros(nF);
rSess = string({R.group}) + "/" + string({R.date});
rCellK = string({R.cell_key});
fovOfR = nan(numel(R),1);
for i = 1:numel(R)
    j = find(fovKey == rSess(i) + "/" + string(R(i).fov), 1);
    if ~isempty(j), fovOfR(i) = j; end
end
% NOTE: locals here are named so they cannot collide with the nested functions
% below -- nested functions SHARE this workspace, so a plain `m`, `f` or `c` here
% would be clobbered by draw()/onKey() and vice versa.
okR = ~isnan(fovOfR) & rCellK(:) ~= "" & ~ismissing(rCellK(:));
% Explicit index loop: `for ck = unique(...)'` binds the WHOLE column at once,
% because a MATLAB for-loop iterates over columns, not elements.
uCK = unique(rCellK(okR));
for ic = 1:numel(uCK)
    ck = uCK(ic);
    fv = unique(fovOfR(okR & rCellK(:) == ck));
    for ja = 1:numel(fv)
        for jb = ja+1:numel(fv)
            nShared(fv(ja),fv(jb)) = nShared(fv(ja),fv(jb)) + 1;
            nShared(fv(jb),fv(ja)) = nShared(fv(jb),fv(ja)) + 1;
        end
    end
end
% A SHARED CELL IS PROOF OF OVERLAP. The same neuron cannot be in two fields that
% do not intersect, so a cell-linked pair is offered whether or not the recorded
% coordinates agree -- coordinates are the thing being corrected here, and using
% them to gate the candidates would let a bad coordinate hide the very pair that
% would fix it. Footprint overlap stays as a SECOND, weaker criterion for pairs
% with no matched cell.
link = (ov > 0) | (nShared > 0);
fprintf('  of those, %d share at least one matched cell (proof of overlap)\n', ...
    nnz(triu(nShared > 0 & ov > 0)));
fprintf('  cell-linked but NOT footprint-overlapping (added by the union): %d\n', ...
    nnz(triu(nShared > 0 & ov == 0)));
fprintf('  total candidate pairs: %d\n\n', nnz(triu(link)));

% ---- build the spanning walk per session ------------------------------------
order = [];  fixedOf = [];
for s = unique(sess)
    idx = find(sess == s & has');
    if numel(idx) < 2, continue; end
    [~,k] = max([F(idx).nroi]);  placed = idx(k);          % anchor = most ROIs
    remaining = setdiff(idx, placed);
    while ~isempty(remaining)
        % Score = overlap fraction, plus a large bonus when the pair shares matched
        % cells. Cell-linked pairs are therefore placed FIRST and become the anchors
        % everything else is aligned against -- they are the only pairs whose
        % overlap is certain rather than inferred from the coordinates being fixed.
        best = 0; bi = []; bf = [];
        for mm = remaining
            sc = ov(mm, placed) + 10*(nShared(mm, placed) > 0);
            [v,p] = max(sc);
            if link(mm, placed(p)) && v > best, best = v; bi = mm; bf = placed(p); end
        end
        if best == 0, break; end                            % rest do not overlap
        order(end+1) = bi;  fixedOf(end+1) = bf; %#ok<AGROW>
        placed(end+1) = bi; %#ok<AGROW>
        remaining = setdiff(remaining, bi);
    end
end
nStep = numel(order);
assert(nStep > 0, 'No overlapping FOV pairs found -- nothing to align.');
fprintf('%d FOVs to align, across %d session(s)\n', nStep, numel(unique(sess(order))));

% ---- GUI --------------------------------------------------------------------
st = struct('k',1,'quit',false);
fig = figure('Color','w','Name','manual FOV alignment','NumberTitle','off', ...
             'Units','normalized','Position',[0.06 0.06 0.86 0.86]);
ax = axes(fig);
set(fig,'KeyPressFcn',@onKey);
draw();
uiwait(fig);

% ---- save -------------------------------------------------------------------
shift = hypot(tx,ty);
moved = hypot(tx-tx0, ty-ty0) > 1e-9;
fovTable = table(fovKey(:), tx, ty, shift, moved | seeded, ones(nF,1), ...
    'VariableNames', {'fov','tx_um','ty_um','shift_um','linked','component'});
save(OUT_MAT,'fovTable','tx','ty','tx0','ty0','order','fixedOf','MIN_OVERLAP');
writetable(fovTable, strrep(OUT_MAT,'.mat','.csv'));
fprintf('\n%d FOVs adjusted by hand (median %.0f um from seed)\n', ...
    sum(moved), median(hypot(tx(moved)-tx0(moved), ty(moved)-ty0(moved))));
fprintf('Saved %s\n  + .csv\n', OUT_MAT);
fprintf('To draw it: point USE_REGISTERED at this file in the cartoon script.\n');

%% ============================ nested ============================
function draw()
    if st.k > nStep, return; end
    m = order(st.k); f = fixedOf(st.k);
    cla(ax);
    % union extent, padded
    xs = [F(f).x_um+tx(f)-EXT(f,1), F(f).x_um+tx(f)+EXT(f,1), ...
          F(m).x_um+tx(m)-EXT(m,1), F(m).x_um+tx(m)+EXT(m,1)];
    ys = [F(f).y_um+ty(f)-EXT(f,2), F(f).y_um+ty(f)+EXT(f,2), ...
          F(m).y_um+ty(m)-EXT(m,2), F(m).y_um+ty(m)+EXT(m,2)];
    pad = 0.08*max(range(xs), range(ys));
    hold(ax,'on');
    show(f, [0 1 0], tx(f), ty(f));      % fixed  = green
    show(m, [1 0 1], tx(m), ty(m));      % moving = magenta
    set(ax,'YDir','normal','Color','k');
    axis(ax,'image'); xlim(ax,[min(xs)-pad max(xs)+pad]); ylim(ax,[min(ys)-pad max(ys)+pad]);
    xlabel(ax,'X (\mum) lateral'); ylabel(ax,'Y (\mum) rostral +');
    d = hypot(tx(m)-tx0(m), ty(m)-ty0(m));
    if nShared(m,f) > 0
        ev = sprintf('%d SHARED CELL(S) -- overlap is proven', nShared(m,f));
    else
        ev = 'NO shared cell -- overlap claimed by coordinates only';
    end
    title(ax, sprintf(['[%d/%d]  %s      %s\nMOVING (magenta): %s     FIXED (green): %s\n' ...
        'offset from seed  %+.1f, %+.1f \\mum  (|%.1f|)   footprint overlap %.0f%%   step %g \\mum\n' ...
        'arrows nudge  |  shift+arrows 10x  |  ENTER accept  |  n skip  |  r reset  |  b back  |  q save+quit'], ...
        st.k, nStep, F(m).sess, ev, ...
        shorten(F(m).fov), shorten(F(f).fov), ...
        tx(m)-tx0(m), ty(m)-ty0(m), d, 100*ov(m,f), STEP_UM), ...
        'FontSize',10,'Interpreter','tex');
end

function show(i, tint, ox, oy)
    aN = IM{i};
    rgb = zeros([size(aN) 3]);
    for c = 1:3, rgb(:,:,c) = tint(c); end
    xe = [F(i).x_um+ox-EXT(i,1), F(i).x_um+ox+EXT(i,1)];
    ye = [F(i).y_um+oy+EXT(i,2), F(i).y_um+oy-EXT(i,2)];
    image(ax, xe, ye, rgb, 'AlphaData', aN);
end

function onKey(~, e)
    if st.k > nStep, return; end
    m = order(st.k);
    s = STEP_UM; if any(strcmp(e.Modifier,'shift')), s = STEP_UM*10; end
    switch e.Key
        case 'leftarrow',  tx(m) = tx(m) - s;
        case 'rightarrow', tx(m) = tx(m) + s;
        case 'uparrow',    ty(m) = ty(m) + s;
        case 'downarrow',  ty(m) = ty(m) - s;
        case 'r',          tx(m) = tx0(m); ty(m) = ty0(m);
        case {'return','y'}, st.k = st.k + 1;
        case {'n','space'},  st.k = st.k + 1;
        case 'b',          st.k = max(1, st.k - 1);
        case 'q',          st.quit = true;
    end
    if st.quit || st.k > nStep
        if isvalid(fig), uiresume(fig); close(fig); end
        return
    end
    draw();
end
end

function s = shorten(s)
if numel(s) > 38, s = [s(1:35) '...']; end
s = strrep(s, '_', '\_');
end
