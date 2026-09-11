%% roi_stitch_check_axes_260801  Which image->stage axis mapping is correct?
% -----------------------------------------------------------------------
% RUN THIS FIRST on any new rig, mounting, or acquisition-coordinate change.
%
% Placing a FOV requires knowing where an image column and an image row go in
% STAGE space. There are EIGHT possibilities (the symmetries of a square): the
% four sign combinations of (col->x, row->y), and those four again TRANSPOSED
% (col->y, row->x) for a 90 degree rotation. Testing only the flips cannot detect
% a rotated frame -- it just returns the best of a wrong set, which is exactly
% how an earlier 4-way version of this test gave a confident wrong answer.
%
% NOTE ON SCOPE. This asks a question about the RIG (how scan mirrors map to
% stage axes), not about the animal. How the prep was mounted is a separate,
% later, stage->ANATOMY transform -- that is what Ventral_surface_ROI_on_cartoon
% handles with AXMAP_ROT90 = [0 1; -1 0] (lateral = +stage y, rostral = -stage x)
% for sessions from 260721 onward. A global rotation cannot change which ROIs
% match, because distances are invariant under it. The mapping tested HERE can,
% because it changes where ROIs sit inside their own FOV.
%
% Test: motorPosition is held fixed and authoritative; only the within-FOV axis
% mapping varies. Overlapping FOVs at similar depth are placed and
% phase-correlated. The correct mapping puts overlapping tissue on top of itself,
% so it shows the smallest residual.
%
% Reading the result: the winner should be clearly separated, and its best pairs
% should reach a few um (the true registration residual). If EVERY mapping is
% large, the stage coordinates are unreliable, not the axis mapping.

cfg = stitch_cfg_260801();
S = load(fullfile(cfg.outDir,'fov_table.mat')); F = S.F; nF = numel(F);

% Work entirely in RAW stage space: motorRaw is the untouched ScanImage value,
% and the half-extents are recomputed from W/H/px rather than read from the
% struct, because the stored ones may already carry the prep rotation (which
% transposes them). Otherwise this test would be validating the very correction
% it is supposed to be measuring.
for i = 1:nF
    if isfield(F,'motorRaw') && ~isempty(F(i).motorRaw), F(i).motor = F(i).motorRaw; end
    F(i).halfW_um = (F(i).W/2) * F(i).px_um;
    F(i).halfH_um = (F(i).H/2) * F(i).px_um;
end

LBL = {'col=+x row=+y  identity', 'col=-x row=+y  fliplr', ...
       'col=+x row=-y  flipud',   'col=-x row=-y  rot180', ...
       'col=+y row=+x  transpose','col=+y row=-x  rot90cw', ...
       'col=-y row=+x  rot90ccw', 'col=-y row=-x  transpose+180'};

P = [];
for a=1:nF, for b=a+1:nF
    w = max(0, min(F(a).motor(1)+F(a).halfW_um,F(b).motor(1)+F(b).halfW_um) - ...
               max(F(a).motor(1)-F(a).halfW_um,F(b).motor(1)-F(b).halfW_um));
    h = max(0, min(F(a).motor(2)+F(a).halfH_um,F(b).motor(2)+F(b).halfH_um) - ...
               max(F(a).motor(2)-F(a).halfH_um,F(b).motor(2)-F(b).halfH_um));
    ovf = (w*h)/max(min(4*F(a).halfW_um*F(a).halfH_um, 4*F(b).halfW_um*F(b).halfH_um),eps);
    if ovf >= 0.5 && abs(F(a).motor(3)-F(b).motor(3)) <= 15, P=[P; a b]; end %#ok<AGROW>
end, end
fprintf('=== roi_stitch_check_axes_260801 ===\n%s\n%d test pairs\n\n', cfg.datasetPath, size(P,1));

R = nan(size(P,1), 8);
for k = 1:size(P,1)
    for v = 1:8, R(k,v) = residual(F, P(k,1), P(k,2), v); end
end

fprintf('%-36s %-36s', 'FOV A', 'FOV B'); fprintf('%8d', 1:8); fprintf('\n');
for k=1:size(P,1)
    fprintf('%-36s %-36s', F(P(k,1)).name(1:min(36,end)), F(P(k,2)).name(1:min(36,end)));
    fprintf('%8.1f', R(k,:)); fprintf('\n');
end

med = median(R,1,'omitnan');
fprintf('\n%-73s','MEDIAN residual (um)'); fprintf('%8.1f', med); fprintf('\n');
fprintf('%-73s','MIN    residual (um)'); fprintf('%8.1f', min(R,[],1,'omitnan')); fprintf('\n');
fprintf('%-73s','frac pairs < 10 um  '); fprintf('%8.2f', mean(R<10,1,'omitnan')); fprintf('\n');

[~,o] = sort(med);
fprintf('\nRanking:\n');
for i = 1:8
    fprintf('  %d. %-30s median %6.1f um   %3.0f%% of pairs under 10 um\n', ...
        i, LBL{o(i)}, med(o(i)), 100*mean(R(:,o(i))<10,'omitnan'));
end
fprintf('\nWINNER: %s\n', LBL{o(1)});
if med(o(2)) < 2*med(o(1))
    fprintf(['CAUTION: runner-up within 2x of the winner -- not a clean separation.\n' ...
             'Add more overlapping pairs before trusting this.\n']);
end

%% ---------------------------------------------------------------------------
function d = residual(F,a,b,v)
% Composite both FOVs onto one canvas at their motor positions under mapping v,
% then phase-correlate. Compositing rather than cropping: after resampling to a
% common pixel size a high-zoom FOV can be only ~240 px, too small for any fixed
% crop window (an earlier version returned all-NaN for exactly this reason).
d = NaN;
try
    A = axmap(norm01(double(imread(F(a).avgPath))), v);
    B = axmap(norm01(double(imread(F(b).avgPath))), v);
    px = max(F(a).px_um, F(b).px_um);
    A = imresize(A, F(a).px_um/px);  B = imresize(B, F(b).px_um/px);
    ax=F(a).motor(1)/px; ay=F(a).motor(2)/px; bx=F(b).motor(1)/px; by=F(b).motor(2)/px;
    ha=[size(A,2) size(A,1)]/2; hb=[size(B,2) size(B,1)]/2;
    x0=min(ax-ha(1),bx-hb(1)); x1=max(ax+ha(1),bx+hb(1));
    y0=min(ay-ha(2),by-hb(2)); y1=max(ay+ha(2),by+hb(2));
    W=round(x1-x0); H=round(y1-y0);
    if W<32||H<32||W>4000||H>4000, return; end
    CA = paste(zeros(H,W), A, round(ay-ha(2)-y0), round(ax-ha(1)-x0));
    CB = paste(zeros(H,W), B, round(by-hb(2)-y0), round(bx-hb(1)-x0));
    if nnz(CA)==0||nnz(CB)==0, return; end
    CA=CA-mean(CA(:)); CB=CB-mean(CB(:));
    w = hann(H)*hann(W)';
    Rr = fft2(CA.*w).*conj(fft2(CB.*w)); Rr = Rr./max(abs(Rr),eps);
    c = fftshift(real(ifft2(Rr)));
    [~,ix]=max(c(:)); [r0,c0]=ind2sub(size(c),ix);
    d = hypot((c0-floor(W/2)-1)*px, (r0-floor(H/2)-1)*px);
catch
end
end

function A = axmap(A, v)
% The eight symmetries of a square, indexed to match LBL above.
switch v
    case 1  % identity
    case 2, A = fliplr(A);
    case 3, A = flipud(A);
    case 4, A = rot90(A,2);
    case 5, A = A.';
    case 6, A = fliplr(A.');
    case 7, A = flipud(A.');
    case 8, A = rot90(A.',2);
end
end
function A = norm01(A), A=A-prctile(A(:),1); A=A/max(prctile(A(:),99.5),eps); A=min(max(A,0),1); end
function C = paste(C, A, r0, c0)
[h,w]=size(A); r0=r0+1; c0=c0+1;
dr=max(r0,1):min(r0+h-1,size(C,1)); dc=max(c0,1):min(c0+w-1,size(C,2));
if isempty(dr)||isempty(dc), return; end
C(dr,dc) = A(dr-r0+1, dc-c0+1);
end
