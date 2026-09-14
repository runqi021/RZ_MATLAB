function varargout = dm_axial_estimator_260913(mode, varargin)
%DM_AXIAL_ESTIMATOR_260913  Estimate axial displacement of a frame against a reference volume.
%
%   E   = dm_axial_estimator_260913('prep', refVol, zAxis, opts)
%   out = dm_axial_estimator_260913('run',  E, frames)
%   z   = dm_axial_estimator_260913('coarse', E, frame)     % argmax tier alone
%
% ---------------------------------------------------------------------------
% WHY NOT ARGMAX OVER SLICES
% The obvious estimator -- correlate the frame against all K reference planes
% and take the best -- is the wrong workhorse. It quantises to the plane step,
% the correlation-vs-z peak is broad and shallow (median peak r was only
% 0.35-0.53 in dm_zstack_demo_260911), and it costs K full-frame correlations
% per frame. It is kept here, but only as a COARSE tier for lock-on and
% recovery.
%
% THE WORKHORSE IS A GRADIENT SOLVE (the axial Lucas-Kanade)
% A closed loop keeps you NEAR z = 0, so linearise about the nearest reference
% plane R0 and solve for the displacement directly:
%
%   I  ~=  (1+g)*R0  +  dz*dR/dz  +  dx*dR/dx  +  dy*dR/dy  +  b
%
% Five unknowns, one weighted least-squares solve, sub-millisecond, sub-micron.
% Because the weights come from the REFERENCE (expected photon noise ~ expected
% intensity) and not from the frame, the whole normal-equation inverse is
% precomputed per plane in 'prep'. Running a frame is then ONE 5-by-Npix
% matrix-vector product.
%
% THE THREE TERMS THAT LOOK OPTIONAL AND ARE NOT
%  * dx, dy -- included even though nothing corrects them in real time. Lateral
%    motion projects onto dR/dz wherever structure is oblique, so leaving them
%    out lets lateral motion masquerade as defocus. Estimating them is how you
%    keep them OUT of dz.
%  * g (gain) -- defocus dims and blurs, so ANY laser/PMT/bleaching drift reads
%    as defocus. This is the single biggest confound for axial estimation and
%    it costs one column to absorb.
%  * b (offset) -- PMT baseline. ScanImage digitises signed with dark near 0.
%
% STRIPS, NOT FRAMES
% A resonant frame takes ~33 ms, which is about a TENTH of a breath cycle. The
% frame is not a snapshot -- it is z-SHEARED top to bottom. Solving per strip
% (default 8 strips of 64 lines) removes that bias and returns estimates at
% ~240 Hz instead of 30, which is what makes a phase model worth fitting.
%
% SIGN CONVENTION -- READ THIS
% dz is the displacement OF THE SAMPLE relative to the reference, in um, on the
% same axis as zAxis. A POSITIVE dz means the frame matches a reference plane
% at HIGHER z than the working plane. The CORRECTION to send is -dz. Getting
% this backwards in a closed loop makes the loop drive the motion instead of
% cancelling it, and it runs away. The lateral helper phase_shift() below
% returns the CORRECTION, not the displacement, and is negated where it is used.
%
% No toolboxes required: no Image Processing, no Signal Processing, no Stats.
%
% Runqi Zhang / 2026-09-13.

switch lower(mode)
    case 'prep',   varargout{1} = prep(varargin{:});
    case 'run',    varargout{1} = run_frames(varargin{:});
    case 'coarse', [varargout{1}, varargout{2}] = coarse_z(varargin{:});
    case 'shift',  varargout{1} = phase_shift(varargin{:});   % exposed so the
                   % analysis shares ONE implementation of the sign convention
                   % rather than keeping a second copy that can drift from it.
    otherwise, error('dm_axial_estimator:mode','unknown mode "%s"', mode);
end
end

%% ============================== PREP ==============================
function E = prep(refVol, zAxis, opts)
%PREP  Precompute everything that does not depend on the incoming frame.
if nargin < 3, opts = struct(); end
E.opts = fill_defaults(opts);
o = E.opts;

refVol = single(refVol);
[ny, nx, K] = size(refVol);
zAxis = zAxis(:).';
assert(numel(zAxis) == K, 'zAxis has %d entries for %d planes', numel(zAxis), K);
assert(all(diff(zAxis) > 0), 'zAxis must be strictly increasing');
assert(K >= 5, 'need at least 5 reference planes, got %d', K);

% Light spatial smoothing WIDENS the linear capture range and suppresses shot
% noise in the gradients. It costs lateral precision, which we do not need --
% x/y are handled post hoc.
if o.smoothPx > 0
    for k = 1:K
        refVol(:,:,k) = blur(refVol(:,:,k), o.smoothPx);
    end
end

% --- axial gradient, central differences on the TRUE z axis -------------
% Non-uniform-safe, so a gated reference with a few dropped planes still works.
Rz = zeros(ny,nx,K,'single');
for k = 1:K
    if k == 1
        Rz(:,:,k) = (refVol(:,:,2) - refVol(:,:,1)) / (zAxis(2)-zAxis(1));
    elseif k == K
        Rz(:,:,k) = (refVol(:,:,K) - refVol(:,:,K-1)) / (zAxis(K)-zAxis(K-1));
    else
        Rz(:,:,k) = (refVol(:,:,k+1) - refVol(:,:,k-1)) / (zAxis(k+1)-zAxis(k-1));
    end
end

% --- strips -------------------------------------------------------------
S = o.nStrips;
edges = round(linspace(0, ny, S+1));
rows  = cell(1,S);
for s = 1:S, rows{s} = (edges(s)+1) : edges(s+1); end

% --- per plane, per strip: pixel subset, design matrix, solved inverse ---
% Only high-|dR/dz| pixels carry axial information. The rest add noise and
% cost. Keeping the top few percent is both faster AND better conditioned.
M    = cell(K,S);   % 5-by-Npix solved operator
idx  = cell(K,S);   % linear indices INTO THE STRIP
colN = cell(K,S);   % column norms used to un-normalise the solution
Askp = cell(K,S);   % scaled design matrix, kept only for the robust refit
wkp  = cell(K,S);
for k = 1:K
    R0k = refVol(:,:,k);
    Rzk = Rz(:,:,k);
    [Rxk, Ryk] = spatial_grad(R0k);
    wFloor = qtile(R0k(:), 0.10);
    for s = 1:S
        r  = rows{s};
        a0 = R0k(r,:);  az = Rzk(r,:);  ax = Rxk(r,:);  ay = Ryk(r,:);

        g   = abs(az(:));
        thr = qtile(g, 1 - o.pixFrac);
        sel = find(g >= thr & isfinite(g));
        if numel(sel) < o.minPix
            [~, ord] = sort(g, 'descend');
            sel = ord(1:min(o.minPix, numel(ord)));
        end
        idx{k,s} = sel;

        A = double([az(sel), ax(sel), ay(sel), a0(sel), ones(numel(sel),1)]);
        % Photon noise: variance ~ intensity, so weight ~ 1/intensity. Use the
        % REFERENCE intensity, not the frame's, so this stays precomputable.
        w = 1 ./ max(double(a0(sel)), double(wFloor));
        w(~isfinite(w) | w <= 0) = 0;

        % Column scaling for conditioning -- dR/dz and R0 differ by orders of
        % magnitude and an unscaled normal-equation solve loses digits.
        cn = sqrt(sum((A.^2) .* w, 1)).';
        cn(cn == 0 | ~isfinite(cn)) = 1;
        As = A ./ cn.';

        AtW = As.' .* w.';
        H   = AtW * As;
        H   = H + o.ridge * trace(H)/5 * eye(5);   % tiny ridge: a flat strip is singular
        M{k,s}    = H \ AtW;
        colN{k,s} = cn;
        if o.robust
            % Kept so a trimmed subset can be re-solved on the fly. Single
            % precision: this is ~14 MB at default settings, and the solve is
            % a 5x5 either way.
            Askp{k,s} = single(As); %#ok<AGROW>
            wkp{k,s}  = single(w);  %#ok<AGROW>
        end
    end
end

% --- coarse tier: decimated, mean-removed, unit-norm planes -------------
dec = o.coarseDecim;
C   = reshape(refVol(1:dec:end, 1:dec:end, :), [], K);
C   = double(C);
C   = C - mean(C,1);
nrm = sqrt(sum(C.^2,1));  nrm(nrm == 0) = 1;
E.coarseRef  = C ./ nrm;
E.coarseSize = size(refVol(1:dec:end, 1:dec:end, 1));

% --- lateral pre-registration reference (middle plane) -------------------
%
% SPEED. This is the loop's bottleneck, not the least squares. Done naively it
% is three 512x512 FFTs per frame and costs more than everything else combined.
% Two fixes, both free:
%   1. the REFERENCE transform is fixed -- compute it once here, not per frame;
%   2. run the correlation DECIMATED. It only has to find the integer shift;
%      the least squares cleans up the remainder, so 2 px precision is plenty
%      and a 2x decimation is 4x less FFT.
E.kMid   = round(K/2);
E.regRef = double(refVol(:,:,E.kMid));
E.win    = hann2(ny, nx);

dr = max(1, round(o.regDecim));
rr = double(refVol(1:dr:end, 1:dr:end, E.kMid));
E.regDecim = dr;
E.regSize  = size(rr);
E.regWin   = hann2(E.regSize(1), E.regSize(2));
E.regRefF  = fft2(rr .* E.regWin);

E.ref   = refVol;  E.Rz = Rz;  E.zAxis = zAxis;
E.rows  = rows;    E.M = M;    E.idx = idx;  E.colN = colN;
E.A     = Askp;    E.w  = wkp;
E.size  = [ny nx]; E.K = K;    E.S = S;
E.dz    = mean(diff(zAxis));
end

%% ============================== ROBUST REFIT ==============================
function p = robust_refit(E, k, s, dvec, p0)
%ROBUST_REFIT  Drop the worst-fitting pixels and solve again.
%
% WHY THIS EXISTS: dynamic fluorescence.
% The model says the frame differs from the reference only by displacement,
% gain and offset. Real preparations break that. A calcium transient is a
% LOCAL brightness change with the shape of one cell -- the global gain term
% does not absorb it, and for any cell sitting off-focus its profile resembles
% dR/dz, so the transient leaks straight into dz. Red blood cells sweeping
% through a vessel do the same thing with the opposite sign.
%
% Both are OUTLIERS: a small fraction of pixels carrying a large residual,
% against a majority that fit well. So fit, throw away the worst trimFrac of
% pixels, and fit again. This is one extra 5x5 solve and it removes most of
% the leak, because activity is sparse in space while displacement is not.
%
% It CANNOT save you if the active pixels are the majority, or if activity is
% spatially locked to the structure you are registering on. That is why the
% guide star should be a structural channel in the first place.
o = E.opts;
As = double(E.A{k,s});
w  = double(E.w{k,s});
r  = As*(p0 .* E.colN{k,s}) - dvec;
thr = qtile(abs(r), 1 - o.trimFrac);
keep = abs(r) <= thr;
if nnz(keep) < 50, p = p0; return; end
Ak = As(keep,:);  wk = w(keep);
AtW = Ak.' .* wk.';
H = AtW * Ak;
H = H + o.ridge * trace(H)/5 * eye(5);
p = (H \ (AtW * dvec(keep))) ./ E.colN{k,s};
if ~all(isfinite(p)), p = p0; end
end

%% ============================== RUN ==============================
function out = run_frames(E, frames, zInit)
%RUN_FRAMES  Estimate dz (and dx, dy, gain) for every strip of every frame.
%
%   out = run_frames(E, frames)         start from the coarse tier
%   out = run_frames(E, frames, zInit)  continue from a previous chunk
%
% A 2 min series at 30 Hz is ~3.8 GB of frames, so the caller reads it in
% chunks. Pass out.zEnd from the previous chunk as zInit or every chunk
% boundary re-runs the coarse tier and can jump.
if nargin < 3, zInit = NaN; end
o = E.opts;
frames = single(frames);
N = size(frames,3);
S = E.S;

out.dz    = nan(S,N);   % um, displacement of the sample (correction is -dz)
out.dx    = nan(S,N);   % px, residual lateral after integer pre-registration
out.dy    = nan(S,N);
out.dxInt = nan(1,N);   % px, integer pre-registration actually applied
out.dyInt = nan(1,N);
out.gain  = nan(S,N);
out.resid = nan(S,N);
out.kUsed = nan(S,N);
out.coarseUsed = false(1,N);
out.stripTime  = ((1:S)' - 0.5) / S;   % fraction of a frame period

zEst = zInit;                            % running estimate, um (NaN = re-lock)
sinceLock = 0;                           % frames since the last coarse check
for n = 1:N
    I = frames(:,:,n);
    if o.smoothPx > 0, I = blur(I, o.smoothPx); end

    % --- integer lateral pre-registration --------------------------------
    % The linear dx/dy terms only cover a pixel or two. Breathing moves the
    % sample further than that laterally, so take the bulk out first and let
    % the least squares clean up the remainder.
    d = phase_shift_cached(E, I);                  % returns the CORRECTION
    if max(abs(d)) > o.maxShiftPx, d = [0 0]; end  % refuse a runaway match
    Ish = circshift(I, d);
    out.dyInt(n) = -d(1);                          % displacement = -correction
    out.dxInt(n) = -d(2);

    % --- coarse tier ------------------------------------------------------
    % NOT on n == 1: a real-time caller runs one frame per call, so n == 1 is
    % EVERY call, and firing coarse every frame quantises dz to the reference
    % plane spacing (seen on the rig 2026-09-13: +1.00 / +0.00 / +2.99 um
    % against a 1 um stack). zInit = NaN is how a caller asks for a re-lock.
    %
    % *** BUT THE OLD RE-LOCK TEST COULD NEVER FIRE. ***
    % It asked whether zEst had drifted more than trustUm from the nearest
    % PLANE -- and zEst is always within half a plane spacing of one, by
    % construction. So if the fine solve ever returned ~0, zEst froze at the
    % first coarse value and the estimator reported a confident CONSTANT
    % forever. That is exactly what happened on 2026-09-13: dz never moved,
    % every lock-in amplitude came out ~0, and nothing flagged it.
    %
    % So coarse is now also run PERIODICALLY as a watchdog. It is cheap
    % (decimated), it only overrides when it disagrees by more than the trust
    % radius, and if the fine tier is working the two agree and it changes
    % nothing. If the fine tier is stuck, this is what breaks it out.
    needCoarse = ~isfinite(zEst);
    if ~needCoarse && o.relockEvery > 0
        sinceLock = sinceLock + 1;
        if sinceLock >= o.relockEvery
            zc = coarse_z(E, Ish);
            sinceLock = 0;
            if abs(zc - zEst) > o.trustUm
                zEst = zc;
                out.coarseUsed(n) = true;
            end
        end
    end
    if needCoarse
        zEst = coarse_z(E, Ish);
        sinceLock = 0;
        out.coarseUsed(n) = true;
    end

    % --- fine tier, per strip --------------------------------------------
    for s = 1:S
        r  = E.rows{s};
        Is = Ish(r,:);
        zs = zEst;
        p  = [0;0;0;0;0];          % defined before any early break
        rs = NaN;
        k  = nearest_k(E, zs);
        for it = 1:o.nIter
            k   = nearest_k(E, zs);
            sel = E.idx{k,s};
            R0  = E.ref(r,:,k);
            dvec = double(Is(sel)) - double(R0(sel));
            p    = (E.M{k,s} * dvec) ./ E.colN{k,s};
            if o.robust, p = robust_refit(E, k, s, dvec, p); end
            rs   = sqrt(mean(dvec.^2));
            zNew = E.zAxis(k) + p(1);
            if ~isfinite(zNew), break; end
            % Refuse a step that leaves the linear zone in one go -- a bad
            % strip should not throw the running estimate across the stack.
            zNew = zs + max(-o.trustUm, min(o.trustUm, zNew - zs));
            if abs(zNew - zs) < 0.02, zs = zNew; break; end
            zs = zNew;
        end
        out.dz(s,n)    = zs;
        out.dx(s,n)    = p(2) + out.dxInt(n);
        out.dy(s,n)    = p(3) + out.dyInt(n);
        out.gain(s,n)  = p(4);
        out.kUsed(s,n) = k;
        % Fit residual, in counts. A strip whose residual jumps has lost the
        % reference (bleaching, a vessel, an out-of-range excursion) and its
        % dz should not be believed -- the analysis gates on this.
        out.resid(s,n) = rs;
    end
    % Carry the frame's median across to the next frame: robust to one bad strip.
    zEst = median(out.dz(:,n), 'omitnan');
end

% Flatten strips into one time series at S x the frame rate. Column-major is
% exactly strip-within-frame order, which is acquisition order.
out.dzSerial = out.dz(:);
out.dxSerial = out.dx(:);
out.dySerial = out.dy(:);
out.zEnd     = zEst;    % hand to the next chunk as zInit
end

%% ============================== COARSE ==============================
function [z, r] = coarse_z(E, I)
%COARSE_Z  Correlate against every reference plane; parabolic peak.
dec = E.opts.coarseDecim;
a = double(reshape(I(1:dec:end, 1:dec:end), [], 1));
a = a - mean(a);
na = norm(a);  if na == 0, z = E.zAxis(E.kMid); r = 0; return; end
rr = (a.' / na) * E.coarseRef;
[r, k] = max(rr);
if k > 1 && k < numel(rr)
    den = rr(k-1) - 2*rr(k) + rr(k+1);
    if den ~= 0, k = k + 0.5*(rr(k-1) - rr(k+1))/den; end
end
z = interp1(1:E.K, E.zAxis, k, 'linear', 'extrap');
end

%% ============================== HELPERS ==============================
function o = fill_defaults(o)
d = struct( ...
    'nStrips',     8, ...    % 8 x 64 lines on a 512-line frame -> ~240 Hz
    'pixFrac',     0.08, ... % keep the top 8% of |dR/dz| pixels
    'minPix',      500, ...
    'smoothPx',    1.0, ...  % sigma of the pre-blur, px
    'coarseDecim', 4, ...
    'regDecim',    2, ...    % decimation for the lateral pre-registration FFT
    'regPeakRatio',8, ...    % phase-corr peak must exceed this x the median
    'trustUm',     2.5, ...  % linear zone ~ half the axial correlation width
    'relockEvery', 30, ...   % watchdog: coarse check every N frames, 0 = never
    'maxShiftPx',  40, ...
    'nIter',       2, ...
    'robust',      true, ... % trim outlier pixels and refit -- see robust_refit
    'trimFrac',    0.10, ... % fraction of pixels dropped per refit
    'ridge',       1e-6);
f = fieldnames(d);
for i = 1:numel(f)
    if ~isfield(o, f{i}) || isempty(o.(f{i})), o.(f{i}) = d.(f{i}); end
end
end

function k = nearest_k(E, z)
[~, k] = min(abs(E.zAxis - z));
end

function [gx, gy] = spatial_grad(A)
% Central differences, replicated at the edges.
gx = zeros(size(A),'single');  gy = zeros(size(A),'single');
gx(:,2:end-1) = (A(:,3:end) - A(:,1:end-2)) / 2;
gx(:,1) = A(:,2)-A(:,1);  gx(:,end) = A(:,end)-A(:,end-1);
gy(2:end-1,:) = (A(3:end,:) - A(1:end-2,:)) / 2;
gy(1,:) = A(2,:)-A(1,:);  gy(end,:) = A(end,:)-A(end-1,:);
end

function B = blur(A, sigma)
% Separable Gaussian without the Image Processing Toolbox.
r = max(1, ceil(3*sigma));
x = -r:r;
g = exp(-(x.^2)/(2*sigma^2));  g = g / sum(g);
B = conv2(g, g, padarray_rep(A, r), 'valid');
B = single(B);
end

function P = padarray_rep(A, r)
P = A([ones(1,r) 1:size(A,1) size(A,1)*ones(1,r)], ...
      [ones(1,r) 1:size(A,2) size(A,2)*ones(1,r)]);
end

function w = hann2(ny, nx)
% Raised cosine in 2-D, written out so the Signal Processing Toolbox is not
% needed (hann() lives there and is not guaranteed on the rig MATLAB).
wy = 0.5 - 0.5*cos(2*pi*(0:ny-1).'/(ny-1));
wx = 0.5 - 0.5*cos(2*pi*(0:nx-1)  /(nx-1));
w  = wy * wx;
end

function d = phase_shift_cached(E, I)
%PHASE_SHIFT_CACHED  Same answer as phase_shift, with the reference transform
% precomputed and the correlation run decimated. Returns the CORRECTION in
% FULL-RESOLUTION pixels: circshift(I, d) aligns I onto the reference.
dr = E.regDecim;
g  = double(I(1:dr:end, 1:dr:end)) .* E.regWin;
R  = E.regRefF .* conj(fft2(g));
R  = R ./ (abs(R) + eps);
c  = real(ifft2(R));
sz = E.regSize;
[pk, ix] = max(c(:));
[iy, jx] = ind2sub(sz, ix);
dy = iy - 1;  if dy > sz(1)/2, dy = dy - sz(1); end
dx = jx - 1;  if dx > sz(2)/2, dx = dx - sz(2); end
d = [dy dx] * dr;

% PEAK QUALITY GATE. Phase correlation is DEGENERATE on a field with no
% texture along one axis -- parallel vessels with nothing crossing them are the
% real case, and the correlation then has a ridge rather than a peak, so the
% argmax lands anywhere along it. Shifting by that garbage corrupts dz far
% worse than not shifting at all. So if the peak does not stand clear of the
% rest of the surface, decline to shift and let the least squares handle the
% lateral term on its own.
%
% Found 2026-09-13: test 9's in-plane vessel volume (pure horizontal stripes)
% silently lost its recovered slope, 1.07 -> 0.09, to exactly this.
med = median(abs(c(:)));
if ~isfinite(pk) || pk < E.opts.regPeakRatio * med
    d = [0 0];
end
end

function d = phase_shift(ref, im)
%PHASE_SHIFT  Integer shift that maps IM back onto REF -- i.e. the CORRECTION.
%
%   circshift(im, d) aligns im onto ref.
%   The DISPLACEMENT of im relative to ref is -d.
%
% Two traps this encodes, both of which have cost real numbers before:
%  1. MATLAB puts the zero-shift peak at index (1,1), not at the centre. An
%     off-by-one here is invisible in a slope and wrong in every displacement.
%  2. The returned value is the correction, not the displacement. Backwards in
%     a corrector, the loop drives the motion instead of cancelling it.
persistent win winSize
sz = size(ref);
if isempty(win) || ~isequal(winSize, sz)
    win = hann2(sz(1), sz(2));  winSize = sz;
end
F = fft2(double(ref) .* win);
G = fft2(double(im)  .* win);
R = F .* conj(G);
R = R ./ (abs(R) + eps);
c = real(ifft2(R));
[~, ix] = max(c(:));
[iy, jx] = ind2sub(sz, ix);
dy = iy - 1;  if dy > sz(1)/2, dy = dy - sz(1); end
dx = jx - 1;  if dx > sz(2)/2, dx = dx - sz(2); end
d = [dy dx];
end

function q = qtile(x, p)
% Quantile without the Statistics Toolbox.
x = x(isfinite(x));
if isempty(x), q = 0; return; end
x = sort(double(x(:)));
n = numel(x);
i = max(1, min(n, round(p * (n-1)) + 1));
q = x(i);
end
