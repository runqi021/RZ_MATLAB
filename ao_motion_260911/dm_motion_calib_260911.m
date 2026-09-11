function dm_motion_calib_260911()
%DM_MOTION_CALIB_260911  Measure how far the sample moves per micron of DM Zernike.
%
% STAGE 0 of the AO motion-correction pipeline.  Produces the only numbers a
% corrector actually needs: the 3x3 matrix G mapping
%
%       [tiltY tipX defocus] coefficients (um)  ->  [dx dy dz] at the sample (um)
%
% so a corrector inverts it:  coeffs = (dr_measured) / G.
%
% It also answers, empirically and against YOUR objective, the two things that
% cannot be read off the Z2C file:
%   - whether the coefficients are um RMS or um PtV (a factor of ~4 in range)
%   - where the response stops being linear and the PSF starts to go
%
% WHAT IT DOES
%   for each of tilt / tip / defocus:
%     save the current DM command, step a coefficient through COEFFS_UM,
%     grab and average NFRAMES at each step, restore the baseline.
%   then: lateral shift by phase correlation against the c = 0 frame;
%         axial shift by matching each frame into a baseline z-stack.
%
% ================================ SAFETY ================================
% Every send is baseline + motion, and the baseline (your system + sample
% aberration correction) is restored in the cleanup handler even on error or
% Ctrl-C.  Commands are checked against the ALPAO +-1 limit BEFORE being sent;
% a step that would clip is skipped, not clamped, because a clamped mirror is
% no longer producing the Zernike you think it is and would silently corrupt
% the fit.  Start with small COEFFS_UM and widen once you have seen it behave.
%
% REQUIRES: ScanImage running (hSI, hAOROIctrl in base), an FOV with structure
% that cross-correlates well (beads, or vasculature), and NO breathing motion
% during the run -- do it on a dead/anaesthetised-still prep or a slide, or the
% motion you are trying to measure gets swamped by the motion you are trying to
% correct.
%
% *** NOT YET RUN ON THE RIG.  Written against the ScanImage 2018b source; the
% *** hardware calls are unexecuted.  Read it before you trust it.
%
% Runqi Zhang / 2026-09-11

%% ========================= USER SETTINGS =========================
Z2C_FILE   = 'D:\RZ_ScanImage_script\SI2018bR1_2018-12-19_4a9264c4fc\BAX331-Z2C.mat';
ROW_TILTY  = 1;        % verified from actuator geometry: row 1 projects 0.988 onto a y-ramp
ROW_TIPX   = 2;        %                                  row 2 projects 0.986 onto an x-ramp
ROW_DEFOC  = 3;        %                                  row 3 projects 0.959 onto r^2

COEFFS_UM  = [-2 -1 -0.5 0 0.5 1 2];   % start small.  Widen after the first run.
NFRAMES    = 10;       % frames averaged per step
SETTLE_S   = 0.3;      % after a send, before grabbing
CHANNEL    = [];       % [] = first channel being displayed

DO_LATERAL = true;
DO_AXIAL   = true;     % needs a z-stack; set false for a quick lateral-only pass
ZSTACK_UM  = 40;       % total z range for the axial reference stack
ZSTEP_UM   = 1;

CMD_LIMIT  = 0.95;     % refuse to send if any actuator would exceed this
OUT_DIR    = fileparts(mfilename('fullpath'));

%% ========================= CONNECT =========================
hSI   = evalin('base','hSI');
hAO   = evalin('base','hAOROIctrl');
hDM   = hAO.hDMctrl;

S = load(Z2C_FILE);
fn = fieldnames(S); Z2C = S.(fn{1});
assert(size(Z2C,2) == 97, 'expected 97 actuators, got %d', size(Z2C,2));

rows  = [ROW_TILTY ROW_TIPX ROW_DEFOC];
names = {'tiltY','tipX','defocus'};

% The baseline is whatever the DM is holding now -- your aberration correction.
% There is no "read current command" call on the ALPAO driver, so it is taken
% from the AO data object and must be re-sent by us at the end.
baseline = get_baseline(hAO);
assert(numel(baseline) == 97, 'baseline command is not 97 long');
restore = onCleanup(@() safe_restore(hDM, baseline));

um_per_px = get_um_per_px(hSI);
fprintf('pixel size %.4f um/px\n', um_per_px);

%% ========================= AXIAL REFERENCE =========================
zs = []; stack = [];
if DO_AXIAL
    fprintf('\n[z-stack] baseline reference, %d um at %d um steps\n', ZSTACK_UM, ZSTEP_UM);
    [stack, zs] = grab_zstack(hSI, hDM, baseline, ZSTACK_UM, ZSTEP_UM, NFRAMES, CHANNEL, SETTLE_S);
end

%% ========================= SWEEP =========================
nC = numel(COEFFS_UM);
IM = cell(3, nC);  applied = nan(3, nC);

for m = 1:3
    fprintf('\n[%s] row %d\n', names{m}, rows(m));
    for k = 1:nC
        c   = COEFFS_UM(k);
        cmd = baseline(:).' + c * Z2C(rows(m),:);
        if max(abs(cmd)) > CMD_LIMIT
            fprintf('   c = %+5.2f um  SKIPPED, would reach |cmd| = %.3f\n', c, max(abs(cmd)));
            continue
        end
        hDM.senddata(cmd);
        pause(SETTLE_S);
        IM{m,k} = grab_avg(hSI, NFRAMES, CHANNEL);
        applied(m,k) = c;
        fprintf('   c = %+5.2f um  max|cmd| %.3f  mean %.0f\n', c, max(abs(cmd)), mean(IM{m,k}(:)));
    end
    hDM.senddata(baseline);  pause(SETTLE_S);
end

%% ========================= ANALYSE =========================
% G(m,:) = [dx dy dz] per micron of coefficient m.  Off-diagonals ARE the
% cross-talk, and the corrector must invert the whole matrix -- three separate
% scalars would leave tip-induced defocus uncorrected.
G = nan(3,3);  fits = struct();
for m = 1:3
    k0 = find(applied(m,:) == 0, 1);
    assert(~isempty(k0) && ~isempty(IM{m,k0}), 'no c = 0 frame for %s', names{m});
    ref = IM{m,k0};
    dx = nan(1,nC); dy = nan(1,nC); dz = nan(1,nC);
    for k = 1:nC
        if isempty(IM{m,k}), continue; end
        if DO_LATERAL
            d = phase_shift(ref, IM{m,k});
            dx(k) = d(2)*um_per_px;  dy(k) = d(1)*um_per_px;   % [row col] -> [y x]
        end
        if DO_AXIAL && ~isempty(stack)
            dz(k) = match_z(stack, zs, IM{m,k});
        end
    end
    ok = ~isnan(applied(m,:));
    G(m,1) = robust_slope(applied(m,ok), dx(ok));
    G(m,2) = robust_slope(applied(m,ok), dy(ok));
    G(m,3) = robust_slope(applied(m,ok), dz(ok));
    fits.(names{m}) = struct('c',applied(m,:),'dx',dx,'dy',dy,'dz',dz);
    fprintf('\n%-8s -> dx %+7.3f  dy %+7.3f  dz %+7.3f   um per um-coefficient\n', ...
            names{m}, G(m,1), G(m,2), G(m,3));
end

%% ========================= REPORT =========================
fprintf('\n================ GAIN MATRIX G ================\n');
fprintf('            dx        dy        dz     (um at sample per um coeff)\n');
for m = 1:3
    fprintf('%-8s %8.3f  %8.3f  %8.3f\n', names{m}, G(m,1), G(m,2), G(m,3));
end
fprintf('\ncondition number %.2f  (large = axes are not separable)\n', cond(G));

% Range: how far can each axis be pushed before an actuator clips, using the
% measured gain.  This is the number that decides whether the DM can cover your
% breathing amplitude at all.
fprintf('\n================ REACHABLE RANGE ================\n');
for m = 1:3
    cmax = (CMD_LIMIT - max(abs(baseline))) / max(abs(Z2C(rows(m),:)));
    fprintf('%-8s coeff +-%.2f um  ->  %+.1f um of motion along its main axis\n', ...
            names{m}, cmax, cmax*max(abs(G(m,:))));
end
fprintf('(baseline aberration correction already uses max|cmd| = %.3f of the +-1 stroke)\n', ...
        max(abs(baseline)));

out = fullfile(OUT_DIR, sprintf('dm_motion_calib_%s.mat', datestr(now,'yymmdd_HHMM')));
save(out, 'G','fits','COEFFS_UM','rows','names','baseline','um_per_px','Z2C_FILE');
fprintf('\nsaved %s\n', out);
plot_fits(fits, names, G, OUT_DIR);
end

%% ========================== HELPERS ==========================
function b = get_baseline(hAO)
% Whatever the DM is holding: system + sample aberration correction.
b = zeros(1,97);
try
    d = hAO.hAOdata.getDMdata();      % [nROIs x 97]
    if ~isempty(d), b = double(d(1,:)); end
catch
    warning(['could not read the current DM command from the AO data object; ' ...
             'assuming a flat mirror.  If your aberration correction was active, ' ...
             'ABORT and fix this -- the calibration would be taken about the ' ...
             'wrong operating point and the baseline would not be restored.']);
end
end

function safe_restore(hDM, baseline)
try, hDM.senddata(baseline); fprintf('[cleanup] DM restored to baseline\n');
catch ME, warning('COULD NOT RESTORE THE DM: %s', ME.message);
end
end

function u = get_um_per_px_mdf(hSI)
% ScanImage's own model: one objectiveResolution (38.9 um/deg) scaling both axes.
% Used only as a sanity check against the measured value -- see stage_ruler.
sf = hSI.hRoiManager.currentRoiGroup.rois(1).get(0);
u  = [sf.sizeXY(1) * hSI.objectiveResolution / sf.pixelResolutionXY(1), ...
      sf.sizeXY(2) * hSI.objectiveResolution / sf.pixelResolutionXY(2)];
end

function u = stage_ruler(hSI, n, ch, stepUm)
% MEASURE um/px with the stage, per axis, instead of trusting the MDF.
%
% Why not the MDF: objectiveResolution is a single number scaling both axes, but
% the stitching work measured x and y disagreeing by ~6% over 144 tile pairs --
% the resonant and galvo axes are not calibrated together.  One scalar is
% therefore ~6% wrong on one axis, and that error would land directly in the DM
% gain.  The stage is the ruler you already trust: 0.7815 um/device-unit, and
% commanded positions land where they are asked.
%
% No calibration grid needed.  Move a known distance, see how far the image went.
p0 = hSI.hMotors.motorPosition;
u  = nan(1,2);
for ax = 1:2
    I = cell(1,2);
    for s = [-1 1]
        p = p0;  p(ax) = p0(ax) + s*stepUm/2;
        hSI.hMotors.motorPosition = p;  pause(0.4);
        I{(s+3)/2} = grab_avg(hSI, n, ch);
    end
    hSI.hMotors.motorPosition = p0;  pause(0.4);
    d = phase_shift(I{1}, I{2});          % [dy dx] in px for a stepUm stage move
    if ax == 1, px = abs(d(2)); else, px = abs(d(1)); end
    u(ax) = stepUm / px;
    fprintf('   axis %d: %.1f um of stage -> %.2f px  =>  %.4f um/px\n', ax, stepUm, px, u(ax));
end
end

function I = grab_avg(hSI, n, ch)
hSI.hStackManager.framesPerSlice = n;
hSI.hStackManager.numSlices = 1;
hSI.startGrab();
while ~strcmp(hSI.acqState,'idle'), pause(0.02); drawnow limitrate; end
f = hSI.hDisplay.lastFrame;
if isempty(ch), idx = 1; else, idx = find(hSI.hDisplay.lastFrameChannels == ch, 1); end
I = single(f{idx});
end

function [stack, zs] = grab_zstack(hSI, hDM, baseline, range_um, step_um, n, ch, settle)
% Baseline z-stack with the STAGE, used only as the axial ruler.
z0 = hSI.hMotors.motorPosition(3);
zs = (-range_um/2 : step_um : range_um/2);
hDM.senddata(baseline); pause(settle);
stack = [];
for k = 1:numel(zs)
    p = hSI.hMotors.motorPosition;  p(3) = z0 + zs(k);
    hSI.hMotors.motorPosition = p;  pause(0.2);
    I = grab_avg(hSI, n, ch);
    if isempty(stack), stack = zeros([size(I) numel(zs)],'single'); end
    stack(:,:,k) = I; %#ok<AGROW>
end
p = hSI.hMotors.motorPosition;  p(3) = z0;  hSI.hMotors.motorPosition = p;
pause(0.3);
end

function d = phase_shift(ref, im)
% Subpixel translation by phase correlation.  Hann-windowed: the frame edges
% move content in and out of view, and an unwindowed FFT reads that as signal.
%
% SIGN, verified against synthetic translations: if `im` is `ref` displaced by
% +5 px, this returns -5.  So d is the shift that maps `im` BACK onto `ref` --
% the correction, not the displacement.  Get this backwards in a corrector and
% it drives the motion instead of cancelling it, doubling the error rather than
% removing it, and the closed loop runs away.  Measured accuracy: exact at zero
% shift, ~0.17 px at +-11 px (the window loses content at large shifts).
w = hann(size(ref,1)) * hann(size(ref,2)).';
A = fft2(double(ref - mean(ref(:))) .* w);
B = fft2(double(im  - mean(im(:)))  .* w);
R = A .* conj(B);  R = R ./ max(abs(R), eps);
c = real(ifft2(R));
[~, i] = max(c(:));
[r, cc] = ind2sub(size(c), i);
% -1 because MATLAB indexes from 1: the zero-shift peak sits at (1,1) and must
% come out as a shift of 0.  Without it every measurement is off by one pixel,
% which is a constant offset in the fit -- invisible in the slope but wrong in
% every reported displacement.
d = [cent(c, r, 1), cent(c, cc, 2)] - 1;
d = d - (d > size(c)/2) .* size(c);        % wrap to signed shift
end

function v = cent(c, i, dim)
% 3-point parabolic peak interpolation along one axis.
n = size(c, dim);
im1 = mod(i-2, n) + 1;  ip1 = mod(i, n) + 1;
if dim == 1, y = [c(im1,1) c(i,1) c(ip1,1)]; else, y = [c(1,im1) c(1,i) c(1,ip1)]; end
den = y(1) - 2*y(2) + y(3);
if den == 0, v = i; else, v = i + 0.5*(y(1) - y(3))/den; end
end

function dz = match_z(stack, zs, im)
% Which z-plane does this frame look like?  Correlation against every slice,
% parabolic interpolation on the peak for sub-step resolution.
n = numel(zs);  r = zeros(1,n);
a = double(im(:)); a = a - mean(a);
for k = 1:n
    b = double(reshape(stack(:,:,k), [], 1)); b = b - mean(b);
    r(k) = (a.'*b) / (norm(a)*norm(b) + eps);
end
[~, k] = max(r);
if k > 1 && k < n
    den = r(k-1) - 2*r(k) + r(k+1);
    if den ~= 0, k = k + 0.5*(r(k-1) - r(k+1))/den; end
end
dz = interp1(1:n, zs, k, 'linear', 'extrap');
end

function s = robust_slope(x, y)
ok = isfinite(x) & isfinite(y);
if nnz(ok) < 3, s = NaN; return; end
p = polyfit(x(ok), y(ok), 1);
s = p(1);
end

function plot_fits(fits, names, G, outDir)
hf = figure('Color','w','Position',[100 100 1100 340]);
for m = 1:3
    ax = subplot(1,3,m); hold(ax,'on');
    f = fits.(names{m});
    plot(ax, f.c, f.dx, 'o-', 'DisplayName','dx');
    plot(ax, f.c, f.dy, 's-', 'DisplayName','dy');
    plot(ax, f.c, f.dz, '^-', 'DisplayName','dz');
    xlabel(ax,'coefficient (um)'); ylabel(ax,'sample motion (um)');
    title(ax, sprintf('%s   [%.2f %.2f %.2f]', names{m}, G(m,1), G(m,2), G(m,3)));
    legend(ax,'Location','best'); box(ax,'off'); grid(ax,'on');
end
exportgraphics(hf, fullfile(outDir,'dm_motion_calib.png'), 'Resolution', 150);
end
