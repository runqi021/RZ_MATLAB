% svd_motion_artifact_260812.m
% -----------------------------------------------------------------------
%  Is the "activity" in this FOV real, or is it motion?
%
%  Companion to Run_SVD_recon_260811.m. Takes the SVD modes it wrote plus the
%  NoRMCorre shifts and the ROI dF/F, and asks the question three ways:
%
%   1. TEMPORAL   correlate each mode's v_k(t) with the rigid shifts sx, sy,
%                 |shift| and frame-to-frame speed. Cheap, but weak on its
%                 own: sx/sy only capture BULK translation, so a mode can be
%                 pure artifact and still correlate poorly (non-rigid warp,
%                 z-motion, scan-phase). Never clear a mode on this alone.
%
%   2. SPATIAL    the decisive one. A rigid translation by (dx,dy) changes the
%                 image by dI = dx*Gx + dy*Gy, so a translation artifact's
%                 eigenimage must lie in the span of the mean image's spatial
%                 gradients. Regress U_k on [Gx Gy] and report R^2:
%                     R^2 -> 1  the mode IS the image shifting
%                     R^2 -> 0  the mode is not explained by translation
%                 This is what tells a red/blue EDGE DIPOLE from a filled soma,
%                 and it needs no shift record at all.
%
%   3. BREATH     if breath_pc1.mat is present, coherence of each mode with
%                 the breathing trace -- a breathing-locked motion mode is the
%                 one that fakes breath-locked calcium.
%
%  Then the same temporal test on every cpSAM ROI's dF/F, so you can see which
%  ROIs inherit the motion.
%
%  Output: <fovPath>\svd_check_260811\motion_artifact.png / .mat
%
%  Runqi Zhang / 2026-08-12

clear; clc; close all;

%% ===================== USER-EDITABLE =====================
fovPath = 'D:\Ventral_surface_summary\Vgat\0730\cell1\roi1_3.2x_x1100y1050_z6_6000f_12lp_00001';
svdDir  = 'svd_check_260811';       % subfolder written by Run_SVD_recon_260811
variants = {'raw','mcmc'};
nMode   = 12;                       % modes to test
maxLagSec = 0.5;                    % allow a small lag in the temporal test
doSave  = true;
% =========================================================

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

pp   = strsplit(regexprep(fovPath,'[\\/]+$',''), {'\','/'});
stem = pp{end};

%% ---- motion signal ----
mcInfoP = fullfile(fovPath, [stem '_ch1_MCinfo.mat']);
assert(isfile(mcInfoP), 'No MCinfo: %s', mcInfoP);
MC = load(mcInfoP, 'sx','sy','r');
sx = MC.sx(:);  sy = MC.sy(:);
mot = struct();
mot.sx    = sx;
mot.sy    = sy;
mot.absSh = hypot(sx, sy);
mot.speed = [0; hypot(diff(sx), diff(sy))];      % frame-to-frame jump
motNames  = {'sx','sy','|shift|','speed'};
motMat    = [mot.sx, mot.sy, mot.absSh, mot.speed];

fprintf('motion: %d frames | sx %.2f..%.2f px, sy %.2f..%.2f px, max speed %.2f px/frame\n', ...
        numel(sx), min(sx), max(sx), min(sy), max(sy), max(mot.speed));

%% ---- optional breath trace ----
brP = fullfile(fovPath,'breath_pc1.mat');
breath = [];
if isfile(brP)
    B = load(brP);
    fn = fieldnames(B);
    for i = 1:numel(fn)
        v = B.(fn{i});
        if isnumeric(v) && isvector(v) && numel(v) > 100
            breath = double(v(:)); break;
        end
    end
end
fprintf('breath trace: %s\n', mat2str(~isempty(breath)));

%% ---- per variant ----
RES = struct();
for iv = 1:numel(variants)
    lab = variants{iv};
    f = fullfile(fovPath, svdDir, lab, 'SVD_result.mat');
    if ~isfile(f); warning('missing %s', f); continue; end
    Sv = load(f); R = Sv.SVD_result;

    Hc = R.Hc; Wc = R.Wc; T = R.T; fps = R.fps;
    K  = min(nMode, size(R.V,2));

    % align motion to the SVD frames (both are the same 2970 acquisition frames)
    m = motMat;
    if size(m,1) ~= T
        n = min(size(m,1), T);
        m = m(end-n+1:end, :);          % front-drop align, same rule as the SVD
        warning('%s: motion %d vs SVD %d frames -- using last %d', ...
                lab, size(motMat,1), T, n);
    end

    % ---- 2. SPATIAL: regress each eigenimage on the mean-image gradients ----
    I = reshape(double(R.mu), [Hc Wc]);
    [Gx, Gy] = gradient(I);
    A = [Gx(:), Gy(:), ones(Hc*Wc,1)];

    % AXIAL (z) motion basis. Defocus convolves the image with a widening PSF,
    % and to first order d/dz of a Gaussian blur is proportional to the
    % Laplacian: sharp bright structures (somata) dim while the surrounding
    % neuropil brightens. So a z-motion mode lives in span{I, lap(I)}, NOT in
    % span{Gx,Gy}. A LOW gradient R^2 therefore means "not in-plane
    % translation" -- it does NOT mean "not an artifact", which is exactly the
    % case in-plane motion correction cannot fix.
    L = del2(I);
    Az = [L(:), I(:), ones(Hc*Wc,1)];

    gradR2 = nan(K,1);
    defocR2 = nan(K,1);
    tcorr  = nan(K, size(m,2));
    for k = 1:K
        u = double(R.U(:,k));
        b = A \ u;                      % least squares on [Gx Gy 1]
        res = u - A*b;
        gradR2(k) = 1 - sum(res.^2)/sum((u-mean(u)).^2);

        bz = Az \ u;
        rz = u - Az*bz;
        defocR2(k) = 1 - sum(rz.^2)/sum((u-mean(u)).^2);

        v = double(R.V(:,k));
        for j = 1:size(m,2)
            tcorr(k,j) = max_abs_xcorr(v, m(:,j), round(maxLagSec*fps));
        end
    end

    % ---- 3. breath coherence-ish: max |xcorr| with the breath trace ----
    bcorr = nan(K,1);
    if ~isempty(breath)
        bb = breath;
        if numel(bb) ~= T
            n = min(numel(bb), T);
            bb = bb(end-n+1:end);
        end
        for k = 1:K
            v = double(R.V(1:numel(bb),k));
            bcorr(k) = max_abs_xcorr(v, bb, round(maxLagSec*fps));
        end
    end

    RES.(lab) = struct('varExp',double(R.varExp(1:K)),'gradR2',gradR2, ...
                       'defocR2',defocR2,'tcorr',tcorr,'bcorr',bcorr, ...
                       'K',K,'fps',fps);

    fprintf('\n=== %s ===\n', lab);
    fprintf('mode  var%%  shiftR2  defocR2   |r| sx    sy   |sh|  speed   breath\n');
    for k = 1:K
        fprintf('%4d %6.2f   %5.2f    %5.2f   %5.2f %5.2f %5.2f %5.2f   %5.2f\n', ...
            k, 100*RES.(lab).varExp(k), gradR2(k), defocR2(k), tcorr(k,1), ...
            tcorr(k,2), tcorr(k,3), tcorr(k,4), bcorr(k));
    end
end

%% ---- ROI dF/F vs motion ----
roi = struct('n',0);
h = dir(fullfile(fovPath,'*_cpSAM_output.mat'));
if ~isempty(h)
    sd = load(fullfile(h(1).folder,h(1).name),'F');
    try
        [fpsR,~] = detect_session_fps(fovPath);
    catch
        fpsR = 30;
    end
    dout = helper.dFF_RZ(double(sd.F),'FPS',fpsR,'BaselineWinSec',20);
    D = dout.dFF;
    n = min(size(D,1), size(motMat,1));
    Dm = D(end-n+1:end,:);  Mm = motMat(end-n+1:end,:);
    rr = nan(size(D,2), size(motMat,2));
    for k = 1:size(D,2)
        for j = 1:size(Mm,2)
            rr(k,j) = max_abs_xcorr(Dm(:,k), Mm(:,j), round(maxLagSec*fpsR));
        end
    end
    roi = struct('n',size(D,2),'r',rr);
    fprintf('\nROI dF/F vs motion (|r|, max over +-%.1fs lag), %d ROIs:\n', maxLagSec, roi.n);
    for j = 1:numel(motNames)
        fprintf('  %-8s median %.2f, max %.2f, n>0.3: %d\n', motNames{j}, ...
                median(rr(:,j)), max(rr(:,j)), nnz(rr(:,j)>0.3));
    end
end

%% ---- figure ----
fh = figure('Color','w','Visible','off','Position',[40 40 1500 900]);
tl = tiledlayout(fh, 3, 2, 'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1,[1 2]); hold(ax,'on');
tt = (0:size(motMat,1)-1)/30;
plot(ax, tt, mot.sx, 'LineWidth',0.6);
plot(ax, tt, mot.sy, 'LineWidth',0.6);
plot(ax, tt, mot.absSh, 'k', 'LineWidth',0.6);
legend(ax, {'sx','sy','|shift|'}, 'Box','off','Location','northeast');
xlabel(ax,'time (s)'); ylabel(ax,'px'); xlim(ax,[tt(1) tt(end)]);
title(ax,'rigid motion estimate (NoRMCorre)'); grid(ax,'on');

for iv = 1:numel(variants)
    lab = variants{iv};
    if ~isfield(RES,lab); continue; end
    Q = RES.(lab);

    ax = nexttile(tl, 2+iv);
    b = bar(ax, [Q.gradR2, Q.defocR2, max(Q.tcorr,[],2)], 'grouped');
    b(1).FaceColor = [0.85 0.33 0.10];
    b(2).FaceColor = [0.47 0.67 0.19];
    b(3).FaceColor = [0.00 0.45 0.74];
    legend(ax, {'in-plane shift: R^2 on [\nabla I]', ...
                'defocus / z: R^2 on [\nabla^2I, I]', ...
                'temporal: max |r| w/ motion'}, ...
           'Box','off','Location','northeast','FontSize',7);
    xlabel(ax,'mode'); ylabel(ax,'score'); ylim(ax,[0 1]);
    title(ax, sprintf('%s -- artifact scores', lab)); grid(ax,'on');

    ax = nexttile(tl, 4+iv);
    scatter(ax, Q.defocR2, Q.bcorr, 60, 100*Q.varExp, 'filled');
    for k = 1:Q.K
        text(ax, Q.defocR2(k), Q.bcorr(k), sprintf(' %d',k), 'FontSize',8);
    end
    xlabel(ax,'defocus R^2  (1 = pure z-motion)');
    ylabel(ax,'|r| with breathing');
    xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,'on');
    cb = colorbar(ax); cb.Label.String = 'var explained (%)';
    title(ax, sprintf('%s -- breathing-locked defocus?', lab));
end

title(tl, sprintf('%s  |  motion vs activity', stem), ...
      'Interpreter','none','FontWeight','bold');

%% ---- save ----
if doSave
    outDir = fullfile(fovPath, svdDir);
    exportgraphics(fh, fullfile(outDir,'motion_artifact.png'), 'Resolution',200);
    save(fullfile(outDir,'motion_artifact.mat'), 'RES','roi','mot','motNames','variants');
    fprintf('\nsaved %s\n', fullfile(outDir,'motion_artifact.{png,mat}'));
end
close(fh);


%% ========================================================================
function v = max_abs_xcorr(a, b, maxLag)
% max |Pearson r| between a and b over lags in +-maxLag, NaN-safe.
%
% Computed as a TRUE per-lag Pearson on the overlapping samples. xcorr(...,
% 'coeff') is not that: it normalises by the full-sequence energy while the
% numerator at lag L sums only N-|L| terms, so it drifts from r as |L| grows.
    a = double(a(:)); b = double(b(:));
    n = min(numel(a), numel(b));
    a = a(1:n); b = b(1:n);
    ok = isfinite(a) & isfinite(b);
    a(~ok) = 0; b(~ok) = 0;
    if std(a)==0 || std(b)==0; v = 0; return; end

    v = 0;
    for L = -maxLag:maxLag
        if L >= 0
            x = a(1+L:end);  y = b(1:end-L);
        else
            x = a(1:end+L);  y = b(1-L:end);
        end
        if numel(x) < 10; continue; end
        r = corr(x, y);
        if isfinite(r); v = max(v, abs(r)); end
    end
end
