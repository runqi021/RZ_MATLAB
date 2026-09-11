%% holo_stim_cgh_sim.m
% Simulation of HOLOGRAPHIC PHOTOSTIMULATION as used in neuroscience
% (computer-generated holography, CGH, on a phase-only SLM).
%
% How this DIFFERS from recording/display holography
% ---------------------------------------------------
%   Recording holography : object  -> record interference  -> reconstruct.
%   Neural stimulation    : desired light pattern (which neurons to hit)
%                           -> COMPUTE a phase-only hologram -> display on SLM
%                           -> the SLM's Fourier transform lights up the targets.
%   It is the INVERSE problem, there is NO reference beam and NO camera, and
%   the SLM can only set PHASE (not amplitude) -- which is why we must SEARCH
%   for the phase with an iterative algorithm (Gerchberg-Saxton).
%
% Shared principle: a lens performs a Fourier transform, so
%       field at the sample (neurons) = FT( field at the SLM ).
%
% PART 1: Gerchberg-Saxton -- design one phase mask to hit many 2D targets.
% PART 2: Prism+lens superposition -- place targets in 3D and verify the
%         focal stack (each spot is sharp only at its assigned depth).
%
% Standalone demo. Safe to run directly in MATLAB.

clear; close all; clc;

%% ===================================================================
%  PART 1 -- GERCHBERG-SAXTON: phase-only hologram for 2D targets
%  ===================================================================

%% -------------------- USER PARAMETERS (Part 1) --------------------
N         = 256;    % SLM / grid size (pixels). Power of 2 -> fast FFT.
nNeurons  = 12;     % number of target cells to stimulate
spotRad   = 2;      % target spot radius (pixels) ~ soma-sized disk
targetSpan= 0.30;   % keep targets within central +/- this fraction (avoid edges)
nIter     = 40;     % Gerchberg-Saxton iterations
rngSeed   = 3;      % reproducible target layout

rng(rngSeed);

% centered SLM pixel-index grids (used throughout)
v = (0:N-1) - floor(N/2);
[Xc, Yc] = meshgrid(v, v);

%% -------------------- Build the TARGET pattern --------------------
% Desired light: bright disks where the neurons are, dark everywhere else.
lo = round(N*(0.5 - targetSpan));
hi = round(N*(0.5 + targetSpan));
tx = randi([lo hi], nNeurons, 1);
ty = randi([lo hi], nNeurons, 1);

targetI = zeros(N);
[gx, gy] = meshgrid(1:N, 1:N);
for k = 1:nNeurons
    targetI( (gx-tx(k)).^2 + (gy-ty(k)).^2 <= spotRad^2 ) = 1;
end
targetAmp = sqrt(targetI);          % work with field amplitude = sqrt(intensity)
targetMask = targetI > 0;           % logical map of where neurons are

% Convenience handles for the (centered) forward/inverse optical FT.
% Forward:  SLM plane -> sample plane.   Inverse: sample -> SLM.
slm2samp = @(A) fftshift(fft2(ifftshift(A)));
samp2slm = @(S) fftshift(ifft2(ifftshift(S)));

%% -------------------- NAIVE attempt (why we need GS) --------------------
% Tempting shortcut: take the inverse FT of the target and just KEEP its phase
% (throw away amplitude, since the SLM is phase-only). It does not work well --
% the reconstruction is full of speckle and the spots are dim/uneven.
A_naive   = exp(1i*angle(samp2slm(targetAmp)));
I_naive   = abs(slm2samp(A_naive)).^2;

%% -------------------- GERCHBERG-SAXTON loop --------------------
% Iterate between the two planes, forcing what we KNOW at each end:
%   * at the SAMPLE plane: amplitude must equal the target (bright spots only)
%   * at the SLM plane   : amplitude must be flat (phase-only device)
% Keep the freely-evolving PHASE each time; that is the degree of freedom GS
% exploits to satisfy both constraints simultaneously.
A = exp(1i*2*pi*rand(N));           % start from random SLM phase
effHist  = zeros(nIter,1);          % energy fraction landing on targets
uniHist  = zeros(nIter,1);          % how even the spots are (1 = perfect)

for it = 1:nIter
    S = slm2samp(A);                          % propagate SLM -> sample
    Iout = abs(S).^2;

    % --- diagnostics on this iteration ---
    onTarget  = sum(Iout(targetMask));
    effHist(it) = onTarget / sum(Iout(:));    % diffraction efficiency
    spotVals  = perSpotIntensity(Iout, tx, ty, spotRad);
    uniHist(it) = 1 - (max(spotVals)-min(spotVals)) / (max(spotVals)+min(spotVals)+eps);

    % --- sample-plane constraint: impose target amplitude, KEEP phase ---
    S = targetAmp .* exp(1i*angle(S));
    A = samp2slm(S);                          % propagate back sample -> SLM

    % --- SLM-plane constraint: phase-only (unit amplitude), KEEP phase ---
    A = exp(1i*angle(A));
end

slmPhase = angle(A);                          % <-- THIS is the hologram we display
I_gs     = abs(slm2samp(A)).^2;               % final reconstruction

%% -------------------- DISPLAY (Part 1) --------------------
figure('Name','Part 1: Gerchberg-Saxton CGH','Color','w', ...
       'Position',[60 60 1200 760]);

subplot(2,3,1);
imagesc(targetI); axis image off; colormap(gca,gray);
title(sprintf('TARGET: %d neurons to stimulate', nNeurons));

subplot(2,3,2);
imagesc(I_naive); axis image off; colormap(gca,gray);
title({'NAIVE phase-only (ifft & keep phase)','-> speckly, uneven, inefficient'});

subplot(2,3,3);
imagesc(slmPhase); axis image off; colormap(gca,hsv); colorbar;
title('SLM PHASE MASK (the computed hologram)');

subplot(2,3,4);
imagesc(I_gs); axis image off; colormap(gca,gray);
title('GS reconstruction -> clean target spots');

subplot(2,3,5);
plot(1:nIter, 100*effHist, 'LineWidth', 1.6); grid on;
xlabel('iteration'); ylabel('on-target energy (%)');
title(sprintf('Efficiency: %.1f%% final', 100*effHist(end)));

subplot(2,3,6);
plot(1:nIter, uniHist, 'LineWidth', 1.6); grid on; ylim([0 1]);
xlabel('iteration'); ylabel('uniformity (1 = perfect)');
title(sprintf('Uniformity: %.2f final', uniHist(end)));

sgtitle('Neural holographic stimulation = SOLVE for a phase-only SLM hologram that focuses light onto chosen cells');

%% ===================================================================
%  PART 2 -- PRISM + LENS SUPERPOSITION: place targets in 3D
%  ===================================================================
% A different, very common CGH method: build the phase as the angle of a SUM
% of simple per-target terms, each a "prism" (blazed grating -> lateral shift)
% plus a "lens" (quadratic phase -> axial shift). One mask, spots in a VOLUME.
%
%   grating(sx,sy): shifts a focal spot to lateral pixel (sx,sy)
%   lens(z)       : moves a spot's best-focus plane to depth z
%
% We VERIFY by computing the focal stack: propagate to several depths and show
% each target is sharp ONLY at its assigned z (defocus blurs it elsewhere).

%% -------------------- USER PARAMETERS (Part 2) --------------------
% Each row: [lateral_x_pixels, lateral_y_pixels, depth_z]  (relative to center)
targets3D = [ -40, -20, -1.0;
               30,  35,  0.0;
              -10,  45,  1.0;
               45, -30,  0.5 ];
defocusScale = 16*pi / N^2;     % maps depth units -> pupil quadratic phase
zPlanes      = [-1.0 -0.5 0.0 0.5 1.0];   % focal-stack depths to render

rho2 = Xc.^2 + Yc.^2;           % pupil radius^2 (for the lens/defocus terms)

%% -------------------- Build the 3D superposition hologram --------------------
M = zeros(N);
for k = 1:size(targets3D,1)
    sx = targets3D(k,1);  sy = targets3D(k,2);  zt = targets3D(k,3);
    grating = 2*pi*(sx*Xc + sy*Yc) / N;     % prism -> lateral position (sx,sy)
    lens    = -defocusScale * zt * rho2;     % lens  -> in-focus at depth zt
    M = M + exp(1i*(grating + lens));        % superpose complex contributions
end
slmPhase3D = angle(M);                       % phase-only mask for all 3D targets
A3D = exp(1i*slmPhase3D);

%% -------------------- Render the FOCAL STACK to verify 3D --------------------
% Defocus to depth z is a quadratic pupil phase D(z). A target's lens term
% cancels D(z) exactly at z = zt -> sharp there, blurred elsewhere.
nz = numel(zPlanes);
stack = cell(nz,1);
for j = 1:nz
    D = exp(1i*defocusScale * zPlanes(j) * rho2);    % propagate to this depth
    stack{j} = abs(slm2samp(A3D .* D)).^2;
end

%% -------------------- DISPLAY (Part 2) --------------------
figure('Name','Part 2: 3D targeting (prism + lens)','Color','w', ...
       'Position',[80 80 1280 520]);

subplot(2, nz, 1);
imagesc(slmPhase3D); axis image off; colormap(gca,hsv);
title('3D SLM phase mask');

% draw the intended targets, color-coded by depth, for reference
subplot(2, nz, 2);
hold on;
for k = 1:size(targets3D,1)
    plot(targets3D(k,1), targets3D(k,2), 'o', 'MarkerSize',8, 'LineWidth',1.5);
    text(targets3D(k,1)+2, targets3D(k,2), sprintf('z=%.1f',targets3D(k,3)), 'FontSize',8);
end
axis equal; grid on; xlim([-N/2 N/2]); ylim([-N/2 N/2]); set(gca,'YDir','reverse');
title('intended targets (x,y,z)'); xlabel('x px'); ylabel('y px');

% focal stack: each spot should snap into focus at its own depth
clim_hi = max(cellfun(@(s) max(s(:)), stack));
for j = 1:nz
    subplot(2, nz, nz + j);
    imagesc(stack{j}, [0 clim_hi]); axis image off; colormap(gca,gray);
    title(sprintf('focal plane z = %.1f', zPlanes(j)));
end

sgtitle('One phase mask, targets in 3D: each cell is in focus only at its assigned depth');

fprintf('Part 1: GS efficiency %.1f%%, uniformity %.2f over %d iters.\n', ...
        100*effHist(end), uniHist(end), nIter);
fprintf('Part 2: %d targets placed across depths %.1f..%.1f.\n', ...
        size(targets3D,1), min(targets3D(:,3)), max(targets3D(:,3)));

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================
function vals = perSpotIntensity(Iout, tx, ty, r)
% Mean intensity within each target spot -- used for the uniformity metric.
    n = numel(tx);
    vals = zeros(n,1);
    [gx, gy] = meshgrid(1:size(Iout,2), 1:size(Iout,1));
    for k = 1:n
        m = (gx-tx(k)).^2 + (gy-ty(k)).^2 <= r^2;
        vals(k) = mean(Iout(m));
    end
end
