%% hologram_fourier_sim.m
% Educational simulation of the *principles* of holography using the FFT.
%
% Core idea
% ---------
% A hologram is NOT just a Fourier transform -- it is the recorded INTENSITY
% of the interference between:
%     O = object wave  (light scattered by the object)
%     R = reference wave (a clean, known beam)
% The film records  H = |O + R|^2  (intensity only -> phase is "lost").
% The reference wave smuggles the lost phase back into the recorded fringes.
%
% Where the Fourier transform enters
% ----------------------------------
%   1. Propagation to the far field (Fraunhofer diffraction) of an object is
%      its 2D spatial Fourier transform:   O = FT(object).
%   2. Reconstruction (re-illuminating the hologram) is ANOTHER Fourier
%      transform of the recorded fringe pattern.
%
% This script:
%   (1) builds an object = random dots on the xy plane
%   (2) forms the object wave O = FT(object)  (far-field / Fourier hologram)
%   (3) adds an OFF-AXIS plane-wave reference R and records H = |O + R|^2
%   (4) reconstructs by FT(H), recovering the dots (+ DC term + twin image)
%
% Runqi Zhang lab -- standalone demo, safe to run directly in MATLAB.

clear; close all; clc;

%% -------------------- USER PARAMETERS --------------------
N           = 512;     % grid size (pixels), use a power of 2 for fast FFT
nDots       = 40;      % number of random object points
objHalfSpan = 0.12;    % object dots confined to central fraction of field (0..0.5)
                       %   (keeps object compact so its FT / twin images separate)
dotSigmaPix = 1.5;     % blur each dot to this Gaussian sigma (pixels), 0 = pure points
refTiltCyc  = 90;      % off-axis reference tilt, in cycles across the full field.
                       %   Larger -> carrier fringes finer, twin images pushed apart.
refAmpFactor = 2.0;    % reference amplitude relative to mean object-wave amplitude.
                       %   >~1 keeps the hologram in the "linear" (good) regime.
rngSeed     = 7;       % random seed for reproducible dot pattern

%% -------------------- 1. OBJECT: random dots --------------------
rng(rngSeed);

% coordinate grids (pixel indices centered at 0)
[xx, yy] = meshgrid(1:N, 1:N);
cx = N/2 + 1;  cy = N/2 + 1;

% random dot positions inside a central box
loBox = round(N*(0.5 - objHalfSpan));
hiBox = round(N*(0.5 + objHalfSpan));
dotX  = randi([loBox, hiBox], nDots, 1);
dotY  = randi([loBox, hiBox], nDots, 1);

object = zeros(N, N);
idx = sub2ind([N N], dotY, dotX);
object(idx) = 1;

% optional Gaussian blur so dots are small blobs (more physical than delta spikes)
if dotSigmaPix > 0
    object = imgaussfilt(object, dotSigmaPix);
end

%% -------------------- 2. OBJECT WAVE: far-field = FT(object) --------------------
% Fraunhofer (far-field) diffraction of the object field is its 2D FT.
% fftshift centers the zero spatial frequency for display & for placing the
% reference. This complex field O carries both amplitude AND phase.
O = fftshift(fft2(ifftshift(object)));

%% -------------------- 3. REFERENCE WAVE + RECORD HOLOGRAM --------------------
% Off-axis reference: a tilted plane wave  R = A * exp(i*2*pi*(fx*x + fy*y)).
% The tilt (carrier frequency) is what separates the real and twin images
% in the reconstruction (this is "off-axis" / Leith-Upatnieks holography).
fx = refTiltCyc / N;          % cycles per pixel along x
fy = refTiltCyc / N;          % cycles per pixel along y  (diagonal tilt)
Aref = refAmpFactor * mean(abs(O(:)));
R = Aref * exp(1i*2*pi*(fx*(xx-cx) + fy*(yy-cy)));

% The hologram is the recorded INTENSITY of the interference. The detector
% sees only |.|^2, discarding the complex phase -- this is the whole point.
H = abs(O + R).^2;
%   Expand: H = |O|^2 + |R|^2 + O.*conj(R) + conj(O).*R
%           \___DC/autocorr___/   \_real img_/  \_twin img_/

%% -------------------- 4. RECONSTRUCTION = FT(hologram) --------------------
% Re-illuminating the developed hologram and propagating to the far field is
% another Fourier transform. The cross terms O.*conj(R) and conj(O).*R are
% plane-wave-modulated, so they land at +/- the carrier offset: the object
% reappears as a real image on one side and a conjugate "twin" on the other,
% with the bright autocorrelation (DC) term in the middle.
recon = fftshift(fft2(ifftshift(H)));

%% -------------------- DISPLAY --------------------
figure('Name','Fourier hologram simulation','Color','w', ...
       'Position',[80 80 1200 800]);

subplot(2,3,1);
imagesc(object); axis image off; colormap(gca,gray);
title(sprintf('(1) Object: %d random dots', nDots));

subplot(2,3,2);
imagesc(log1p(abs(O))); axis image off; colormap(gca,parula);
title('(2) Object wave |O| = |FT(object)|  (log)');

subplot(2,3,3);
imagesc(angle(O)); axis image off; colormap(gca,hsv);
title('Object-wave PHASE (lost on recording)');

subplot(2,3,4);
imagesc(H); axis image off; colormap(gca,gray);
title('(3) Recorded hologram H = |O+R|^2');

% zoom on a patch to reveal the interference carrier fringes
subplot(2,3,5);
zc = cx-30:cx+30;
imagesc(H(zc,zc)); axis image off; colormap(gca,gray);
title('Zoom: interference fringes (carrier)');

subplot(2,3,6);
imagesc(log1p(abs(recon))); axis image off; colormap(gca,gray);
title('(4) Reconstruction = FT(H): real + twin + DC');

sgtitle('Holography principles via FFT:  object \rightarrow FT \rightarrow interfere with reference \rightarrow record |.|^2 \rightarrow FT back');

%% -------------------- WHAT THE METRICS MEAN --------------------
% - |O|  (panel 2): the diffraction / far-field amplitude. This is the part a
%        naive "just FFT the dots" view gives you -- but it has NO phase recorded.
% - angle(O) (panel 3): the phase a real detector CANNOT measure directly.
%        Holography's job is to preserve this.
% - H (panel 4): the actual hologram -- real, non-negative intensity with
%        interference fringes that encode O's phase relative to R.
% - recon (panel 6): re-Fourier-transforming H gives back the object (real
%        image), its mirror-conjugate (twin image, diagonally opposite), and a
%        bright central autocorrelation blob (the |O|^2+|R|^2 / DC term).
%        The off-axis reference tilt is what physically separates these.

fprintf('Done. Off-axis carrier = %d cycles/field; reference/object amp ratio = %.2f\n', ...
        refTiltCyc, refAmpFactor);
