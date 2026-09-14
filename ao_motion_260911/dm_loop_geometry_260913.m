function S = dm_loop_geometry_260913()
%DM_LOOP_GEOMETRY_260913  Set the scan geometry that sets the loop's frame rate.
%
% The correction loop leaves 2*sin(pi*f*tau) of the motion, so DELAY IS THE ONLY
% LEVER and the frame period is most of it. This is the script that changes it.
%
%   dm_axial_closeloop / dm_loop_validate INHERIT geometry and never set it --
%   deliberately, since a reference stack is only valid for the geometry it was
%   taken at. So set the geometry here first, THEN take a reference, THEN run.
%
% ============================ THE LADDER ============================
% Predicted median removal across the 150-recording ventral archive, from each
% recording's own breath spectrum through the loop's transfer function:
%
%   512x256, 58.3 Hz, delay 1.46 frames (as run 2026-09-13) ....... 63%
%   + extrapolate the disturbance (gx fix, already patched) ....... 76%
%   + 128 lines -> 110 Hz ......................................... 86%
%   + flyback at its 126 us floor -> 122 Hz ....................... 88%
%
% The last rung is worth 2 points and costs nothing. The 128-line rung is worth
% 10 and costs either pixel shape or field height -- that is the real decision,
% and PIXEL_MODE below is where you make it.
%
% ======================= CROP OR SUBSAMPLE =======================
% Frame rate depends ONLY on line count, never on scan amplitude, so both modes
% give the same rate. They differ in what you keep. At 2x zoom (0.9462 um/px
% measured, NOT the 1.4815 the header claims -- see project_autostitch_ncc):
%
%   'crop'       484 x 121 um field, 0.946 x 0.946 um pixels
%                forceSquarePixels TRUE: ScanImage shrinks the y amplitude with
%                the lines, so pixels stay square and you see less.
%
%   'subsample'  484 x 242 um field, 0.946 x 1.892 um pixels (2:1)
%                forceSquarePixels FALSE and the y amplitude left alone, so the
%                field is kept and y is sampled half as finely.
%
% WHICH IS BETTER IS NOT SETTLED. Subsampling keeps more vasculature in view,
% which a five-parameter gradient fit wants; cropping keeps the high spatial
% frequencies that carry dR/dz, which is where the axial signal lives. THE
% ESTIMATOR GAIN DECIDES IT: it was 0.932 at 512x256. Take a reference each way
% and run the sweep; whichever holds closest to 0.932 is right for your sample.
%
% *** A GEOMETRY CHANGE INVALIDATES YOUR REFERENCE STACK. *** A 256-line
% reference cannot register 128-line frames. Re-take it -- under a minute.
%
% ============================== SAFETY ==============================
% Nothing here moves the stage, the DM or the laser. It only writes scan
% settings, and only while idle. The previous values are printed so you can put
% them back, and returned in S.before.
%
% *** NOT YET RUN ON THE RIG. *** Written against the ScanImage 2018b source.
%
% Runqi Zhang / 2026-09-13

%% ========================= USER SETTINGS =========================
APPLY        = false;    % <<< false reports only. Set true to actually change it.

TARGET_LINES = 128;      % 256 = now (58 Hz) | 128 = 110-122 Hz | 116 = 120 Hz at 1 ms flyback
PIXEL_MODE   = 'subsample';   % 'crop' | 'subsample'  -- see above
FLYBACK_US   = 126.4;    % 126.4 = the floor (1 resonant cycle). 1011 = the default.
                         % ResScan rounds UP to whole cycles, so only multiples
                         % of 1/scannerFrequency are reachable.
%% =================================================================

hSI = evalin('base','hSI');
if ~strcmpi(hSI.acqState,'idle')
    error('dm_geom:notIdle','ScanImage is not idle -- stop the acquisition first.');
end
hR = hSI.hRoiManager;  hS = hSI.hScan2D;

S.before = snap(hSI);
fprintf('\n================ CURRENT ================\n');
show(S.before);

if ~APPLY
    fprintf('\n[report only] set APPLY = true to change it.\n');
    fprintf('predicted at the target geometry:\n');
    predict(hS.scannerFrequency, TARGET_LINES, FLYBACK_US*1e-6, S.before);
    return
end

%% ========================= APPLY =========================
% ORDER MATTERS. linesPerFrame cannot differ from pixelsPerLine while
% forceSquarePixelation is on, so that flag comes off first. forceSquarePixels
% is what decides crop vs subsample, so it is set before the line count -- with
% it TRUE, ScanImage recomputes the y amplitude for us; with it FALSE we put the
% old amplitude back afterwards, because setting linesPerFrame nudges it.
samSlow0 = hR.scanAngleMultiplierSlow;

hR.forceSquarePixelation = false;
switch lower(PIXEL_MODE)
    case 'crop',      hR.forceSquarePixels = true;
    case 'subsample', hR.forceSquarePixels = false;
    otherwise, error('dm_geom:mode','PIXEL_MODE must be ''crop'' or ''subsample''');
end
hR.linesPerFrame = TARGET_LINES;
if strcmpi(PIXEL_MODE,'subsample')
    hR.scanAngleMultiplierSlow = samSlow0;      % keep the field
end
hS.flybackTimePerFrame = FLYBACK_US*1e-6;

S.after = snap(hSI);
fprintf('\n================ APPLIED ================\n');
show(S.after);

%% ========================= VERIFY =========================
% Read back rather than assume: the flyback quantises, and forceSquarePixels can
% overrule an amplitude you thought you set.
fprintf('\n---- checks ----\n');
ok = true;
if S.after.lines ~= TARGET_LINES
    fprintf('  !! linesPerFrame is %d, asked for %d\n', S.after.lines, TARGET_LINES); ok = false;
end
fbWant = FLYBACK_US*1e-6;
nCyc   = ceil(fbWant*hS.scannerFrequency);
fbGot  = S.after.flyback;
fprintf('  flyback asked %.1f us -> %.1f us (%d resonant cycle%s, rounded up)\n', ...
        FLYBACK_US, fbGot*1e6, nCyc, repmat('s',1,nCyc>1));
if strcmpi(PIXEL_MODE,'subsample') && abs(S.after.samSlow - samSlow0) > 1e-6
    fprintf('  !! y amplitude moved %.4f -> %.4f despite subsample mode\n', samSlow0, S.after.samSlow);
    ok = false;
end
if strcmpi(PIXEL_MODE,'crop') && abs(S.after.aspect - 1) > 0.02
    fprintf('  !! pixels are %.2f:1, expected square in crop mode\n', S.after.aspect); ok = false;
end
fprintf('  frame rate %.1f Hz (was %.1f)\n', S.after.rate, S.before.rate);
fprintf('  loop delay at 1.46 frames: %.1f ms (was %.1f)\n', ...
        1460/S.after.rate, 1460/S.before.rate);
if ok, fprintf('  geometry applied as requested\n'); end

fprintf('\n*** RE-TAKE THE REFERENCE STACK. *** A %d-line reference cannot\n', S.before.lines);
fprintf('    register %d-line frames. dm_loop_validate takes one inline.\n', S.after.lines);
fprintf('\nto put it back:\n');
fprintf('    hSI.hRoiManager.forceSquarePixelation = %d;\n', S.before.fsPixelation);
fprintf('    hSI.hRoiManager.forceSquarePixels     = %d;\n', S.before.fsPixels);
fprintf('    hSI.hRoiManager.linesPerFrame         = %d;\n', S.before.lines);
fprintf('    hSI.hRoiManager.scanAngleMultiplierSlow = %.6f;\n', S.before.samSlow);
fprintf('    hSI.hScan2D.flybackTimePerFrame       = %.9f;\n\n', S.before.flyback);
end

%% ========================== HELPERS ==========================
function s = snap(hSI)
hR = hSI.hRoiManager; hS = hSI.hScan2D;
s.px          = hR.pixelsPerLine;
s.lines       = hR.linesPerFrame;
s.samSlow     = hR.scanAngleMultiplierSlow;
s.samFast     = hR.scanAngleMultiplierFast;
s.fsPixelation= hR.forceSquarePixelation;
s.fsPixels    = hR.forceSquarePixels;
s.zoom        = hR.scanZoomFactor;
s.flyback     = hS.flybackTimePerFrame;
s.linePeriod  = hR.linePeriod;
s.rate        = hR.scanFrameRate;
% Pixel size from the SCAN FIELD, not from objectiveResolution alone -- and note
% that number is ONE value for both axes, while the axes were measured to differ
% by ~6% (resonant vs galvo). Treat the y figure as nominal.
sf = hR.currentRoiGroup.rois(1).get(0);
s.umPxX = sf.sizeXY(1)*hSI.objectiveResolution / sf.pixelResolutionXY(1);
s.umPxY = sf.sizeXY(2)*hSI.objectiveResolution / sf.pixelResolutionXY(2);
s.aspect = s.umPxY / s.umPxX;
end

function show(s)
fprintf('  frame        %d x %d px   zoom %.2f\n', s.px, s.lines, s.zoom);
fprintf('  pixel        %.3f x %.3f um   (%.2f:1)\n', s.umPxX, s.umPxY, s.aspect);
fprintf('  field        %.0f x %.0f um\n', s.px*s.umPxX, s.lines*s.umPxY);
fprintf('  amplitude    fast %.4f  slow %.4f\n', s.samFast, s.samSlow);
fprintf('  square flags pixelation %d  pixels %d\n', s.fsPixelation, s.fsPixels);
fprintf('  line period  %.2f us   flyback %.1f us\n', s.linePeriod*1e6, s.flyback*1e6);
fprintf('  FRAME RATE   %.1f Hz  -> loop delay %.1f ms at 1.46 frames\n', s.rate, 1460/s.rate);
end

function predict(fRes, nLines, fb, before)
lp   = 1/(2*fRes);                       % bidirectional: two lines per cycle
fbQ  = max(ceil(fb*fRes),1)/fRes;        % as ResScan will actually round it
rate = 1/(nLines*lp + fbQ);
fprintf('  %d lines, flyback %.1f us -> %.1f Hz (delay %.1f ms at 1.46 frames)\n', ...
        nLines, fbQ*1e6, rate, 1460/rate);
fprintf('  vs now %.1f Hz: frame period %.2f -> %.2f ms\n', ...
        before.rate, 1000/before.rate, 1000/rate);
end
