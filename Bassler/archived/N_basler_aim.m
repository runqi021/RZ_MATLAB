function N_basler_aim()
close all; clc;
% Preview both Basler cameras and draw crop ROI.
% Outputs ROI coordinates to paste into N_Basler_trig2P.m
%
% For each camera:
%   - Live feed window opens
%   - Press 'c' to draw crop rectangle, then Enter to confirm
%   - Press 'r' to reset to full frame
%   - Press 'q' to accept and move to next camera
%
% At the end, prints the ROI strings ready to paste.

%% =================== SETTINGS ===========================================
cam1_serial = "22730662";  % acA2040-90umNIR (breathing)
cam1_binH   = 1;
cam1_binV   = 1;
cam1_pixel  = "Mono8";

cam2_serial = "22607866";  % acA1300-200um (whisking)
cam2_binH   = 1;
cam2_binV   = 1;
cam2_pixel  = "Mono8";

%% =================== PATHS ==============================================
pyExe     = "C:\Program Files\Python314\python.exe";
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "basler_dual_acq.py");

%% =================== PREVIEW CAM 1 ======================================
fprintf("=== Camera 1: breathing (SN=%s) ===\n", cam1_serial);
fprintf("  'c' = draw crop    'r' = reset    'q' = done\n\n");
cmd = sprintf('"%s" "%s" preview --serial "%s" --fps 30 --binH %d --binV %d --pixel %s', ...
    pyExe, pyScript, cam1_serial, cam1_binH, cam1_binV, cam1_pixel);
[~, out] = system(cmd);
fprintf("%s\n", out);

cam1_roi = "";
roiMatch = regexp(out, 'ROI=(\S+)', 'tokens');
if ~isempty(roiMatch) && ~strcmp(roiMatch{1}{1}, 'none')
    cam1_roi = roiMatch{1}{1};
end

%% =================== PREVIEW CAM 2 ======================================
fprintf("\n=== Camera 2: whisking (SN=%s) ===\n", cam2_serial);
fprintf("  'c' = draw crop    'r' = reset    'q' = done\n\n");
cmd = sprintf('"%s" "%s" preview --serial "%s" --fps 150 --binH %d --binV %d --pixel %s', ...
    pyExe, pyScript, cam2_serial, cam2_binH, cam2_binV, cam2_pixel);
[~, out] = system(cmd);
fprintf("%s\n", out);

cam2_roi = "";
roiMatch = regexp(out, 'ROI=(\S+)', 'tokens');
if ~isempty(roiMatch) && ~strcmp(roiMatch{1}{1}, 'none')
    cam2_roi = roiMatch{1}{1};
end

%% =================== OUTPUT =============================================
fprintf("\n");
fprintf("========================================\n");
fprintf("  Copy these into N_Basler_trig2P.m:\n");
fprintf("========================================\n");
fprintf('cam1_roi = "%s";\n', cam1_roi);
fprintf('cam2_roi = "%s";\n', cam2_roi);
fprintf("========================================\n");
end
