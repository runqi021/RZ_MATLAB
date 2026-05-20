function N_basler_trig2P()
close all; clc;
% Single-run dual Basler acquisition triggered by 2P scanner.
%
%   1. Arms both cameras (waiting for Line3 triggers)
%   2. Records when triggers arrive (live preview windows)
%   3. Saves and EXITS when triggers stop (no trigger for splitGap_s)
%   Run this script again for the next acquisition.
%
% Use N_Basler_aim.m to preview cameras and get ROI coordinates.

%% =================== USER SETTINGS ======================================
cam1_roi = "525,369,904,850";
cam2_roi = "528,350,512,532";

saveDir    = "C:\test_basler";

splitGap_s = 10;          % seconds of no triggers → done, save and exit

% ---- Camera 1: breathing (1 frame per 2P trigger) ----------------------
cam1_serial   = "22730662"; % acA2040-90umNIR
cam1_name     = "breathing";
cam1_trigger  = "hardware"; % Line3 falling edge → 1 frame
cam1_fps      = 30;         % expected 2P rate
cam1_exposure = 5000;       % microseconds (5ms)
cam1_pixel    = "Mono8";
cam1_binH     = 1;
cam1_binV     = 1;

% ---- Camera 2: whisking (burst of N frames per 2P trigger) -------------
cam2_serial   = "22607866"; % acA1300-200um
cam2_name     = "whisking";
cam2_trigger  = "burst";    % Line3 falling edge → burstN frames at cam2_fps
cam2_fps      = 150;        % Hz, internal burst rate
cam2_burstN   = 5;          % frames per trigger
cam2_exposure = 3000;       % microseconds (3ms — must be short for 150 Hz)
cam2_pixel    = "Mono8";
cam2_binH     = 1;
cam2_binV     = 1;

%% =================== PATHS ==============================================
pyExe     = "C:\Program Files\Python314\python.exe";
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "basler_dual_acq.py");

%% =================== ACQUIRE ============================================
% Build per-camera args
cam1Args = sprintf(['--cam1-serial "%s" --cam1-name %s --cam1-trigger %s ' ...
    '--cam1-fps %g --cam1-exposure %g --cam1-binH %d --cam1-binV %d --cam1-pixel %s'], ...
    cam1_serial, cam1_name, cam1_trigger, cam1_fps, cam1_exposure, cam1_binH, cam1_binV, cam1_pixel);
if strlength(cam1_roi) > 0
    cam1Args = sprintf('%s --cam1-roi %s', cam1Args, cam1_roi);
end

cam2Args = sprintf(['--cam2-serial "%s" --cam2-name %s --cam2-trigger %s ' ...
    '--cam2-fps %g --cam2-burst %d --cam2-exposure %g --cam2-binH %d --cam2-binV %d --cam2-pixel %s'], ...
    cam2_serial, cam2_name, cam2_trigger, cam2_fps, cam2_burstN, cam2_exposure, ...
    cam2_binH, cam2_binV, cam2_pixel);
if strlength(cam2_roi) > 0
    cam2Args = sprintf('%s --cam2-roi %s', cam2Args, cam2_roi);
end

cmd = sprintf('"%s" "%s" acquire --saveDir "%s" --split-gap %g --single-run %s %s', ...
    pyExe, pyScript, saveDir, splitGap_s, cam1Args, cam2Args);

fprintf("Arming cameras... (waiting for 2P triggers on Line3)\n");
fprintf(">> %s\n", cmd);
status = system(cmd);

if status ~= 0
    fprintf("\n*** Error (code %d). Check output above. ***\n", status);
else
    fprintf("\nDone. Run again for next acquisition.\n");
end
end
