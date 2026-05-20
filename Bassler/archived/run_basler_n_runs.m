function run_basler_n_runs()
close all; clc;
% Run Basler acquisitions via pypylon (Python).
% Requires: pip install pypylon opencv-python (in Python 3.14)

%% ------------------- USER SETTINGS --------------------------------------
filePath      = "C:\260407_KA_electro_dorsal_whisking";  % <-- paste your save path here
saveDir       = filePath;
baseName      = "basler";

fps           = 90;                        % Hz
dur_s         = 60;                        % seconds per run
period_s      = 5*60;                      % seconds between run STARTS
nRuns         = 60;                        % TOTAL NUMBER OF ACQUISITIONS

pixelFormat   = "Mono8";                   % "Mono8", "Mono12", "Mono12Packed"
binH          = 3;                         % horizontal binning (1=none, 2, 4)
binV          = 3;                         % vertical   binning (1=none, 2, 4)

%% ------------------- PATHS ----------------------------------------------
pyExe   = "C:\Program Files\Python314\python.exe";   % Python with pypylon
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "basler_acq.py");

commonArgs = sprintf('--fps %g --pixel %s --binH %d --binV %d', ...
    fps, pixelFormat, binH, binV);

%% ------------------- PREVIEW --------------------------------------------
fprintf("=== Live Preview ===\n");
fprintf("  'c' = draw crop ROI,  'r' = reset to full frame,  'q' = proceed\n");
cmd = sprintf('"%s" "%s" %s preview', pyExe, pyScript, commonArgs);
fprintf(">> %s\n", cmd);
[~, previewOut] = system(cmd);
fprintf("%s", previewOut);

% Parse ROI from preview output
roiArg = "";
roiMatch = regexp(previewOut, 'ROI=(\S+)', 'tokens');
if ~isempty(roiMatch)
    roiStr = roiMatch{1}{1};
    if ~strcmp(roiStr, 'none')
        roiArg = sprintf('--roi %s', roiStr);
        fprintf("Crop ROI: %s\n", roiStr);
    end
end

%% ------------------- ACQUIRE --------------------------------------------
resp = input("Start acquisition? [Y/n]: ", "s");
if ~isempty(resp) && ~strcmpi(resp, 'y')
    fprintf("Aborted.\n");
    return
end

cmd = sprintf('"%s" "%s" %s acquire --saveDir "%s" --baseName %s --dur %g --period %g --nRuns %d %s', ...
    pyExe, pyScript, commonArgs, saveDir, baseName, dur_s, period_s, nRuns, roiArg);
fprintf(">> %s\n", cmd);
system(cmd);

fprintf("\nDone.\n");
end
