function run_basler_n_runs_n_cam()
close all; clc;
% Multi-camera Basler N-runs free-run acquisition (untriggered).
% Calls basler_dual_acq.py with the acquire-nruns subcommand.
% RAM-first lossless pipeline: pylon transport buffer -> queue.Queue
% -> background writer_thread -> cv2.VideoWriter FFV1 AVI (lossless).
%
% Each enabled camera runs free-run at its requested fps. For each run:
%   open AVI -> grab dur_s seconds -> close AVI + save timestamps.csv ->
%   idle until period_s has elapsed since run start -> next run.
% Total nRuns runs.
%
% Output layout (matches basler_gui.m):
%   <saveDir>/<camName>/<camName>_<YYYYMMDD_HHMMSS>_runNNN/
%       <camName>_<YYYYMMDD_HHMMSS>_runNNN.avi
%       timestamps.csv

%% ------------------- USER SETTINGS --------------------------------------
saveDir       = "C:\Users\Admin\Desktop\260505_breathing_wt";

dur_s         = 60;            % seconds per run
period_s      = 5*60;          % seconds between RUN STARTS
nRuns         = 200;            % total number of runs

% Per-camera config (set enable=false to skip a slot).
% serial: "" -> auto-pick first available; or paste exact SN string.
% exposure: 0 -> auto-cap to ~90% of frame period.
% roi: leave "" -- the preview step writes it after you press 'q'.
cams(1).enable    = true;
cams(1).serial    = "";
cams(1).name      = "cam1";
cams(1).fps       = 100;
cams(1).pixel     = "Mono8";
cams(1).binH      = 3;
cams(1).binV      = 3;
cams(1).exposure  = 0;
cams(1).roi       = "";

cams(2).enable    = true;
cams(2).serial    = "";
cams(2).name      = "cam2";
cams(2).fps       = 100;
cams(2).pixel     = "Mono8";
cams(2).binH      = 2;
cams(2).binV      = 2;
cams(2).exposure  = 0;
cams(2).roi       = "";

cams(3).enable    = false;
cams(3).serial    = "";
cams(3).name      = "cam3";
cams(3).fps       = 90;
cams(3).pixel     = "Mono8";
cams(3).binH      = 1;
cams(3).binV      = 1;
cams(3).exposure  = 0;
cams(3).roi       = "";

%% ------------------- PATHS ----------------------------------------------
pyExe     = "C:\Program Files\Python314\python.exe";
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "basler_dual_acq.py");

%% ------------------- ENUMERATE CAMERAS ----------------------------------
fprintf("=== Enumerating Basler cameras ===\n");
[~, listOut] = system(sprintf('"%s" "%s" list', pyExe, pyScript));
fprintf("%s\n", listOut);

% Parse all serial numbers from the list output (one per "SN=..." line).
snTokens = regexp(string(listOut), 'SN=(\S+)', 'tokens');
availSerials = string.empty;
for s = 1:numel(snTokens)
    availSerials(end+1) = string(snTokens{s}{1}); %#ok<SAGROW>
end
fprintf("[enum] Parsed %d serial(s): %s\n", numel(availSerials), strjoin(availSerials, ", "));

% Auto-assign serials to any enabled cam with empty serial, in order.
% Without this, every preview call would pick devices[0] and you'd see the
% same camera previewed N times (acquire is fine — it tracks 'used').
usedSerials = strings(0);
for k = 1:numel(cams)
    if ~cams(k).enable, continue; end
    if strlength(cams(k).serial) > 0
        usedSerials(end+1) = cams(k).serial; %#ok<SAGROW>
    end
end
for k = 1:numel(cams)
    if ~cams(k).enable || strlength(cams(k).serial) > 0, continue; end
    picked = "";
    for s = 1:numel(availSerials)
        if ~any(availSerials(s) == usedSerials)
            picked = availSerials(s);
            break;
        end
    end
    if picked == ""
        fprintf(2, "WARNING: no free camera to auto-assign to %s — preview/acquire may fail.\n", cams(k).name);
        continue;
    end
    cams(k).serial = picked;
    usedSerials(end+1) = picked; %#ok<SAGROW>
    fprintf("[enum] Auto-assigned %s -> SN=%s\n", cams(k).name, picked);
end

%% ------------------- PREVIEW EACH ENABLED CAMERA ------------------------
fprintf("=== Live Preview ===\n");
fprintf("  'c' = draw crop ROI    'r' = reset to full    'q' = accept\n");

for k = 1:numel(cams)
    if ~cams(k).enable, continue; end
    fprintf("\n--- Preview %s (SN=%s) ---\n", cams(k).name, cams(k).serial);
    serialArg = "";
    if strlength(cams(k).serial) > 0
        serialArg = sprintf('--serial "%s"', cams(k).serial);
    end
    cmd = sprintf('"%s" "%s" preview %s --fps %g --binH %d --binV %d --pixel %s', ...
        pyExe, pyScript, serialArg, cams(k).fps, ...
        cams(k).binH, cams(k).binV, cams(k).pixel);
    fprintf(">> %s\n", cmd);
    [~, out] = system(cmd);
    fprintf("%s", out);

    % Parse ROI from preview output
    roiMatch = regexp(out, 'ROI=(\S+)', 'tokens');
    if ~isempty(roiMatch)
        roiStr = roiMatch{1}{1};
        if ~strcmp(roiStr, 'none')
            cams(k).roi = string(roiStr);
            fprintf("  %s ROI: %s\n", cams(k).name, roiStr);
        else
            cams(k).roi = "";
        end
    end
end

%% ------------------- CONFIRM PLAN ---------------------------------------
fprintf("\n=== Acquisition Plan ===\n");
fprintf("  saveDir : %s\n", saveDir);
fprintf("  dur     : %g s  (%g min)\n", dur_s, dur_s/60);
fprintf("  period  : %g s  (%g min, start-to-start)\n", period_s, period_s/60);
totalSec = period_s * (nRuns - 1) + dur_s;
fprintf("  nRuns   : %d   (total session ~%.1f min)\n", nRuns, totalSec/60);
fprintf("  Cameras :\n");
nEnabled = 0;
for k = 1:numel(cams)
    if ~cams(k).enable, continue; end
    nEnabled = nEnabled + 1;
    roiTxt = "full";
    if strlength(cams(k).roi) > 0, roiTxt = cams(k).roi; end
    snTxt = cams(k).serial;
    if strlength(snTxt) == 0, snTxt = "AUTO"; end
    fprintf("    [%d] %-12s  SN=%s  %g Hz  bin %dx%d  pix=%s  ROI=%s\n", ...
        nEnabled, cams(k).name, snTxt, cams(k).fps, ...
        cams(k).binH, cams(k).binV, cams(k).pixel, roiTxt);
end
if nEnabled == 0
    fprintf("\nERROR: no cameras enabled. Exiting.\n"); return
end

resp = input(sprintf('\nStart acquisition? [Y/n]: '), "s");
if ~isempty(resp) && ~strcmpi(resp, 'y')
    fprintf("Aborted.\n"); return
end

%% ------------------- ACQUIRE --------------------------------------------
camArgs = "";
slot = 0;
for k = 1:numel(cams)
    if ~cams(k).enable, continue; end
    slot = slot + 1;
    c = cams(k);
    serial = c.serial;
    if strlength(serial) == 0, serial = "AUTO"; end
    camArgs = sprintf(['%s --cam%d-serial "%s" --cam%d-name %s ' ...
        '--cam%d-fps %g --cam%d-pixel %s ' ...
        '--cam%d-binH %d --cam%d-binV %d --cam%d-exposure %g'], ...
        camArgs, slot, serial, slot, c.name, ...
        slot, c.fps, slot, c.pixel, ...
        slot, c.binH, slot, c.binV, slot, c.exposure);
    if strlength(c.roi) > 0
        camArgs = sprintf('%s --cam%d-roi %s', camArgs, slot, c.roi);
    end
end

cmd = sprintf('"%s" "%s" acquire-nruns --saveDir "%s" --dur %g --period %g --nRuns %d %s', ...
    pyExe, pyScript, saveDir, dur_s, period_s, nRuns, camArgs);
fprintf("\n>> %s\n", cmd);
system(cmd);

fprintf("\nDone.\n");
end
