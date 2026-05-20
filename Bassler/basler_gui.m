function basler_gui()
% Dual Basler camera GUI — controls acquisition via basler_dual_acq.py.
% Preview pops up as OpenCV windows. All parameters editable in GUI.

%% ─── Python path ────────────────────────────────────────────────────────
pyExe    = "C:\Program Files\Python314\python.exe";
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "basler_dual_acq.py");

%% ─── Enumerate cameras ─────────────────────────────────────────────────
[~, out] = system(sprintf('"%s" "%s" list', pyExe, pyScript));
serials = regexp(out, 'SN=(\d+)', 'tokens');
models  = regexp(out, 'Model=(\S+)', 'tokens');
sn1 = ""; sn2 = ""; sn3 = ""; mdl1 = ""; mdl2 = ""; mdl3 = "";
if numel(serials) >= 1, sn1 = string(serials{1}{1}); mdl1 = string(models{1}{1}); end
if numel(serials) >= 2, sn2 = string(serials{2}{1}); mdl2 = string(models{2}{1}); end
if numel(serials) >= 3, sn3 = string(serials{3}{1}); mdl3 = string(models{3}{1}); end

%% ─── State ──────────────────────────────────────────────────────────────
state.runIdx = 0;
state.acquiring = false;
state.process = [];

%% ─── Main figure ────────────────────────────────────────────────────────
fig = uifigure('Name', 'Basler Multi-Camera Acquisition', ...
    'Position', [100 100 840 410], 'CloseRequestFcn', @onClose);

%% ─── Camera 1 panel ────────────────────────────────────────────────────
p1 = uipanel(fig, 'Title', sprintf('Cam1  SN:%s (%s)', sn1, mdl1), ...
    'Position', [5 195 275 210]);
c1 = createCameraPanel(p1, "cam1", sn1, ...
    "burst", 10, 300, 2000, 1, 1, "", sn1 ~= "");

%% ─── Camera 2 panel ────────────────────────────────────────────────────
p2 = uipanel(fig, 'Title', sprintf('Cam2  SN:%s (%s)', sn2, mdl2), ...
    'Position', [282 195 275 210]);
c2 = createCameraPanel(p2, "cam2", sn2, ...
    "hardware", 1, 30, 5000, 1, 1, "", sn2 ~= "");

%% ─── Camera 3 panel ────────────────────────────────────────────────────
p3 = uipanel(fig, 'Title', sprintf('Cam3  SN:%s (%s)', sn3, mdl3), ...
    'Position', [559 195 275 210]);
c3 = createCameraPanel(p3, "cam3", sn3, ...
    "hardware", 1, 30, 5000, 1, 1, "", sn3 ~= "");

%% ─── Bottom controls ───────────────────────────────────────────────────
bp = uipanel(fig, 'Position', [5 5 829 185]);

uilabel(bp, 'Text', 'Save dir:', 'Position', [5 150 55 20]);
saveDirField = uieditfield(bp, 'text', 'Value', 'C:\test_basler', ...
    'Position', [60 150 700 20]);
uibutton(bp, 'Text', '...', 'Position', [765 150 30 20], ...
    'ButtonPushedFcn', @(~,~) browseSaveDir());

uilabel(bp, 'Text', 'Gap:', 'Position', [5 125 30 20]);
splitGapField = uieditfield(bp, 'numeric', 'Value', 5, ...
    'Position', [38 125 40 20]);
uilabel(bp, 'Text', 's', 'Position', [80 125 15 20]);

runLabel = uilabel(bp, 'Text', 'Run: 000', ...
    'Position', [735 80 70 22], 'FontWeight', 'bold');

previewBtn = uibutton(bp, 'Text', 'PREVIEW', ...
    'Position', [10 75 120 35], 'FontSize', 13, ...
    'ButtonPushedFcn', @(~,~) onPreview());
recordBtn = uibutton(bp, 'Text', 'ARM & RECORD', ...
    'Position', [140 75 180 35], 'FontSize', 13, ...
    'BackgroundColor', [0.3 0.8 0.3], ...
    'ButtonPushedFcn', @(~,~) onRecord());
stopBtn = uibutton(bp, 'Text', 'STOP', ...
    'Position', [330 75 120 35], 'FontSize', 13, ...
    'BackgroundColor', [0.9 0.3 0.3], 'Enable', 'off', ...
    'ButtonPushedFcn', @(~,~) onStop());

statusLabel = uilabel(bp, 'Text', 'Status: Idle', ...
    'Position', [5 10 820 25], 'FontSize', 12);

%% ─── Timer for monitoring background acquisition ───────────────────────
pollTimer = timer('ExecutionMode', 'fixedRate', 'Period', 1, ...
    'TimerFcn', @(~,~) pollAcquisition());

%% ═══════════════════════════════════════════════════════════════════════
%  NESTED FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════

function controls = createCameraPanel(parent, tag, serial, ...
        trigDefault, burstDefault, fpsDefault, expDefault, ...
        binHDefault, binVDefault, roiDefault, enableDefault)
    y = 170;
    dy = 24;

    controls.enable = uicheckbox(parent, 'Text', 'Enable', ...
        'Value', enableDefault, 'Position', [5 y 65 20]);
    uilabel(parent, 'Text', 'Name:', 'Position', [75 y 38 20]);
    controls.name = uieditfield(parent, 'text', 'Value', tag, ...
        'Position', [115 y 145 20]);
    y = y - dy;

    uilabel(parent, 'Text', 'Trigger:', 'Position', [5 y 45 20]);
    controls.trigger = uidropdown(parent, ...
        'Items', {'hardware','burst'}, ...
        'Value', trigDefault, ...
        'Position', [52 y 80 20]);
    uilabel(parent, 'Text', 'Edge:', 'Position', [137 y 32 20]);
    controls.edge = uidropdown(parent, ...
        'Items', {'falling','rising'}, ...
        'Value', 'falling', ...
        'Position', [170 y 80 20]);
    y = y - dy;

    uilabel(parent, 'Text', 'Burst N:', 'Position', [5 y 48 20]);
    controls.burstN = uieditfield(parent, 'numeric', 'Value', burstDefault, ...
        'Position', [55 y 45 20]);
    uilabel(parent, 'Text', 'FPS:', 'Position', [110 y 30 20]);
    controls.fps = uieditfield(parent, 'numeric', 'Value', fpsDefault, ...
        'Position', [142 y 50 20]);
    y = y - dy;

    uilabel(parent, 'Text', 'Exp:', 'Position', [5 y 30 20]);
    controls.exposure = uieditfield(parent, 'numeric', 'Value', expDefault, ...
        'Position', [37 y 55 20]);
    uilabel(parent, 'Text', 'µs', 'Position', [94 y 18 20]);
    uilabel(parent, 'Text', 'binH:', 'Position', [118 y 30 20]);
    controls.binH = uieditfield(parent, 'numeric', 'Value', binHDefault, ...
        'Position', [150 y 30 20]);
    uilabel(parent, 'Text', 'binV:', 'Position', [185 y 30 20]);
    controls.binV = uieditfield(parent, 'numeric', 'Value', binVDefault, ...
        'Position', [217 y 30 20]);
    y = y - dy;

    uilabel(parent, 'Text', 'ROI:', 'Position', [5 y 28 20]);
    controls.roi = uieditfield(parent, 'text', 'Value', roiDefault, ...
        'Position', [35 y 100 20]);
    controls.resetRoiBtn = uibutton(parent, 'Text', 'Full', ...
        'Position', [138 y 32 20], ...
        'ButtonPushedFcn', @(~,~) set(controls.roi, 'Value', ''));
    controls.previewBtn = uibutton(parent, 'Text', 'Preview', ...
        'Position', [175 y 70 20], ...
        'ButtonPushedFcn', @(~,~) previewOne(serial, controls));
    y = y - dy;

    controls.frameLabel = uilabel(parent, 'Text', 'Frames: 0', ...
        'Position', [5 y 100 20], 'FontWeight', 'bold');
    controls.bufferLabel = uilabel(parent, 'Text', 'Buffer: 0', ...
        'Position', [110 y 100 20]);
    y = y - 20;

    controls.statusLabel = uilabel(parent, 'Text', 'Idle', ...
        'Position', [5 y 240 20], 'FontColor', [0.4 0.4 0.4]);

    controls.serial = serial;

    % Callbacks must be wired AFTER all fields exist — anonymous-function
    % closures capture `controls` by value, so a callback defined earlier
    % would see a stale snapshot missing later fields.
    controls.trigger.ValueChangedFcn = @(dd,~) onTriggerChange(dd, controls);
    controls.enable.ValueChangedFcn  = @(cb,~) onEnableChange(cb, controls);

    % Apply initial enable + trigger-dependent graying
    onEnableChange(controls.enable, controls);
end

function onTriggerChange(dd, controls)
    if ~controls.enable.Value
        return;
    end
    % FPS is editable in both modes:
    %   hardware -> used to cap exposure + stamp AVI metadata at 2P rate
    %   burst    -> internal frame rate of each burst
    controls.fps.Enable = 'on';
    if strcmp(dd.Value, 'burst')
        controls.burstN.Enable = 'on';
    else % hardware
        controls.burstN.Enable = 'off';
    end
end

function onEnableChange(cb, controls)
    tf = cb.Value;
    v = 'off';
    if tf, v = 'on'; end
    controls.name.Enable       = v;
    controls.trigger.Enable    = v;
    controls.edge.Enable       = v;
    controls.exposure.Enable   = v;
    controls.binH.Enable       = v;
    controls.binV.Enable       = v;
    controls.roi.Enable        = v;
    controls.previewBtn.Enable = v;
    if tf
        onTriggerChange(controls.trigger, controls);
    else
        controls.burstN.Enable = 'off';
        controls.fps.Enable    = 'off';
    end
end

function previewOne(serial, controls)
    fps = controls.fps.Value;
    if strcmp(controls.trigger.Value, 'hardware')
        fps = 30; % just for preview
    end
    cmd = sprintf('"%s" "%s" preview --serial "%s" --fps %g --binH %d --binV %d --pixel Mono8', ...
        pyExe, pyScript, serial, fps, controls.binH.Value, controls.binV.Value);
    statusLabel.Text = sprintf('Status: Previewing %s ...', controls.name.Value);
    drawnow;
    [~, out] = system(cmd);

    % Parse ROI from output
    roiMatch = regexp(out, 'ROI=(\S+)', 'tokens');
    if ~isempty(roiMatch) && ~strcmp(roiMatch{1}{1}, 'none')
        controls.roi.Value = roiMatch{1}{1};
    end
    statusLabel.Text = 'Status: Idle';
end

function browseSaveDir()
    d = uigetdir(saveDirField.Value, 'Select save directory');
    if d ~= 0
        saveDirField.Value = d;
    end
end

function args = buildCamArgs(controls, camNum)
    c = controls;
    args = sprintf(['--cam%d-serial "%s" --cam%d-name %s --cam%d-trigger %s ' ...
        '--cam%d-edge %s --cam%d-fps %g --cam%d-exposure %g ' ...
        '--cam%d-binH %d --cam%d-binV %d --cam%d-pixel Mono8'], ...
        camNum, c.serial, camNum, c.name.Value, camNum, c.trigger.Value, ...
        camNum, c.edge.Value, camNum, c.fps.Value, camNum, c.exposure.Value, ...
        camNum, c.binH.Value, camNum, c.binV.Value, camNum);
    if strcmp(c.trigger.Value, 'burst')
        args = sprintf('%s --cam%d-burst %d', args, camNum, c.burstN.Value);
    end
    if strlength(c.roi.Value) > 0
        args = sprintf('%s --cam%d-roi %s', args, camNum, c.roi.Value);
    end
end

function onPreview()
    statusLabel.Text = 'Status: Preview — enabled cams in order. Press q in each window.';
    drawnow;
    if c1.enable.Value && strlength(c1.serial) > 0
        previewOne(c1.serial, c1);
    end
    if c2.enable.Value && strlength(c2.serial) > 0
        previewOne(c2.serial, c2);
    end
    if c3.enable.Value && strlength(c3.serial) > 0
        previewOne(c3.serial, c3);
    end
    statusLabel.Text = 'Status: Idle';
end

function onRecord()
    if state.acquiring
        return;
    end

    % Collect enabled cameras; pass first enabled as --cam1-*, second as --cam2-*, etc.
    selected = {};
    if c1.enable.Value && strlength(c1.serial) > 0, selected{end+1} = c1; end
    if c2.enable.Value && strlength(c2.serial) > 0, selected{end+1} = c2; end
    if c3.enable.Value && strlength(c3.serial) > 0, selected{end+1} = c3; end

    if isempty(selected)
        statusLabel.Text = 'Status: No cameras enabled — tick Enable on at least one panel.';
        return;
    end

    state.runIdx = state.runIdx + 1;
    runLabel.Text = sprintf('Run: %03d', state.runIdx);

    camArgs = "";
    for k = 1:numel(selected)
        camArgs = sprintf('%s %s', camArgs, buildCamArgs(selected{k}, k));
    end

    cmd = sprintf('"%s" "%s" acquire --saveDir "%s" --split-gap %g --single-run %s', ...
        pyExe, pyScript, saveDirField.Value, splitGapField.Value, camArgs);

    % Lock controls
    setControlsEnabled(false);
    state.acquiring = true;
    state.pidFile = fullfile(tempdir, 'basler_acq.pid');
    state.pidSeen = false;
    state.pid = "";
    if isfile(state.pidFile), delete(state.pidFile); end
    statusLabel.Text = 'Status: Armed — waiting for triggers on Line3 ...';
    if c1.enable.Value, c1.statusLabel.Text = 'Armed'; end
    if c2.enable.Value, c2.statusLabel.Text = 'Armed'; end
    if c3.enable.Value, c3.statusLabel.Text = 'Armed'; end
    drawnow;

    % Launch Python in a new console window (returns immediately).
    % The console shows all diagnostic output (Line3 state, frame counts).
    system(sprintf('start "Basler Acquisition" %s', cmd));
    start(pollTimer);
end

function onStop()
    if state.acquiring
        % Kill Python process by PID
        if isfile(state.pidFile)
            pid = strip(fileread(state.pidFile));
            system(sprintf('taskkill /F /PID %s 2>nul', pid));
        end
        stopAcquisition('Stopped by user');
    end
end

function pollAcquisition()
    if ~state.acquiring, return; end
    if ~state.pidSeen
        % Wait for PID file to appear (Python writes it at startup)
        if ~isfile(state.pidFile), return; end
        state.pid = strip(fileread(state.pidFile));
        state.pidSeen = true;
        return;
    end
    % PID was seen before — check if process is still alive
    [~, out] = system(sprintf('tasklist /FI "PID eq %s" /NH 2>nul', state.pid));
    if ~contains(out, 'python', 'IgnoreCase', true)
        stopAcquisition('Done');
    end
end

function parseOutput(output)
    % Extract frame counts and fps from Python output
    % e.g. "[breathing] Run 1 saved: 6000 frames  120 MB  30.0 Hz"
    matches = regexp(output, '\[(\w+)\] Run \d+ saved: (\d+) frames\s+\d+ MB\s+([\d.]+) Hz', 'tokens');
    for i = 1:numel(matches)
        camName = matches{i}{1};
        nFrames = matches{i}{2};
        fps = matches{i}{3};
        if strcmp(camName, c1.name.Value)
            c1.frameLabel.Text = sprintf('Frames: %s', nFrames);
            c1.statusLabel.Text = sprintf('Saved @ %s Hz', fps);
        elseif strcmp(camName, c2.name.Value)
            c2.frameLabel.Text = sprintf('Frames: %s', nFrames);
            c2.statusLabel.Text = sprintf('Saved @ %s Hz', fps);
        elseif strcmp(camName, c3.name.Value)
            c3.frameLabel.Text = sprintf('Frames: %s', nFrames);
            c3.statusLabel.Text = sprintf('Saved @ %s Hz', fps);
        end
    end
end

function stopAcquisition(msg)
    stop(pollTimer);
    state.acquiring = false;
    if isfile(state.pidFile), delete(state.pidFile); end
    setControlsEnabled(true);
    statusLabel.Text = sprintf('Status: %s  (Run %03d)', msg, state.runIdx);
end

function setControlsEnabled(tf)
    v = 'on';
    if ~tf, v = 'off'; end
    previewBtn.Enable = v;
    recordBtn.Enable = v;
    saveDirField.Enable = v;
    splitGapField.Enable = v;
    if tf
        stopBtn.Enable = 'off';
    else
        stopBtn.Enable = 'on';
    end
end

function onClose(~, ~)
    if state.acquiring
        onStop();
    end
    try stop(pollTimer); catch, end
    try delete(pollTimer); catch, end
    delete(fig);
end

end
