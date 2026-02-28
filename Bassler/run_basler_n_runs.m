function run_basler_n_runs()
% Run Basler acquisitions for a fixed number of runs:
%   - record dur_s seconds at fps
%   - start a new run every period_s seconds
%   - stop after nRuns runs total
% Saves a new AVI per run with a timestamped filename.

%% ------------------- USER SETTINGS --------------------------------------
gentlDeviceID = 1;                         % videoinput('gentl', ID)
saveDir       = fullfile(pwd, "basler_acq");
baseName      = "basler";

fps           = 100;                       % Hz
dur_s         = 60;                         % seconds per run (set to 60 for 1 min)
period_s      = 5*60;                    % seconds between run STARTS (set to 5*60 for 5 min)
nRuns         = 50;                         % TOTAL NUMBER OF ACQUISITIONS

roi           = [];                        % [] or [x y w h]
useDisk       = true;

stopDelay_s   = dur_s + 0.5;               % safety margin

%% ------------------- SETUP ----------------------------------------------
if ~exist(saveDir, "dir"), mkdir(saveDir); end

vid = videoinput('gentl', gentlDeviceID);
ctl = baslerController(vid);

triggerconfig(ctl.vid, 'immediate');
ctl.vid.TriggerRepeat = 0;

ctl.params.AcquisitionFrameRate = fps;

% Optional exposure sanity (Basler typically in us)
try
    ctl.params.ExposureTime = min(ctl.params.ExposureTime, 1e6/fps);
catch
end

ctl.params.FramesPerTrigger = round(fps * dur_s);

if ~isempty(roi)
    ctl.params.ROIPosition = roi;
end

bwMBs = ctl.params.BytesPerFrame * ctl.params.AcquisitionFrameRate / 1e6;
fprintf("Estimated bandwidth (Mono8): %.1f MB/s\n", bwMBs);

% Cleanup handler
masterT = [];
stopT   = [];
c = onCleanup(@()cleanup_all()); %#ok<NASGU>

% Run counter
runIdx = 0;

%% ------------------- TIMERS ---------------------------------------------
stopT = timer( ...
    'ExecutionMode','singleShot', ...
    'StartDelay', stopDelay_s, ...
    'TimerFcn', @stop_acq);

masterT = timer( ...
    'ExecutionMode','fixedRate', ...
    'Period', period_s, ...
    'BusyMode','drop', ...
    'TimerFcn', @start_acq);

% Start periodic acquisitions (first one happens immediately)
start_acq();
start(masterT);

fprintf("Running %d runs. Ctrl+C to stop early.\n", nRuns);

%% ------------------- NESTED CALLBACKS -----------------------------------
    function start_acq(~,~)
        if runIdx >= nRuns
            end_sequence();
            return;
        end

        if isrunning(ctl.vid)
            warning("Camera still running; skipping this cycle.");
            return;
        end

        runIdx = runIdx + 1;

        ts = string(datetime("now","Format","yyyyMMdd_HHmmss"));
        fname = fullfile(saveDir, sprintf("%s_%s_run%03d.avi", baseName, ts, runIdx));

        if useDisk
            vw = VideoWriter(fname, 'Grayscale AVI');   % DO NOT open(vw)
            ctl.logger = vw;                            % sets LoggingMode='disk'
        else
            ctl.vid.LoggingMode = 'memory';
        end

        fprintf("[Run %d/%d] START %s\n", runIdx, nRuns, fname);

        start(ctl.vid);

        try, stop(stopT); catch, end
        start(stopT);
    end

    function stop_acq(~,~)
        if ~isrunning(ctl.vid)
            return;
        end

        stop(ctl.vid);

        % Close logger if possible (toolbox often closes it; keep try/catch)
        try
            if isa(ctl.logger,'VideoWriter')
                close(ctl.logger);
            end
        catch
        end

        fprintf("[Run %d/%d] STOP. FramesAcquired=%d, DiskLoggerFrameCount=%d\n", ...
            runIdx, nRuns, ctl.vid.FramesAcquired, ctl.vid.DiskLoggerFrameCount);

        if runIdx >= nRuns
            end_sequence();
        end
    end

    function end_sequence()
        fprintf("\n=== Completed %d runs — stopping. ===\n", nRuns);
        try, if ~isempty(masterT) && isvalid(masterT), stop(masterT); end, catch, end
        try, if ~isempty(stopT)   && isvalid(stopT),   stop(stopT);   end, catch, end
        cleanup_all();
    end

    function cleanup_all()
        fprintf("Cleaning up timers + camera...\n");

        try, if ~isempty(masterT) && isvalid(masterT), delete(masterT); end, end
        try, if ~isempty(stopT)   && isvalid(stopT),   delete(stopT);   end, end

        try
            if isa(ctl.logger,'VideoWriter')
                close(ctl.logger);
            end
        catch
        end

        try
            if isrunning(ctl.vid), stop(ctl.vid); end
        catch
        end

        try
            delete(ctl);
        catch
        end
    end
end
