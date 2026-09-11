function breath_fixedmetric_to_gui(rootDir, metric)
% breath_fixedmetric_to_gui  Expose the fixed-metric traces to breathing_peak_gui_pc1.
%
% Writes one breath_pc1.mat into each run folder so the existing peak GUI can
% open these traces with no change to the GUI itself.  breathing_peak_gui_pc1.m
% scans recursively for breath_pc1.mat and needs only the breathTrace field
% (t / fps / eigImg / varExp / pc are optional and are filled in here).
%
% WHAT THE GUI DOES TO THE TRACE
%   breathing_peak_gui_pc1.m detrends and z-scores whatever it loads (its lines
%   208-209) so that one set of peak thresholds works everywhere.  That is the
%   right behaviour for picking peaks and the wrong basis for comparing breath
%   depth -- z-scoring is exactly the normalisation the fixed-metric pipeline
%   exists to avoid.  So: use the GUI to check DETECTION and TIMING run by run,
%   and use breath_fixedmetric_analyze.m for amplitude.  The numbers you see on
%   the GUI's y-axis are not comparable across runs; the ones in the analyze
%   figures are.
%
% USAGE
%   breath_fixedmetric_to_gui                                  % disp, default root
%   breath_fixedmetric_to_gui('D:\breath_tracking_test','fb')  % frozen-basis PC1
%
%   then in MATLAB:  breathing_peak_gui_pc1
%   and point "Browse Master Folder" at ROOT_DIR -- it will list all runs.

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR  = 'D:\260728_vglut2_soma-g8s\phys';
METRIC    = 'disp';    % 'disp' (rigid displacement, px) | 'fb' | 'pv'
OVERWRITE = false;     % false = refuse to clobber an existing breath_pc1.mat
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT_DIR = rootDir; end
if nargin >= 2 && ~isempty(metric),  METRIC   = metric;  end

matFile = fullfile(ROOT_DIR, 'breath_fixedmetric.mat');
assert(isfile(matFile), ...
    'breath_fixedmetric.mat not found in %s -- run breath_fixedmetric_extract.py first.', ROOT_DIR);
S = load(matFile);
assert(isfield(S, METRIC), 'metric "%s" not present in breath_fixedmetric.mat', METRIC);

% Prefer the FULL-LENGTH trace when the extractor provides it.
%
% disp/fb/pv are all truncated to T = the shortest run, which is what makes their
% amplitudes comparable BETWEEN runs -- correct for the fixed-metric comparison,
% wrong for the trace exported here. Every downstream use of breath_pc1.mat is
% temporal (peak/trough GUIs, breath phase), and the phase pipelines take
% T = min(numel(breath), nCaFrames), so a clipped trace silently truncates the
% whole analysis. On 260728_vglut2 that cost 150 of 200 s -- 25% of each recording.
%
% disp_full / fb_full are NaN-padded to Tmax; run_len marks each run's true end.
fullName = [METRIC '_full'];
useFull  = isfield(S, fullName) && isfield(S, 'run_len');
if useFull
    X       = double(S.(fullName));
    runLen  = double(S.run_len(:));
elseif isfield(S, [METRIC '_full'])
    error('%s present but run_len missing -- re-run breath_fixedmetric_extract.py', fullName);
else
    X       = double(S.(METRIC));
    runLen  = repmat(size(X,1), size(X,2), 1);
    warning(['no %s in this .mat -- falling back to the run-truncated trace.\n' ...
             '         Re-run breath_fixedmetric_extract.py to get full-length traces.'], fullName);
end
% PER-RUN fps when the extractor provides it. The breath cam is 2P-triggered, so
% runs with different scan configs run at different rates (30-50 Hz on
% 260728_vglut2); writing the session median into every breath_pc1.mat would
% rescale most of them. fps_run is the effective (post-stride) rate per run.
fpsSession = double(S.fps);           % median; header prints only
if isfield(S, 'fps_run')
    fpsRun = double(S.fps_run(:));
else
    fpsRun = [];
end
fps = fpsSession;                     % overwritten per run inside the loop, and
                                      % it is `fps` that gets SAVED into each file
runNames = cellstr(string(S.run_names(:)));
isBase   = logical(S.is_baseline(:));
nR       = numel(runNames);

unitStr = 'a.u.';
if strcmp(METRIC, 'disp'), unitStr = 'px'; end
fprintf('breath_fixedmetric_to_gui: metric "%s" (%s), %d runs, %.4f fps\n', ...
    METRIC, unitStr, nR, fpsSession);
if useFull
    fprintf('  full-length traces: %d-%d frames (%.1f-%.1f s) per run\n', ...
        min(runLen), max(runLen), min(runLen)/fpsSession, max(runLen)/fpsSession);
    if isfield(S,'T')
        fprintf('  (the run-truncated arrays used for cross-run comparison are %d frames)\n', ...
                round(double(S.T)));
    end
end

nWrote = 0; nSkip = 0;
for j = 1:nR
    folder = fullfile(ROOT_DIR, runNames{j});
    if ~isfolder(folder)
        fprintf(2, '  missing folder, skipped: %s\n', runNames{j});
        nSkip = nSkip + 1;
        continue;
    end
    outFile = fullfile(folder, 'breath_pc1.mat');
    % Refuse only files this pipeline did NOT write -- i.e. genuine
    % breath_svd_pc1.m output, which uses the same filename and would be
    % destroyed.  Our own files carry src_metric and are safe to regenerate,
    % otherwise every re-run or metric switch would need OVERWRITE set by hand.
    if isfile(outFile) && ~OVERWRITE
        try
            ours = isfield(load(outFile, 'src_metric'), 'src_metric');
        catch
            ours = false;
        end
        if ~ours
            fprintf(2, ['  refusing to overwrite a file this pipeline did not write\n' ...
                '    (looks like breath_svd_pc1.m output; set OVERWRITE=true to replace): %s\n'], outFile);
            nSkip = nSkip + 1;
            continue;
        end
    end

    % Trim the NaN pad back to this run's own length, and use THIS run's rate.
    breathTrace = X(1:min(runLen(j), size(X,1)), j);
    breathTrace = breathTrace(~isnan(breathTrace));
    if isempty(fpsRun), fps = fpsSession; else, fps = fpsRun(j); end
    t           = (0:numel(breathTrace)-1)' / fps;
    eigImg      = double(S.u1_img);
    % Inspiration-minus-expiration map. Without this the GUI's image panel is
    % blank: compute_diff_map wants either diffImg or the full U/V/sv factors,
    % and this pipeline stores neither the factors nor the frames in the .mat.
    if isfield(S, 'diff_imgs')
        diffImg = double(S.diff_imgs(:, :, j));
    else
        diffImg = [];
    end
    varExp      = double(S.axis_var(:))';
    pc          = 1;
    src_metric  = METRIC;
    src_file    = matFile;
    is_baseline = isBase(j);

    save(outFile, 'breathTrace', 't', 'fps', 'eigImg', 'diffImg', 'varExp', 'pc', ...
        'src_metric', 'src_file', 'is_baseline');
    nWrote = nWrote + 1;
end

fprintf('  wrote %d breath_pc1.mat, skipped %d\n', nWrote, nSkip);
fprintf('\nnext:  breathing_peak_gui_pc1\n');
fprintf('       Browse Master Folder -> %s\n', ROOT_DIR);
fprintf('\nNOTE: the GUI z-scores each trace, so its y-axis is NOT comparable\n');
fprintf('      across runs. Use it for peak/timing QC; use\n');
fprintf('      breath_fixedmetric_analyze.m for amplitude.\n');

end
