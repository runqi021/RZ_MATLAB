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
OVERWRITE = true;     % false = refuse to clobber an existing breath_pc1.mat
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT_DIR = rootDir; end
if nargin >= 2 && ~isempty(metric),  METRIC   = metric;  end

matFile = fullfile(ROOT_DIR, 'breath_fixedmetric.mat');
assert(isfile(matFile), ...
    'breath_fixedmetric.mat not found in %s -- run breath_fixedmetric_extract.py first.', ROOT_DIR);
S = load(matFile);
assert(isfield(S, METRIC), 'metric "%s" not present in breath_fixedmetric.mat', METRIC);

X        = double(S.(METRIC));
fps      = double(S.fps);
runNames = cellstr(string(S.run_names(:)));
isBase   = logical(S.is_baseline(:));
nR       = numel(runNames);
t        = (0:size(X, 1) - 1)' / fps;

unitStr = 'a.u.';
if strcmp(METRIC, 'disp'), unitStr = 'px'; end
fprintf('breath_fixedmetric_to_gui: metric "%s" (%s), %d runs, %.4f fps\n', ...
    METRIC, unitStr, nR, fps);

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

    breathTrace = X(:, j);
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
