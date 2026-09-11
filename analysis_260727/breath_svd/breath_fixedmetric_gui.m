function breath_fixedmetric_gui(rootDir)
% breath_fixedmetric_gui  Pick the breathing ROI, then run the whole pipeline.
%
% One call does the lot, the way breath_svd_pc1.m does:
%   1. pops up a preview, you drag the ROI box, ENTER confirms
%   2. extracts the three metrics for every run in the folder
%   3. runs the analysis and opens the QC figures
%
% The preview is OVERLAID WITH BREATHING-BAND MOTION POWER (hot = moving), so
% you can see the motion field while choosing instead of guessing from anatomy.
% That matters here: chest wall motion is not spatially uniform, and for
% breathing that is real physiology -- ribcage and abdominal compartments move
% with different amplitudes and phases.  On the 260723 session one region held
% 0.68x of its amplitude over 75 min while another fell to 0.10x, and their
% motion axes differed by ~70 degrees.  Where you put the box IS the definition
% of "breathing amplitude", so put it on one named compartment.
%
% ONE ROI PER SESSION, not per video -- a per-video box is a per-video spatial
% weighting, i.e. the very non-comparability this pipeline removes.  Runs are
% registered to a common reference, so the single box follows the anatomy as
% the FOV drifts (up to 9 px over a session on this rig).
%
% USAGE
%   breath_fixedmetric_gui
%   breath_fixedmetric_gui('D:\breath_tracking_test')
%
% To reuse an ROI you already picked and skip straight to the compute, set
% PICK_ROI = false below.
%
% Video decoding is FFV1, which MATLAB's VideoReader cannot handle, so the
% heavy lifting runs in Python (dlc310).  This function is the launcher.

%% ------C:\Users\Admin\Desktop\260826_ChAT-soma-g8s_vagotomized----------------------- USER PARAMETERS -----------------------------
ROOT_DIR = 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\phys';
PYEXE    = 'C:\Users\Admin\.conda\envs\dlc310\python.exe';   % NOT cellpose-gpu
PICK_ROI = true;     % false = keep the existing breath_roi.mat and just compute
DO_EXTRACT = true;   % false = skip extraction (reuse breath_fixedmetric.mat)
DO_ANALYZE = true;
DO_BRIDGE  = true;   % write per-run breath_pc1.mat for the peak/trough GUIs
% Which trace the peak/trough GUIs receive.  'fb' by default for TIMING
% PRECISION, not because the others fail to find breaths.  On the 260723
% session all three found the same breaths (119-121 per run, and pv agrees with
% fb at F1 0.946 once the match tolerance is half a breath cycle).  They differ
% in how precisely the peak is placed: median peak-time error vs fb is ~100 ms
% for pv (20% of a cycle, 90th pct 233 ms) versus 2-3x tighter for disp, and fb
% has the steadiest intervals (IBI CV 0.073 vs 0.091 pv, 0.132 disp).
% Per-video PC1 blends chest components with different phases, and the blend
% shifts run to run, which is what smears its peak times.
% If you only need breath COUNT or RATE, any of the three is fine.
%
% This does NOT make fb the amplitude metric.  The GUIs save event INDICES, and
% breath_fixedmetric_analyze applies those same indices to every metric, so
% amplitude is still read off `disp` in physical pixels.  Detect on the trace
% that finds breaths best; measure on the trace whose units mean something.
GUI_METRIC = 'fb';   % 'disp' | 'fb' | 'pv'
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT_DIR = rootDir; end

% This file sits TWO levels below the repo root (analysis_260727\breath_svd),
% so repoRoot needs two fileparts -- one would resolve to analysis_260727 and
% detect_session_fps, Chronux and +helper would not be found. Same correction
% the rest of analysis_260727 already carries.
here     = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(here));
addpath(here); addpath(repoRoot);

assert(isfolder(ROOT_DIR), 'ROOT_DIR not found: %s', ROOT_DIR);
assert(isfile(PYEXE), ['Python not found: %s\n' ...
    'Must be an env with working BLAS -- cellpose-gpu segfaults on matmul.'], PYEXE);

% The two .py workers are siblings of this file, not of the repo root.
roiPy     = fullfile(here, 'breath_fixedmetric_roi.py');
extractPy = fullfile(here, 'breath_fixedmetric_extract.py');
assert(isfile(roiPy),     'missing %s', roiPy);
assert(isfile(extractPy), 'missing %s', extractPy);

%% --- 1. ROI -----------------------------------------------------------------
if PICK_ROI
    fprintf('=== 1/3  ROI selection ===\n');
    fprintf('A preview window will open. Drag a box, then press ENTER.\n');
    fprintf('Hot colour = strong breathing-band motion.\n\n');
    cmd = sprintf('"%s" "%s" "%s"', PYEXE, roiPy, ROOT_DIR);
    status = system(cmd);
    if status ~= 0
        error('ROI selection failed (status %d).', status);
    end
else
    fprintf('=== 1/3  ROI selection skipped (PICK_ROI = false) ===\n');
end

roiFile = fullfile(ROOT_DIR, 'breath_roi.mat');
if isfile(roiFile)
    Rr = load(roiFile);
    fprintf('ROI in use: x=%d y=%d w=%d h=%d', Rr.roi_xywh(1), Rr.roi_xywh(2), ...
        Rr.roi_xywh(3), Rr.roi_xywh(4));
    if isfield(Rr, 'power_frac')
        fprintf('  (%.1f%% of frame breathing-band power)', 100 * Rr.power_frac);
    end
    fprintf('\n');
else
    fprintf('No breath_roi.mat -- the extractor will use the full frame.\n');
end

%% --- 2. extract -------------------------------------------------------------
if DO_EXTRACT
    fprintf('\n=== 2/3  extraction ===\n');
    fprintf('Full decode of every run; several minutes. Cached runs are skipped,\n');
    fprintf('but changing the ROI clears the cache and forces a re-decode.\n\n');
    t0 = tic;
    cmd = sprintf('"%s" -u "%s" "%s"', PYEXE, extractPy, ROOT_DIR);
    status = system(cmd);
    if status ~= 0
        error('extraction failed (status %d).', status);
    end
    fprintf('extraction took %.1f min\n', toc(t0) / 60);
else
    fprintf('\n=== 2/3  extraction skipped (DO_EXTRACT = false) ===\n');
end

%% --- 3. analyze -------------------------------------------------------------
if DO_ANALYZE
    fprintf('\n=== 3/4  analysis ===\n');
    breath_fixedmetric_analyze(ROOT_DIR);
end

%% --- 4. hand the traces to the peak/trough GUIs -----------------------------
% The extractor writes ONE consolidated breath_fixedmetric.mat per session, but
% breathing_peak_gui_pc1.m expects a breath_pc1.mat inside each run folder.  That
% translation runs here rather than being a step to remember.
if DO_BRIDGE
    fprintf('\n=== 4/4  handing traces to the peak/trough GUIs ===\n');
    breath_fixedmetric_to_gui(ROOT_DIR, GUI_METRIC);
end

fprintf('\ndone. Detection workflow:\n');
fprintf('  breathing_peak_gui_pc1      %% Browse Master Folder -> %s\n', ROOT_DIR);
fprintf('  breathing_trough_gui_pc1    %% after saving peaks\n');
fprintf('Then re-run breath_fixedmetric_analyze to use your curated events:\n');
fprintf('  breath_fixedmetric_analyze(''%s'')\n', ROOT_DIR);

end
