function falloff_run_all_260805(rootDir, which)
% FALLOFF_RUN_ALL_260805  Batch driver for the whole optical-penetration project.
%
% THIS is the batch script.  falloff_260804.m and falloff_projections_260805.m
% each do ONE FOV; this table says which FOVs exist and what settings each
% session needs, and runs them all.
%
%   falloff_run_all_260805                 % everything
%   falloff_run_all_260805([], 'falloff')  % falloff figures only
%   falloff_run_all_260805([], 'proj')     % projection figures only
%   falloff_run_all_260805([], 'spfft')    % spatial Fourier figures only
%   falloff_run_all_260805([], 'scatter')  % just the cross-session scatter
%
% Adding a session = adding rows to SESSIONS below.  Nothing else changes.
%
% WHY A TABLE AND NOT AUTO-DETECTION
%   Four things differ per session and NONE of them can be safely inferred from
%   the files:
%     group      wildtype / shiverer -- also sets the fit range
%     lut        which %-to-mW calibration was in force (boundary 2026-07-23)
%     zSurface   motor z of the pia; found from where (top-bot) PEAKS, and it is
%                NOT the top of the stack (260114: stack starts 50 um above it;
%                251104: pia sits 20 um BELOW motor 0)
%     depthMode  'filename' when the filename z IS depth, 'motor' when it is
%                motor z
%     ch         SI channel of the VESSEL dye when the session saved more than
%                one ([] for single-channel sessions).  260806 saved
%                channelSave = [1;3] -- 1 is the SST cell label, 3 the vessels --
%                and since channels interleave page by page, reading the wrong
%                one measures the other label instead of failing.
%   The mount (which side view is coronal) IS derived, from the session date,
%   because that rule is temporal and documented.
%
% Cached per FOV in the data folder, so re-runs are seconds:
%   falloff_<tag>_metrics.mat    per-slice metrics + fixed-bin histograms
%   falloff_proj_<tag>_vol.mat   raw isotropic volume + per-plane setpoint
% Delete those to force a re-read of the TIFFs.

%% --------------------------- USER PARAMETERS -------------------------------
ROOT = 'C:\fall-off';

FIT_WT   = [10 150];        % fit depth range, wildtype (was 10-200 until 260806)
FIT_SHIV = [10 500];        % fit depth range, shiverer

% Depths at which single-plane curves are drawn, [first, step].  Wildtype stacks
% stop near 200 um so they get a 50 um step; shiverer runs to 600 and gets 100.
DIST_WT   = [25  50];
DIST_SHIV = [50 100];

%  session folder            group       lut             zSurf  depthMode   FOV tags                              ch
SESSIONS = { ...
 '251104_wt_fitc',          'wildtype', 'pre_260723',    -20, 'motor',    {'*col01_row02*','*col04_row02*'},      []
 '260728_vglut2_vessel',    'wildtype', 'pre_260723',      0, 'filename', {'roi1'},                               []
 '260114_shiver_vessel',    'shiverer', 'pre_260723',      0, 'motor',    {'tile20','tile23','tile27'},           []
 '260804_shiver_dbh_vessel','shiverer', 'post_260723',     0, 'filename', {'roi1','roi2'},                        []
 '260806_sst_vessel',       'wildtype', 'post_260723',     0, 'filename', {'roi1'},                                3 };
%% ---------------------------------------------------------------------------

if nargin < 1 || isempty(rootDir), rootDir = ROOT; end
if nargin < 2 || isempty(which),   which   = 'all'; end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'spatial-fourier'));
addpath(fileparts(thisDir));

doFall = any(strcmpi(which, {'all','falloff'}));
doProj = any(strcmpi(which, {'all','proj'}));
doFFT  = any(strcmpi(which, {'all','spfft'}));
doScat = any(strcmpi(which, {'all','scatter'}));

nFov = sum(cellfun(@numel, SESSIONS(:,6)));

fprintf('\n=== falloff_run_all: %d sessions, %d FOVs, mode ''%s'' ===\n', ...
    size(SESSIONS,1), nFov, which);

t0 = tic;
for s = 1:size(SESSIONS,1)
    [sess, grp, lut, zSurf, dMode, tags, ch] = SESSIONS{s,:};
    d = fullfile(rootDir, sess);
    if ~isfolder(d), warning('missing session folder: %s -- skipped', d); continue; end
    if strcmpi(grp, 'shiverer')
        fitD = FIT_SHIV;  distS = DIST_SHIV;
    else
        fitD = FIT_WT;    distS = DIST_WT;
    end

    for k = 1:numel(tags)
        fprintf('\n---- %s / %s  (%s, fit %g-%g um, planes %g:%g:) ----\n', ...
            sess, tags{k}, grp, fitD, distS);
        if doFall
            falloff_260804(d, tags{k}, zSurf, lut, fitD, dMode, distS, ch);
        end
        if doProj
            % mount left on 'auto': derived from the session date
            falloff_projections_260805(d, tags{k}, lut, [], zSurf, dMode, [], ch);
        end
        if doFFT
            spatial_fourier_260806(d, tags{k}, zSurf, lut, dMode, distS, ch);
        end
    end
end

if doScat
    fprintf('\n---- cross-session scatter ----\n');
    falloff_summary_scatter_260805(rootDir);
end

fprintf('\n=== done in %.1f min ===\n', toc(t0)/60);
end
