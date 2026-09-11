% manual_F_override_260817.m
% -----------------------------------------------------------------------
%  Replace ONE ROI's Cellpose fluorescence with a hand-drawn ImageJ measurement,
%  in place, in that recording's *_ch1_dFF.mat.
%
%  WHY THIS IS SAFE TO DO PER-COLUMN. dFF is computed by
%      dFF = (F - movmedian(F, win, 1)) ./ movmedian(F, win, 1)
%  and movmedian along dim 1 is independent per column, so recomputing a single
%  column gives exactly what recomputing the whole matrix would give. There is
%  no cross-ROI term anywhere in the dF/F path -- no neuropil subtraction, no
%  common-mode removal -- so no other ROI's numbers change.
%
%  IT PROVES ITS OWN PARAMETERS BEFORE WRITING. The stored `params` records fps
%  and baselineWinSec but not the full call, so the script first recomputes the
%  EXISTING dFF column from the EXISTING F column and requires it to reproduce
%  what is already in the file to ~1e-12. If it does not, the reconstruction is
%  wrong and the script aborts having changed nothing. Without that check a
%  silently different baseline window would put a trace on the page that no
%  other ROI in the archive is comparable with.
%
%  THE ORIGINAL IS BACKED UP FIRST, to *_ch1_dFF_preManualF_<stamp>.mat, and the
%  script refuses to run if that name already exists. Provenance (source csv,
%  roi, cell id, when, and the fact that this column is not Cellpose's) is
%  written into the mat as `manualF`, so a file that has been overridden always
%  says so from the inside.
%
%  WHAT THIS DOES NOT DO: ca_spike_data.mat is untouched. Events were detected
%  from the OLD trace, so every event-based number for this cell -- PSTH, logZ,
%  phase, the polar dot, its class -- still describes the Cellpose trace. Only
%  the dF/F and raw-F trace panels change. Re-run detection separately if the
%  events are meant to follow the new trace.
%
%  Runqi Zhang / 2026-08-17
% -----------------------------------------------------------------------

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);

%% ===================== USER-EDITABLE =====================
fd     = 'D:\Ventral_surface_summary\Vgat\0730\deep\roi1_2.4x_x1300y900_z250_3000f_30lp_00001';
csvF   = fullfile(fd, 'cell279_manual_F.csv');
roi    = 32;      % ROI index INSIDE this recording's dFF matrix
cellId = 277;     % registry cell id, for the provenance record only
note   = ['hand-drawn ImageJ ROI replacing the cpSAM mask for cell 277. ' ...
          'NOTE the csv is named cell279 but sits in the z250 recording, whose ' ...
          'registered cell is 277 (= ROI 32); cell 279 is z265/ROI29. ' ...
          'Confirmed by RZ 2026-08-17 that 277 is meant.'];
% =========================================================

assert(isfolder(fd),  'folder not found: %s', fd);
assert(isfile(csvF),  'csv not found: %s', csvF);

d = dir(fullfile(fd,'*_ch1_dFF.mat'));
assert(numel(d)==1, 'expected exactly one *_ch1_dFF.mat in %s (found %d)', fd, numel(d));
matF = fullfile(d.folder, d.name);

%% ---- load ----
S = load(matF);
for f = {'dFF','F_roi','F_roi_raw','params','dFFout'}
    assert(isfield(S,f{1}), 'missing variable %s in %s', f{1}, d.name);
end
[T, N] = size(S.dFF);
assert(roi>=1 && roi<=N, 'roi %d out of range 1..%d', roi, N);

M = readmatrix(csvF);
assert(size(M,2)>=2, 'expected an ImageJ "Slice,Mean" csv: %s', csvF);
Fm = M(:,2);
assert(numel(Fm)==T, ['csv has %d samples but the recording has %d frames. ' ...
    'A manual trace must be measured on the same movie the dFF came from.'], numel(Fm), T);
assert(all(isfinite(Fm)), 'csv contains non-finite values');

fps = S.params.fps;
bws = S.params.baselineWinSec;
fprintf('recording : %s\n', d.name);
fprintf('roi %d (cell %d) | T=%d N=%d | fps=%g baselineWinSec=%g tossFrames=%g\n', ...
        roi, cellId, T, N, fps, bws, S.params.tossFrames);

%% ---- prove the reconstruction on the EXISTING column, before writing ----
% dFF_RZ is called by the pipeline as dFF_RZ(F_roi,'FPS',fps,'BaselineWinSec',bws),
% which leaves DropFirstSec at its default of 0. Recompute the old column and
% require it to match what is stored.
chk = helper.dFF_RZ(S.F_roi(:,roi), 'FPS', fps, 'BaselineWinSec', bws);
err = max(abs(chk.dFF - S.dFF(:,roi)));
fprintf('self-check : recomputed old dFF matches stored to %.3g\n', err);
assert(err < 1e-10, ['ABORT -- could not reproduce the stored dFF from the ' ...
    'stored F with fps=%g, BaselineWinSec=%g (max err %.3g). The parameters ' ...
    'in `params` do not describe how this file was made; nothing was changed.'], ...
    fps, bws, err);

%% ---- back up ----
stamp = datestr(now,'yymmdd_HHMMSS'); %#ok<TNOW1,DATST>
bakF  = regexprep(matF, '_dFF\.mat$', sprintf('_dFF_preManualF_%s.mat', stamp));
assert(~isfile(bakF), 'backup already exists: %s', bakF);
copyfile(matF, bakF);
fprintf('backup    : %s\n', bakF);

%% ---- patch the column ----
oldF   = S.F_roi_raw(:,roi);
oldDff = S.dFF(:,roi);
newDff = helper.dFF_RZ(Fm, 'FPS', fps, 'BaselineWinSec', bws);

S.F_roi_raw(:,roi) = Fm;
S.F_roi(:,roi)     = Fm;          % identical here: tossFrames = 0
S.dFF(:,roi)       = newDff.dFF;
if isfield(S.dFFout,'F_dff') && size(S.dFFout.F_dff,1)==T, S.dFFout.F_dff(:,roi) = Fm;            end
if isfield(S.dFFout,'dFF')   && size(S.dFFout.dFF,1)==T,   S.dFFout.dFF(:,roi)   = newDff.dFF;    end

%% ---- provenance, written INTO the file ----
rec = struct('roi',roi, 'cellId',cellId, 'sourceCsv',csvF, ...
             'appliedOn',datestr(now,'yyyy-mm-dd HH:MM:SS'), ... %#ok<TNOW1,DATST>
             'appliedBy',mfilename, 'backup',bakF, 'note',note, ...
             'fps',fps, 'baselineWinSec',bws);
if isfield(S,'manualF'), S.manualF(end+1) = rec; else, S.manualF = rec; end

save(matF, '-struct', 'S', '-v7.3');

%% ---- report ----
fprintf('\nROI %d replaced with the manual trace:\n', roi);
fprintf('  raw F  mean   %8.0f -> %8.0f counts\n', mean(oldF), mean(Fm));
fprintf('  raw F  CV     %8.4f -> %8.4f\n', std(oldF)/mean(oldF), std(Fm)/mean(Fm));
fprintf('  dFF    max    %8.3f -> %8.3f\n', max(oldDff), max(newDff.dFF));
fprintf('  dFF    std    %8.4f -> %8.4f\n', std(oldDff), std(newDff.dFF));
fprintf('  corr(old raw F, manual)  = %.3f\n', corr(oldF, Fm));
fprintf('\nwrote %s\n', matF);
fprintf('ca_spike_data.mat NOT touched -- events still come from the old trace.\n');
