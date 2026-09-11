% make_cell_link_260824.m
% -----------------------------------------------------------------------
%  Build a cell_link.mat for the 260824 vagotomised Vglut2 session from its
%  roi_match curation, so that every tool keyed to cell_link can read this
%  session's cell identity -- above all temporal_phase_cell_gui_260812.m,
%  which looks for exactly this file and this schema.
%
%  WHY THIS EXISTS
%    This session was matched with the newer roi_match pipeline, whose output is
%    roi_match_out_260824\roi_match_curated.{csv,mat}. cell_link_260727.m was
%    never run on it, so there is no cell_link.mat -- and the GUI's two identity
%    sources are the archive registry (which does not contain this session) and
%    cell_link (absent). With neither, the GUI silently falls back to plain ROI
%    indices, and ROI 17 is NOT cell 17.
%
%    The cell NUMBERS here are the curation's own: cell 17 in this file is the
%    cell 17 of gui_renders_260824\..._cell017_pooled*.png, because both read
%    curatedCellMembers / roi_match_curated.csv. Nothing is renumbered.
%
%  THIS IS A RECONSTRUCTION, NOT A cell_link_260727 RUN. The schema is copied
%  field-for-field from D:\Ventral_surface_summary\Vglut2\0810\cell_pooled\
%  cell_link.mat so that consumers cannot tell the difference, but the grouping
%  comes from the curation file rather than from the matcher being re-run.
%  link.grouping_source records that.
%
%  TWO DELIBERATE DIFFERENCES FROM THE ARCHIVED cell_link FILES
%
%   1. rec_path points into the ARCHIVE (D:\Ventral_surface_summary\Vglut2\0824\
%      <site>\<recording>), not at the acquisition drive. The archived 0810/0730
%      links still point at C:\260810_...\phys, which breaks as soon as that
%      session folder is cleaned off the acquisition drive. This session's data
%      now lives in the archive, so the link points there and stays valid.
%      Consequence: open the ARCHIVE copy of a recording in the GUI, not the C:
%      copy, or the GUI will report "cell N is NOT in the folder you typed".
%
%   2. fov_name holds the FULL recording name. cell_link_260727.m ran fileparts
%      on a dotted folder name and truncated it ("roi1_3x_12" for
%      "roi1_3x_12.5lp_..."), which is why the GUI's header tells you never to
%      read link.fov_name. Nothing here depends on the bug being reproduced.
%
%  TOSSED ROIs ARE INCLUDED, with cell_id = NaN. That is not cosmetic: every
%  scanner (Ventral_surface_polar_coh_vs_rayleigh_260808, polar_coh_rayleigh_260824,
%  ventral_map_*) reads
%        if isnan(Tk.cell_id(i)), tossed_set(kk) = true;
%  to EXCLUDE a rejected ROI. Dropping those rows instead would make each one
%  fall through to the "no registry entry -> one cell per ROI" path and put the
%  curation's rejects back into the analysis.
%
%  OUTPUT  D:\Ventral_surface_summary\Vglut2\0824\cell_pooled\cell_link.mat
%
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; clc;

%% ===================== USER-EDITABLE =====================
archiveSession = 'D:\Ventral_surface_summary\Vglut2\0824';
curatedMat     = fullfile(archiveSession, 'roi_match_out_260824', 'roi_match_curated.mat');
outFile        = fullfile(archiveSession, 'cell_pooled', 'cell_link.mat');
genotype       = 'Vglut2';
dateStr        = '0824';
overwrite      = false;      % refuse to clobber an existing cell_link by default
% =========================================================

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);

assert(isfile(curatedMat), 'curation not found: %s', curatedMat);
if isfile(outFile) && ~overwrite
    error('make_cell_link_260824:exists', ...
          '%s already exists. Set overwrite = true to replace it.', outFile);
end

M = load(curatedMat, 'curated');
C = M.curated;
R = C.roi;
nObs  = numel(R.roi);
fovs  = string(C.fov_name(:));
nCell = numel(C.curatedCellMembers);

fprintf('\n=========== make_cell_link_260824 ===========\n');
fprintf('curation : %s\n', curatedMat);
fprintf('          %d ROI observations, %d FOVs, %d curated cells\n', nObs, numel(fovs), nCell);

%% ---- locate each recording inside the archive -------------------------
% archive_session_260811 clusters recordings into <site> folders by
% motorPosition, so the site a recording landed in is discovered here rather
% than recomputed -- recomputing risks disagreeing with what is on disk.
recPathOf = strings(numel(fovs),1);
siteOf    = strings(numel(fovs),1);
sites = dir(fullfile(archiveSession, 'cell*'));
sites = sites([sites.isdir]);
for i = 1:numel(fovs)
    for s = 1:numel(sites)
        cand = fullfile(archiveSession, sites(s).name, char(fovs(i)));
        if isfolder(cand)
            recPathOf(i) = string(cand);
            siteOf(i)    = string(sites(s).name);
            break;
        end
    end
end
missing = find(recPathOf == "");
if ~isempty(missing)
    fprintf(2, 'NOT FOUND in the archive:\n');
    fprintf(2, '   %s\n', fovs(missing));
    error('make_cell_link_260824:missingRec', ...
          '%d recording(s) are in the curation but not under %s -- archive them first.', ...
          numel(missing), archiveSession);
end
fprintf('archive  : all %d recordings located under %d site folder(s)\n', ...
        numel(fovs), numel(unique(siteOf)));

%% ---- per-recording facts ---------------------------------------------
nRoiDff  = nan(numel(fovs),1);   nRoiSpk = nan(numel(fovs),1);
nFrames  = nan(numel(fovs),1);   fpsRec  = nan(numel(fovs),1);
hasBrth  = false(numel(fovs),1); spkCount = cell(numel(fovs),1);
for i = 1:numel(fovs)
    fp = char(recPathOf(i));
    dd = dir(fullfile(fp, '*_ch1_dFF.mat'));
    if ~isempty(dd)
        D = load(fullfile(dd(1).folder, dd(1).name), 'dFF');
        nRoiDff(i) = size(D.dFF,2);  nFrames(i) = size(D.dFF,1);
    end
    fpsRec(i) = detect_session_fps(fp, NaN);
    hasBrth(i) = isfile(fullfile(fp,'breath_peak_pc1.mat')) && ...
                 isfile(fullfile(fp,'breath_insp_start_pc1.mat'));
    csf = fullfile(fp,'ca_spike_data.mat');
    if isfile(csf)
        CA = load(csf, 'roi_spikes');
        nRoiSpk(i)  = numel(CA.roi_spikes);
        spkCount{i} = arrayfun(@(r) double(sum(r.spike_train(:)>0)), CA.roi_spikes(:));
    end
end
nNoSpk = nnz(isnan(nRoiSpk));
if nNoSpk > 0
    fprintf(2, 'note: %d recording(s) have no ca_spike_data.mat -- has_spikes = false there\n', nNoSpk);
end

%% ---- cell id per observation -----------------------------------------
% Observations in a curated cell take that cell's number. Everything else --
% tossed (grpOf < 0) and anything the curation left unassigned -- gets NaN,
% which is the archive's "exclude me" marker. That matches the CSV, which also
% lists only the observations belonging to a cell, so this file and every figure
% already rendered from the CSV describe exactly the same set of cells.
cellIdOf = nan(nObs,1);
for c = 1:nCell
    mem = C.curatedCellMembers{c}(:)';
    cellIdOf(mem) = c;
end
isTossed = false(nObs,1);
if isfield(C,'grpOf'), isTossed = C.grpOf(:) < 0; end
nAssigned  = nnz(~isnan(cellIdOf));
nTossedRow = nnz(isTossed & isnan(cellIdOf));
nOther     = nnz(isnan(cellIdOf) & ~isTossed);
fprintf('identity : %d obs in cells, %d tossed, %d unassigned (all NaN -> excluded)\n', ...
        nAssigned, nTossedRow, nOther);

cellSizeOf = nan(nObs,1);
for c = 1:nCell
    mem = C.curatedCellMembers{c}(:)';
    cellSizeOf(mem) = numel(mem);
end

% cell centroid = mean over the cell's observations, as in cell_link_260727
cellX = nan(nObs,1); cellY = nan(nObs,1); cellZ = nan(nObs,1);
for c = 1:nCell
    mem = C.curatedCellMembers{c}(:)';
    cellX(mem) = mean(R.x(mem));
    cellY(mem) = mean(R.y(mem));
    cellZ(mem) = mean(R.z(mem));
end

%% ---- build obsT -------------------------------------------------------
obs        = (1:nObs)';
recNameOf  = fovs(R.fov(:));
recPathCol = recPathOf(R.fov(:));
cellKeyCol = strings(nObs,1);
statusCol  = strings(nObs,1);
nSpkCol    = nan(nObs,1);
hasSpkCol  = false(nObs,1);
for i = 1:nObs
    cellKeyCol(i) = sprintf('%s#%d', recNameOf(i), R.roi(i));
    if isnan(cellIdOf(i))
        if isTossed(i), statusCol(i) = "tossed"; else, statusCol(i) = "unassigned"; end
    elseif cellSizeOf(i) > 1
        statusCol(i) = "grouped";
    else
        statusCol(i) = "ungrouped";
    end
    sc = spkCount{R.fov(i)};
    if ~isempty(sc) && R.roi(i) >= 1 && R.roi(i) <= numel(sc)
        nSpkCol(i)   = sc(R.roi(i));
        hasSpkCol(i) = true;
    end
end

obsT = table(obs, cellIdOf, cellKeyCol, cellSizeOf, statusCol, recNameOf, ...
             double(R.roi(:)), uint16(R.lab(:)), R.x(:), R.y(:), R.z(:), ...
             R.cx_px(:), R.cy_px(:), cellX, cellY, cellZ, ...
             hasSpkCol, nSpkCol, recPathCol, ...
    'VariableNames', {'obs','cell_id','cell_key','cell_size','status', ...
                      'rec_name','roi_index','maskL_label','x_um','y_um','z_um', ...
                      'cx_px','cy_px','cell_x_um','cell_y_um','cell_z_um', ...
                      'has_spikes','n_spikes','rec_path'});

%% ---- build recT -------------------------------------------------------
nRoiMatcher = accumarray(R.fov(:), 1, [numel(fovs) 1]);
maxRoiIdx   = accumarray(R.fov(:), double(R.roi(:)), [numel(fovs) 1], @max, 0);
labContig   = false(numel(fovs),1);
for i = 1:numel(fovs)
    rr = sort(double(R.roi(R.fov(:) == i)));
    labContig(i) = ~isempty(rr) && isequal(rr(:)', 1:numel(rr));
end
recT = table(fovs, recPathOf, true(numel(fovs),1), nRoiMatcher, nRoiSpk, nRoiDff, ...
             maxRoiIdx, labContig, nFrames, fpsRec, hasBrth, ...
    'VariableNames', {'rec_name','rec_path','in_matcher','n_roi_matcher', ...
                      'n_roi_spikes','n_roi_dff','max_roi_index', ...
                      'labels_contiguous','n_frames_ca','fps','has_breath'});

%% ---- assemble link ----------------------------------------------------
link = struct();
link.obsT     = obsT;
link.recT     = recT;
link.cells    = C.curatedCellMembers(:);
link.cellKey  = strings(nCell,1);
link.cellSize = zeros(nCell,1);
for c = 1:nCell
    mem = C.curatedCellMembers{c}(:)';
    link.cellSize(c) = numel(mem);
    link.cellKey(c)  = cellKeyCol(mem(1));
end
link.cellOf = cellIdOf;
if isfield(C,'grpOf'), link.grpOf = C.grpOf(:); else, link.grpOf = nan(nObs,1); end
link.grouping_source = 'roi_match_curated_260824 (reconstructed by make_cell_link_260824)';
link.roi        = R;
link.fov_name   = fovs;              % FULL names, not fileparts-truncated
link.fov_folder = recPathOf;

link.audit = struct( ...
    'grouping_source',   link.grouping_source, ...
    'n_obs_total',       nObs, ...
    'n_obs_in_cells',    nAssigned, ...
    'n_obs_tossed',      nTossedRow, ...
    'n_obs_unassigned',  nOther, ...
    'n_cells',           nCell, ...
    'max_cell_id',       nCell, ...
    'n_cells_multi',     nnz(link.cellSize > 1), ...
    'n_cells_single',    nnz(link.cellSize == 1), ...
    'n_rec_scanned',     numel(fovs), ...
    'rec_missing_spikes',fovs(isnan(nRoiSpk)), ...
    'rec_missing_breath',fovs(~hasBrth), ...
    'n_obs_with_spikes', nnz(hasSpkCol));

link.cfg = struct( ...
    'rootPath',    archiveSession, ...
    'genotype',    genotype, ...
    'dateStr',     dateStr, ...
    'curatedFile', curatedMat, ...
    'linkFile',    outFile, ...
    'builtBy',     'make_cell_link_260824.m', ...
    'builtOn',     '2026-08-25', ...
    'note',        ['rec_path points into the ARCHIVE, not the acquisition drive; ' ...
                    'fov_name is the full recording name, not fileparts-truncated']);

%% ---- consistency checks BEFORE writing --------------------------------
% A cell_link that addresses the wrong column of dFF would silently mis-assign
% every trace downstream, so the index range is checked against the real data
% rather than trusted.
bad = {};
for i = 1:numel(fovs)
    rr = double(R.roi(R.fov(:) == i));
    if ~isnan(nRoiDff(i)) && (any(rr < 1) || any(rr > nRoiDff(i)))
        bad{end+1} = sprintf('%s: roi index outside dFF (has %d ROIs, saw %s)', ...
            fovs(i), nRoiDff(i), mat2str(unique(rr(rr<1 | rr>nRoiDff(i)))')); %#ok<SAGROW>
    end
    if ~isnan(nRoiSpk(i)) && ~isnan(nRoiDff(i)) && nRoiSpk(i) ~= nRoiDff(i)
        bad{end+1} = sprintf('%s: ca_spike_data has %d ROIs but dFF has %d', ...
            fovs(i), nRoiSpk(i), nRoiDff(i)); %#ok<SAGROW>
    end
end
if ~isempty(bad)
    fprintf(2,'CONSISTENCY CHECK FAILED:\n'); fprintf(2,'   %s\n', bad{:});
    error('make_cell_link_260824:consistency','refusing to write cell_link');
end
fprintf('checks   : OK -- every ROI index addresses a real dFF column\n');

%% ---- write ------------------------------------------------------------
outDir = fileparts(outFile);
if ~isfolder(outDir), mkdir(outDir); end
save(outFile, 'link');
fprintf('\nwrote %s\n', outFile);
fprintf('  obsT %d rows (%d with a cell id), recT %d rows, %d cells\n', ...
        height(obsT), nAssigned, height(recT), nCell);
fprintf('  cells spanning >1 recording: %d\n', nnz(link.cellSize > 1));
fprintf('Done.\n');
