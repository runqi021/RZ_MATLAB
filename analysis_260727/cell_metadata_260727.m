function cell_metadata_260727()
%% cell_metadata_260727  The one metadata table everything else should be read from.
% -----------------------------------------------------------------------
% Until now the answers were scattered: identity in cell_link.csv, geometry in
% roi_match_curated.csv, coherence in sig_rois.csv, modulation in
% breath_time_peth_cells.csv. Nothing tied an ROI to its cell, its place, and its
% physiology in one row, so every figure had to re-join them by hand and could
% quietly join them wrongly.
%
% This builds that join once, from whatever has been run, and writes it at TWO
% levels because both are genuinely needed:
%
%   roi_metadata.csv   ONE ROW PER ROI OBSERVATION (recording x roi index)
%       what it is      : rec_name, roi_index, maskL_label
%       where it is     : x/y/z_um (stage, comparable across recordings)
%                         cx/cy_px (within its own FOV, for drawing on that image)
%       which cell      : cell_id, cell_size, status (grouped/ungrouped/tossed)
%       its own activity: n_events, event rate, recording duration
%       its CELL's score: modulation, latency, p, q -- repeated on every member row
%                         so an ROI-level spatial map needs no join at all
%
%   cell_metadata.csv  ONE ROW PER CELL
%       identity, mean coordinates, how many recordings and events it pools,
%       and every score, once.
%
% The ROI-level table is the one to use for anatomy (each ROI has a real position
% in its own recording); the cell-level table is the one to use for statistics
% (each cell is one independent unit, an ROI is not).
%
% SOURCES, all optional -- whatever is missing is simply left as NaN and reported:
%   cell_pool.mat              identity, geometry, per-observation event counts
%   breath_time_peth_data.mat  PRIMARY modulation score + latency
%   coherence_polar_data.mat   secondary coherence measure, matched per observation
%   cell_coherence_pooled.mat  pooled per-cell coherence, if it was run
%
% Output: <dataset>\analysis_260727\roi_metadata.csv, cell_metadata.csv,
%         metadata.mat
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'coh_ca_breath')); addpath(fileparts(here));
cfg = coh_cfg_260727();

fprintf('\n============ cell_metadata_260727 ============\n');
fprintf('dataset: %s\n', cfg.rootPath);
pool = ensure_pool_260727();   % builds it if missing
assert(isfield(pool,'obsT'), ['This cell_pool.mat predates the obsT field. ' ...
       'Re-run cell_pool_260727.m so geometry is carried through.']);
obsT = pool.obsT;
nObs = height(obsT);
fprintf('observations: %d | cells: %d\n', nObs, numel(pool.cells));

%% ---- per-observation activity from the pool ----
obs = pool.obs;
o_nEv = nan(nObs,1); o_dur = nan(nObs,1); o_usable = false(nObs,1); o_fps = nan(nObs,1);
for i = 1:numel(obs)
    j = obs(i).obs;
    if isnan(j) || j < 1 || j > nObs, continue; end
    o_usable(j) = obs(i).usable;
    if ~isnan(obs(i).n_spikes_used), o_nEv(j) = obs(i).n_spikes_used; end
    k = obs(i).rec;
    if ~isnan(k) && k >= 1 && k <= numel(pool.rec) && pool.rec(k).usable
        o_dur(j) = pool.rec(k).T / pool.rec(k).fps;
        o_fps(j) = pool.rec(k).fps;
    end
end
o_rate = o_nEv ./ max(o_dur,eps) * 60;

%% ---- PRIMARY: time-domain modulation, per cell ----
nCell = numel(pool.cells);
[c_mod, c_lat, c_p, c_q, c_peak, c_base, c_nEvC, c_durC, c_nBr, c_tested, ...
 c_recr, c_prec, c_precz, c_precq, c_mad, c_supz, c_supq, c_peakOverMean] = ...
    deal(nan(nCell,1));
[c_sig, c_sigPre, c_sigSup] = deal(false(nCell,1));
TRIGGER  = 'onset';   % 'onset' | 'peak' -- which alignment's scores go in the metadata
pethFile = fullfile(cfg.outRoot,'breath_time',lower(TRIGGER),'breath_time_peth_data.mat');
assert(isfile(pethFile), ['No PETH results for TRIGGER=''%s'' at\n  %s\n' ...
    'Run breath_time_peth_260727.m with that TRIGGER first.'], TRIGGER, pethFile);
if isfile(pethFile)
    B = load(pethFile,'R','sig_exc','q_exc','sig_sup','q_sup','sig_pre','q_pre','params');
    for c = 1:min(numel(B.R), nCell)
        R = B.R(c);
        if isempty(R.cell_id) || isnan(R.cell_id), continue; end
        c_tested(c) = R.tested;  c_nEvC(c) = R.n_events;  c_durC(c) = R.duration_s;
        c_nBr(c)  = R.n_accepted_breaths;
        c_mod(c)  = R.mod_exc_z;            % PRIMARY score
        c_lat(c)  = R.preferred_latency_s;
        c_peak(c) = R.peak_rate;  c_base(c) = R.mean_peth_rate;
        c_peakOverMean(c) = R.peak_over_mean;
        c_p(c)    = R.p_exc;
        c_recr(c) = R.recruitment;   c_prec(c) = R.precision_fraction;
        c_precz(c)= R.precision_z;   c_mad(c)  = R.latency_mad_s;
        c_supz(c) = R.mod_sup_z;
        if c <= numel(B.q_exc),   c_q(c)      = B.q_exc(c);   end
        if c <= numel(B.sig_exc), c_sig(c)    = B.sig_exc(c); end
        if c <= numel(B.q_pre),   c_precq(c)  = B.q_pre(c);   end
        if c <= numel(B.sig_pre), c_sigPre(c) = B.sig_pre(c); end
        if c <= numel(B.q_sup),   c_supq(c)   = B.q_sup(c);   end
        if c <= numel(B.sig_sup), c_sigSup(c) = B.sig_sup(c); end
    end
    fprintf('joined time-domain modulation for %d cells\n', nnz(~isnan(c_mod)));
else
    fprintf('breath_time_peth_data.mat not found -- modulation columns left NaN\n');
end

%% ---- SECONDARY: per-observation coherence, matched on (rec_name, roi_index) ----
o_cohR = nan(nObs,1); o_cohTh = nan(nObs,1); o_cohSig = false(nObs,1);
if isfile(cfg.cohData)
    C = load(cfg.cohData);
    if isfield(C,'recNames') && isfield(C,'roiIdx')
        key  = string(C.recNames(:)) + "|" + string(C.roiIdx(:));
        okey = string(obsT.rec_name) + "|" + string(obsT.roi_index);
        [tf,loc] = ismember(okey, key);
        o_cohR(tf)   = C.PP.r(loc(tf));
        o_cohTh(tf)  = C.PP.th(loc(tf));
        o_cohSig(tf) = C.PP.r(loc(tf)) >= C.confC;
        fprintf('joined coherence for %d of %d observations\n', nnz(tf), nObs);
    else
        fprintf(['coherence_polar_data.mat has no recNames/roiIdx (it predates the\n' ...
                 '  explicit join key) -- coherence columns left NaN. Re-run the foundation.\n']);
    end
else
    fprintf('coherence_polar_data.mat not found -- coherence columns left NaN\n');
end

%% ---- ROI-LEVEL TABLE: cell scores repeated on every member row ----
cid = obsT.cell_id;
val = ~isnan(cid) & cid >= 1 & cid <= nCell;
g = @(v) local_expand(v, cid, val, nObs);

roiT = table( ...
    obsT.obs, obsT.rec_name, obsT.roi_index, obsT.maskL_label, ...
    obsT.cell_id, obsT.cell_key, obsT.cell_size, obsT.status, ...
    obsT.x_um, obsT.y_um, obsT.z_um, obsT.cx_px, obsT.cy_px, ...
    obsT.cell_x_um, obsT.cell_y_um, obsT.cell_z_um, ...
    o_nEv, o_dur, o_rate, o_fps, o_usable, ...
    g(c_mod), g(c_lat)*1000, g(c_p), g(c_q), g(c_sig), g(c_tested), ...
    g(c_peakOverMean), g(c_recr), g(c_prec), g(c_precz), g(c_sigPre), ...
    g(c_mad)*1000, g(c_supz), g(c_sigSup), ...
    o_cohR, o_cohTh, rad2deg(o_cohTh), o_cohSig, obsT.rec_path, ...
    'VariableNames', { ...
    'obs','rec_name','roi_index','maskL_label', ...
    'cell_id','cell_key','cell_size','status', ...
    'x_um','y_um','z_um','cx_px','cy_px', ...
    'cell_x_um','cell_y_um','cell_z_um', ...
    'roi_n_events','roi_duration_s','roi_rate_per_min','fps','roi_usable', ...
    'cell_mod_exc_z','cell_latency_ms','cell_p_exc','cell_q_exc','cell_sig_exc','cell_tested', ...
    'cell_peak_over_mean','cell_recruitment','cell_precision','cell_precision_z','cell_sig_precision', ...
    'cell_latency_mad_ms','cell_mod_sup_z','cell_sig_sup', ...
    'roi_coherence_r','roi_coherence_th_rad','roi_coherence_th_deg','roi_coherence_sig','rec_path'});

%% ---- CELL-LEVEL TABLE ----
cx = nan(nCell,1); cy = nan(nCell,1); cz = nan(nCell,1);
csz = zeros(nCell,1); crecs = strings(nCell,1); crois = strings(nCell,1); ckey = strings(nCell,1);
cCohMax = nan(nCell,1);
for c = 1:nCell
    m = find(cid == c);
    if isempty(m), continue; end
    csz(c) = numel(m);
    cx(c) = mean(obsT.x_um(m),'omitnan');
    cy(c) = mean(obsT.y_um(m),'omitnan');
    cz(c) = mean(obsT.z_um(m),'omitnan');
    crecs(c) = strjoin(cellstr(string(obsT.rec_name(m))), '|');
    crois(c) = strjoin(string(obsT.roi_index(m)), '|');
    ckey(c)  = obsT.cell_key(m(1));
    if any(~isnan(o_cohR(m))), cCohMax(c) = max(o_cohR(m)); end
end
keep = csz > 0;
cellT = table((1:nCell)', ckey, csz, cx, cy, cz, c_nEvC, c_durC, c_nBr, ...
              c_nEvC./max(c_durC,eps)*60, c_tested, ...
              c_mod, c_lat*1000, c_peak, c_base, c_peakOverMean, c_p, c_q, c_sig, ...
              c_recr, c_prec, c_precz, c_precq, c_sigPre, c_mad*1000, ...
              c_supz, c_supq, c_sigSup, ...
              cCohMax, crecs, crois, ...
    'VariableNames', {'cell_id','cell_key','n_obs','x_um','y_um','z_um', ...
                      'n_events','duration_s','n_accepted_breaths','rate_per_min','tested', ...
                      'mod_exc_z','latency_ms','peak_ev_s','mean_ev_s','peak_over_mean', ...
                      'p_exc','q_exc','sig_exc', ...
                      'recruitment','precision','precision_z','q_precision','sig_precision', ...
                      'latency_mad_ms','mod_sup_z','q_sup','sig_sup', ...
                      'best_roi_coherence_r','recordings','roi_indices'});
cellT = cellT(keep,:);

%% ---- report ----
fprintf('\n---- metadata ----\n');
fprintf('  ROI rows  : %d  (%d grouped into multi-obs cells)\n', height(roiT), nnz(obsT.cell_size>1));
fprintf('  cell rows : %d  (%d multi-observation)\n', height(cellT), nnz(cellT.n_obs>1));
fprintf('  with a modulation score : %d cells\n', nnz(~isnan(cellT.mod_exc_z)));
fprintf('  significant (excitation): %d cells | precision %d | suppression %d\n', nnz(cellT.sig_exc), nnz(cellT.sig_precision), nnz(cellT.sig_sup));
fprintf('  with coherence          : %d ROIs\n', nnz(~isnan(roiT.roi_coherence_r)));
fprintf('  stage coords present    : %d of %d ROIs\n', nnz(~isnan(roiT.x_um)), height(roiT));
if nnz(~isnan(roiT.x_um)) < height(roiT)
    fprintf('    (NaN coords mean identity grouping -- no matcher run for this dataset)\n');
end

%% ---- save ----
if ~isfolder(cfg.outRoot), mkdir(cfg.outRoot); end
writetable(roiT,  fullfile(cfg.outRoot,'roi_metadata.csv'));
writetable(cellT, fullfile(cfg.outRoot,'cell_metadata.csv'));
save(fullfile(cfg.outRoot,'metadata.mat'),'roiT','cellT','cfg','-v7.3');
fprintf('\nSaved roi_metadata.csv (%d rows) + cell_metadata.csv (%d rows) + metadata.mat to\n  %s\n', ...
        height(roiT), height(cellT), cfg.outRoot);
end

%% ========================= helpers =========================
function out = local_expand(v, cid, val, nObs)
% Repeat a per-cell value onto every ROI row belonging to that cell.
if islogical(v), out = false(nObs,1); else, out = nan(nObs,1); end
out(val) = v(cid(val));
end
