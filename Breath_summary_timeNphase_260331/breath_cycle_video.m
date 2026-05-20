clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% breath_cycle_video.m
%  Batch: collects ALL spiking ROIs across ALL FOVs under masterFolder.
%  Reads actual MC TIFF frames, averages across breath onsets (+-window),
%  crops each ROI, stitches into montage, writes video.

%% ========================= USER PARAMETERS =========================
masterFolder   = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx        = [];          % session indices to skip
fps_img        = 30;          % imaging frame rate (Hz) — fallback
nDrop          = 30;          % frames to drop from start
PixelSizeBase  = 1.7778;      % um/px at 1x zoom (ScanImage default)
crop_um        = 25;          % ROI crop box size (um)
gamma_val      = 0.6;         % gamma correction
clip_lo        = 45;           % brightness clip percentile (low)
clip_hi        = 96;        % brightness clip percentile (high)
vidWindow_s    = 1;           % total window around each breath onset (seconds)
vidQuality     = 95;          % MJPEG quality (1-100)
% ====================================================================

%% 1 -- Discover sessions + load coherence significance map
allMat = dir(fullfile(masterFolder, '**', '*cpSAM_output.mat'));
keep = true(numel(allMat),1);
for ii = flip(skipIdx(:)')
    if ii >= 1 && ii <= numel(allMat), keep(ii) = false; end
end
allMat = allMat(keep);
nSess = numel(allMat);
fprintf('Found %d sessions\n', nSess);
assert(nSess > 0, 'No cpSAM_output.mat found under %s', masterFolder);

% Load coherence significance map (from breath_combined_summary)
fov_map_path = fullfile(masterFolder, 'breath_combined_summary', 'fov_map.mat');
assert(isfile(fov_map_path), 'fov_map.mat not found: %s', fov_map_path);
tmp_fm = load(fov_map_path);
fov_map = tmp_fm.fov_map;   % {sessName, [sig_roi_ids]; ...}
fprintf('Loaded fov_map: %d FOVs\n', size(fov_map, 1));

nPhaseBins = [];  % auto: 2 * mean_cycle_frames (set after data collection)

%% 2 -- Collect all ROIs across sessions (sig from coherence)
all_frames  = {};   % {tsz x tsz x nWinFr} per ROI — time domain
all_pframes = {};   % {tsz x tsz x nPhaseBins} per ROI — phase domain
all_perims  = {};
all_labels  = {};
all_tsz     = [];
all_ncyc    = [];
all_cycle_s = [];
all_is_sig  = [];   % 1 = coherence-significant, 0 = not
all_breath_trig = {};

for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    fprintf('\n[%d/%d] %s\n', kk, nSess, sessName);

    try
        % -- Required files --
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
        mc_all = dir(fullfile(folderPath, '*MC*.tif'));
        % Filter out AVG, cp_masks, cpSAM — keep only actual MC video TIFFs
        keep_mc = true(numel(mc_all), 1);
        for im = 1:numel(mc_all)
            nm = mc_all(im).name;
            if contains(nm, {'AVG','cp_masks','cpSAM','mask'}, 'IgnoreCase', true)
                keep_mc(im) = false;
            end
        end
        mc_all = mc_all(keep_mc);
        % Prefer MC_MC over MC
        is_mc2 = contains({mc_all.name}, 'MC_MC');
        mc_hits = [mc_all(is_mc2); mc_all(~is_mc2)];
        if isempty(sam_hits) || isempty(bp_file) || isempty(csv_hits) || isempty(mc_hits)
            fprintf('  Missing files, skipping.\n');
            continue;
        end

        % -- Load SAM for maskL only --
        SAMload = load(fullfile(folderPath, sam_hits(1).name), 'maskL');
        maskL = SAMload.maskL;

        % -- MC TIFF info --
        mcTifPath = fullfile(folderPath, mc_hits(1).name);
        tif_info  = imfinfo(mcTifPath);
        nTifFrames = numel(tif_info);
        T = nTifFrames - nDrop;
        fs = detect_session_fps(folderPath, fps_img);
        fprintf('  MC TIFF: %s (%d frames, fs=%d)\n', mc_hits(end).name, nTifFrames, fs);

        % -- Zoom -> um_per_px (from _meta.mat or any TIFF with metadata) --
        um_per_px = PixelSizeBase;
        meta_hits = dir(fullfile(folderPath, '*_meta.mat'));
        if ~isempty(meta_hits)
            mt = load(fullfile(meta_hits(1).folder, meta_hits(1).name));
            if isfield(mt,'pixelSize_um') && isfinite(mt.pixelSize_um) && mt.pixelSize_um > 0
                um_per_px = mt.pixelSize_um;
            elseif isfield(mt,'zoomFactor') && isfinite(mt.zoomFactor) && mt.zoomFactor > 0
                um_per_px = PixelSizeBase / mt.zoomFactor;
            end
        end
        fprintf('  um_per_px = %.4f\n', um_per_px);

        % -- Spiking ROIs + per-ROI spike times --
        ca_file = fullfile(folderPath, 'ca_spike_data.mat');
        if ~isfile(ca_file), fprintf('  No ca_spike_data.mat, skipping.\n'); continue; end
        in = load(ca_file);
        roiAll_id = find(in.ifSpike == 1);
        nSpk = numel(roiAll_id);
        if nSpk == 0, fprintf('  No spiking ROIs.\n'); continue; end
        % Store spike frame indices per spiking ROI
        roi_spike_idx = cell(nSpk, 1);
        for rr = 1:nSpk
            roi_spike_idx{rr} = in.roi_spikes(roiAll_id(rr)).spike_idx(:);
        end

        % Look up coherence-sig ROIs for this session from fov_map
        sig_roi_ids = [];
        for fm = 1:size(fov_map, 1)
            if strcmp(fov_map{fm,1}, sessName)
                sig_roi_ids = fov_map{fm, 2};
                break;
            end
        end
        roiAll_sig = ismember(roiAll_id, sig_roi_ids);
        fprintf('  %d spiking ROIs (%d coherence-sig, %d non-sig)\n', ...
            nSpk, sum(roiAll_sig), sum(~roiAll_sig));

        % -- Breath onsets --
        bp = load(fullfile(bp_file(1).folder, bp_file(1).name));
        onsets = bp.insp_onset_idx;
        onsets(onsets < nDrop) = [];
        onsets = onsets(:) - nDrop;
        onsets = onsets(onsets > nDrop & onsets <= T - nDrop);
        if numel(onsets) < 2, fprintf('  Too few breath peaks.\n'); continue; end

        b_frames = sort(onsets);
        one_cycle = mean(diff(b_frames) / fs);

        % -- Reconstruct breathing signal --
        dlc_raw  = readmatrix(fullfile(folderPath, csv_hits(1).name), 'NumHeaderLines', 3);
        data_dlc = dlc_raw;
        data_dlc(1:min(nDrop,size(data_dlc,1)), :) = [];

        dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
        dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
        fp = bp.findpeak_params;
        nActiveDots = numel(fp.dot_selection);
        traces_b = NaN(size(data_dlc,1), nActiveDots);
        for d = 1:nActiveDots
            di  = dot_idx.(fp.dot_selection{d});
            xc  = data_dlc(:, dot_cols(di,1));
            yc  = data_dlc(:, dot_cols(di,2));
            pc  = data_dlc(:, dot_cols(di,3));
            bad = pc < fp.likelihood_thr;
            switch fp.coord_types{di}
                case 'x',     sig = xc;
                case 'y',     sig = yc;
                case '-x',    sig = -xc;
                case '-y',    sig = -yc;
                case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
                otherwise,    sig = xc;
            end
            sig(bad) = NaN;
            if sum(~isnan(sig)) >= 2
                sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');
            end
            traces_b(:,d) = sig;
        end
        switch fp.combine_method
            case 'sum', breath = sum(traces_b, 2, 'omitnan');
            otherwise,  breath = mean(traces_b, 2, 'omitnan');
        end
        if fp.inverted, breath = -breath; end
        breath = detrend(breath);

        % -- Breath-triggered average breathing (z-scored, raw snippets) --
        half_win = round(vidWindow_s / 2 * fs);
        nWinFr   = 2 * half_win + 1;
        breath = (breath - mean(breath)) / std(breath);  % z-score once
        breath_snips = nan(numel(b_frames), nWinFr);
        for bb = 1:numel(b_frames)
            i1 = b_frames(bb) - half_win;
            i2 = b_frames(bb) + half_win;
            if i1 >= 1 && i2 <= size(breath,1)
                breath_snips(bb,:) = breath(i1:i2);
            end
        end
        ok_b = ~all(isnan(breath_snips),2);
        if any(ok_b)
            all_breath_trig{end+1} = mean(breath_snips(ok_b,:),1,'omitnan');
        end

        % -- Precompute ROI crop coords --
        [imgH, imgW] = deal(tif_info(1).Height, tif_info(1).Width);
        cpx = round(crop_um / um_per_px / 2);
        tsz = 2*cpx + 1;

        roi_coords = struct('rid',{},'sig',{},'spk_rr',{},'cy',{},'cx',{},'r1',{},'r2',{},'c1',{},'c2',{},'dy',{},'dx',{},'ph',{},'pw',{});
        for rr = 1:nSpk
            rid = roiAll_id(rr);
            [ry, rx] = find(maskL == rid);
            if isempty(ry), continue; end
            cy = round(mean(ry)); cx = round(mean(rx));
            r1 = max(1,cy-cpx); r2 = min(imgH,cy+cpx);
            c1 = max(1,cx-cpx); c2 = min(imgW,cx+cpx);
            dy = cpx-(cy-r1); dx = cpx-(cx-c1);
            roi_coords(end+1) = struct('rid',rid,'sig',roiAll_sig(rr),'spk_rr',rr,'cy',cy,'cx',cx, ...
                'r1',r1,'r2',r2,'c1',c1,'c2',c2, ...
                'dy',dy,'dx',dx,'ph',r2-r1+1,'pw',c2-c1+1);
        end
        nCrop = numel(roi_coords);
        if nCrop == 0
            fprintf('  No valid ROI crops.\n');
            continue;
        end

        % -- Per-ROI: which breath onsets have a spike in that cycle? --
        nOnsets = numel(b_frames);
        % onset_has_spike(bb, rr) = true if ROI rr spiked in cycle bb
        onset_has_spike = false(nOnsets, nCrop);
        for rr = 1:nCrop
            spk = roi_spike_idx{roi_coords(rr).spk_rr};  % spike frames in trimmed space
            for bb = 1:nOnsets
                cyc_start = b_frames(bb);
                cyc_end = cyc_start + round(one_cycle * fs);
                if bb < nOnsets
                    cyc_end = b_frames(bb+1);
                end
                if any(spk >= cyc_start & spk < cyc_end)
                    onset_has_spike(bb, rr) = true;
                end
            end
        end

        % -- Init per-ROI accumulators --
        crop_sum = zeros(tsz, tsz, nWinFr, nCrop);
        crop_cnt = zeros(nWinFr, nCrop);  % per-ROI counts

        % -- Determine needed TIFF frames --
        needed_idx = false(nTifFrames, 1);
        for bb = 1:nOnsets
            if ~any(onset_has_spike(bb,:)), continue; end
            for ff = 1:nWinFr
                tidx = b_frames(bb) - half_win + ff - 1 + nDrop;
                if tidx >= 1 && tidx <= nTifFrames
                    needed_idx(tidx) = true;
                end
            end
        end
        fprintf('  Reading %d / %d MC frames for %d ROIs...\n', ...
            sum(needed_idx), nTifFrames, nCrop);

        % -- Read frames, crop each ROI, accumulate (only spike-active cycles) --
        for tidx = find(needed_idx)'
            fr = double(imread(mcTifPath, tidx, 'Info', tif_info));
            trimmed = tidx - nDrop;
            for bb = 1:nOnsets
                ff = trimmed - (b_frames(bb) - half_win) + 1;
                if ff >= 1 && ff <= nWinFr
                    for rr = 1:nCrop
                        if ~onset_has_spike(bb, rr), continue; end
                        rc = roi_coords(rr);
                        patch = fr(rc.r1:rc.r2, rc.c1:rc.c2);
                        crop_sum(rc.dy+1:rc.dy+rc.ph, rc.dx+1:rc.dx+rc.pw, ff, rr) = ...
                            crop_sum(rc.dy+1:rc.dy+rc.ph, rc.dx+1:rc.dx+rc.pw, ff, rr) + patch;
                        crop_cnt(ff, rr) = crop_cnt(ff, rr) + 1;
                    end
                end
            end
        end

        % -- Average + normalize + store each ROI --
        for rr = 1:nCrop
            rc = roi_coords(rr);
            crop_avg = crop_sum(:,:,:,rr);
            for ff = 1:nWinFr
                if crop_cnt(ff, rr) > 0
                    crop_avg(:,:,ff) = crop_avg(:,:,ff) / crop_cnt(ff, rr);
                end
            end

            % Global normalization across all frames
            all_px = crop_avg(:);
            vlo = prctile(all_px, clip_lo);
            vhi = prctile(all_px, clip_hi);
            if vhi <= vlo, vhi = vlo + eps; end
            crop_avg = max(0, min(1, (crop_avg - vlo)/(vhi - vlo))) .^ gamma_val;

            % Perimeter
            mc_v = maskL(rc.r1:rc.r2, rc.c1:rc.c2);
            pr = bwperim(mc_v == rc.rid);
            pt = false(tsz, tsz);
            pt(rc.dy+1:rc.dy+rc.ph, rc.dx+1:rc.dx+rc.pw) = pr;

            all_frames{end+1}  = crop_avg;
            all_perims{end+1}  = pt;
            all_labels{end+1}  = sprintf('S%d #%d', kk, rc.rid);
            all_tsz(end+1)     = tsz;
            all_ncyc(end+1)    = nWinFr;
            all_cycle_s(end+1) = one_cycle;
            all_is_sig(end+1)  = rc.sig;
        end

        % -- Phase-domain accumulation (0 to 4pi = 2 cycles per onset) --
        sess_nPhaseBins = 2 * round(one_cycle * fs);  % match data resolution
        sess_nPhaseBins = max(sess_nPhaseBins, 4);
        phase_sum = zeros(tsz, tsz, sess_nPhaseBins, nCrop);
        phase_cnt = zeros(sess_nPhaseBins, nCrop);
        phase_edges = linspace(0, 4*pi, sess_nPhaseBins+1);

        nOnsets_p = numel(b_frames);
        for bb = 1:nOnsets_p-2   % need bb+2 for 0-4pi
            if ~any(onset_has_spike(bb, :)), continue; end
            % Get cycle lengths
            len1 = b_frames(bb+1) - b_frames(bb);
            len2 = b_frames(bb+2) - b_frames(bb+1);
            if len1 < 2 || len2 < 2, continue; end

            % Frames spanning 2 cycles: b_frames(bb) to b_frames(bb+2)-1
            for f = b_frames(bb):b_frames(bb+2)-1
                tidx = f + nDrop;
                if tidx < 1 || tidx > nTifFrames, continue; end

                % Compute phase (0 to 4pi)
                if f < b_frames(bb+1)
                    ph = 2*pi * (f - b_frames(bb)) / len1;
                else
                    ph = 2*pi + 2*pi * (f - b_frames(bb+1)) / len2;
                end

                % Find bin
                bin = find(ph >= phase_edges(1:end-1) & ph < phase_edges(2:end), 1);
                if isempty(bin), continue; end

                fr = double(imread(mcTifPath, tidx, 'Info', tif_info));
                for rr = 1:nCrop
                    if ~onset_has_spike(bb, rr), continue; end
                    rc = roi_coords(rr);
                    patch = fr(rc.r1:rc.r2, rc.c1:rc.c2);
                    phase_sum(rc.dy+1:rc.dy+rc.ph, rc.dx+1:rc.dx+rc.pw, bin, rr) = ...
                        phase_sum(rc.dy+1:rc.dy+rc.ph, rc.dx+1:rc.dx+rc.pw, bin, rr) + patch;
                    phase_cnt(bin, rr) = phase_cnt(bin, rr) + 1;
                end
            end
        end

        % Average + normalize phase crops per ROI
        for rr = 1:nCrop
            pavg = phase_sum(:,:,:,rr);
            for bb2 = 1:sess_nPhaseBins
                if phase_cnt(bb2, rr) > 0
                    pavg(:,:,bb2) = pavg(:,:,bb2) / phase_cnt(bb2, rr);
                end
            end
            all_px = pavg(:);
            vlo = prctile(all_px, clip_lo);
            vhi = prctile(all_px, clip_hi);
            if vhi <= vlo, vhi = vlo + eps; end
            pavg = max(0, min(1, (pavg - vlo)/(vhi - vlo))) .^ gamma_val;
            all_pframes{end+1} = pavg;  % [tsz x tsz x sess_nPhaseBins]
        end

        nSig = sum([roi_coords.sig]);
        fprintf('  %d ROIs collected (%d sig, %d non-sig, cycle=%.3fs)\n', ...
            nCrop, nSig, nCrop-nSig, one_cycle);

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
    end
end

nTotal = numel(all_frames);
fprintf('\n=== Total ROIs collected: %d ===\n', nTotal);
assert(nTotal > 0, 'No ROIs found across any session.');
fprintf('  %d sig (coherence), %d non-sig\n', sum(all_is_sig), sum(~all_is_sig));

% Video frame count: fixed window at 30 Hz
mean_cycle_s = mean(all_cycle_s);
vidFps       = 30;
nVidFrames   = 2 * round(vidWindow_s / 2 * vidFps) + 1;
fprintf('Mean cycle = %.3f s | Video: %d frames, +/-%.2f s at %d Hz\n', ...
    mean_cycle_s, nVidFrames, vidWindow_s/2, vidFps);

%% 3 -- Resample to common tile size and frame count
tsz_common = max(all_tsz);

frames_u = zeros(tsz_common, tsz_common, nVidFrames, nTotal);
perims_u = false(tsz_common, tsz_common, nTotal);

for ii = 1:nTotal
    fr = all_frames{ii};   % [tsz x tsz x nWinFr_native]
    nOrig = size(fr, 3);

    % Resize spatial if needed
    if all_tsz(ii) ~= tsz_common
        fr2 = zeros(tsz_common, tsz_common, nOrig);
        for ff = 1:nOrig
            fr2(:,:,ff) = imresize(fr(:,:,ff), [tsz_common tsz_common], 'bilinear');
        end
        fr = fr2;
        perims_u(:,:,ii) = imresize(double(all_perims{ii}), ...
            [tsz_common tsz_common], 'nearest') > 0.5;
    else
        perims_u(:,:,ii) = all_perims{ii};
    end

    % Resample temporal if needed
    if nOrig == nVidFrames
        frames_u(:,:,:,ii) = fr;
    else
        t_orig = linspace(0, 1, nOrig);
        t_new  = linspace(0, 1, nVidFrames);
        for r = 1:size(fr,1)
            for c = 1:size(fr,2)
                frames_u(r,c,:,ii) = interp1(t_orig, squeeze(fr(r,c,:)), t_new, 'pchip');
            end
        end
    end
end

%% -- Grand-mean breath trace --
nBsess = numel(all_breath_trig);
breath_all = nan(nBsess, nVidFrames);
for ss = 1:nBsess
    bt = all_breath_trig{ss};
    nOrig = numel(bt);
    if nOrig == nVidFrames
        breath_all(ss,:) = bt;
    else
        breath_all(ss,:) = interp1(linspace(0,1,nOrig), bt, ...
            linspace(0,1,nVidFrames), 'pchip');
    end
end
breath_mean = mean(breath_all, 1, 'omitnan');
bmn = min(breath_mean); bmx = max(breath_mean);
if bmx > bmn
    breath_mean = (breath_mean - bmn) / (bmx - bmn);
else
    breath_mean = zeros(size(breath_mean));
end

%% 4 -- Montage grid layout: sig LEFT | gap | non-sig RIGHT
idx_sig  = find(all_is_sig == 1);
idx_nsig = find(all_is_sig == 0);
nSig  = numel(idx_sig);
nNsig = numel(idx_nsig);

% Each group gets its own column count
nC_sig  = max(1, ceil(sqrt(nSig)));
nR_sig  = ceil(max(nSig,1) / nC_sig);
nC_nsig = max(1, ceil(sqrt(nNsig)));
nR_nsig = ceil(max(nNsig,1) / nC_nsig);

nR = max(nR_sig, nR_nsig);
gp = 2;
gap_col = tsz_common;  % one empty tile-width column between groups

w_sig  = nC_sig  * (tsz_common+gp) + gp;
w_nsig = nC_nsig * (tsz_common+gp) + gp;
mw = w_sig + gap_col + w_nsig;
mh = nR * (tsz_common+gp) + gp;

% Compute tile origins
torg = zeros(nTotal, 2);

% Sig group (left)
for ii = 1:nSig
    [gr, gc] = ind2sub([nR, nC_sig], ii);
    torg(idx_sig(ii),1) = gp + (gr-1)*(tsz_common+gp) + 1;
    torg(idx_sig(ii),2) = gp + (gc-1)*(tsz_common+gp) + 1;
end

% Non-sig group (right, offset by w_sig + gap_col)
x_off = w_sig + gap_col;
for ii = 1:nNsig
    [gr, gc] = ind2sub([nR, nC_nsig], ii);
    torg(idx_nsig(ii),1) = gp + (gr-1)*(tsz_common+gp) + 1;
    torg(idx_nsig(ii),2) = x_off + gp + (gc-1)*(tsz_common+gp) + 1;
end

fprintf('Layout: %d sig (%dx%d) | gap | %d non-sig (%dx%d) | %d rows\n', ...
    nSig, nR_sig, nC_sig, nNsig, nR_nsig, nC_nsig, nR);

% Static perimeter mask
pmask = false(mh, mw);
for ii = 1:nTotal
    ro = torg(ii,1); co = torg(ii,2);
    if ro == 0, continue; end  % skip if no valid crop
    pmask(ro:ro+tsz_common-1, co:co+tsz_common-1) = ...
        pmask(ro:ro+tsz_common-1, co:co+tsz_common-1) | perims_u(:,:,ii);
end

%% 5 -- Scale + pre-render labels + breath trace strip
vsc = max(1, ceil(100 / tsz_common));
trace_h = round(tsz_common * vsc / 2);
pad_h   = 4;
bar_h   = trace_h + pad_h + 42;
frame_h = mh*vsc + bar_h;
frame_w = mw*vsc;

% Pre-render breath trace (white on black)
trace_strip = zeros(trace_h, frame_w, 3, 'uint8');
margin_x = 10;
plot_w = frame_w - 2*margin_x;
for pp = 1:nVidFrames-1
    x1 = margin_x + round((pp-1)/(nVidFrames-1) * plot_w);
    x2 = margin_x + round(pp/(nVidFrames-1) * plot_w);
    y1 = trace_h - 2 - round(breath_mean(pp) * (trace_h - 4));
    y2 = trace_h - 2 - round(breath_mean(pp+1) * (trace_h - 4));
    nPts = max(abs(x2-x1), abs(y2-y1)) + 1;
    xs = round(linspace(x1, x2, nPts));
    ys = round(linspace(y1, y2, nPts));
    xs = max(1, min(frame_w, xs));
    ys = max(1, min(trace_h, ys));
    for kp = 1:numel(xs)
        trace_strip(ys(kp), xs(kp), :) = 255;
    end
end

fprintf('Montage: %d sig + %d nsig, tile=%dpx, scale=%dx -> %dx%d video\n', ...
    nSig, nNsig, tsz_common, vsc, frame_w, frame_h);

% Pre-render ROI labels
lbl_frame = zeros(frame_h, frame_w, 3, 'uint8');
has_lbl = false;
try
    lpos = zeros(nTotal, 2);
    for ii = 1:nTotal
        lpos(ii,:) = [(torg(ii,2)-1)*vsc+2, (torg(ii,1)-1)*vsc+2];
    end
    lbl_frame = insertText(lbl_frame, lpos, all_labels, ...
        'FontSize', 10, 'BoxOpacity', 0, ...
        'TextColor', 'yellow', 'AnchorPoint', 'LeftTop');
    has_lbl = true;
catch
end
lbl_mask = any(lbl_frame > 0, 3);

%% 6 -- Write video
vidPath = fullfile(masterFolder, 'breathCycleF_allFOV.avi');
vw = VideoWriter(vidPath, 'Motion JPEG AVI');
vw.FrameRate = vidFps;
vw.Quality   = vidQuality;
open(vw);

for tt = 1:nVidFrames
    % Build montage from real averaged pixel data
    fg = zeros(mh, mw);
    for ii = 1:nTotal
        if torg(ii,1) == 0, continue; end
        tile = frames_u(:,:,tt,ii);
        tile = max(0, min(1, tile));
        ro = torg(ii,1); co = torg(ii,2);
        fg(ro:ro+tsz_common-1, co:co+tsz_common-1) = tile;
    end

    % Grayscale -> RGB, yellow perimeter
    frgb = repmat(uint8(fg*255), [1 1 3]);
    for cc = 1:2
        p = frgb(:,:,cc); p(pmask) = 255; frgb(:,:,cc) = p;
    end
    p = frgb(:,:,3); p(pmask) = 0; frgb(:,:,3) = p;

    % Scale up
    frgb = imresize(frgb, vsc, 'nearest');

    % Bottom strip: breath trace + moving white cursor
    strip = trace_strip;
    cursor_x = margin_x + round((tt-1)/(nVidFrames-1) * plot_w);
    cursor_x = max(1, min(frame_w, cursor_x));
    strip(:, max(1,cursor_x-1):min(frame_w,cursor_x+1), :) = 255;

    bottom = zeros(bar_h, frame_w, 3, 'uint8');
    bottom(pad_h+1:pad_h+trace_h, :, :) = strip;
    frgb = cat(1, frgb, bottom);

    % Composite ROI labels
    if has_lbl
        for cc = 1:3
            p = frgb(:,:,cc); lp = lbl_frame(:,:,cc);
            p(lbl_mask) = lp(lbl_mask);
            frgb(:,:,cc) = p;
        end
    end

    % Time text
    t_sec = -vidWindow_s/2 + (tt-1)/vidFps;
    try
        frgb = insertText(frgb, [4, frame_h-38], ...
            sprintf('time = %+.3f s', t_sec), ...
            'FontSize', 30, 'BoxOpacity', 0, ...
            'TextColor', 'white', 'AnchorPoint', 'LeftTop');
    catch
    end

    writeVideo(vw, frgb);
end

close(vw);
fprintf('\nTime video saved: %s\n', vidPath);
fprintf('  %d ROIs, %d frames, %d fps playback\n', nTotal, nVidFrames, vidFps);

%% 7 -- Phase video (0 to 4pi)
% Common phase bin count = 2 * mean cycle frames
nPhaseBins = 2 * round(mean_cycle_s * vidFps);
nPhaseBins = max(nPhaseBins, 4);
fprintf('Phase video: %d bins (0-4pi)\n', nPhaseBins);

% Resample phase frames to common tile size + common bin count
pframes_u = zeros(tsz_common, tsz_common, nPhaseBins, nTotal);
for ii = 1:nTotal
    pfr = all_pframes{ii};   % [tsz x tsz x sess_nPhaseBins]
    nOrig = size(pfr, 3);

    % Spatial resize
    if all_tsz(ii) ~= tsz_common
        pfr2 = zeros(tsz_common, tsz_common, nOrig);
        for ff = 1:nOrig
            pfr2(:,:,ff) = imresize(pfr(:,:,ff), [tsz_common tsz_common], 'bilinear');
        end
        pfr = pfr2;
    end

    % Temporal resample to common nPhaseBins
    if nOrig == nPhaseBins
        pframes_u(:,:,:,ii) = pfr;
    else
        t_orig = linspace(0, 1, nOrig);
        t_new  = linspace(0, 1, nPhaseBins);
        for r = 1:tsz_common
            for c = 1:tsz_common
                pframes_u(r,c,:,ii) = interp1(t_orig, squeeze(pfr(r,c,:)), t_new, 'pchip');
            end
        end
    end
end

% Pre-render phase breath trace (idealized: 2 cosine cycles, 0-4pi)
phase_trace = zeros(trace_h, frame_w, 3, 'uint8');
for pp = 1:nPhaseBins-1
    ph1 = 4*pi * (pp-1) / (nPhaseBins-1);
    ph2 = 4*pi * pp / (nPhaseBins-1);
    x1 = margin_x + round((pp-1)/(nPhaseBins-1) * plot_w);
    x2 = margin_x + round(pp/(nPhaseBins-1) * plot_w);
    % Cosine wave normalized to [0,1]: peak at 0, trough at pi
    y1 = trace_h - 2 - round((0.5 + 0.5*cos(ph1)) * (trace_h - 4));
    y2 = trace_h - 2 - round((0.5 + 0.5*cos(ph2)) * (trace_h - 4));
    nPts = max(abs(x2-x1), abs(y2-y1)) + 1;
    xs = round(linspace(x1, x2, nPts));
    ys = round(linspace(y1, y2, nPts));
    xs = max(1, min(frame_w, xs));
    ys = max(1, min(trace_h, ys));
    for kp = 1:numel(xs)
        phase_trace(ys(kp), xs(kp), :) = 255;
    end
end

% Write phase video
vidPath_ph = fullfile(masterFolder, 'breathCycleF_allFOV_phase.avi');
vw2 = VideoWriter(vidPath_ph, 'Motion JPEG AVI');
vw2.FrameRate = vidFps;
vw2.Quality   = vidQuality;
open(vw2);

for tt = 1:nPhaseBins
    fg = zeros(mh, mw);
    for ii = 1:nTotal
        if torg(ii,1) == 0, continue; end
        tile = pframes_u(:,:,tt,ii);
        tile = max(0, min(1, tile));
        ro = torg(ii,1); co = torg(ii,2);
        fg(ro:ro+tsz_common-1, co:co+tsz_common-1) = tile;
    end

    frgb = repmat(uint8(fg*255), [1 1 3]);
    for cc = 1:2
        p = frgb(:,:,cc); p(pmask) = 255; frgb(:,:,cc) = p;
    end
    p = frgb(:,:,3); p(pmask) = 0; frgb(:,:,3) = p;

    frgb = imresize(frgb, vsc, 'nearest');

    % Bottom strip: phase trace + cursor
    strip = phase_trace;
    cursor_x = margin_x + round((tt-1)/(nPhaseBins-1) * plot_w);
    cursor_x = max(1, min(frame_w, cursor_x));
    strip(:, max(1,cursor_x-1):min(frame_w,cursor_x+1), :) = 255;

    bottom = zeros(bar_h, frame_w, 3, 'uint8');
    bottom(pad_h+1:pad_h+trace_h, :, :) = strip;
    frgb = cat(1, frgb, bottom);

    if has_lbl
        for cc = 1:3
            p = frgb(:,:,cc); lp = lbl_frame(:,:,cc);
            p(lbl_mask) = lp(lbl_mask);
            frgb(:,:,cc) = p;
        end
    end

    % Phase text
    ph_val = 4*pi * (tt-1) / (nPhaseBins-1);
    try
        frgb = insertText(frgb, [4, frame_h-38], ...
            sprintf('phase = %.2f rad (%.0f%s)', ph_val, ph_val/pi*180, char(176)), ...
            'FontSize', 30, 'BoxOpacity', 0, ...
            'TextColor', 'white', 'AnchorPoint', 'LeftTop');
    catch
    end

    writeVideo(vw2, frgb);
end

close(vw2);
fprintf('\nPhase video saved: %s\n', vidPath_ph);
fprintf('  %d ROIs, %d bins (0-4pi), %d fps playback\n', nTotal, nPhaseBins, vidFps);
