%% calcium_breath_combined.m
% Per-ROI combined figure: Col 1 = breath-triggered time domain,
%                          Col 2 = multitaper spectrum / coherence / polar.
%
% Requires (all in folderPath):
%   *_dFF.mat, ca_spike_data.mat, breath_peak_data.mat, *snapshot_200.csv
%
% Layout per ROI (3 × 2, all square):
%   (1,1) Breath-trig dF/F        (1,2) PSD  (ROI + breathing)
%   (2,1) Ca²⁺ spike raster       (2,2) Coherence vs frequency
%   (3,1) dt_last vs dt_next       (3,2) Polar  (phase / coherence)
close all; clc; clear all;
%% ============================== USER ==============================
folderPath      = "C:\Users\zhang\Desktop\roi5_1400-1230-0_x4.4_15lp_6000f";
FPS             = 30;       % fallback imaging fps
TW              = 3;        % multitaper time-bandwidth  (K = 2*TW-1 tapers)
fmin            = 0.01;     % Hz  (PSD lower bound)
fmax            = 10;       % Hz  (PSD upper bound; [] → fs/2)
f_breath_search = [0.2, 0.5]; % Hz  (search band for breathing peak + coherence)
saveFigs        = true;

%% ============================== LOAD ==============================

% --- *_dFF.mat ---
dffEntry = dir(fullfile(folderPath, '*_dFF.mat'));
assert(numel(dffEntry) == 1, "Expected 1 *_dFF.mat, found %d.", numel(dffEntry));
dffMatPath = fullfile(dffEntry.folder, dffEntry.name);
[~, stem]  = fileparts(dffMatPath);
D   = load(dffMatPath);
dFF = double(D.dFF);   % [T × N]
[T, N] = size(dFF);
if isfield(D,'params') && isfield(D.params,'FPS')
    fs = double(D.params.FPS);
else
    fs = FPS;
    warning("params.FPS not found — using fallback FPS=%g", FPS);
end
if isempty(fmax), fmax = fs / 2; end

% --- ca_spike_data.mat ---
spikeFile = fullfile(folderPath, 'ca_spike_data.mat');
assert(isfile(spikeFile), "ca_spike_data.mat not found — run calcium_spike_gui first.");
Ks      = load(spikeFile);
roiList = find(logical(Ks.ifSpike(:)'));
caData  = Ks.roi_spikes;
assert(~isempty(roiList), "No spiking ROIs found.");

% --- breath_peak_data.mat ---
bpFile = fullfile(folderPath, 'breath_peak_data.mat');
assert(isfile(bpFile), "breath_peak_data.mat not found — run breathing_peak_gui first.");
BP          = load(bpFile);
insp_onsets = sort(BP.insp_onset_idx(:));
fp          = BP.findpeak_params;
if isfield(fp,'camera_fps'), fs_breath = double(fp.camera_fps);
else,                        fs_breath = FPS; end

% --- DLC CSV → breath signal [T × 1] aligned with dFF ---
csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv'));
assert(~isempty(csv_hits), "No *snapshot_200.csv found.");
dlc_raw  = readmatrix(fullfile(folderPath, csv_hits(1).name), 'NumHeaderLines', 3);
dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
dot_map  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
nDots    = numel(fp.dot_selection);
traces   = NaN(size(dlc_raw,1), nDots);
for d = 1:nDots
    di  = dot_map.(fp.dot_selection{d});
    xc  = dlc_raw(:, dot_cols(di,1));
    yc  = dlc_raw(:, dot_cols(di,2));
    pc  = dlc_raw(:, dot_cols(di,3));
    bad = pc < fp.likelihood_thr;
    switch fp.coord_types{di}
        case 'x',         sig =  xc;
        case 'y',         sig =  yc;
        case '-x',        sig = -xc;
        case '-y',        sig = -yc;
        case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
        otherwise,        sig =  xc;
    end
    sig(bad) = NaN;
    if sum(~isnan(sig)) >= 2
        sig = fillmissing(sig, 'linear', 'EndValues','nearest');
    end
    traces(:,d) = sig;
end
switch fp.combine_method
    case 'sum',  breath_dlc = sum(traces,  2, 'omitnan');
    otherwise,   breath_dlc = mean(traces, 2, 'omitnan');
end
if fp.inverted, breath_dlc = -breath_dlc; end
breath_dlc = detrend(breath_dlc);
breath_dlc = (breath_dlc - mean(breath_dlc)) / std(breath_dlc);
nB = numel(breath_dlc);
if nB >= T, breath_dlc = breath_dlc(1:T);
else,       breath_dlc(end+1:T) = 0; end

% --- chronux path ---
chronuxDir = fullfile(fileparts(mfilename('fullpath')), '..', 'chronux_2_12');
if isfolder(chronuxDir) && ~contains(path, 'chronux')
    addpath(genpath(chronuxDir));
end

fprintf("T=%d  N=%d  %d spiking ROIs  %d insp onsets\n", ...
    T, N, numel(roiList), numel(insp_onsets));

%% ============================== BREATHING PSD + PEAK + FWHM ==============================
params_breath.Fs     = fs_breath;
params_breath.tapers = [TW, 2*TW - 1];
params_breath.pad    = 0;
params_breath.fpass  = [fmin, min(fmax, fs_breath/2)];
params_breath.err    = [2, 0.05];
[Sb, fb, Serr_b] = mtspectrumc(breath_dlc, params_breath);
fbC = fb(:);  Sb = Sb(:);

mask_bp      = fbC >= f_breath_search(1) & fbC <= f_breath_search(2);
[~, rel_idx] = max(Sb(mask_bp));
idx_bp       = find(mask_bp, 1) + rel_idx - 1;
f_breath_peak = fbC(idx_bp);
half_max      = Sb(idx_bp) / 2;

idx_lo = idx_bp;
while idx_lo > 1       && Sb(idx_lo) > half_max, idx_lo = idx_lo - 1; end
idx_hi = idx_bp;
while idx_hi < numel(fbC) && Sb(idx_hi) > half_max, idx_hi = idx_hi + 1; end
f_fwhm = [fbC(idx_lo), fbC(idx_hi)];
fprintf("Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f] Hz\n", ...
    f_breath_peak, f_fwhm(1), f_fwhm(2));

%% ============================== INSPIRATION TIMING ==============================
Delta_k = diff(insp_onsets) / fs;
win_sec = mean(Delta_k) * 2.5;
win     = round(win_sec * fs);
t_w     = (-win:win)' / fs;
nW      = 2*win + 1;

insp_valid = insp_onsets(insp_onsets > win & insp_onsets <= T - win);
nInsp      = numel(insp_valid);
fprintf("Breath-trig window: ±%.3f s  (%d valid onsets)\n", win_sec, nInsp);

% pre-compute breath-triggered dFF for all ROIs
dff_trig = cell(N, 1);
for i = 1:N
    seg = zeros(nInsp, nW);
    for e = 1:nInsp
        c = insp_valid(e);
        seg(e,:) = dFF(c-win : c+win, i);
    end
    dff_trig{i} = seg;
end

% breath-triggered autocorrelogram (same window applied to breathing signal)
breath_trig_seg = zeros(nInsp, nW);
for e = 1:nInsp
    c = insp_valid(e);
    breath_trig_seg(e,:) = breath_dlc(c-win : c+win);
end
breath_autocorr = mean(breath_trig_seg, 1);
if max(abs(breath_autocorr)) > 0
    breath_autocorr = breath_autocorr / max(abs(breath_autocorr));
end

%% ============================== CHRONUX PARAMS ==============================
params_mt.Fs     = fs;
params_mt.tapers = [TW, 2*TW - 1];
params_mt.pad    = 0;
params_mt.fpass  = [fmin, fmax];
params_mt.err    = [2, 0.05];

params_coh.Fs     = fs;
params_coh.tapers = [TW, 2*TW - 1];
params_coh.pad    = 0;
params_coh.fpass  = f_breath_search;
params_coh.err    = [2, 0.01];   % jackknife → phistd + Cerr

% resample breath to calcium fs if needed
if abs(fs_breath - fs) > 0.1
    breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
else
    breath_rs = breath_dlc;
end
T_coh     = min(T, numel(breath_rs));
breath_rs = breath_rs(1:T_coh);

% confC: same for every ROI (depends only on T_coh and tapers)
sig0 = dFF(1:T_coh, roiList(1)) - mean(dFF(1:T_coh, roiList(1)));
[~, ~, ~, ~, ~, ~, ~, confC_val] = coherencyc(breath_rs, sig0, params_coh);

cmap = lines(numel(roiList));

%% ============================== PER-ROI FIGURES ==============================
for k = 1:numel(roiList)
    roi_id = roiList(k);
    col    = cmap(k,:);
    sig_c  = dFF(:, roi_id) - mean(dFF(:, roi_id));

    % --- PSD ---
    [Smt, fmt, Serr_mt] = mtspectrumc(sig_c, params_mt);
    fmt = fmt(:);  Smt = Smt(:);

    % --- coherence ---
    [~, C_coh, phi_coh, ~, ~, ~, f_coh, ~, phistd_coh, Cerr_coh] = ...
        coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
    f_coh = f_coh(:);  C_coh = C_coh(:);  phi_coh = phi_coh(:);

    mask_fwhm = f_coh >= f_fwhm(1) & f_coh <= f_fwhm(2);
    if ~any(mask_fwhm), mask_fwhm = true(size(f_coh)); end
    r_k  = mean(C_coh(mask_fwhm));
    th_k = angle(mean(exp(1i * phi_coh(mask_fwhm))));   % circular mean
    dphi = mean(phistd_coh(mask_fwhm));
    rlo  = max(0,   mean(Cerr_coh(1, mask_fwhm)));
    rhi  = min(1,   mean(Cerr_coh(2, mask_fwhm)));

    % --- breath-trig quantities ---
    ca_ev    = caData(roi_id).spike_idx(:);
    seg_dff  = dff_trig{roi_id};
    dff_mean = mean(seg_dff, 1);

    % spike raster (sorted by proximity to onset)
    ca_rel  = cell(nInsp, 1);
    dt_near = nan(nInsp, 1);
    for e = 1:nInsp
        c    = insp_valid(e);
        spks = ca_ev(ca_ev >= c-win & ca_ev <= c+win);
        trel = (spks - c) / fs;
        ca_rel{e} = trel;
        if ~isempty(trel), dt_near(e) = min(abs(trel)); end
    end
    [~, sort_idx] = sort(dt_near, 'ascend', 'MissingPlacement','last');

    % phase scatter
    dt_last_spk = nan(numel(ca_ev), 1);
    dt_next_spk = nan(numel(ca_ev), 1);
    for e = 1:numel(ca_ev)
        c   = ca_ev(e);
        bef = insp_onsets(insp_onsets <= c);
        aft = insp_onsets(insp_onsets >  c);
        if ~isempty(bef), dt_last_spk(e) = (c - bef(end)) / fs; end
        if ~isempty(aft), dt_next_spk(e) = (aft(1) - c)   / fs; end
    end
    keep3 = ~isnan(dt_last_spk) & ~isnan(dt_next_spk);
    xl3   = dt_last_spk(keep3);
    xn3   = dt_next_spk(keep3);

    % ----------------------------------------------------------------
    % Figure
    % ----------------------------------------------------------------
    fig = figure('Color','w', 'Units','normalized', 'Position',[0.03 0.02 0.48 0.76], ...
        'Name', sprintf('ROI#%02d', roi_id));
    sgtitle(sprintf('%s  —  ROI #%d', strrep(stem,'_','\_'), roi_id), 'FontSize',10);

    % (1,1) breath-triggered dF/F
    ax11 = subplot(3,2,1);
    hold(ax11,'on');
    plot(ax11, t_w, seg_dff', 'Color',[0.82 0.82 0.82], 'LineWidth',0.6);
    plot(ax11, t_w, dff_mean, 'Color', col, 'LineWidth',2);
    xline(ax11, 0, '--', 'Color',[0.4 0.4 0.4]);
    xlim(ax11, [-win_sec, win_sec]);
    xlabel(ax11, 'Time from inspiration (s)');
    ylabel(ax11, '\DeltaF/F');
    % overlay breath autocorrelogram on right axis
    yyaxis(ax11, 'right');
    plot(ax11, t_w, breath_autocorr, 'k--', 'LineWidth',1.5);
    ylabel(ax11, 'Autocorrelation');
    ylim(ax11, [-0.1, 1.1]);
    ax11.YAxis(2).Color = [0 0 0];
    yyaxis(ax11, 'left');
    ax11.YAxis(1).Color = col;
    hold(ax11,'off');  box(ax11,'on');

    % (2,1) spike raster
    ax21 = subplot(3,2,3);
    hold(ax21,'on');
    for row = 1:nInsp
        trel = ca_rel{sort_idx(row)};
        for tt = 1:numel(trel)
            plot(ax21, [trel(tt) trel(tt)], [row-0.4 row+0.4], 'k', 'LineWidth',0.8);
        end
    end
    xline(ax21, 0, 'r--', 'LineWidth',1.2);
    T_breath = 1 / f_breath_peak;   % expected breath period (s)
    xline(ax21, -T_breath, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',1);
    xline(ax21,  T_breath, ':', 'Color',[0.6 0.6 0.6], 'LineWidth',1);
    xlim(ax21, [-win_sec, win_sec]);  ylim(ax21, [0, nInsp+1]);
    set(ax21, 'YDir','reverse');
    xlabel(ax21, 'Time from inspiration (s)');
    ylabel(ax21, 'Breath #');
    hold(ax21,'off');  box(ax21,'on');

    % (3,1) dt_last vs dt_next phase scatter
    ax31 = subplot(3,2,5);
    hold(ax31,'on');
    scatter(ax31, xl3, xn3, 18, 'k', 'filled', 'MarkerFaceAlpha',0.4);
    for dv = mean(Delta_k) + [-std(Delta_k), 0, std(Delta_k)]
        xv = linspace(0, dv, 50);
        plot(ax31, xv, dv - xv, '--', 'Color',[0.72 0.72 0.72], 'LineWidth',0.8);
    end
    all_vals = [xl3(:); xn3(:)];
    ax_lim   = max(mean(Delta_k) + std(Delta_k), max(all_vals(~isnan(all_vals)), [], 'omitnan'));
    plot(ax31, [0 ax_lim], [0 ax_lim], '-', 'Color',[0.2 0.6 0.9], 'LineWidth',1.2);
    mx = mean(xl3);  sx = std(xl3);
    my = mean(xn3);  sy = std(xn3);
    plot(ax31, [mx-sx mx+sx], [my my],     '-', 'Color',[0.85 0.33 0.1], 'LineWidth',2);
    plot(ax31, [mx mx],       [my-sy my+sy],'-', 'Color',[0.85 0.33 0.1], 'LineWidth',2);
    plot(ax31, mx, my, '+', 'Color',[0.85 0.33 0.1], 'MarkerSize',10, 'LineWidth',2);
    xlabel(ax31, 'dt_{last} (s)');  ylabel(ax31, 'dt_{next} (s)');
    hold(ax31,'off');  box(ax31,'on');

    % (1,2) PSD — current ROI + breathing
    ax12 = subplot(3,2,2);
    hold(ax12,'on');
    yLo = Serr_mt(1,:)';  yHi = Serr_mt(2,:)';
    patch(ax12, [fmt; flipud(fmt)], [yLo; flipud(yHi)], col, ...
          'EdgeColor','none', 'FaceAlpha',0.2);
    plot(ax12, fmt, Smt, 'Color', col, 'LineWidth',1.5);
    yLo_b = Serr_b(1,:)';  yHi_b = Serr_b(2,:)';
    patch(ax12, [fbC; flipud(fbC)], [yLo_b; flipud(yHi_b)], [0.3 0.3 0.3], ...
          'EdgeColor','none', 'FaceAlpha',0.15);
    plot(ax12, fbC, Sb, 'k', 'LineWidth',2);
    ax12.XScale = 'log';  ax12.YScale = 'log';
    ax12.XLim   = [fmin, fmax];
    grid(ax12,'on');
    xlabel(ax12, 'Frequency (Hz)');  ylabel(ax12, 'PSD');
    hold(ax12,'off');  box(ax12,'on');

    % (2,2) spike phase histogram + cosine fit
    phi_spikes = 2*pi * xl3 ./ (xl3 + xn3) - pi;   % ∈ [–π, π], –π = at inspiration

    nBins    = 24;
    edges    = linspace(-pi, pi, nBins+1);
    ctrs     = (edges(1:end-1) + edges(2:end)) / 2;
    counts   = histcounts(phi_spikes, edges);

    % wrap to [–3π/2, 3π/2] for display
    left_mask  = ctrs >= pi/2;
    right_mask = ctrs <= -pi/2;
    ctrs_ext   = [ctrs(left_mask)-2*pi,  ctrs,  ctrs(right_mask)+2*pi];
    cnts_ext   = [counts(left_mask),     counts, counts(right_mask)];

    % cosine fit via circular mean
    r_vec = mean(exp(1i * phi_spikes));
    phi0  = angle(r_vec);
    Rlen  = abs(r_vec);
    B     = numel(phi_spikes) / nBins;
    A     = 2 * B * Rlen;
    th_fit  = linspace(-1.5*pi, 1.5*pi, 300);
    fit_ext = B + A * cos(th_fit - phi0);

    ax22 = subplot(3,2,4);
    hold(ax22,'on');
    bar(ax22, ctrs_ext, cnts_ext, 1, 'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
    plot(ax22, th_fit, fit_ext, 'Color', col, 'LineWidth',2);
    xline(ax22, 0, 'r--', 'LineWidth',1.2);
    xline(ax22, -pi, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    xline(ax22,  pi, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    ax22.XLim   = [-1.5*pi, 1.5*pi];
    ax22.XTick  = [-pi -pi/2 0 pi/2 pi];
    ax22.XTickLabel = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    xlabel(ax22, 'Breathing phase (rad)');
    ylabel(ax22, 'Spike count');
    hold(ax22,'off');  box(ax22,'on');

    % (3,2) polar — phase / coherence for this ROI
    t6  = subplot(3,2,6);
    pos = t6.Position;
    delete(t6);
    ax32 = polaraxes(fig, 'Position', pos);
    hold(ax32,'on');
    % significance circle
    th_circ = linspace(0, 2*pi, 360);
    polarplot(ax32, th_circ, repmat(confC_val, 1, 360), 'k--', 'LineWidth',1.2);
    % radial CI bar
    polarplot(ax32, [th_k th_k], [rlo rhi], '-', 'Color', col, 'LineWidth',2);
    % phase arc (±phistd at constant radius)
    phi_arc = linspace(th_k - dphi, th_k + dphi, 40);
    polarplot(ax32, phi_arc, r_k * ones(1,40), '-', 'Color', col, 'LineWidth',1.5);
    % dot
    ph = polarplot(ax32, th_k, r_k, 'o', ...
                   'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'MarkerSize',8);
    ph.UserData = struct('roi', roi_id, 'phistd', dphi, 'rlo', rlo, 'rhi', rhi);
    ax32.RLim              = [0, 1];
    ax32.ThetaZeroLocation = 'top';
    ax32.ThetaDir          = 'counterclockwise';
    hold(ax32,'off');

    % data cursor on polar dot
    dcm = datacursormode(fig);
    dcm.UpdateFcn = @(~, info) polarTip(info);

    % Manually square all cartesian axes by cropping Position to min(w,h)
    % (axis square leaves Position non-square; this guarantees ax.Position IS the square box)
    drawnow;
    for axi = {ax11, ax21, ax31, ax12, ax22}
        a = axi{1};
        p = a.Position;
        s = min(p(3), p(4));
        a.Position = [p(1)+(p(3)-s)/2, p(2)+(p(4)-s)/2, s, s];
    end
    drawnow;

    % Coherence inset in bottom-left corner of PSD tile (ax12)
    rp  = ax12.Position;   % [left bottom s s] — now truly square
    inW = 0.38 * rp(3);    inH = 0.38 * rp(4);
    inX = rp(1) + 0.03 * rp(3);
    inY = rp(2) + 0.03 * rp(4);
    ax_ins = axes(fig, 'Position',[inX inY inW inH]);
    hold(ax_ins,'on');
    plot(ax_ins, f_coh, C_coh, 'Color', col, 'LineWidth',1.2);
    patch(ax_ins, [f_fwhm(1) f_fwhm(2) f_fwhm(2) f_fwhm(1)],[0 0 1 1], ...
          [1 0.6 0.6],'EdgeColor','none','FaceAlpha',0.25);
    xline(ax_ins, f_breath_peak,'r:','LineWidth',1.2);
    yline(ax_ins, confC_val,'k--','LineWidth',1);
    % center xlim on peak so FWHM patch appears centered
    hw = diff(f_breath_search) * 0.6;
    ax_ins.XLim  = [f_breath_peak - hw, f_breath_peak + hw];
    ax_ins.YLim  = [0, 1];
    ax_ins.XScale = 'linear';
    ax_ins.FontSize = 7;
    ax_ins.Color = [1 1 1 0.85];
    grid(ax_ins,'on');
    xlabel(ax_ins,'Hz','FontSize',7);  ylabel(ax_ins,'Coh','FontSize',7);
    hold(ax_ins,'off');  box(ax_ins,'on');
    axis(ax_ins,'square');

    if saveFigs
        outPdf = fullfile(folderPath, sprintf('%s_ROI%02d_combined.pdf', stem, roi_id));
        exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');
        fprintf("Saved: %s\n", outPdf);
    end
end

fprintf("Done.\n");

%% --- datacursor callback ---
function txt = polarTip(info)
    src = info.Target;
    ud  = [];
    if isprop(src,'UserData'), ud = src.UserData; end
    th = info.Position(1);
    r  = info.Position(2);
    if isstruct(ud) && isfield(ud,'roi')
        txt = sprintf('ROI %d\nPhase = %.3g rad (%.1f\xB0) \xB1 %.3g rad\nCoh = %.3g  [%.3g, %.3g]', ...
                      ud.roi, th, rad2deg(th), ud.phistd, r, ud.rlo, ud.rhi);
    else
        txt = sprintf('Phase = %.3g rad (%.1f\xB0)\nCoh = %.3g', th, rad2deg(th), r);
    end
end
