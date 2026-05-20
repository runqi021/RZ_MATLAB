clear; close all; clc;
%% Batch_breath_trigAVG_260325.m  –  Batch breathing-triggered analysis
%  Loops through all subfolders of masterFolder that contain ca_spike_data.mat
%  (i.e. sessions with responsive ROIs). Runs the full breathing-triggered
%  analysis pipeline per session: spectrum, stacked dFF, triggered averages,
%  dt_last raster, phase analysis, avg projection overlay.
%
%  ── USER PARAMETERS ──────────────────────────────────────────────────────
masterFolder = 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing';
skipIdx      = [];          % session indices to skip
fps_img      = 30;          % imaging frame rate (Hz)
nDrop        = 30;          % frames to drop from start
dFF_scale    = 0.20;        % dF/F scale bar value (20%)
% ──────────────────────────────────────────────────────────────────────────

%% 1 ── Discover sessions with ca_spike_data.mat
allMat = dir(fullfile(masterFolder, '**', 'ca_spike_data.mat'));

skipIdx = flip(skipIdx);
keep = true(numel(allMat),1);
keep(skipIdx) = false;
allMat = allMat(keep);

nSess = numel(allMat);
fprintf('Found %d sessions with ca_spike_data.mat\n', nSess);
assert(nSess > 0, 'No ca_spike_data.mat found under %s', masterFolder);

nProcessed = 0;
nSkipped   = 0;
nFailed    = 0;
failedSessions = {};

%% 2 ── Main loop
for kk_sess = 1:nSess
    folderPath = allMat(kk_sess).folder;
    [~, sessName] = fileparts(folderPath);
    fprintf('\n========== [%d/%d] %s ==========\n', kk_sess, nSess, sessName);
    close all;

    try

    %% ---- Output subfolder ----
    outDir = fullfile(folderPath, 'breath_trigger_analysis');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    %% ---- Auto-find required files ----
    csv_hits = dir(fullfile(folderPath, '*snapshot_best-10.csv'));
    sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
    bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

    if isempty(csv_hits)
        warning('No DLC CSV found in %s, skipping.', sessName);
        nSkipped = nSkipped + 1; continue;
    end
    if isempty(sam_hits)
        warning('No SAM output found in %s, skipping.', sessName);
        nSkipped = nSkipped + 1; continue;
    end
    if isempty(bp_file)
        warning('No breath_peak_data.mat found in %s, skipping.', sessName);
        nSkipped = nSkipped + 1; continue;
    end

    dlc_csv = fullfile(folderPath, csv_hits(1).name);
    SAM     = fullfile(folderPath, sam_hits(1).name);

    fprintf('DLC CSV : %s\n', csv_hits(1).name);
    fprintf('SAM     : %s\n', sam_hits(1).name);

    %% ---- Load SAM → dFF ----
    SAMload = load(SAM);
    F_raw = SAMload.F;
    F = F_raw;
    %F(1:nDrop, :) = [];
    dFFout = helper.dFF_RZ(F);
    dFF = dFFout.dFF;

    %% ---- Basic params ----
    [T, N_roi] = size(dFF);
    fps_breath = fps_img;
    fs = fps_img;
    t_img = (0:T-1)' / fps_img;
    roi_ids = arrayfun(@(k) sprintf('%02d', k), 1:N_roi, 'UniformOutput', false);

    %% ---- Load calcium spike ----
    in = load(fullfile(folderPath, 'ca_spike_data.mat'));
    ifSpk = in.ifSpike; roiSpk_id = find(ifSpk == 1);
    ca_spk_data = in.roi_spikes;

    nSpkROI = numel(roiSpk_id);

    ca_spk_train = zeros(T, nSpkROI, 'double');
    ca_spk_id = {}; t_spk = {};

    for k = 1:nSpkROI
        i = roiSpk_id(k);
        ca_spk_train(:, k) = ca_spk_data(i).spike_train;
        ca_spk_id{k} = ca_spk_data(i).spike_idx;
        t_spk{k} = ca_spk_data(i).spike_t;
    end

    %% ---- Load breath peak ----
    breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));

    insp_onsets = breath_peak.insp_onset_idx;
    Sb = breath_peak.insp_onsets_train;

    insp_onsets(insp_onsets<nDrop) = [];
    insp_onsets = insp_onsets - nDrop;
    Sb(1:nDrop) = [];

    %% ---- DLC breathing input ----
    dlc_raw = readmatrix(dlc_csv, 'NumHeaderLines', 3);

    data_dlc = dlc_raw;
    data_dlc(1:nDrop, :) = [];
    frame_idx = data_dlc(:,1);
    t_breath  = frame_idx / fps_breath;

    dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
    dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);

    fp          = breath_peak.findpeak_params;
    nF_dlc      = size(data_dlc, 1);
    nActiveDots = numel(fp.dot_selection);
    traces      = NaN(nF_dlc, nActiveDots);

    fprintf('Breath signal: %d dots, combine=%s, invert=%d, lk_thr=%.2f\n', ...
        nActiveDots, fp.combine_method, fp.inverted, fp.likelihood_thr);

    for d = 1:nActiveDots
        di  = dot_idx.(fp.dot_selection{d});
        xc  = data_dlc(:, dot_cols(di,1));
        yc  = data_dlc(:, dot_cols(di,2));
        pc  = data_dlc(:, dot_cols(di,3));
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
            sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');
        end
        traces(:,d) = sig;
        fprintf('  %s -> %s\n', fp.dot_selection{d}, fp.coord_types{di});
    end

    switch fp.combine_method
        case 'sum',  breath = sum(traces,  2, 'omitnan');
        otherwise,   breath = mean(traces, 2, 'omitnan');
    end

    if fp.inverted, breath = -breath; end

    breath = detrend(breath);
    breath = (breath - mean(breath)) / std(breath);

    fprintf('Breathing trace: %d samples after dropping first %d frames\n', ...
            numel(breath), nDrop);

    %% ---- Compute Delta_k early (needed by Figs 2,3,4,5) ----
    b_frames = sort(insp_onsets(:));
    Delta_k  = diff(b_frames) / fs;       % breath interval durations (s)

    if numel(insp_onsets) < 2
        warning('Fewer than 2 inspiration onsets in %s, skipping.', sessName);
        nSkipped = nSkipped + 1; continue;
    end

    %% ---- Chronux spectrum ----
    TW_spec = 4;
    params_spec.Fs     = fps_img;
    params_spec.tapers = [TW_spec, 2*TW_spec - 1];
    params_spec.pad    = 0;
    params_spec.fpass  = [0.1 15];
    params_spec.err    = [2 0.05];

    [Sk_raw, fk_raw] = mtspectrumc(breath, params_spec);
    figure('Color','w');
    plot(fk_raw, Sk_raw, 'k', 'LineWidth', 1.2);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    xlim([0 10]);
    grid('on');
    axis square;

    outFig = fullfile(outDir, 'breathing_spectrum_raw.pdf');
    exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
    fprintf('Saved breathing spectrum:\n  %s\n', outFig);

    %% ---- Figure 1: Breath + stacked dFF (all ROIs, FULL recording) ----
    T_img = numel(breath);

    amp = max(dFF(:)) - min(dFF(:));
    if amp == 0, amp = 1; end; gap = 0.8 * amp;

    Y = zeros(T_img, N_roi);
    for k = 1:N_roi
        Y(:,k) = dFF(:,k) + gap*(N_roi - k);
    end

    yt  = 0:gap:gap*(N_roi-1);
    ytl = roi_ids(end:-1:1);

    figure('Position',[200 200 900 700], 'Color','w');
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    % TOP PANEL: breathing + breathing events (Sb)
    ax1 = nexttile(1);
    plot(ax1, t_img, -breath, 'k', 'LineWidth', 1); hold(ax1,'on');

    idx_b = find(Sb(:) ~= 0);
    if ~isempty(idx_b)
        t_b = t_img(idx_b);
        y_b = breath(idx_b);
        plot(ax1, t_b, -y_b, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 3);
    end
    xlim(ax1, [t_img(1) t_img(end)]);
    xlabel(ax1, 'Time (s)');
    ylabel(ax1, 'Breathing (z)');
    title(ax1, 'Breathing trace (DLC) with detected inspirations');
    grid(ax1,'off');
    hold(ax1,'off');

    % BOTTOM PANEL: dFF
    ax2 = nexttile(2, [2 1]);
    hold(ax2,'on');

    for kk = 1:N_roi
        if ismember(kk, roiSpk_id)
            trace_color = [1 0 0];
        else
            trace_color = [0 0 0];
        end
        plot(ax2, t_img, Y(:,kk), 'Color', trace_color, 'LineWidth', 0.5);
    end

    xlim(ax2, [t_img(1) t_img(end)]);
    ylim(ax2, [-gap, gap*(N_roi-1)+gap]);

    yticks(ax2, yt);
    yticklabels(ax2, ytl);

    xlabel(ax2, 'Time (s)');
    ylabel(ax2, 'ROI#');
    title(ax2, 'dF/F');
    grid(ax2,'off');

    % SCALE BAR
    xl_ax2 = xlim(ax2);
    x_sb = xl_ax2(2) - 2;
    y_sb = -0.3 * gap;

    plot(ax2, [x_sb x_sb], [y_sb y_sb+dFF_scale], 'k', 'LineWidth', 2);
    text(ax2, x_sb + 0.5, y_sb + dFF_scale/2, ...
        sprintf('%.1f dF/F', dFF_scale), ...
        'VerticalAlignment','middle','HorizontalAlignment','left', ...
        'FontSize', 9);
    hold(ax2,'off');

    outFig = fullfile(outDir, 'breath_N_dFF.pdf');
    exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
    fprintf('Saved breathing and dFF traces:\n  %s\n', outFig);

    %% ---- Figure 2: Breath + spiking ROIs + spike ticks (30 breaths) ----
    xl2 = min(30 * mean(Delta_k), t_img(end));

    Y2 = zeros(T_img, nSpkROI);
    for k = 1:nSpkROI
        i = roiSpk_id(k);
        Y2(:,k) = dFF(:,i) + 1.5*gap*(nSpkROI - k);
    end

    yt_spk  = fliplr(1.5*gap*(nSpkROI - (1:nSpkROI)));
    ytl_spk = arrayfun(@(id) sprintf('ROI#%02d', id), roiSpk_id(end:-1:1), 'UniformOutput', false);

    figure('Position',[200 200 900 700], 'Color','w');
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    ax1 = nexttile(1);
    plot(ax1, t_img, -breath, 'Color', 'k', 'LineWidth',1.2); xlim([0 xl2]);
    xline(ax1, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

    ax2 = nexttile(2, [2 1]);
    plot(ax2, t_img, Y2, 'Color', 'k', 'LineWidth',0.8); xlim([0 xl2]);
    xline(ax2, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

    off = 1.5*gap*(nSpkROI - (1:nSpkROI));

    tickCenter = off - 0.50*gap;
    tickHalf   = 0.12*gap;

    hold on
    for k = 1:nSpkROI
        spk = find(ca_spk_train(:,k));
        if isempty(spk), continue; end

        x = spk;
        xx = [x(:) x(:)]';
        yy = [(tickCenter(k)-tickHalf)*ones(numel(spk),1), ...
              (tickCenter(k)+tickHalf)*ones(numel(spk),1)]';

        plot(xx/fs, yy, 'k', 'LineWidth', 1);
    end

    % SCALE BAR
    x_sb2 = xl2 - 1;
    y_sb2 = -0.5 * 1.5*gap;
    plot(ax2, [x_sb2 x_sb2], [y_sb2 y_sb2+dFF_scale], 'k', 'LineWidth', 2);
    text(ax2, x_sb2 + 0.3, y_sb2 + dFF_scale/2, ...
        sprintf('%.1f dF/F', dFF_scale), ...
        'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize', 9);

    yticks(ax2, yt_spk);
    yticklabels(ax2, ytl_spk);
    ylabel(ax2, 'ROI#');
    xlabel('Time (s)')

    outFig = fullfile(outDir, sprintf('breath_ca_spk_%ds.pdf', round(xl2)));
    exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
    fprintf('Saved breath+ca spk (%ds):\n  %s\n', round(xl2), outFig);

    %% ---- Figure 3: Inspiration-triggered dFF (win = 1 mean breath cycle) ----
    win  = round(mean(Delta_k) * fs);     % +/- 1 mean breath cycle (frames)
    t_evt = (-win:win) / fs;

    T_img = size(dFF,1);
    insp_ev = insp_onsets(insp_onsets > win & insp_onsets < T_img - win);
    nEv = numel(insp_ev);

    fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);

    dff_evt_seg_all = cell(N_roi,1);

    for i = 1:N_roi
        roi_trace = dFF(:, i);

        seg = zeros(nEv, 2*win+1);
        for e = 1:nEv
            c = insp_ev(e);
            seg(e,:) = roi_trace(c-win : c+win);
        end

        dff_evt_seg_all{i} = seg;
    end

    % Per-ROI event-triggered overlays (all trials + mean)
    figure('Name', (sprintf('Event-triggered dFF (n=%d breaths)', nEv)), 'Color','w');

    for i = 1:N_roi
        seg = dff_evt_seg_all{i};
        if isempty(seg), continue; end

        dff_avg = mean(seg, 1);

        subplot(ceil(N_roi/10), 10, i);
        hold on;

        plot(t_evt, seg', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.7);
        plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);

        xline(0, '--', 'Color', [0.4 0.4 0.4]);
        yline(0, '--', 'Color', [0.6 0.6 0.6]);

        xlabel('Time from inspiration onset (s)');
        ylabel('dF/F');
        title(sprintf('ROI#%d', i));

        hold off;
        axis square;
    end

    outFig = fullfile(outDir, 'triggered_avg.pdf');
    exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
    fprintf('Saved triggered average:\n  %s\n', outFig);

    %% ---- Figure 4: Per-ROI calcium-triggered breathing overlays ----
    win_ca = round(mean(Delta_k) * fs);   % +/- 1 mean breath cycle (frames)
    t_ca_evt = (-win_ca:win_ca) / fs;

    for k = 1:nSpkROI
        roi_id = roiSpk_id(k);
        ca_frames = ca_spk_id{k};
        nCaEv = numel(ca_frames);

        if nCaEv == 0, continue; end

        % LEFT: breath-triggered dFF overlay
        seg_breath = dff_evt_seg_all{roi_id};
        dff_avg = mean(seg_breath, 1);

        % RIGHT: calcium-triggered breathing segments
        ca_valid = ca_frames(ca_frames > win_ca & ca_frames < T_img - win_ca);
        nCaValid = numel(ca_valid);
        seg_ca = zeros(nCaValid, 2*win_ca+1);
        for e = 1:nCaValid
            c = ca_valid(e);
            seg_ca(e,:) = breath(c-win_ca : c+win_ca);
        end
        breath_avg = mean(seg_ca, 1);

        fig4 = figure('Color','w', 'Name', sprintf('ROI#%02d triggered avg', roi_id), ...
                      'Position', [100 100 900 400]);

        % LEFT subplot: breath-triggered dFF
        subplot(1,2,1);
        hold on;
        plot(t_evt, seg_breath', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
        plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);
        xline(0, '--', 'Color', [0.4 0.4 0.4]);
        yline(0, '--', 'Color', [0.6 0.6 0.6]);
        xlabel('Time from inspiration onset (s)');
        ylabel('dF/F');
        title(sprintf('ROI#%02d  breath-trig dFF (n=%d)', roi_id, nEv));
        axis square;
        hold off;

        % RIGHT subplot: calcium-triggered breathing
        subplot(1,2,2);
        hold on;
        if nCaValid > 0
            plot(t_ca_evt, -seg_ca', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
            plot(t_ca_evt, -breath_avg, 'Color', 'k', 'LineWidth', 2);
        end
        xline(0, '--', 'Color', [0.4 0.4 0.4]);
        yline(0, '--', 'Color', [0.6 0.6 0.6]);
        xlabel('Time from Ca^{2+} spike (s)');
        ylabel('Breathing (z-score)');
        title(sprintf('ROI#%02d  Ca-trig breath (n=%d)', roi_id, nCaValid));
        axis square;
        hold off;

        outPdf = fullfile(outDir, sprintf('ROI%02d_triggered_avg.pdf', roi_id));
        exportgraphics(fig4, outPdf, 'ContentType','vector', 'BackgroundColor','none');
        outFigFile = fullfile(outDir, sprintf('ROI%02d_triggered_avg.fig', roi_id));
        savefig(fig4, outFigFile);
        fprintf('Saved triggered avg for ROI#%02d:\n  %s\n  %s\n', roi_id, outPdf, outFigFile);
    end

    %% ---- Figure 5: dt_last analysis (raster + histogram) ----
    win_raster = mean(Delta_k)*2;       % seconds
    win5       = round(win_raster * fs); % frames

    dtlast_results = struct();

    for k = 1:nSpkROI
        roi_id   = roiSpk_id(k);
        ca_ev    = ca_spk_id{k}(:);
        nCaTotal = numel(ca_ev);

        if nCaTotal == 0, continue; end

        % Assign each Ca event to its breath interval
        [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
        keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
        ca_in   = ca_ev(keep_ca);
        bin_in  = bin_idx(keep_ca);
        nEv_in  = numel(ca_in);

        if nEv_in == 0
            fprintf('ROI#%02d : no calcium events within breath intervals, skipping.\n', roi_id);
            continue;
        end

        % dt_last: time from last inspiration (at or before)
        dt_last = (ca_in - b_frames(bin_in)) / fs;

        % Build insp_rel for raster
        insp_rel_in = cell(nEv_in, 1);
        for e = 1:nEv_in
            c      = ca_in(e);
            nearby = insp_onsets(insp_onsets >= (c-win5) & insp_onsets <= (c+win5));
            insp_rel_in{e} = (nearby - c) / fs;
        end

        % Histogram bins (one-sided, 0 to max interval)
        half_max   = max(Delta_k);
        bin_number = round(half_max / (1/fs));
        if bin_number < 10, bin_number = 10; end
        if mod(bin_number,2), bin_number = bin_number + 1; end

        bin_edges_dt = linspace(0, half_max, bin_number);
        H_obs        = histcounts(dt_last, bin_edges_dt, 'Normalization','probability') * 100;

        fprintf('ROI#%02d : %d Ca events in %d breath intervals\n', roi_id, nEv_in, numel(Delta_k));

        % Figure: raster (2 cols) + histogram (1 col)
        fig5 = figure('Color', 'w', ...
            'Name', sprintf('dt_last: ROI#%02d', roi_id), ...
            'Position', [0 0 1100 400]);

        % Subplot 1: Raster sorted by dt_last ascending
        ax_rast = subplot(1,3,[1 2]);
        hold(ax_rast, 'on');

        [~, sort_idx] = sort(dt_last, 'ascend');

        for row = 1:nEv_in
            e       = sort_idx(row);
            t_ticks = insp_rel_in{e};
            for tt = 1:numel(t_ticks)
                plot(ax_rast, [t_ticks(tt) t_ticks(tt)], [row-0.4  row+0.4], ...
                    'k', 'LineWidth', 1);
            end
        end

        xline(ax_rast, 0, 'k--', 'LineWidth', 1);
        xlim(ax_rast, [-win_raster, win_raster]);
        ylim(ax_rast, [0 nEv_in+1]);
        set(ax_rast, 'YDir', 'reverse');
        xlabel(ax_rast, 'Time from Ca^{2+} event (s)');
        ylabel(ax_rast, sprintf('Ca^{2+} event # \n (sorted by dt_{last})'), 'HorizontalAlignment', 'center');
        title(ax_rast, sprintf('ROI#%02d  (%d events in %d breaths)', roi_id, nEv_in, sum(Sb)));
        hold(ax_rast, 'off');

        % Subplot 2: dt_last histogram
        ax_hist = subplot(1,3,3);
        hold(ax_hist, 'on');

        stairs(ax_hist, bin_edges_dt, [H_obs H_obs(end)], 'Color', 'k', 'LineWidth', 2);
        histogram(ax_hist, 'BinEdges', bin_edges_dt, 'BinCounts', H_obs, ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        xlabel(ax_hist, 'Time from last inspiration (s)');
        ylabel(ax_hist, 'Proportion of Ca^{2+} events (%)');
        xlim(ax_hist, [0 half_max + 0.1]);
        ylim(ax_hist, [0 45]);
        axis(ax_hist, 'square');
        hold(ax_hist, 'off');

        % Save figure
        outFig_fig = fullfile(outDir, sprintf('calcium_breath_dtlast_ROI%02d.fig', roi_id));
        savefig(fig5, outFig_fig);
        outFig_pdf = fullfile(outDir, sprintf('calcium_breath_dtlast_ROI%02d.pdf', roi_id));
        exportgraphics(fig5, outFig_pdf, 'ContentType', 'vector', 'BackgroundColor','none');
        fprintf('Saved dt_last figure for ROI#%02d:\n  %s\n  %s\n', roi_id, outFig_fig, outFig_pdf);

        % Store per-ROI results
        dtlast_results(k).roi_id    = roi_id;
        dtlast_results(k).dt_last   = dt_last;
        dtlast_results(k).ca_in     = ca_in;
        dtlast_results(k).bin_in    = bin_in;
        dtlast_results(k).nEv_in    = nEv_in;
        dtlast_results(k).H_obs     = H_obs;
        dtlast_results(k).bin_edges = bin_edges_dt;
        dtlast_results(k).half_max  = half_max;
    end

    %% ---- Filter breathing signal for Hilbert ----
    d = designfilt("lowpassfir", ...
        PassbandFrequency=0.2,StopbandFrequency=0.25, ...
        PassbandRipple=1,StopbandAttenuation=60, ...
        DesignMethod="equiripple");
    breath_filt = filtfilt(d,breath);

    %% ---- Hilbert transform → phase ----
    y = hilbert(breath_filt);
    phi = angle(y);

    %% ---- Figure 6: Phase histogram per spiking ROI (KDE + cosine fit) ----
    nBins_phase = 40;
    bin_edges_phase = linspace(-pi, pi, nBins_phase+1);
    bin_ctrs_phase  = (bin_edges_phase(1:end-1) + bin_edges_phase(2:end)) / 2;

    nKDE = 200;
    xi_kde = linspace(-pi, pi, nKDE);

    for k = 1:nSpkROI
        roi_id = roiSpk_id(k);
        spk_frames = ca_spk_id{k};
        spk_frames = spk_frames(spk_frames >= 1 & spk_frames <= numel(phi));
        if isempty(spk_frames), continue; end

        ca_phi = phi(spk_frames);
        nPh    = numel(ca_phi);

        % Circular statistics
        mu_vec  = mean(exp(1i * ca_phi));
        theta0  = angle(mu_vec);
        R       = abs(mu_vec);

        % Histogram (percentage)
        H_pct = histcounts(ca_phi, bin_edges_phase, 'Normalization','percentage');

        % KDE (circular)
        binW   = 2*pi/nBins_phase;
        xi_kde = linspace(-pi, pi, nKDE).';

        [~,~,bw0] = ksdensity(ca_phi(:));
        bw0 = min(max(bw0, 0.05), 1.0);

        ca_phi_wrap = [ca_phi(:); ca_phi(:)+2*pi; ca_phi(:)-2*pi];
        f_kde = ksdensity(ca_phi_wrap, xi_kde, 'Bandwidth', bw0);

        f_kde_pct = f_kde / trapz(xi_kde, f_kde) * 100 * binW;

        % Cosine fit on histogram
        X_cos = [cos(bin_ctrs_phase(:) - theta0), ones(nBins_phase, 1)];
        coeff = X_cos \ H_pct(:);
        A_fit = coeff(1);
        B_fit = coeff(2);

        th_fine = linspace(-pi, pi, 200);
        y_fine  = A_fit * cos(th_fine - theta0) + B_fit;

        % Figure: 2 panels side by side
        fig6 = figure('Color','w', 'Position', [100 100 900 420]);

        % LEFT: KDE
        ax_kde = subplot(1,2,1);
        hold(ax_kde, 'on');

        histogram(ax_kde, 'BinEdges', bin_edges_phase, 'BinCounts', H_pct, ...
            'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

        plot(ax_kde, xi_kde, f_kde_pct, 'k', 'LineWidth', 2);

        [kde_pks, kde_locs] = findpeaks(f_kde_pct, xi_kde, 'MinPeakProminence', 0.3);
        xline(ax_kde, 0, '--', 'Color', [0.6 0.6 0.6]);
        xlim(ax_kde, [-pi pi]);
        ylim([0 25]);
        xticks(ax_kde, [-pi, -pi/2, 0, pi/2, pi]);
        xticklabels(ax_kde, {'-pi', '-pi/2', '0', 'pi/2', 'pi'});
        xlabel(ax_kde, 'Phase (rad)');
        ylabel(ax_kde, '% events per bin');
        nKDE_peaks = numel(kde_pks);
        if nKDE_peaks > 1
            title(ax_kde, sprintf('ROI#%02d  KDE (%d peaks)', roi_id, nKDE_peaks));
        else
            title(ax_kde, sprintf('ROI#%02d  KDE', roi_id));
        end
        axis(ax_kde, 'square');
        hold(ax_kde, 'off');

        % RIGHT: Cosine fit
        ax_cos = subplot(1,2,2);
        hold(ax_cos, 'on');

        histogram(ax_cos, 'BinEdges', bin_edges_phase, 'BinCounts', H_pct, ...
            'FaceColor', [0.7 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

        plot(ax_cos, th_fine, y_fine, 'r', 'LineWidth', 2);

        xline(ax_cos, 0, '--', 'Color', [0.6 0.6 0.6]);
        xlim(ax_cos, [-pi pi]);
        ylim([0 25]);
        xticks(ax_cos, [-pi, -pi/2, 0, pi/2, pi]);
        xticklabels(ax_cos, {'-pi', '-pi/2', '0', 'pi/2', 'pi'});
        xlabel(ax_cos, 'Phase (rad)');
        ylabel(ax_cos, '% events per bin');
        title(ax_cos, sprintf('ROI#%02d  Cosine fit', roi_id));
        axis(ax_cos, 'square');
        hold(ax_cos, 'off');

        % Pack parameters (per ROI)
        [gaussian_peak_pct, idxMax] = max(f_kde_pct);
        gaussian_peak_phase = xi_kde(idxMax);

        cos_baseline  = B_fit;
        cos_amplitude = A_fit;
        cos_peak      = A_fit + B_fit;

        phase_params = struct();
        phase_params.roi_id = roi_id;
        phase_params.n      = nPh;
        phase_params.theta0 = theta0;
        phase_params.R      = R;

        phase_params.gaussian_peak_phase = gaussian_peak_phase;
        phase_params.gaussian_peak_pct   = gaussian_peak_pct;

        phase_params.kde_pks  = kde_pks(:);
        phase_params.kde_locs = kde_locs(:);

        phase_params.cos_baseline  = cos_baseline;
        phase_params.cos_amplitude = cos_amplitude;
        phase_params.cos_peak      = cos_peak;

        % Save
        outPdf6 = fullfile(outDir, sprintf('phase_hist_ROI%02d.pdf', roi_id));
        exportgraphics(fig6, outPdf6, 'ContentType','vector', 'BackgroundColor','none');
        outFig6 = fullfile(outDir, sprintf('phase_hist_ROI%02d.fig', roi_id));
        savefig(fig6, outFig6);
        outMat6 = fullfile(outDir, sprintf('phase_params_ROI%02d.mat', roi_id));
        save(outMat6, 'phase_params');
        fprintf('ROI#%02d : theta0=%.2f rad, R=%.3f, %d KDE peak(s)\n', ...
            roi_id, theta0, R, nKDE_peaks);
    end

    %% ---- Figure 7: Average projection with ROI mask overlay ----
    avg_hits = dir(fullfile(folderPath, '*_minusDark_MC_AVG_for_CP.tif'));
    if ~isempty(avg_hits)
        avgImg   = double(imread(fullfile(folderPath, avg_hits(1).name)));

        lo = prctile(avgImg(:), 1);
        hi = prctile(avgImg(:), 99);
        avgImg = (avgImg - lo) / (hi - lo);
        avgImg(avgImg < 0) = 0;
        avgImg(avgImg > 1) = 1;

        avgRGB = repmat(avgImg, [1 1 3]);

        maskL    = SAMload.maskL;
        outlines = bwperim(maskL > 0);
        yellow   = [1 1 0];
        for ch = 1:3
            plane = avgRGB(:,:,ch);
            plane(outlines) = yellow(ch);
            avgRGB(:,:,ch) = plane;
        end

        % Save 1: no labels
        outNoLabel = fullfile(outDir, 'avg_projection_ROImask.png');
        imwrite(avgRGB, outNoLabel);
        fprintf('Saved avg projection (no labels):\n  %s\n', outNoLabel);

        % Save 2: with ROI labels
        [imgH, imgW, ~] = size(avgRGB);
        dpi = 150;
        fig7 = figure('Color','k', ...
            'Units','inches', 'Position',[1 1 imgW/dpi imgH/dpi]);
        ax7 = axes(fig7, 'Units','normalized', 'Position',[0 0 1 1]);
        imshow(avgRGB, 'Parent', ax7, 'Border','tight');
        hold(ax7, 'on');

        offset_px  = 35;
        dirs       = [0 -1; 1 -1; 1 0; 1 1;
                      0  1;-1  1;-1 0;-1 -1];
        dirs       = dirs ./ vecnorm(dirs, 2, 2);
        txt_hw     = [20 12];
        placed_boxes = zeros(0, 4);

        for k = 1:nSpkROI
            roi_id = roiSpk_id(k);
            [r, c] = find(maskL == roi_id);
            if isempty(r), continue; end
            cy = mean(r);  cx = mean(c);

            outward = [cx - imgW/2, cy - imgH/2];
            outward = outward / (norm(outward) + eps);

            best_pos   = [cx, cy];
            best_score = -Inf;

            for dd = 1:size(dirs, 1)
                tx = cx + offset_px * dirs(dd, 1);
                ty = cy + offset_px * dirs(dd, 2);
                if tx < 1 || tx > imgW || ty < 1 || ty > imgH, continue; end

                bx = [tx - txt_hw(1)/2, ty - txt_hw(2)/2, ...
                      tx + txt_hw(1)/2, ty + txt_hw(2)/2];

                ov = 0;
                for p = 1:size(placed_boxes, 1)
                    pb = placed_boxes(p, :);
                    if bx(1)<pb(3) && bx(3)>pb(1) && bx(2)<pb(4) && bx(4)>pb(2)
                        ov = ov + 1;
                    end
                end

                txi = round(max(1, min(imgW, tx)));
                tyi = round(max(1, min(imgH, ty)));
                on_roi = maskL(tyi, txi) > 0;

                score = -ov*100 - on_roi*50 + dot(dirs(dd,:), outward)*10;
                if score > best_score
                    best_score = score;
                    best_pos   = [tx, ty];
                end
            end

            text(ax7, best_pos(1), best_pos(2), sprintf('%d', roi_id), ...
                'Color','y', 'FontSize', 8, 'FontWeight','bold', ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle');

            placed_boxes(end+1, :) = [best_pos(1)-txt_hw(1)/2, best_pos(2)-txt_hw(2)/2, ...
                                      best_pos(1)+txt_hw(1)/2, best_pos(2)+txt_hw(2)/2]; %#ok<AGROW>
        end
        hold(ax7, 'off');

        outLabeled = fullfile(outDir, 'avg_projection_ROImask_labeled.png');
        exportgraphics(fig7, outLabeled, 'ContentType','image', 'Resolution', 600);
        fprintf('Saved avg projection (labeled):\n  %s\n', outLabeled);
    else
        warning('No avg projection TIF found in %s, skipping Figure 7.', sessName);
    end

    %% ---- Save master file ----
    outMat_master = fullfile(outDir, 'breath_master_working.mat');
    save(outMat_master, ...
        'folderPath', ...
        'fps_img', 'fps_breath', 'fs', ...
        'breath_peak', ...
        'roiSpk_id', 'nSpkROI', 'ca_spk_data', 'ca_spk_id', 't_spk', ...
        'win', 't_evt', 'win5', ...
        'breath', 't_img', 'dFF', 'Sb', ...
        'Delta_k', 'b_frames', ...
        'dtlast_results', ...
        'dff_evt_seg_all');
    fprintf('Saved master file:\n  %s\n', outMat_master);

    nProcessed = nProcessed + 1;

    catch ME
        nFailed = nFailed + 1;
        failedSessions{end+1} = sessName; %#ok<AGROW>
        fprintf(2, '  ** FAILED: %s\n', ME.message);
        continue;
    end
end

%% 3 ── Summary
fprintf('\n==============================\n');
fprintf('Batch complete: %d processed, %d skipped, %d failed out of %d sessions\n', ...
    nProcessed, nSkipped, nFailed, nSess);
if nFailed > 0
    fprintf('Failed sessions:\n');
    for j = 1:numel(failedSessions)
        fprintf('  - %s\n', failedSessions{j});
    end
end

%% ---- Local functions ----

function S = breath_scan_stat(t, w)
% BREATH_SCAN_STAT  Max number of events in any sliding window of width w.
t = sort(t(:));
n = numel(t);
if n == 0,  S = 0;  return;  end

S    = 1;
left = 1;
for right = 2:n
    while t(right) - t(left) >= w
        left = left + 1;
    end
    S = max(S, right - left + 1);
end
end
