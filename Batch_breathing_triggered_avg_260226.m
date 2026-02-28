% Batch
folderPath = 'D:\RUNQI\phys\processed\breathing';
list = dir(folderPath);

ct = 0;
for k = 1:numel(list)
    dlc_csv = []; SAM =[];
    if strcmp(list(k).name,'.') || strcmp(list(k).name,'..'); continue;
    else
        currentPath = fullfile(folderPath, list(k).name);
        sublist = dir(currentPath);
        for i = 1:numel(sublist)
            if endsWith(sublist(i).name, 'snapshot_best-20.csv')
                dlc_csv = fullfile(currentPath, sublist(i).name);
            end
            if endsWith(sublist(i).name, 'cpSAM_output.mat')
                SAM = fullfile(currentPath, sublist(i).name); 
            end
        end

        if isempty(dlc_csv) || isempty(SAM); continue; end
        ct = ct + 1;
        %% Run
        % %% load
        % SAM = "D:\RUNQI\phys\processed\breathing\pFN_roi4_z0_512x512_10x_2000f_00001\pFN_roi4_z0_512x512_10x_2000f_00001_ch1_minusDark_MC_cpSAM_output.mat";
        % dlc_csv = "D:\RUNQI\phys\processed\breathing\pFN_roi4_z0_512x512_10x_2000f_00001\Basler_acA1300-200um__22607866__20260224_022725977DLC_Resnet101_260224_vglut2_breathingFeb25shuffle1_snapshot_best-20.csv";
                
        % Load F
        SAMload = load(SAM);
        F_raw = SAMload.F;
        % toss first second (system steady state）
        F = F_raw;
        F(1:30, :) = [];
        
        % Inverted scanner signal to trigger bassler acq (Falling edge, 2P acq start, bassler acq start, synched). 
        % dFF
        dFFout = helper.dFF_RZ(F);
        dFF = dFFout.dFF;
        
        % ------------------- BASIC PARAMS ---------------------------
        [T, N_roi] = size(dFF);
        
        fps_img  = 30;   % FPS
        fps_breath = fps_img; % stim trigger by flyback
        nDrop = 30;      % for breathing video
        
        t_img = (0:T-1)' / fps_img;
        roi_ids = arrayfun(@(k) sprintf('%02d', k), 1:N_roi, 'UniformOutput', false);
        
        % ------------------- DLC BREATHING INPUT --------------------
        % DLC CSV
        use_dot = 'mean';    % 'dot1', 'dot2', or 'mean'
        
        % DLC CSV has 3 header lines: scorer / bodyparts / coords
        dlc_raw = readmatrix(dlc_csv, 'NumHeaderLines', 3);
        
        data_dlc = dlc_raw;
        data_dlc(1:30, :) = [];
        frame_idx = data_dlc(:,1);
        
        x1 = data_dlc(:,2);  y1 = data_dlc(:,3);  p1 = data_dlc(:,4); s1 = hypot(x1, y1);
        x2 = data_dlc(:,5);  y2 = data_dlc(:,6);  p2 = data_dlc(:,7); s2 = hypot(x2, y2);
        x3 = data_dlc(:,8);  y3 = data_dlc(:,9);  p3 = data_dlc(:,10); s3 = hypot(x3, y3);
        
        % Raw time for breathing video (before alignment)
        t_breath = frame_idx / fps_breath;
        
        % % Choose breathing signal
        % switch lower(use_dot)
        %     case 'dot1'
        %         breath = s1;
        %     case 'dot2'
        %         breath = s2;
        %     otherwise
        %         breath = (s1 + s2)/2;   % average of both dots
        % end
        
        breath = (x1 + x3)/2;
        
        % Optional: use only high-likelihood samples, interpolate others
        good = (p1 > 0.6) & (p2 > 0.6);
        if any(good)
            breath(~good) = interp1(t_breath(good), breath(good), ...
                                         t_breath(~good), 'linear', 'extrap');
        end
        
        % detrend + z-score breathing
        breath = detrend(breath);
        breath = (breath - mean(breath)) / std(breath);
        
        fprintf('Breathing trace: %d samples after dropping first %d frames\n', ...
                numel(breath), nDrop);
        
        % ---- Breathing event detection: inspiration onsets from derivative ----
        
        fs = fps_img;              % imaging / breathing sampling rateke
        T_img = numel(breath);     % length of breathing trace (after resample)
        
        f_band = [0.1 2.5];
        
        % Detect strong negative-going segments:
        %    inspiration onset ≈ large NEGATIVE slope, so we look for peaks in -breath_vel
        
        peak = breath;
        % Choose thresholds
        minProm = 2 * std(peak);      % tweak factor 0.3–1 if too many/few
        minDist = round(fs / f_band(2));         % at least one cycle at the max breathing freq
        
        [insp_amp, insp_onsets] = findpeaks(peak, ...
            'MinPeakProminence', minProm);
        
        %[pk_vel, loc_vel] = findpeaks(peak, ...
        %    'MinPeakProminence', minProm, ...
        %    'MinPeakDistance',   minDist);
        
        % --- binary spike/event trains S(t, ROI) ---
        Sb = zeros(T_img,1);   % 0/1 events
        
        Sb(insp_onsets) = 1;
        
        % frame shift (real in_onset is at peak rise)
        preFrames = 0;
        insp_onsets = insp_onsets - preFrames;
        insp_onsets(insp_onsets < 1) = 1;        % clamp to valid range
        
        % 5) Remove events too close to edges for windowing later
        win = round(1 * fs);                   % e.g. ±2.5 s window later
        valid = (insp_onsets > win) & (insp_onsets < (T_img - win));
        insp_onsets = insp_onsets(valid);
        insp_amp = insp_amp(valid);
        
        t_insp = t_img(insp_onsets);             % event times in seconds
        
        fprintf('Detected %d inspiration-onset events.\n', numel(insp_onsets));
        
        % 6) QC plot
        figure('Name','Breathing event detection','Position',[200 200 800 600]);
        
        plot(t_img, breath, 'k'); hold on;
        plot(t_img(insp_onsets), breath(insp_onsets), 'ro', 'MarkerFaceColor','r');
        xline(t_insp, 'r:');
        xlabel('Time (s)');
        ylabel('d/dt (breath)');
        title('Breathing derivative with detected negative-slope peaks');
        grid on;
        
        % Chronux spectrum
        TW_spec = 4;
        params_spec.Fs     = fps_img;
        params_spec.tapers = [TW_spec, 2*TW_spec - 1];
        params_spec.pad    = 0;
        params_spec.fpass  = [0.1 15];
        params_spec.err    = [2 0.05];
        
        [Sk_raw, fk_raw, Sconfk_raw] = mtspectrumc(breath, params_spec);
        figure('Color','White');
        plot(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
        %hold on;
        %plot(fk_raw, Sconfk_raw(1,:), 'k--','LineWidth', 0.6);
        %hold on;
        %plot(fk_raw, Sconfk_raw(2,:), 'k--','LineWidth', 0.6);
        %hold on;
        xlabel('Frequency (Hz)');
        ylabel('Power');
        xlim([0 10]);  
        grid('on');
        axis square;
        
        % Save breathing spectrum
        outFig = fullfile(currentPath, 'breathing_spectrum_raw.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        fprintf('Saved breathing spectrum:\n  %s\n', outFig);
        
        % %%
        % [Sk_Sb, fk_Sb, Sconfk_Sb] = mtspectrumc(Sb, params_spec);
        % 
        % figure;
        % plot(fk_Sb, Sk_Sb,                'k',  'LineWidth', 1.2);
        % %hold on;
        % %plot(fk_Sb, Sconfk_Sb(1,:), 'k--','LineWidth', 0.6);
        % %hold on;
        % %plot(fk_Sb, Sconfk_Sb(2,:), 'k--','LineWidth', 0.6);
        % %hold on;
        % xlabel('Frequency (Hz)');
        % ylabel('Power');
        % xlim([0 10]); 
        % grid('on');
        
        % Select ROI to keep
        keepIdx = linspace(1, N_roi, N_roi);
        
        %keepIdx = [3 4 5 7 11];
        F_keep = F(:, keepIdx);
        dFF_keep = dFF(:, keepIdx);
        nKeep = numel(keepIdx);
        % ------------------- FIGURE 1: BREATH + STACKED dFF (Top 10 ROIs) --------
        
        % plv_vals(isnan(plv_vals)) = 0;
        % % sort PLV in DESCENDING order (highest first)
        % [~, idx_sort] = sort(plv_vals, 'descend');
        % % ---- choose only top 10 ROIs ----
        % nTop = min(10, N);
        % idx_top = idx_sort(1:nTop);
        % % extract only these ROIs
        % X_top = dFF_keep(:, idx_top);     % [T_img x nTop]
        
        % global amplitude
        amp = max(dFF_keep(:)) - min(dFF_keep(:));
        if amp == 0, amp = 1; end; gap = 0.8 * amp;
        
        % build stacked matrix (highest at top)
        Y = zeros(T_img, nKeep);
        for k = 1:nKeep
            Y(:,k) = dFF_keep(:,k) + gap*(nKeep - k);
        end
        
        % y-ticks and labels
        yt  = 0:gap:gap*(nKeep-1);                    % positions (bottom → top)
        ytl = roi_ids(keepIdx(end:-1:1));            % labels (bottom = lowest PLV of top 10)
        
        figure('Position',[200 200 900 700]);
        tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
        
        % ---- TOP PANEL: breathing + breathing events (Sb) ----
        ax1 = nexttile(1);
        plot(ax1, t_img, -breath, 'k', 'LineWidth', 1); hold(ax1,'on');
        
        % --- breathing spike overlay (Sb) ---
        % assume Sb is [T_img x 1] or [1 x T_img], binary 0/1
        idx_b = find(Sb(:) ~= 0);
        if ~isempty(idx_b)
            t_b = t_img(idx_b);
            y_b = breath(idx_b);
            plot(ax1, t_b, -y_b, 'ro', 'MarkerFaceColor','r', ...
                 'MarkerSize', 3);   % breathing event markers
        end
        xlim(ax1, [t_img(1) t_img(end)])
        xlabel(ax1, 'Time (s)');
        ylabel(ax1, 'Breathing (z)');
        title(ax1, 'Breathing trace (DLC) with detected inspirations');
        grid(ax1,'off');
        hold(ax1,'off');
        
        % BOTTOM PANEL: dFF
        ax2 = nexttile(2, [2 1]); 
        hold(ax2,'on');
        
        for kk = 1:nKeep
            roi_idx = keepIdx(kk);
        
            % main trace
            plot(ax2, t_img, Y(:,kk), 'Color', 'k', 'LineWidth', 0.5);
        
            % % --- ROI spike overlay using S ---
            % idx_evt = find(S(:, roi_idx) ~= 0);   % S is [T_img x N]
            % if ~isempty(idx_evt)
            %     t_ev = t_img(idx_evt);
            %     y_ev = Y(idx_evt, kk)+0.03;           % already offset version
            %     plot(ax2, t_ev, y_ev, '.', ...
            %          'MarkerFaceColor', 'w', ...
            %          'MarkerEdgeColor', 'r', ...
            %          'MarkerSize', 8);
            % end
        end
        
        xlim(ax2, [t_img(1) t_img(end)]);
        ylim(ax2, [-gap, gap*(nKeep-1)+gap]);
        
        yticks(ax2, yt);
        yticklabels(ax2, ytl);
        
        xlabel(ax2, 'Time (s)');
        ylabel(ax2, 'ROI#');
        title(ax2, '\DeltaF/F');
        grid(ax2,'off');
        
        % SCALE BAR
        dFF_scale = 0.2;     % 0.2 ΔF/F vertical bar
        
        % position vertical bar near right
        x0 = 1*(t_img(end));
        y0 = gap*(nKeep-1) + 0.2*gap;
        
        plot(ax2, [x0 x0], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
        text(ax2, x0+0.5, dFF_scale + 0.05*dFF_scale, ...
            sprintf('%.2f \\DeltaF/F', dFF_scale), ...
            'VerticalAlignment','bottom','HorizontalAlignment','center', ...
            'FontSize', 10);
        
        hold(ax2,'off');
        
        % Save breathing spectrum
        outFig = fullfile(currentPath, 'breath_N_dFF.png');
        exportgraphics(gcf, outFig, 'ContentType','vector');
        fprintf('Saved breathing and dFF traces:\n  %s\n', outFig);
        
        % ---- Figure 2. Event-triggered segments per significant ROI (RAW dFF) ----
        
        fs   = fps_img;
        win  = round(0.5 * fs);          % ±2.5 s
        t_evt = (-win:win) / fs;
        
        % Make sure insp_onsets are valid for this window
        T_img = size(dFF_keep,1);
        valid = (insp_onsets > win) & (insp_onsets < (T_img - win));
        insp_onsets = insp_onsets(valid);
        nEv = numel(insp_onsets);
        
        fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);
        
        dff_evt_seg_all = cell(nKeep,1);
        
        for i = 1:nKeep
            roi_trace = dFF_keep(:, i);          % RAW ΔF/F
        
            seg = zeros(nEv, 2*win+1);
            for e = 1:nEv
                c = insp_onsets(e);
                seg(e,:) = roi_trace(c-win : c+win);
            end
        
            dff_evt_seg_all{i} = seg;
        end
        
        % ---------- Per-ROI event-triggered overlays (all trials + mean, RAW dFF) ----------
        
        figure('Name', (sprintf('Event-triggered dFF (n=%d breaths)', size(seg,1))));
        
        for i = 1:nKeep
            seg = dff_evt_seg_all{i};    % nEv x (2*win+1)
            if isempty(seg)
                continue;
            end
        
            % average for this ROI (RAW ΔF/F)
            dff_avg = mean(seg, 1);
        
            hold on;
        
            subplot(ceil(nKeep/10), 10, i);
            % all individual trials
            plot(t_evt, seg', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.7); hold on;
            % overlay mean
            plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);
        
            % reference lines
            xline(0, '--', 'Color', [0.4 0.4 0.4]);  % inspiration onset at 0
            yline(0, '--',  'Color', [0.6 0.6 0.6]);
        
            xlabel('Time from inspiration onset (s)');
            ylabel('\DeltaF/F');
            title(sprintf('ROI#%d', keepIdx(i)));
        
            %grid on;
            hold off;
            axis square;
        
        end
        
        outFig = fullfile(currentPath, 'triggered_avg.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        fprintf('Saved triggered average:\n  %s\n', outFig);
%%
    end
    close all;
end
fprintf(sprintf('%d files processed \n',ct));

