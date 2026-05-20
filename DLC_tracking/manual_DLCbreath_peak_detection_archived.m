%% Manual DLC breath peak detection
%%
% %% ---- Breathing event detection: fine-tune parameters ----
%
% % =========================================================
% % USER-EDITABLE PEAK DETECTION PARAMETERS
% % =========================================================
% minProm_factor = 1.5;      % Prominence = factor × std(breath).
%                             %   Increase → fewer peaks.  Decrease → more peaks.
% minDist_sec    = 0.3;      % Min seconds between consecutive peaks.
% minHeight      = -Inf;     % Min absolute peak height (z-score). -Inf = disabled.
% minWidth_sec   = 0;        % Min peak width in seconds. 0 = disabled.
% snap_win_sec   = 0.2;     % Manual-add snaps to nearest local max within this window (s).
% preFrames      = 0;        % Shift all detected peaks earlier by this many frames.
% % =========================================================
%
% fs    = fps_img;
% T_img = numel(breath);
%
% peak     = breath;
% minProm  = minProm_factor * std(peak);
% minDist  = round(minDist_sec * fs);
% minWidth = round(minWidth_sec * fs);
% snap_win = round(snap_win_sec * fs);
%
% % Build findpeaks argument list dynamically
% fp_args = {'MinPeakProminence', minProm, 'MinPeakDistance', minDist};
% if isfinite(minHeight)
%     fp_args = [fp_args, {'MinPeakHeight', minHeight}];
% end
% if minWidth > 0
%     fp_args = [fp_args, {'MinPeakWidth', minWidth}];
% end
%
% [insp_amp, insp_onsets] = findpeaks(peak, fp_args{:});
%
% % Frame shift
% insp_onsets = insp_onsets - preFrames;
% insp_onsets(insp_onsets < 1) = 1;
%
% fprintf('Auto-detected %d inspiration-onset events.\n', numel(insp_onsets));
%
% % Freeze auto result (used for saving and audit trail)
% insp_onsets_auto = sort(insp_onsets(:));
%
% %% ---- Initial auto-detection QC plot ----
%
% figure('Name', 'Breathing event detection (auto)', 'Position', [200 200 1000 400]);
% plot(t_img, breath, 'k'); hold on;
% plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), 'ro', 'MarkerFaceColor', 'r');
% xline(t_img(insp_onsets_auto), 'r:');
% xlabel('Time (s)'); ylabel('Breath (z)');
% title(sprintf('Auto-detected peaks: n=%d  |  minProm=%.2f  minDist=%.2f s', ...
%     numel(insp_onsets_auto), minProm, minDist_sec));
% grid on; hold off;
%
% %% ---- Interactive peak editor ----
% % Left-click  = add a peak (snaps to nearest local max within snap_win_sec)
% % Right-click = delete nearest existing peak
% % Press Enter = finish editing
%
% insp_onsets_added   = [];            % frame indices manually added
% insp_onsets_deleted = [];            % frame indices removed (from auto set)
% insp_onsets_edit    = insp_onsets_auto;   % working copy
%
% fig_edit = figure('Name', 'Peak Editor | Left-click=Add  Right-click=Delete  Enter=Done', ...
%     'Position', [100 80 1100 480]);
% ax_edit = axes('Parent', fig_edit);
%
% fprintf('\n=== INTERACTIVE PEAK EDITOR ===\n');
% fprintf('  Left-click  : add peak (snaps to local max within %.2f s)\n', snap_win_sec);
% fprintf('  Right-click : delete nearest existing peak\n');
% fprintf('  Press Enter : finish editing\n\n');
%
% while ishandle(fig_edit)
%     % Refresh plot
%     cla(ax_edit);
%     plot(ax_edit, t_img, breath, 'k', 'LineWidth', 1);
%     hold(ax_edit, 'on');
%     if ~isempty(insp_onsets_edit)
%         plot(ax_edit, t_img(insp_onsets_edit), breath(insp_onsets_edit), ...
%             'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
%     end
%     xlabel(ax_edit, 'Time (s)'); ylabel(ax_edit, 'Breath (z)');
%     title(ax_edit, sprintf('Peak Editor  |  n=%d  (+%d / -%d from auto)  |  L=Add  R=Delete  Enter=Done', ...
%         numel(insp_onsets_edit), numel(insp_onsets_added), numel(insp_onsets_deleted)));
%     grid(ax_edit, 'on'); hold(ax_edit, 'off');
%     drawnow;
%
%     % Wait for user click
%     try
%         [x_click, ~, button] = ginput(1);
%     catch
%         break;
%     end
%     if isempty(x_click) || ~ishandle(fig_edit)
%         break;   % Enter pressed → done
%     end
%
%     % Convert clicked x (time in s) → nearest frame index
%     [~, fc] = min(abs(t_img - x_click));
%
%     if button == 1   % Left-click: add peak
%         lo = max(1, fc - snap_win);
%         hi = min(numel(breath), fc + snap_win);
%         [~, rel_idx] = max(breath(lo:hi));
%         new_frame = lo + rel_idx - 1;
%         if ~ismember(new_frame, insp_onsets_edit)
%             insp_onsets_edit = sort([insp_onsets_edit(:); new_frame]);
%             insp_onsets_added(end+1) = new_frame;
%             insp_onsets_deleted(insp_onsets_deleted == new_frame) = [];
%             fprintf('Added   frame %5d  (t = %.2f s)\n', new_frame, t_img(new_frame));
%         else
%             fprintf('Peak at frame %d already exists — skipped.\n', new_frame);
%         end
%
%     elseif button == 3   % Right-click: delete nearest
%         if ~isempty(insp_onsets_edit)
%             [~, ni] = min(abs(t_img(insp_onsets_edit) - x_click));
%             del_frame = insp_onsets_edit(ni);
%             insp_onsets_edit(ni) = [];
%             if ismember(del_frame, insp_onsets_auto)
%                 insp_onsets_deleted(end+1) = del_frame;
%             end
%             insp_onsets_added(insp_onsets_added == del_frame) = [];
%             fprintf('Deleted frame %5d  (t = %.2f s)\n', del_frame, t_img(del_frame));
%         end
%     end
% end
%
% % Finalise
% insp_onsets_final   = sort(insp_onsets_edit(:));
% insp_onsets_added   = sort(insp_onsets_added(:));
% insp_onsets_deleted = sort(insp_onsets_deleted(:));
%
% fprintf('\nEditing complete: %d final peaks  (+%d added, -%d deleted from %d auto)\n', ...
%     numel(insp_onsets_final), numel(insp_onsets_added), ...
%     numel(insp_onsets_deleted), numel(insp_onsets_auto));
%
% % Overwrite insp_onsets so all downstream code is unaffected
% insp_onsets = insp_onsets_final;
% t_insp      = t_img(insp_onsets);
%
% % Recompute binary event train from final peaks
% Sb = zeros(T_img, 1);
% Sb(insp_onsets) = 1;
%
% %% ---- Save: final peak detection figure + .mat ----
%
% fig_qc_final = figure('Name', 'Breathing peak detection (final)', ...
%     'Position', [200 200 1000 400]);
% plot(t_img, breath, 'k', 'LineWidth', 1); hold on;
% if ~isempty(insp_onsets_auto)
%     plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), ...
%         'o', 'Color', [0.65 0.65 0.65], 'MarkerFaceColor', [0.65 0.65 0.65], ...
%         'MarkerSize', 4, 'DisplayName', 'Auto peaks');
% end
% if ~isempty(insp_onsets_added)
%     plot(t_img(insp_onsets_added), breath(insp_onsets_added), ...
%         'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Manually added');
% end
% if ~isempty(insp_onsets_final)
%     plot(t_img(insp_onsets_final), breath(insp_onsets_final), ...
%         'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Final accepted');
%     xline(t_img(insp_onsets_final), 'r:');
% end
% legend('Location', 'best');
% xlabel('Time (s)'); ylabel('Breath (z)');
% title(sprintf('Final: %d peaks  (auto=%d, +%d added, -%d deleted)', ...
%     numel(insp_onsets_final), numel(insp_onsets_auto), ...
%     numel(insp_onsets_added), numel(insp_onsets_deleted)));
% grid on; hold off;
%
% outFig_peaks = fullfile(folderPath, 'breath_peak_detection_final.png');
% exportgraphics(fig_qc_final, outFig_peaks, 'ContentType', 'vector');
% fprintf('Saved final peak detection figure:\n  %s\n', outFig_peaks);
%
% % findpeaks parameter struct
% findpeak_params = struct( ...
%     'minProm_factor',  minProm_factor, ...
%     'minProm',         minProm, ...
%     'minDist_sec',     minDist_sec, ...
%     'minDist_frames',  minDist, ...
%     'minHeight',       minHeight, ...
%     'minWidth_sec',    minWidth_sec, ...
%     'minWidth_frames', minWidth, ...
%     'snap_win_sec',    snap_win_sec, ...
%     'preFrames',       preFrames, ...
%     'fps',             fs ...
% );
%
% outMat_peaks = fullfile(folderPath, 'breath_peak_data.mat');
% save(outMat_peaks, ...
%     'findpeak_params', ...
%     'insp_onsets_auto', ...
%     'insp_onsets_added', ...
%     'insp_onsets_deleted', ...
%     'insp_onsets_final');
% fprintf('Saved peak detection data:\n  %s\n', outMat_peaks);
