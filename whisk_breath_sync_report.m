function whisk_breath_sync_report(syncMat)
% whisk_breath_sync_report  Assess thermal<->Basler synchronisation for one run.
%
% Both cameras are driven by the SAME 400 Hz WFG trigger, but their absolute
% clocks are unrelated, so frame-level alignment is measured from a shared
% physical signal (breathing). This:
%   - reports per-camera fps / dropped frames (inter-frame-interval histograms)
%   - band-passes each camera's breathing proxy and cross-correlates them to
%     recover the integer-frame offset:   basler_frame = thermal_frame + LAG
%   - writes <run>_sync_report.png and appends the result into the _sync.mat
%
% Input: the *_sync.mat written by whisk_breath_sync_extract.py.
%   whisk_breath_sync_report('C:\...\cam1_..._run001_sync.mat')

if nargin < 1 || isempty(syncMat)
    error('Pass the *_sync.mat path from whisk_breath_sync_extract.py');
end
S = load(syncMat);
[outDir, stem] = fileparts(syncMat);

FS        = 400;            % shared WFG trigger rate (Hz) -- both cameras
BP        = [1 12];        % breathing band (Hz) for the cross-correlation
INTEG_MS  = 0.97648;       % FLIR thermal integration time, ms (frame center = trigger + INTEG/2); Basler exposure 0.25 ms
MAXLAG_S  = 6;             % search +/- this many seconds for the offset

% ---------- per-camera timing ----------
cam_s = double(S.basler_cam_ns - S.basler_cam_ns(1)) / 1e9;   % Basler camera clock
[fps_b, drop_b, isi_b] = timing_stats(cam_s, FS);
ther_t = double(S.thermal_t(:));
[fps_t, drop_t, isi_t] = timing_stats(ther_t, FS);
Nb = numel(cam_s); Nt = numel(ther_t);

% ---------- breathing cross-correlation ----------
bdecim = double(S.basler_bdecim);
fs_b = FS / bdecim;
% put both proxies on a common FS grid on their overlapping support
tb = (0:numel(S.basler_breath)-1)' / fs_b;
xb = resample_to(double(S.basler_breath(:)), tb, FS);
tt = (0:Nt-1)' / FS;
xt = resample_to(double(S.thermal_nostril(:)), tt, FS);
xb = bandpass_simple(xb, FS, BP);
xt = bandpass_simple(xt, FS, BP);
n  = min(numel(xb), numel(xt));
xb = zscore(xb(1:n)); xt = zscore(xt(1:n));

maxlag = round(MAXLAG_S * FS);
[c, lags] = xcorr(xb, xt, maxlag, 'coeff');     % positive lag => basler later than thermal
[~, im] = max(abs(c));
lag = lags(im);                                  % frames: basler_frame ~= thermal_frame + lag
rho = c(im);

% ---------- figure ----------
f = figure('Color','w','Position',[80 80 1180 760],'Visible','off');
tl = tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');

nexttile; histogram(isi_b*1e3,80); xline(1e3/FS,'r'); title('Basler inter-frame interval (ms)'); xlabel('ms'); ylabel('count');
nexttile; histogram(isi_t*1e3,80); xline(1e3/FS,'r'); title('Thermal inter-frame interval (ms)'); xlabel('ms');
nexttile; imagesc(S.thermal_var_img); axis image off; colormap(gca,hot); title('Thermal per-pixel variance (nostril)');

nexttile([1 3]);
shift = lag;                                     % shift thermal by +lag to overlay on basler
tsec = (0:n-1)/FS;
plot(tsec, xb, 'k'); hold on;
xt_sh = circshift(xt, shift);
plot(tsec, xt_sh, 'r');
xlim([min(20,tsec(end)*0.3) min(28,tsec(end)*0.3+8)]);   % ~8 s window
legend({'Basler nostril','Thermal nostril (shifted)'},'Location','northeast');
title(sprintf('Breathing proxies overlaid at LAG = %d frames (%.2f ms)', lag, lag*1e3/FS));
xlabel('s (Basler frame time)');

nexttile([1 2]);
plot(lags*1e3/FS, c, 'b'); hold on; xline(lag*1e3/FS,'r');
xlabel('lag (ms): basler - thermal'); ylabel('corr'); title(sprintf('Cross-correlation (peak \\rho=%.2f)', rho));

ax = nexttile; axis(ax,'off');
txt = {
  sprintf('FILE: %s', stem)
  ''
  sprintf('Basler : %d frames, %.3f Hz, dropped %d', Nb, fps_b, drop_b)
  sprintf('Thermal: %d frames, %.3f Hz, dropped %d', Nt, fps_t, drop_t)
  sprintf('Thermal is %d frames (%.2f s) shorter', Nb-Nt, (Nb-Nt)/FS)
  ''
  sprintf('LAG = %d frames  (basler = thermal + LAG)', lag)
  sprintf('     = %.2f ms   peak corr = %.2f', lag*1e3/FS, rho)
  ''
  sprintf('integ %.4f ms -> frame center +%.3f ms', INTEG_MS, INTEG_MS/2)
  '(constant sub-frame offset, ignore for behavior)'
};
text(ax,0,1,txt,'VerticalAlignment','top','FontName','Consolas','FontSize',9,'Interpreter','none');
title(tl, sprintf('Thermal <-> Basler sync   |   %s', stem),'Interpreter','none');

outPng = fullfile(outDir, [stem '_report.png']);
exportgraphics(f, outPng, 'Resolution', 150);
close(f);

% ---------- persist result ----------
sync_result = struct('lag_frames',lag,'lag_ms',lag*1e3/FS,'corr',rho, ...
    'fps_basler',fps_b,'fps_thermal',fps_t,'drop_basler',drop_b, ...
    'drop_thermal',drop_t,'Nb',Nb,'Nt',Nt,'integ_ms',INTEG_MS,'fs',FS); %#ok<NASGU>
save(syncMat, 'sync_result', '-append');

fprintf('\n=== SYNC RESULT ===\n');
fprintf('Basler  %d frames @ %.3f Hz, dropped %d\n', Nb, fps_b, drop_b);
fprintf('Thermal %d frames @ %.3f Hz, dropped %d\n', Nt, fps_t, drop_t);
fprintf('LAG = %d frames (%.2f ms), corr %.2f  =>  basler = thermal + LAG\n', lag, lag*1e3/FS, rho);
fprintf('saved %s\n', outPng);
end

% ---------------- helpers ----------------
function [fps, ndrop, isi] = timing_stats(t, fs_nom)
isi = diff(t(:));
med = median(isi);
fps = (numel(t)-1) / (t(end)-t(1));
ndrop = round(t(end)*fs_nom) + 1 - numel(t);
end

function y = resample_to(x, t, fs)
tu = (0:floor(t(end)*fs))' / fs;
y = interp1(t, x, tu, 'linear', 'extrap');
end

function y = bandpass_simple(x, fs, band)
% zero-phase Butterworth band-pass (Signal Processing Toolbox)
[b,a] = butter(2, band/(fs/2), 'bandpass');
y = filtfilt(b, a, x);
end
