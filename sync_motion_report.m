function sync_motion_report(motMat)
% sync_motion_report  Lock the thermal<->Basler integer frame offset using shared
% MOTION energy (abrupt movements appear on the same WFG pulse in both cameras),
% which is robust to the cameras measuring breathing differently.
%
% Input: *_syncmot.mat from sync_motion_extract.py
%   sync_motion_report('C:\...\cam1_..._run001_syncmot.mat')

S = load(motMat);
[outDir, stem] = fileparts(motMat);
FS = 400;                  % shared WFG trigger (Hz)
MAXLAG_S = 6;              % search +/- this for the offset
HP = 0.5;                  % high-pass (Hz) to emphasise transients over slow drift

bm = double(S.basler_motion(:));  bm(1) = bm(2);    % first diff is 0
tm = double(S.thermal_motion(:)); tm(1) = tm(2);
hp = @(x) filtfilt(butter(2, HP/(FS/2), 'high'), 1, x - mean(x));
zb = zscore(hp(bm));  zt = zscore(hp(tm));

maxlag = round(MAXLAG_S*FS);
n = min(numel(zb), numel(zt));
zb = zb(1:n); zt = zt(1:n);
[c, lags] = xcorr(zb, zt, maxlag, 'coeff');             % +lag => basler later than thermal
[rho, im] = max(c);  lag = lags(im);                    % basler_frame = thermal_frame + lag

% breathing proxies (secondary, nostril) for the overlay sanity-check
bp = @(x) filtfilt(butter(2,[1 12]/(FS/2),'bandpass'), 1, x - mean(x));
bn = zscore(bp(double(S.basler_nostril(:))));
tn = zscore(bp(double(S.thermal_nostril(:))));

% timing stats
cam_s = double(S.basler_cam_ns - S.basler_cam_ns(1))/1e9;
fps_b = (numel(cam_s)-1)/cam_s(end);
ther_t = double(S.thermal_t(:));
fps_t = (numel(ther_t)-1)/ther_t(end);
Nb = double(S.Nb); Nt = double(S.Nt);

f = figure('Color','w','Position',[60 60 1200 780],'Visible','off');
tl = tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');

nexttile([1 2]);
tsec=(0:n-1)/FS; plot(tsec,zb,'k'); hold on; plot(tsec,circshift(zt,lag),'r');
xlim([0 min(tsec(end),130)]); legend({'Basler motion','Thermal motion (shifted)'});
title(sprintf('Whole-recording motion energy at LAG=%d frames (%.2f ms), \\rho=%.2f',lag,lag*1e3/FS,rho));
xlabel('s');

ax=nexttile; axis(ax,'off');
txt = {sprintf('FILE: %s',stem); '';
  sprintf('Basler : %d fr, %.3f Hz',Nb,fps_b);
  sprintf('Thermal: %d fr, %.3f Hz',Nt,fps_t);
  sprintf('Thermal %d fr (%.2f s) shorter',Nb-Nt,(Nb-Nt)/FS); '';
  sprintf('LAG = %d frames',lag);
  sprintf('    = %.2f ms,  rho = %.2f',lag*1e3/FS,rho); '';
  'basler_frame = thermal_frame + LAG'};
text(ax,0,1,txt,'VerticalAlignment','top','FontName','Consolas','FontSize',10,'Interpreter','none');

nexttile; plot(lags*1e3/FS,c,'b'); hold on; xline(lag*1e3/FS,'r');
xlabel('lag (ms): basler - thermal'); ylabel('corr'); title('Motion cross-correlation'); xlim([-200 200]);

nexttile([1 2]);
plot(tsec, bn(1:n), 'k'); hold on;
plot(tsec, circshift(tn(1:n), lag), 'r');
xlim([60 68]); legend({'Basler nostril','Thermal nostril (shifted)'});
title('Breathing proxies at the motion-derived LAG (8 s)'); xlabel('s');

nexttile; histogram(diff(cam_s)*1e3,80); xline(1e3/FS,'r'); title('Basler ISI (ms)'); xlabel('ms');

title(tl, sprintf('Thermal <-> Basler MOTION sync   |   %s', stem),'Interpreter','none');
outPng = fullfile(outDir,[stem '_motion_report.png']);
exportgraphics(f,outPng,'Resolution',150); close(f);

sync_result = struct('lag_frames',lag,'lag_ms',lag*1e3/FS,'corr',rho, ...
    'fps_basler',fps_b,'fps_thermal',fps_t,'Nb',Nb,'Nt',Nt,'fs',FS,'method','motion'); %#ok<NASGU>
save(motMat,'sync_result','-append');

fprintf('\n=== MOTION SYNC ===\n');
fprintf('LAG = %d frames (%.2f ms), corr %.2f  =>  basler = thermal + LAG\n',lag,lag*1e3/FS,rho);
fprintf('saved %s\n', outPng);
end
