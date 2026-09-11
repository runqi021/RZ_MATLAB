function dlc_angle_bilateral(cleanMat)
% dlc_angle_bilateral  Bilateral L-R whisker-ANGLE coordination by cross-
% correlation (xcorr only -- no coherence yet).
%
% Loads the cleaned angle from dlc_clean_whisker, orients each whisker so
% protraction is positive (removes the left/right mirror sign), band-passes to
% the whisking band, z-scores, and cross-correlates within +/-60 ms (< half a
% whisk cycle). Positive lag => R leads L.
%
%   dlc_angle_bilateral('C:\...\..._clean.mat')

if nargin < 1
  cleanMat = ['C:\260613_breathing_thermalNbasler\5916297\archived\' ...
    'cam1_20260613_192624_run001\' ...
    'cam1_20260613_192624_run001DLC_Resnet50_5916297_experimentalJun14shuffle1_snapshot_best-10_clean.mat'];
end

BP       = [4 25];   % whisking band (Hz)
maxLagS  = 0.06;     % xcorr search window (s)

[folder, stem] = fileparts(cleanMat);
Sd = load(cleanMat); out = Sd.out; t = Sd.t; fps = Sd.fps; T = numel(t);

% (1) orient each side protraction-positive (flip the mirror sign)
aL = orient(out.L.angleFilled);
aR = orient(out.R.angleFilled);

% (2) band-pass + z-score, then xcorr within +/-60 ms
[b,a] = butter(2, BP/(fps/2), 'bandpass');
zL = zscore(filtfilt(b,a,aL));
zR = zscore(filtfilt(b,a,aR));
ml = round(maxLagS*fps);
[c,lags] = xcorr(zL, zR, ml, 'coeff');
[rho,im] = max(c); lag = lags(im); lag_ms = lag*1e3/fps;

% ---- figure ----
f = figure('Color','w','Position',[60 60 1200 760],'Visible','off');
subplot(3,1,1);
plot(t, aL, 'Color',[0 .55 0]); hold on; plot(t, aR, 'Color',[.8 .2 0]);
legend({'L','R'}); ylabel('angle (deg)'); xlabel('s'); grid on;
title('oriented whisker angle (protraction +)   L vs R');

win = round(4*fps); v = movvar(zL,win)+movvar(zR,win); [~,jc] = max(v);
i0 = max(1,jc-round(win/2)); i1 = min(T,i0+win);
subplot(3,1,2);
plot(t(i0:i1), aL(i0:i1), 'Color',[0 .55 0]); hold on;
plot(t(i0:i1), aR(i0:i1), 'Color',[.8 .2 0]);
legend({'L','R'}); ylabel('angle (deg)'); xlabel('s'); grid on;
title(sprintf('busiest bout  (%.0f-%.0f s)', t(i0), t(i1)));

subplot(3,1,3);
plot(lags*1e3/fps, c, 'k', 'LineWidth',1); hold on; xline(lag_ms,'r');
xlabel('lag ms   (L vs R)'); ylabel('corr (coeff)'); grid on;
title(sprintf('L-R angle xcorr:  peak %.2f @ %+0.0f ms   (%s)', rho, lag_ms, leadstr(lag)));

sgtitle(stem, 'Interpreter','none');
png = fullfile(folder, [stem '_anglexcorr.png']);
exportgraphics(f, png, 'Resolution',150); close(f);
save(fullfile(folder, [stem '_anglexcorr.mat']), 'aL','aR','zL','zR','c','lags','rho','lag_ms','fps','BP');
fprintf('L-R angle xcorr: peak %.2f @ %+0.0f ms (%s)\nsaved %s\n', rho, lag_ms, leadstr(lag), png);
end

% ===== helpers =====
function y = orient(x)
% flip so the dominant (top-1%) deflections are positive -> protraction +
[~,ix] = sort(abs(x),'descend');
k = ix(1:max(1, round(0.01*numel(x))));
y = x * sign(median(x(k)));
end

function s = leadstr(lag)
if lag > 0, s = 'R leads L'; elseif lag < 0, s = 'L leads R'; else, s = 'in phase'; end
end
