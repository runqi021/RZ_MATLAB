function dlc_clean_whisker(csvPath)
% dlc_clean_whisker  Physics-based cleaning of DLC whisker tracking.
%
% Flags bad-position frames by KINEMATICS (Hampel on the tip x,y) and GEOMETRY
% (base->tip distance), NOT by likelihood -- because likelihood is the model's
% uncertainty, not a per-frame correctness flag, and it drops across a whole
% hard epoch even on frames whose cross is actually placed correctly. So:
%   * isolated bad-position frames (short runs) are interpolated,
%   * long bad runs are left as gaps and marked INVALID (never interpolated),
%   * sustained low-likelihood stretches are flagged as a separate lowConf mask
%     for review / for pulling into the training set -- they are NOT auto-removed
%     if their position is geometrically/kinematically plausible.
%
% Outputs (saved to <stem>_clean.mat as struct `out` with fields L and R):
%   tip      [T x 2]  cleaned tip xy (isolated bad frames interpolated; long NaN)
%   sweep    [T x 1]  tip projected on its sweep axis (axis fit on trusted frames)
%   valid    [T x 1]  position is trustworthy  -> use for xcorr / coherence
%   lowConf  [T x 1]  sustained low-likelihood epoch -> review / add to training
%   p        [T x 1]  raw tip likelihood
%
%   dlc_clean_whisker('C:\...\xxxDLC_Resnet50_...csv')

if nargin < 1
  csvPath = ['C:\260613_breathing_thermalNbasler\5916297\archived\' ...
    'cam1_20260613_192624_run001\' ...
    'cam1_20260613_192624_run001DLC_Resnet50_5916297_experimentalJun14shuffle1_snapshot_best-10.csv'];
end

% ---- user params ----
kHampel     = 5;     % MAD multiplier for kinematic/geometry outlier flagging
winSec      = 0.05;  % half-window for the kinematic local median (s)
pLow        = 0.5;   % likelihood considered "uncertain"
lowEpochSec = 0.25;  % sustained low-likelihood run length flagged as lowConf (s)
maxFillSec  = 0.10;  % bad runs shorter than this -> interpolate; longer -> invalid

[folder, stem] = fileparts(csvPath);

% ---- fps from the run's camera timestamps, else 400 ----
fps = 400; ts = fullfile(folder, 'timestamps.csv');
if isfile(ts), A = readmatrix(ts); fps = (size(A,1)-1)/((A(end,2)-A(1,2))/1e9); end

M = readmatrix(csvPath, 'NumHeaderLines', 3);
T = size(M,1); t = (0:T-1)'/fps;
win     = max(1, round(winSec*fps));
lowRun  = max(1, round(lowEpochSec*fps));
maxFill = max(1, round(maxFillSec*fps));

% column layout: frame, vL0(x y l), vL1, vR0, vR1
sides = struct('name',{'L','R'}, 'base',{[2 3],[8 9]}, 'tip',{[5 6],[11 12]}, 'pcol',{7,13});
out = struct();

for s = 1:2
  base = M(:, sides(s).base);
  tip  = M(:, sides(s).tip);
  p    = M(:, sides(s).pcol);

  % (1) GEOMETRY: base->tip distance is a near-rigid segment; a tip that jumps
  %     onto the wrong thing breaks it, regardless of likelihood.
  d = vecnorm(tip - base, 2, 2);
  geomBad = abs(d - median(d,'omitnan')) > kHampel * 1.4826 * mad(d,1);

  % (2) KINEMATICS: neither dot can teleport -> Hampel outliers on each coord,
  %     for the TIP and the BASE (a bad base frame swings the angle too).
  [~, kxt] = hampel(tip(:,1),  win, kHampel);  [~, kyt] = hampel(tip(:,2),  win, kHampel);
  [~, kxb] = hampel(base(:,1), win, kHampel);  [~, kyb] = hampel(base(:,2), win, kHampel);
  bad = geomBad | kxt | kyt | kxb | kyb;  % physically-wrong-position frames

  % long bad runs -> invalid (never interpolated); short runs -> interpolate
  [rs, re] = runs(bad);
  invalidLong = false(T,1);
  for r = 1:numel(rs)
    if (re(r)-rs(r)+1) > maxFill, invalidLong(rs(r):re(r)) = true; end
  end
  [tipC, tipF]   = clean_dot(tip,  bad, invalidLong);   % gapped + fully-filled
  [baseC, baseF] = clean_dot(base, bad, invalidLong);
  valid = ~invalidLong & all(isfinite(tipC),2) & all(isfinite(baseC),2);

  % (3) sustained low-likelihood epochs (independent of position) -> review mask
  lowConf = false(T,1);
  [ls, le] = runs(p < pLow);
  for r = 1:numel(ls)
    if (le(r)-ls(r)+1) >= lowRun, lowConf(ls(r):le(r)) = true; end
  end

  % WHISKER ANGLE = base->tip direction, measured RELATIVE to its own mean
  % direction (over valid frames). Referencing to the mean centers the range at
  % 0, so there is no +/-180 wrap to unwrap -- unwrap on raw atan2 jumps a full
  % turn whenever the whisker sits near the +/-180 boundary (the left one does).
  % Now safe because the bad frames that used to flip it have been removed.
  v = tipF - baseF; m = mean(v(valid,:),1); m = m/(norm(m)+eps);
  ang = rad2deg(atan2(m(1)*v(:,2) - m(2)*v(:,1), m(1)*v(:,1) + m(2)*v(:,2)));
  ang = ang - median(ang(valid));
  angGap = ang; angGap(~valid) = NaN;

  % sweep (tip projected on its PC axis) kept for reference
  vf = valid & ~lowConf;
  if nnz(vf) < 10, vf = valid; end
  mu = mean(tipC(vf,:), 1);
  [V,~] = eig(cov(tipC(vf,:) - mu));
  sw = (tipC - mu) * V(:,end);  swF = (tipF - mu) * V(:,end);

  out.(sides(s).name) = struct('tip',tipC,'tipFilled',tipF,'base',baseC, ...
    'baseFilled',baseF,'angle',angGap,'angleFilled',ang,'sweep',sw,'sweepFilled',swF, ...
    'valid',valid,'lowConf',lowConf,'p',p,'badFrac',mean(bad), ...
    'lowConfFrac',mean(lowConf),'validFrac',mean(valid));
  fprintf(['%s tip: %.1f%% bad-position (geom|kin), %.1f%% in low-conf epochs ' ...
    '(>%.0f ms), %.1f%% usable\n'], sides(s).name, 100*mean(bad), ...
    100*mean(lowConf), 1000*lowEpochSec, 100*mean(valid));
end

% ---- diagnostic figure ----
f = figure('Color','w','Position',[60 60 1200 720],'Visible','off');
for s = 1:2
  S  = out.(sides(s).name);
  ax = subplot(2,1,s);
  yyaxis(ax,'left'); hold(ax,'on');
  plot(ax, t, S.angleFilled, 'k-', 'LineWidth',0.5);   % cleaned whisker angle (deg)
  fin = isfinite(S.angleFilled);
  yl = prctile(S.angleFilled(fin), [0.5 99.5]);
  pad = 0.10*max(diff(yl),eps); yl = [yl(1)-pad yl(2)+pad];
  hp = gobjects(0);
  hp = shade(ax, t, ~S.valid,  yl, [1 .85 .85], hp);   % invalid (long bad) - pink
  hp = shade(ax, t,  S.lowConf,yl, [.85 .9 1],  hp);   % low-conf epoch  - blue
  uistack(hp, 'bottom');
  ylim(ax, yl); ylabel(ax,'angle (deg)');
  yyaxis(ax,'right'); plot(ax, t, S.p, 'Color',[0 .5 0]); ylim(ax,[0 1]);
  ylabel(ax,'likelihood');
  title(ax, sprintf(['%s whisker angle  (black=clean, pink=invalid, ' ...
    'blue=low-conf epoch, green=likelihood)'], sides(s).name));
  xlabel(ax,'s'); grid(ax,'on');
end
sgtitle(stem, 'Interpreter','none');
png = fullfile(folder, [stem '_clean.png']);
exportgraphics(f, png, 'Resolution',150); close(f);
save(fullfile(folder, [stem '_clean.mat']), 'out', 't', 'fps');
fprintf('saved %s\n', png);
end

% ===== helpers =====
function [s, e] = runs(mask)        % start/end indices of true-runs
m = [false; mask(:); false]; dm = diff(m); s = find(dm==1); e = find(dm==-1)-1;
end

function y = fill_interp(y, fillable)
idx = (1:numel(y))'; good = isfinite(y) & ~fillable;
if nnz(good) >= 2
  y(fillable) = interp1(idx(good), y(good), idx(fillable), 'linear', NaN);
end
end

function [cl, clF] = clean_dot(xy, bad, invalidLong)
% gapped (long bad runs left NaN) + fully linearly-interpolated continuous copy
cl = xy; cl(bad,:) = NaN;
fillable = isnan(cl(:,1)) & ~invalidLong;
cl(:,1) = fill_interp(cl(:,1), fillable);
cl(:,2) = fill_interp(cl(:,2), fillable);
clF = cl;
clF(:,1) = fill_interp(clF(:,1), isnan(clF(:,1)));
clF(:,2) = fill_interp(clF(:,2), isnan(clF(:,2)));
clF = fillmissing(clF, 'nearest');
end

function sw = local_sweep(xy)       % raw sweep on finite frames (for display)
g = all(isfinite(xy),2); mu = mean(xy(g,:),1);
[V,~] = eig(cov(xy(g,:) - mu)); sw = (xy - mu) * V(:,end);
end

function hp = shade(ax, t, mask, yl, col, hp)
[s, e] = runs(mask);
for r = 1:numel(s)
  hp(end+1) = patch(ax, [t(s(r)) t(e(r)) t(e(r)) t(s(r))], ...
    [yl(1) yl(1) yl(2) yl(2)], col, 'EdgeColor','none'); %#ok<AGROW>
end
end
