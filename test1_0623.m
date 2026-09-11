% ================= ZERO-PHASE FILTER La =================
% Assumes:
%   La   = whisker angle vector
%   fpsW = sampling rate, e.g. 400 Hz

fs = fpsW;

x = La(:);

% filtfilt cannot handle NaNs, so interpolate first
x = fillmissing(x, 'linear', 'EndValues', 'nearest');

% Optional: remove huge DC offset for cleaner plotting
% Comment this out if you want absolute angle preserved in LP signal.
x0 = x - mean(x, 'omitnan');
% figure; plot(t, x); hold on; plot(t,x0);

% ---------------- Low-pass < 2 Hz ----------------
ord = 4;
fc_lp = 5;   % Hz

[b_lp, a_lp] = butter(ord, fc_lp/(fs/2), 'low');
La_lp5 = filtfilt(b_lp, a_lp, x0);

% [b_lp, a_lp] = butter(ord, 1/(fs/2), 'low');
% La_lp1 = filtfilt(b_lp, a_lp, x0);
% [b_lp, a_lp] = butter(ord, 3/(fs/2), 'low');
% La_lp3 = filtfilt(b_lp, a_lp, x0);
% figure; plot(t, x0-La_lp3, Color='k'); hold on; plot(t, x0-La_lp1); hold on; plot(t, La_bp5_50, 'r'); hold on;

% [b_lp, a_lp] = butter(ord, 1/(fs/2), 'low');
% La_lp1 = filtfilt(b_lp, a_lp, x0);
% [b_lp, a_lp] = butter(ord, 5/(fs/2), 'low');
% La_lp5 = filtfilt(b_lp, a_lp, x0);
% figure; plot(t, La_lp5, Color='k'); hold on; plot(t, La_lp1); hold on; plot(t, x0, 'g', 'LineWidth',0.1); hold on;
%%
% ---------------- Band-pass 5-50 Hz ----------------
bp = [5 50]; % Hz

[b_bp, a_bp] = butter(ord, bp/(fs/2), 'bandpass');
La_bp5_50 = filtfilt(b_bp, a_bp, x0);

% ---------------- Quick plot ----------------
t = (0:numel(La)-1)' / fs;

figure('Color','w'); hold on;
plot(t, x0, 'Color', [0.7 0.7 0.7]);
plot(t, La_lp5, 'k', 'LineWidth', 2);
plot(t, La_bp5_50, 'r');
xlabel('time (s)');
ylabel('La angle, demeaned deg');
legend({'raw La', 'LP < 2 Hz', 'BP 5-50 Hz'});
box off;
%%


%% ================= WHISK-EPOCH DETECTION =================

% ---------------- envelope detection ----------------
thr = 5;
%THR_FRAC   = 0.2;  % threshold = THR_FRAC * 95th percentile envelope
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;     % s
MERGE_GAP  = 0.2;   % s

env_raw = abs(hilbert(La_bp5_50));
env = movmean(env_raw, max(1, round(ENV_SMOOTH * fs)));

%thr = THR_FRAC * prctile(env, 95);

is_whisk = env > thr;

d = diff([false; is_whisk; false]);
s = find(d == 1);
e = find(d == -1) - 1;

ep = [t(s), t(e)];

% merge nearby epochs
if ~isempty(ep)
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= MERGE_GAP
            ep2(end,2) = ep(i,2);
        else
            ep2(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end
    ep = ep2;
end

% remove short epochs
if ~isempty(ep)
    ep = ep(ep(:,2) - ep(:,1) >= MIN_DUR, :);
end

fprintf('%d whisk epochs detected, %.1f%% of time\n', ...
    size(ep,1), 100 * sum(ep(:,2)-ep(:,1)) / t(end));

% ---------------- plot ----------------
figure('Color','w'); hold on;

yl = [min([-env; La_bp5_50]), max([env; La_bp5_50])];

for q = 1:size(ep,1)
    patch(ep(q,[1 2 2 1]), yl([1 1 2 2]), ...
        [0.3 0.75 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end

plot(t, La_bp5_50, 'k');
plot(t, env, 'r', 'LineWidth', 1.5);
plot(t, -env, 'r', 'LineWidth', 1.5);
yline(thr, 'r--');
yline(-thr, 'r--');

xlabel('time (s)');
ylabel('La BP angle / envelope');
title('La whisk-epoch detection');
legend({'whisk epoch', 'La BP 5-60 Hz', '+/- envelope', '', 'threshold'}, ...
    'Location','best');
box off;







