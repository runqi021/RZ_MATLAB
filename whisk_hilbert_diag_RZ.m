function whisk_hilbert_diag_RZ()
% whisk_hilbert_diag_RZ  Hilbert-transform inspection for ONE whisker, ONE
% session (default 5916300 n4). One figure, four stacked x-linked panels:
%   (1) RAW signal           : demeaned DLC angle (NaN-filled, NOT band-passed)
%   (2) BP signal            : band-pass 5-50 Hz, zero-phase filtfilt
%   (3) Hilbert AMPLITUDE    : |hilbert(BP)|  (instantaneous amplitude, raw)
%   (4) Hilbert ENVELOPE     : movmean(|hilbert(BP)|, ENV_SMOOTH)  (smoothed)
%
% (amplitude and envelope are the same quantity; "envelope" is just the
%  smoothed amplitude used for whisk-epoch detection.)

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
animal    = "5916300";
kRun      = 4;
SIDE      = "L";        % "L" or "R"

fpsW      = 400;
BP        = [5 50];     % whisk bandpass (Hz)
ENV_SMOOTH= 0.05;       % s, envelope smoothing (movmean)
ZOOM      = [];         % [t0 t1] s; [] = full trace
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
if strcmpi(SIDE,"L"), ang = La; else, ang = Ra; end
t  = (0:numel(ang)-1)'/fpsW;

raw  = fillmissing(ang-mean(ang,'omitnan'),'linear','EndValues','nearest'); % raw demeaned
bp   = filtfilt(bw,aw, raw);                                                % BP 5-50
amp  = abs(hilbert(bp));                                                    % instantaneous amplitude
env  = movmean(amp, max(1,round(ENV_SMOOTH*fpsW)));                         % smoothed envelope

fprintf('%s n%d (%s): dur=%.0fs | amp median=%.2f, env median=%.2f deg\n', ...
    animal, kRun, SIDE, t(end), median(amp), median(env));

zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

fig = figure('Color','w','Position',[60 60 1240 860]);

ax1 = subplot(4,1,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, raw, '-', 'Color',[0.3 0.3 0.3]);
ylabel(ax1,'raw (deg)');
title(ax1, sprintf('%s n%d  whisker %s  -  Hilbert diagnostic', animal,kRun,SIDE), 'Interpreter','none');

ax2 = subplot(4,1,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, bp, '-', 'Color',[0 0.45 0.75]);
yline(ax2,0,'k-');
ylabel(ax2,sprintf('BP %g-%g Hz (deg)',BP(1),BP(2)));

ax3 = subplot(4,1,3); hold(ax3,'on'); grid(ax3,'on');
plot(ax3, t, bp, '-', 'Color',[0.75 0.85 0.95]);          % BP faint, for context
plot(ax3, t, amp, '-', 'Color',[0.85 0.2 0.2], 'LineWidth',1.2);
plot(ax3, t,-amp, '-', 'Color',[0.85 0.2 0.2], 'LineWidth',1.2);
ylabel(ax3,'Hilbert amplitude (deg)');
title(ax3,'|hilbert(BP)|  (instantaneous amplitude, \pm)', 'Interpreter','tex');

ax4 = subplot(4,1,4); hold(ax4,'on'); grid(ax4,'on');
plot(ax4, t, amp, '-', 'Color',[0.9 0.7 0.7]);            % raw amplitude faint
plot(ax4, t, env, '-', 'Color',[0.6 0.1 0.6], 'LineWidth',1.6);
plot(ax4, t,-env, '-', 'Color',[0.6 0.1 0.6], 'LineWidth',1.6);
ylabel(ax4,'Hilbert envelope (deg)');
xlabel(ax4,'time (s)');
title(ax4, sprintf('movmean(|hilbert(BP)|, %.0f ms)  (smoothed envelope, \\pm)', 1000*ENV_SMOOTH), 'Interpreter','tex');

linkaxes([ax1 ax2 ax3 ax4],'x'); xlim(ax1, zoomwin);

% ===== FIGURE 2: raw overlaid with |+hilbert| + |-hilbert| =====
pos_env = amp;                 % positive Hilbert envelope (+|hilbert(BP)|)
neg_env = -amp;                % negative Hilbert envelope (-|hilbert(BP)|)
hilbSum = abs(pos_env) + abs(neg_env);   % sum of the two -> full envelope band (= 2*amp)

figure('Color','w','Position',[120 200 1240 460]); hold on; grid on;
plot(t, raw,     '-', 'Color',[0.5 0.5 0.5]);
plot(t, hilbSum, '-', 'Color',[0.6 0.1 0.6], 'LineWidth',1.6);
yline(0,'k-');
xlim(zoomwin); xlabel('time (s)'); ylabel('deg');
legend({'raw','|+hilbert| + |-hilbert|'}, 'Location','northeastoutside');
title(sprintf('%s n%d  whisker %s  -  raw (grey) + (|+hilbert| + |-hilbert|) (purple)', ...
      animal,kRun,SIDE), 'Interpreter','none');
end

% ================= helpers =================
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
