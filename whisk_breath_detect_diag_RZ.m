function whisk_breath_detect_diag_RZ()
% whisk_breath_detect_diag_RZ  Diagnostic for the inspiration-onset detection
% used by whisk_breath_freq_hist_RZ.m / whisk_breath_raster_RZ.m.
% Shows, for ONE session, the breath signal with detected inspiration onsets
% (troughs of the band-passed inhale-up breath) marked, plus per-breath
% instantaneous frequency. Use it to sanity-check the counting.
%   TOP : raw (grey) + band-passed (black) breath, red v = detected onsets
%   MID : zoom of TOP over ZOOM window
%   BOT : per-breath frequency (1/IBI) at each onset, with 3 & 5 Hz lines

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916296";
kRun      = 3;
BR_BP     = [2 15];     % breath bandpass before onset detection (Hz)
BR_PROM   = 0.5;        % inspiration trough prominence (x std of bandpassed breath)
ZOOM      = [20 30];    % s, zoom window for the middle panel
BASAL_HZ  = 3; SNIFF_HZ = 5;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
assert(isfile(Pn.breath),'no breath .mat: %s', Pn.breath);
Bs=load(Pn.breath); br=Bs.breath(:); fbB=double(Bs.fps); t=(0:numel(br)-1)'/fbB;

% --- method A: band-passed troughs (what the hist/raster use now) ---
[b2,a2] = butter(2, BR_BP/(fbB/2), 'bandpass');
brf = filtfilt(b2, a2, fillmissing(br,'linear'));
[trough,iloc] = findpeaks(-brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fbB/BR_BP(2)));
tInsp = (iloc-1)/fbB; yInsp = -trough;
f = 1./diff(tInsp);
% --- method B: RAW troughs (the original, no bandpass; over-counts?) ---
brd = br - mean(br,'omitnan');
[~,ilocR] = findpeaks(-brd, 'MinPeakProminence', BR_PROM*std(brd,'omitnan'), 'MinPeakDistance', round(0.06*fbB));
tInspR = (ilocR-1)/fbB; fR = 1./diff(tInspR);

fprintf('%s n%d: fps=%.1f, dur=%.0fs\n', animal, kRun, fbB, t(end));
fprintf('  BP  : %d onsets, %.2f/s, median %.2f Hz\n', numel(tInsp),  numel(tInsp)/t(end),  median(f,'omitnan'));
fprintf('  RAW : %d onsets, %.2f/s, median %.2f Hz\n', numel(tInspR), numel(tInspR)/t(end), median(fR,'omitnan'));

figure('Color','w','Position',[60 60 1200 760]);

ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, br-mean(br,'omitnan'), '-', 'Color',[0.7 0.7 0.7]);
plot(ax1, t, brf, 'k-');
plot(ax1, tInsp, yInsp, 'rv', 'MarkerFaceColor','r', 'MarkerSize',4);
xlim(ax1,[t(1) t(end)]); ylabel(ax1,'breath (inhale up)');
title(ax1, sprintf('%s n%d  —  %d inspiration onsets (red), %.2f breaths/s', ...
    animal, kRun, numel(tInsp), numel(tInsp)/t(end)));
legend(ax1, {'raw (demeaned)', sprintf('BP %g-%g Hz',BR_BP(1),BR_BP(2)), 'detected onset'}, ...
    'Location','northeastoutside');

ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, brd, '-', 'Color',[0.7 0.7 0.7]);          % raw demeaned
plot(ax2, t, brf, 'k-');                                  % bandpassed
plot(ax2, tInsp,  interp1(t,brf,tInsp),  'rv', 'MarkerFaceColor','r', 'MarkerSize',7);    % BP onsets
plot(ax2, tInspR, interp1(t,brd,tInspR), 'bo', 'MarkerSize',7, 'LineWidth',1);            % RAW onsets
xlim(ax2, ZOOM); ylabel(ax2,'breath'); xlabel(ax2,'time (s)');
legend(ax2, {'raw','BP', sprintf('BP onset (%.2f/s)',numel(tInsp)/t(end)), ...
    sprintf('RAW onset (%.2f/s)',numel(tInspR)/t(end))}, 'Location','northeastoutside');
title(ax2, sprintf('zoom %g-%g s — count real breaths vs red(BP)/blue(RAW) marks', ZOOM(1), ZOOM(2)));

ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
plot(ax3, tInsp(2:end), f, '.-', 'Color',[0.1 0.4 0.8], 'MarkerSize',8);
yline(ax3, BASAL_HZ,'k--'); yline(ax3, SNIFF_HZ,'r--');
xlim(ax3,[t(1) t(end)]); ylim(ax3,[0 max(15,prctile(f,99))]);
xlabel(ax3,'time (s)'); ylabel(ax3,'breath freq (Hz)');
title(ax3, sprintf('per-breath frequency (1/IBI)  —  median %.2f Hz', median(f,'omitnan')));
linkaxes([ax1 ax3],'x');
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
