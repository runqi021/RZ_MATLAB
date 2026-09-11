function breath_hilbert_peak_diag_RZ()
% breath_hilbert_peak_diag_RZ  Overlay Hilbert "peaks" on the DETRENDED breath
% trace to see how well they land on the real maxima.
% Hilbert peak = where angle(hilbert(x)) crosses 0 upward (= analytic-signal max).
% No bandpass -- just detrend, so you can judge whether raw Hilbert is usable.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916296";
kRun      = 3;
BP        = [2 15];     % breath bandpass before Hilbert (Hz) -- TUNE THIS
ZOOM      = [20 30];    % s, zoom window
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
assert(isfile(Pn.breath),'no breath .mat: %s', Pn.breath);
Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); t=(0:numel(br)-1)'/fb;

x    = detrend(fillmissing(br,'linear'));         % DETRENDED breath (no bandpass)
phi  = angle(hilbert(x));  pk  = find(phi(1:end-1)<0 & phi(2:end)>=0);   % raw Hilbert peaks
[b2,a2] = butter(2, BP/(fb/2), 'bandpass');
xb   = filtfilt(b2,a2,x);                          % BANDPASSED breath
phib = angle(hilbert(xb)); pkb = find(phib(1:end-1)<0 & phib(2:end)>=0); % BP Hilbert peaks

fprintf('%s n%d: fps=%.1f, dur=%.0fs\n', animal, kRun, fb, t(end));
fprintf('  RAW Hilbert : %d peaks (%.2f/s)\n', numel(pk),  numel(pk)/t(end));
fprintf('  BP  Hilbert : %d peaks (%.2f/s)  [%g-%g Hz]\n', numel(pkb), numel(pkb)/t(end), BP(1),BP(2));

figure('Color','w','Position',[60 80 1200 620]);
ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
plot(ax1, t, x, '-', 'Color',[0.7 0.7 0.7]);
plot(ax1, t, xb, 'k-');
plot(ax1, t(pk),  x(pk),   'r.', 'MarkerSize',6);
plot(ax1, t(pkb), xb(pkb), 'gv', 'MarkerFaceColor','g', 'MarkerSize',4);
xlim(ax1,[t(1) t(end)]); ylabel(ax1,'breath');
legend(ax1, {'detrended','BP', sprintf('RAW peaks (%.1f/s)',numel(pk)/t(end)), ...
    sprintf('BP peaks (%.1f/s)',numel(pkb)/t(end))}, 'Location','northeastoutside');
title(ax1, sprintf('%s n%d  —  raw Hilbert (red) vs band-passed Hilbert (green)', animal, kRun));

ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t, x, '-', 'Color',[0.7 0.7 0.7]);
plot(ax2, t, xb, 'k-');
plot(ax2, t(pk),  x(pk),   'r.', 'MarkerSize',10);
plot(ax2, t(pkb), xb(pkb), 'gv', 'MarkerFaceColor','g', 'MarkerSize',8);
xlim(ax2, ZOOM); xlabel(ax2,'time (s)'); ylabel(ax2,'breath');
title(ax2, sprintf('zoom %g-%g s  (BP %g-%g Hz — green should land one-per-breath)', ZOOM(1),ZOOM(2),BP(1),BP(2)));
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
