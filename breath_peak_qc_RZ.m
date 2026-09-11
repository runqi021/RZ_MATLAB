% breath_peak_qc_RZ.m  (run as a script)
% QC the inspiration-onset detection on ALL regenerated *_breath.mat at once.
% Inspiration onset = trough of the band-passed inhale-up breath (the event used
% by whisk_breath_*). One montage panel per session: a WINSEC window of the
% band-passed breath with detected onsets (red v) + per-session median rate.
% Console prints a table with QC flags (too-short IBIs = double counts, long gaps
% = misses, out-of-range median rate).

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
BR_BP    = [2 10];     % breathing band (Hz) before onset detection
BR_PROM  = 0.5;        % inspiration trough prominence (x std of bandpassed breath)
WINSEC   = 12;         % s, display window per panel
T0_FRAC  = 0.40;       % start the display window this far into each recording
NCOL     = 3;          % montage columns
SHORT_HZ = 12;         % IBI faster than this -> likely double-detection (QC flag)
LONG_S   = 1.0;        % IBI longer than this -> likely a miss (QC flag)
% ======================================================================

d = dir(fullfile(char(dataRoot), '**', '*_breath.mat'));
assert(~isempty(d), 'no *_breath.mat under %s', dataRoot);
[~,o] = sort({d.name}); d = d(o);
N = numel(d);
fprintf('QC: %d breath files\n\n', N);
fprintf('%-14s %5s %6s %6s %7s %7s %7s\n','session','dur','nOns','medHz','IQRHz','%fast','%long');

figure('Color','w','Position',[40 40 1500 900]);
tl = tiledlayout(ceil(N/NCOL), NCOL, 'TileSpacing','compact','Padding','compact');
for i = 1:N
    B = load(fullfile(d(i).folder, d(i).name));
    br = B.breath(:); fps = double(B.fps); t = (0:numel(br)-1)'/fps;
    if isfield(B,'animal'), tagn = sprintf('%s n%g', char(string(B.animal)), B.run);
    else, tk = regexp(d(i).name,'(\d+)_.*n(\d+)','tokens','once'); tagn = strjoin(tk,' n'); end
    [b2,a2] = butter(2, BR_BP/(fps/2), 'bandpass');
    brf = filtfilt(b2, a2, fillmissing(br,'linear'));
    [~,iloc] = findpeaks(-brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fps/BR_BP(2)));
    tIns = (iloc-1)/fps; ibi = diff(tIns); f = 1./ibi;
    medHz = median(f,'omitnan'); iqrHz = iqr(f(~isnan(f)));
    pFast = 100*mean(f > SHORT_HZ);            % suspicious double-counts
    pLong = 100*mean(ibi > LONG_S);            % suspicious misses
    flag = '';
    if medHz<1 || medHz>SHORT_HZ, flag=[flag ' RATE?']; end
    if pFast>5, flag=[flag ' DBL?']; end
    if pLong>5, flag=[flag ' MISS?']; end %#ok<AGROW>
    fprintf('%-14s %5.0f %6d %6.2f %7.2f %6.0f%% %6.0f%%%s\n', tagn, t(end), numel(iloc), medHz, iqrHz, pFast, pLong, flag);

    % --- montage panel: WINSEC window ---
    nexttile; hold on; grid on;
    t0 = max(0, T0_FRAC*t(end)); t1 = min(t(end), t0+WINSEC);
    w = t>=t0 & t<=t1; mk = tIns>=t0 & tIns<=t1;
    plot(t(w), brf(w), 'k-', 'LineWidth',0.7);
    plot(tIns(mk), interp1(t,brf,tIns(mk)), 'rv', 'MarkerFaceColor','r', 'MarkerSize',5);
    xlim([t0 t1]); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    ttl = sprintf('%s  %.2f Hz (n%d)', tagn, medHz, numel(iloc));
    if ~isempty(flag), title(ttl,'Color',[0.85 0 0]); else, title(ttl); end
end
title(tl, sprintf('Inspiration-onset QC  —  BP %g-%g Hz, prom %.1f, %g s window @ %.0f%%', ...
    BR_BP(1),BR_BP(2),BR_PROM,WINSEC,100*T0_FRAC));
fprintf('\nflags: RATE? median outside 1-%g Hz | DBL? >5%% IBI faster than %g Hz | MISS? >5%% IBI longer than %gs\n', ...
    SHORT_HZ, SHORT_HZ, LONG_S);
