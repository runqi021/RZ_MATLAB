function whisk_epoch_single_overlay_RZ()
% whisk_epoch_single_overlay_RZ  Detailed overlay for ONE auto-detected whisking
% epoch of ONE session.
%   TOP : band-passed signals, breathing (left y-axis) + whisker L & R (right
%         y-axis, both whiskers share that axis).
%   BOT : Hilbert instantaneous phase overlay of the three (breath, L, R).
% Epoch index uses the same auto-detection as whisk_epoch_analyze_auto_RZ.m.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916296";
kRun      = 3;          % nK
EP        = 7;          % which detected whisk epoch
THR_FRAC  = 0.30;       % whisk-epoch detection (must match the analysis script)
MIN_DUR   = 1.0;
MERGE_GAP = 0.20;
BP        = [6 30];     % whisk bandpass (Hz)
BR_BP     = [2 15];     % breathing bandpass (Hz)
PAD       = 0.3;        % s context shown either side of the epoch
PSD_FPASS = [0 40];     % Hz, PSD display band
PSD_SEC   = 1.5;        % s, FFT segment length for averaged periodogram (50% overlap)
fpsW      = 400;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot,'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

% ---- whisker ----
wcsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun));
M = dlc_gate_interp(wcsv, 0.6);    % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
t  = (0:numel(La)-1)'/fpsW;
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
Laf = fillmissing(La,'linear'); Raf = fillmissing(Ra,'linear');   % raw angle (PSD uses detrended raw)
env = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;
phiL = angle(hilbert(xL)); phiR = angle(hilbert(xR));

% ---- breathing (resampled to whisk grid, BP) ----
brw = nan(size(t)); phiB = nan(size(t));
Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kRun)), dataRoot);
assert(isfile(Pn.breath), 'no breath .mat: %s', Pn.breath);
Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
[b2,a2]=butter(2,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
brw = interp1(tBr, brf, t, 'linear', NaN);
braw = interp1(tBr, fillmissing(br,'linear'), t, 'linear', NaN);   % raw breath resampled (PSD)
ok = ~isnan(brw); phiB(ok) = angle(hilbert(brw(ok)));

% ---- epoch window ----
ep = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
assert(EP>=1 && EP<=size(ep,1), 'epoch %d out of range (1..%d)', EP, size(ep,1));
w  = t>=ep(EP,1)-PAD & t<=ep(EP,2)+PAD;     % epoch + context (time series)
we = t>=ep(EP,1)    & t<=ep(EP,2);          % epoch only (scatter)
xl = [ep(EP,1)-PAD ep(EP,2)+PAD];
colL=[0 0.5 0]; colR=[0 0.4 0.85]; colB=[0 0 0];

figure('Color','w','Position',[60 60 1150 860]);

% ---- TOP: band-passed signals, breath (left) + whiskers (right) ----
ax1 = subplot(3,2,[1 2]); hold(ax1,'on'); grid(ax1,'on');
yyaxis(ax1,'left');
hB = plot(ax1, t(w), brw(w), '-', 'Color', colB, 'LineWidth', 1.6);
ylabel(ax1, sprintf('breathing (BP %g-%g Hz)', BR_BP(1), BR_BP(2)));
ax1.YColor = colB;
yyaxis(ax1,'right');
hL = plot(ax1, t(w), xL(w), '-', 'Color', colL, 'LineWidth', 1.2);
hR = plot(ax1, t(w), xR(w), '-', 'Color', colR, 'LineWidth', 1.2);
ylabel(ax1, sprintf('whisker angle (deg, BP %g-%g Hz)', BP(1), BP(2)));
ax1.YColor = [0 0 0];
xlim(ax1, xl);
xline(ax1, ep(EP,1), 'r--'); xline(ax1, ep(EP,2), 'r--');
legend(ax1, [hB hL hR], {'breath','whisk L','whisk R'}, 'Location','northeastoutside');
title(ax1, sprintf('%s n%d  epoch %d/%d  [%.2f-%.2f s, %.2f s]', ...
    animal, kRun, EP, size(ep,1), ep(EP,1), ep(EP,2), ep(EP,2)-ep(EP,1)));

% ---- MIDDLE: Hilbert phase overlay ----
ax2 = subplot(3,2,[3 4]); hold(ax2,'on'); grid(ax2,'on');
plot(ax2, t(w), phiB(w), '.', 'Color', colB, 'MarkerSize', 5);
plot(ax2, t(w), phiL(w), '.', 'Color', colL, 'MarkerSize', 5);
plot(ax2, t(w), phiR(w), '.', 'Color', colR, 'MarkerSize', 5);
xline(ax2, ep(EP,1), 'r--'); xline(ax2, ep(EP,2), 'r--');
xlim(ax2, xl); ylim(ax2,[-pi pi]); yticks(ax2,[-pi 0 pi]); yticklabels(ax2,{'-\pi','0','\pi'});
xlabel(ax2,'time (s)'); ylabel(ax2,'Hilbert phase');
legend(ax2, {'breath','whisk L','whisk R'}, 'Location','northeastoutside');
title(ax2, 'Hilbert instantaneous phase (0 = insp peak / protraction peak)');

% ---- BOTTOM-LEFT: xL vs xR scatter for THIS epoch ----
ax3 = subplot(3,2,5); hold(ax3,'on'); grid(ax3,'on');
scatter(ax3, xL(we), xR(we), 1, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
axis(ax3,'equal');
pf = polyfit(xL(we), xR(we), 1); xs = [min(xL(we)) max(xL(we))];
plot(ax3, xs, polyval(pf,xs), 'r-', 'LineWidth', 1.6);            % linear fit
xlabel(ax3,'xL (deg, BP)'); ylabel(ax3,'xR (deg, BP)');
[rr,pp] = corr(xL(we), xR(we));
title(ax3, sprintf('xL vs xR  (epoch %d, r=%.2f  p=%.1e, slope=%.2f)', EP, rr, pp, pf(1)));

% ---- BOTTOM-RIGHT: averaged-FFT periodogram PSD for THIS epoch ----
% |FFT|^2 of overlapping segments, averaged; band = +/-1 STD in LOG space.
ax4 = subplot(3,2,6); hold(ax4,'on'); grid(ax4,'on');
seg = min(round(PSD_SEC*fpsW), floor(sum(we)/3)); nov = round(seg/2);
[fS,mB,sB] = psd_fft_avg(fillmissing(braw(we),'linear'), fpsW, seg, nov);
[~ ,mW,sW] = psd_fft_avg([fillmissing(Laf(we),'linear') fillmissing(Raf(we),'linear')], fpsW, seg, nov);
psdband(ax4, fS, mB, sB, [0 0 0]);
psdband(ax4, fS, mW, sW, [0 0.4 0.85]);
hb=plot(ax4, fS, 10.^mB, 'k-', 'LineWidth',1.5);
hw=plot(ax4, fS, 10.^mW, '-', 'Color',[0 0.4 0.85], 'LineWidth',1.5);
set(ax4,'YScale','log'); xlim(ax4, PSD_FPASS);
ib = fS>=1 & fS<=PSD_FPASS(2); ff=fS(ib);
[~,ipb]=max(mB(ib)); [~,ipw]=max(mW(ib));
xlabel(ax4,'freq (Hz)'); ylabel(ax4,'power');
legend(ax4,[hb hw], {sprintf('breath (pk %.1f Hz)',ff(ipb)), ...
    sprintf('whisk (pk %.1f Hz)',ff(ipw))}, 'Location','northeast');
title(ax4, sprintf('FFT avg periodogram, %g s seg (epoch %d, %d segs)', PSD_SEC, EP, max(1,floor((sum(we)-seg)/nov)+1)));
end

function [f, mlog, slog] = psd_fft_avg(X, fps, seg, nov)
% averaged periodogram: |fft|^2 of detrended, Hann-windowed overlapping segments
% (pooled across columns of X). mlog/slog = mean/STD of log10(power) across segments.
    seg = max(round(seg),8); nfft = 2^nextpow2(seg); w = hann(seg); step = max(seg-nov,1);
    P = [];
    for c = 1:size(X,2)
        x = X(:,c); st = 1:step:(numel(x)-seg+1);
        for k = 1:numel(st)
            s = detrend(x(st(k):st(k)+seg-1)) .* w;
            P(end+1,:) = abs(fft(s,nfft)).^2; %#ok<AGROW>
        end
    end
    f = (0:nfft-1)*(fps/nfft); half = f <= fps/2; f = f(half);
    if isempty(P), mlog = nan(1,nnz(half)); slog = zeros(1,nnz(half)); return; end
    Plog = log10(P(:,half) + eps);
    mlog = mean(Plog,1);
    if size(Plog,1) < 2, slog = zeros(size(mlog)); else, slog = std(Plog,0,1); end
end
function psdband(ax, f, mlog, slog, col)
    if all(isnan(mlog)), return; end
    patch(ax, [f fliplr(f)], [10.^(mlog+slog) fliplr(10.^(mlog-slog))], col, ...
        'FaceAlpha',0.15, 'EdgeColor','none');
end

% ================= helpers =================
function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac*prctile(env,95);
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1;
    ep = [t(s) t(e)];
    if ~isempty(ep)
        m = ep(1,:);
        for i=2:size(ep,1)
            if ep(i,1)-m(end,2) <= mergeGap, m(end,2)=ep(i,2); else, m(end+1,:)=ep(i,:); end %#ok<AGROW>
        end
        ep = m(m(:,2)-m(:,1) >= minDur, :);
    end
end
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
