% thermal_nostril_breath_view.m
% View the single-video thermal nostril extraction from
% thermal_nostril_breath_single.py (<ats_stem>_nostrilC.mat).
%
% Default disk-ROI QC ("is it good?" pass before the interactive ROI draw).
% FIGURE 1 (RAW): avg projection (deg C, blue-white-red) + raw ROI temp + raw PSD.
% FIGURE 2 (BP) : bandpassed + inverted (inhale up) trace + PSD.
% (Uses subplot, not tiledlayout: R2021b has a colorbar+tiledlayout+imagesc bug.)

matPath = "D:\260615_thermalNbasler\5840027\cam1_20260615_203627_run001\Rec-A6753sc_00122-167_03_35_47_124-0111_nostrilC.mat";
BP      = [1.5 12.5];   % bandpass (Hz)
INVERT  = true;     % inhale cools the nostril -> invert so inhale reads as a rise

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'mod','bluewhitered'));
S = load(matPath);
fps   = double(S.fps_trace);
sides = {'L','R'};

% ---- FIGURE 1: RAW (avg proj + raw trace + raw PSD) ----
figure('Color','w','Position',[60 60 1250 720]);
for i = 1:2
    sd  = sides{i};
    avg = double(S.(sd+"_avg"));
    tr  = double(S.(sd+"_trace_C"));
    sidelbl = char(S.(sd+"_side"));
    t   = (0:numel(tr)-1)/fps;

    ax = subplot(2,3,(i-1)*3 + 1);
    imagesc(ax, avg); axis(ax,'image','off'); colormap(ax, bluewhitered(256)); colorbar(ax);
    title(ax, sprintf('%s avg proj (\\circC)', sidelbl));

    ax = subplot(2,3,(i-1)*3 + 2);
    plot(ax, t, tr, '-'); grid(ax,'on'); xlabel(ax,'s'); ylabel(ax,'\circC');
    title(ax, sprintf('%s RAW ROI temp (mean %.2f \\circC)', sidelbl, mean(tr,'omitnan')));

    ax = subplot(2,3,(i-1)*3 + 3);
    x = tr - mean(tr,'omitnan'); nfft = 2^nextpow2(numel(x));
    P = abs(fft(x, nfft)).^2; fr = (0:nfft-1)*(fps/nfft);
    keep = fr <= min(40, fps/2);
    plot(ax, fr(keep), P(keep), '-'); grid(ax,'on'); hold(ax,'on');
    xline(ax, BP(1), 'r--'); xline(ax, BP(2), 'r--');
    inb = fr>=BP(1) & fr<=BP(2); [~,ip] = max(P(inb)); frb = fr(inb);
    xlabel(ax,'Hz'); ylabel(ax,'power');
    title(ax, sprintf('%s raw PSD (in-band peak %.2f Hz)', sidelbl, frb(ip)));
end
sgtitle(sprintf('Nostril thermal breathing — RAW  (win %dpx, disk r=%.1f)', ...
    double(S.win), double(S.disk_radius)));

% ---- FIGURE 2: BANDPASSED + inverted ----
figure('Color','w','Position',[140 100 1150 600]);
for i = 1:2
    sd  = sides{i};
    tr  = double(S.(sd+"_trace_C"));
    sidelbl = char(S.(sd+"_side"));
    [bb,aa] = butter(2, BP/(fps/2), 'bandpass');
    trbp = filtfilt(bb, aa, tr - mean(tr,'omitnan'));
    if INVERT, trbp = -trbp; end
    t = (0:numel(trbp)-1)/fps;

    ax = subplot(2,2,i);
    plot(ax, t, trbp, '-'); grid(ax,'on'); xlabel(ax,'s'); ylabel(ax,'\circC (BP, inhale up)');
    title(ax, sprintf('%s breathing %g-%g Hz', sidelbl, BP(1), BP(2)));

    ax = subplot(2,2,i+2);
    nfft = 2^nextpow2(numel(trbp));
    P = abs(fft(trbp, nfft)).^2; fr = (0:nfft-1)*(fps/nfft);
    keep = fr <= min(40, fps/2);
    plot(ax, fr(keep), P(keep), '-'); grid(ax,'on'); hold(ax,'on');
    xline(ax, BP(1), 'r--'); xline(ax, BP(2), 'r--');
    inb = fr>=BP(1) & fr<=BP(2); [~,ip] = max(P(inb)); frb = fr(inb);
    xlabel(ax,'Hz'); ylabel(ax,'power');
    title(ax, sprintf('%s PSD (peak %.2f Hz)', sidelbl, frb(ip)));
end
invstr = ''; if INVERT, invstr = ' (inverted: inhale up)'; end
sgtitle(sprintf('Nostril thermal breathing — BANDPASS %g-%g Hz%s', BP(1), BP(2), invstr));
