% drawROI_N_lpsub_powermap.m
% SIMPLE batch ROI labeller (original flow) but you draw the ellipse on the
% 2-10 Hz BREATH-SNR map (so you can see where breathing is). For each video:
%   1) two windows pop up (LEFT then RIGHT) showing the breath-SNR map ->
%      draw an ellipse (or circle), DOUBLE-CLICK to confirm.
%   2) a traces figure pops out (LEFT / RIGHT / AVERAGE detrended).
%   3) CLOSE that figure -> choose Save & Next / Redo / Stop.
% Saves <stem>_nostrilROI.mat (ROI) + <stem>_breath.mat (avg breathing).
% LP-subtraction detrend, zero-phase, no flips. RUN IN MATLAB (interactive).
close all; clc; clear;

% ============================ USER-EDITABLE ============================
VIDEOS_DIR = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
dataRoot   = "D:\260615_thermalNbasler";
ROI_SHAPE  = 'ellipse';   % 'ellipse' | 'circle'
STAT       = 'mean';      % per-frame ROI readout: 'mean' | 'median' | 'max' | 'min'
LP_CUT     = 1;           % Hz, low-pass baseline subtracted off
INVERT     = true;        % inhale cools the nostril -> invert so inhale reads as a rise
BAND       = [2 10];      % breathing band (Hz) for the SNR map you draw on
FINE_BP    = [2 10];      % finer band saved as breath_bp ([] = skip)
START_AT   = 1;           % resume from this index in the (sorted) csv list
% ======================================================================

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here),'mod','bluewhitered'));

d = dir(fullfile(char(VIDEOS_DIR), '*DLC*.csv'));
files = sort(string({d.name}));
assert(~isempty(files), 'no DLC csvs in %s', VIDEOS_DIR);
fprintf('BATCH: %d videos in %s\n', numel(files), VIDEOS_DIR);

saved = strings(0); stopped = false;
for i = START_AT:numel(files)
    csv = fullfile(VIDEOS_DIR, files(i));
    fprintf('\n[%d/%d] %s\n', i, numel(files), files(i));
    P = thermal_resolve_paths(csv, dataRoot);
    S = load(P.nostrilC);
    fps = double(S.L_stack_fps);
    [lb,la] = butter(2, LP_CUT/(fps/2), 'low');

    while true   % redo loop
        close all;
        R = draw_and_detrend(S, ROI_SHAPE, STAT, lb, la, INVERT, BAND, fps, ...
                             sprintf('[%d/%d] %s', i, numel(files), files(i)));
        detrA = (R.L.detr + R.R.detr) / 2;
        t = (0:numel(detrA)-1)/fps;
        fh = plot_three(R, detrA, t, fps, LP_CUT, INVERT, sprintf('[%d/%d] %s', i, numel(files), files(i)));
        drawnow; waitfor(fh);          % examine the traces, CLOSE the figure to get the prompt

        choice = questdlg(sprintf('%s  (%d/%d)\nKeep this ROI?', files(i), i, numel(files)), ...
            'ROI decision', 'Save & Next', 'Redo', 'Stop', 'Save & Next');
        if isempty(choice), choice = 'Stop'; end

        if strcmp(choice, 'Redo'),  continue;  end
        if strcmp(choice, 'Save & Next')
            save_outputs(P, R, detrA, t, fps, LP_CUT, FINE_BP, INVERT, STAT, ROI_SHAPE, csv, BAND);
            saved(end+1) = files(i); %#ok<SAGROW>
        end
        if strcmp(choice, 'Stop'), stopped = true; end
        break;
    end
    if stopped, break; end
end
fprintf('\n==== %d saved%s ====\n', numel(saved), ternch(stopped, ' (stopped early)', ''));

% ============================ local functions ============================
function R = draw_and_detrend(S, ROI_SHAPE, STAT, lb, la, INVERT, BAND, fps, tag)
    sides = {'L','R'}; R = struct();
    for i = 1:2
        sd  = sides{i};
        stk = double(S.(sd+"_stack")); avg = double(S.(sd+"_avg"));
        sz  = size(avg); sidelbl = char(S.(sd+"_side")); cc = (sz(1)+1)/2;
        Smap = snrmap(stk, fps, sz, BAND);              % 2-10 Hz breath-SNR map (draw on this)
        figure('Color','w','Position',[120+ (i-1)*640, 120, 600, 560], 'Name',sprintf('%s nostril SNR',sidelbl));
        ax = axes; imagesc(ax, Smap); axis(ax,'image'); colormap(ax, hot); colorbar(ax);
        ax.CLim = [0 max(Smap(:))+eps]; hold(ax,'on'); plot(ax, cc, cc, 'c+','MarkerSize',10,'LineWidth',1.0);
        title(ax, sprintf('%s  —  breath SNR (%g-%g Hz)\nDRAW ROI on the hot-spot, double-click to confirm', ...
              sidelbl, BAND(1), BAND(2)));
        if strcmpi(ROI_SHAPE,'circle'), h = drawcircle(ax); else, h = drawellipse(ax); end
        wait(h); mask = createMask(h);
        flat = reshape(stk, size(stk,1), []); roipix = flat(:, mask(:));
        switch STAT
            case 'mean',   tr = mean(roipix, 2, 'omitnan');
            case 'median', tr = median(roipix, 2, 'omitnan');
            case 'max',    tr = max(roipix, [], 2, 'omitnan');
            case 'min',    tr = min(roipix, [], 2, 'omitnan');
            otherwise,     error('STAT must be mean | median | max | min');
        end
        tr = fillmissing(tr, 'linear', 'EndValues','nearest'); tr(~isfinite(tr)) = 0;
        base = filtfilt(lb, la, tr); detr = tr - base;
        if INVERT, detr = -detr; end
        R.(sd).side = sidelbl; R.(sd).mask = mask; R.(sd).trace = tr; R.(sd).detr = detr;
        R.(sd).fps = double(S.(sd+"_stack_fps")); R.(sd).npix = nnz(mask);
        R.(sd).shape = ROI_SHAPE; R.(sd).stat = STAT;
        if strcmpi(ROI_SHAPE,'circle')
            R.(sd).center = h.Center; R.(sd).radius = h.Radius;
        else
            R.(sd).center = h.Center; R.(sd).semiaxes = h.SemiAxes; R.(sd).angle = h.RotationAngle;
        end
        fprintf('  %s: ROI %d px (peak SNR %.1f)\n', sidelbl, nnz(mask), max(Smap(:)));
    end
    fprintf('  %s\n', tag);
end

function Smap = snrmap(stk, fps, sz, band)
% per-pixel breath SNR = power(band) / power(band(2)..min(40,Nyq)), floored
    F = reshape(stk, size(stk,1), []); nyq = fps/2;
    [bb,ab] = butter(2, band/nyq, 'bandpass');
    lo = band(2); hi = min(40, 0.95*nyq);
    Fc = F; Fc(~isfinite(Fc)) = 0;
    Pin = var(filtfilt(bb,ab,Fc), 0, 1);
    if hi > lo*1.05
        [bn,an] = butter(2, [lo hi]/nyq, 'bandpass'); Pno = var(filtfilt(bn,an,Fc), 0, 1);
    else
        Pno = max(var(Fc,0,1) - Pin, 0);
    end
    snr = Pin ./ (Pno + eps); snr(Pin < 0.10*max(Pin)) = 0; snr(all(isnan(F),1)) = 0;
    Smap = reshape(snr, sz);
end

function fh = plot_three(R, detrA, t, fps, LP_CUT, INVERT, tag)
    fh = figure('Color','w','Position',[200 120 1100 640]);
    ax1 = subplot(3,1,1); plot(ax1, t, R.L.detr,'-'); grid(ax1,'on'); ylabel(ax1,'\circC');
    title(ax1, sprintf('LEFT detrended (raw - LP %g Hz)', LP_CUT));
    ax2 = subplot(3,1,2); plot(ax2, t, R.R.detr,'-'); grid(ax2,'on'); ylabel(ax2,'\circC');
    title(ax2, 'RIGHT detrended');
    ax3 = subplot(3,1,3); plot(ax3, t, detrA,'-'); grid(ax3,'on'); ylabel(ax3,'\circC'); xlabel(ax3,'s');
    x = detrA - mean(detrA,'omitnan'); nf = 2^nextpow2(numel(x));
    Pw = abs(fft(x,nf)).^2; fr = (0:nf-1)*(fps/nf); inb = fr>=LP_CUT & fr<=15;
    [~,ip] = max(Pw(inb)); frb = fr(inb);
    title(ax3, sprintf('L+R AVERAGE (breathing, peak %.2f Hz)', frb(ip)));
    linkaxes([ax1 ax2 ax3],'x');
    invstr = ''; if INVERT, invstr = ' (inhale up)'; end
    sgtitle(sprintf('%s  —  LP-sub %g Hz%s    [CLOSE this figure to choose]', tag, LP_CUT, invstr), 'Interpreter','none');
end

function save_outputs(P, R, detrA, t, fps, LP_CUT, FINE_BP, INVERT, STAT, ROI_SHAPE, csv, BAND)
    R.src = char(csv); R.nostrilC = P.nostrilC;
    save(P.nostrilROI, '-struct', 'R');
    if ~isempty(FINE_BP)
        [fbb,faa] = butter(2, FINE_BP/(fps/2), 'bandpass'); detrA_bp = filtfilt(fbb, faa, detrA);
    else
        detrA_bp = [];
    end
    B = struct();
    B.breath = detrA(:); B.breath_bp = detrA_bp(:); B.fps = fps; B.t = t(:);
    B.method = 'lpsub'; B.lp_cut = LP_CUT; B.fine_bp = FINE_BP; B.inverted = INVERT;
    B.roi_stat = STAT; B.roi_shape = ROI_SHAPE; B.snr_band = BAND; B.animal = P.animal; B.run = P.k;
    B.src_csv = char(csv); B.src_ats = P.ats;
    save(P.breath, '-struct', 'B');
    fprintf('  saved ROI + breath for %s n%d\n', P.animal, P.k);
end

function s = ternch(c, a, b)
    if c, s = a; else, s = b; end
end
