function ca_spike_kernel_fit_260722()
%CA_SPIKE_KERNEL_FIT_260722  Measure the calcium kernel from supervised labels.
%
% Step 1 of the calibrated-deconvolution spike detector. Uses the hand-curated
% events in ca_spike_data.mat as ground truth to measure, per cell type:
%
%   - rise time (10-90%), per event (immune to click jitter, unlike an average)
%   - decay tau from a single-exponential fit
%   - the implied AR(1) coefficient  g = exp(-1/(fps*tau))
%   - whether AR(1) (instantaneous rise) is even an adequate model
%
% These replace the hardcoded g = 0.93 in calcium_spike_gui.m:761.
%
% Read-only w.r.t. session folders. Writes a summary + figure to
%   <ROOT_DIR>\_spike_kernel_260722\
%
% Events are used only if they are (a) high-confidence, SNR >= SNR_MIN in
% robust noise units, and (b) isolated, so a neighbouring event cannot
% contaminate the decay -- the bias that inflated the first pass.

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = 'D:\Ventral_surface_summary';
OUT_SUB    = '_spike_kernel_260722';
SNR_MIN    = 4;        % keep events >= this many robust sigma (high-confidence)
ISO_PRE    = 1.0;      % s: no other labelled event this far BEFORE
ISO_POST   = 3.0;      % s: ...or this far AFTER (>= 3*tau for the slow groups)
BASE_WIN   = [-1.0 -0.3];  % s: local pre-event baseline window
FIT_WIN    = 3.0;      % s: decay fit window after the peak
FIT_FLOOR  = 0.15;     % fit decay only while y > this fraction of amplitude
FALLBACK_FPS = 30;
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

outDir = fullfile(ROOT_DIR, OUT_SUB);
if ~isfolder(outDir), mkdir(outDir); end

d = dir(fullfile(ROOT_DIR, '**', 'ca_spike_data.mat'));
assert(~isempty(d), 'No ca_spike_data.mat under %s', ROOT_DIR);
fprintf('ca_spike_kernel_fit: %d labelled session(s)\n\n', numel(d));

E = struct('group',{},'session',{},'roi',{},'snr',{},'amp',{}, ...
           'rise',{},'tau',{},'r2',{},'wave',{});
tGrid = [];   % common time grid for the stacked waveforms

for i = 1:numel(d)
    fo  = d(i).folder;
    rel = fo(numel(ROOT_DIR)+2:end);
    grp = strtok(rel, filesep);
    if strcmp(grp, 'Vglut2_test'), grp = 'Vglut2'; end   % same animal line

    dh = dir(fullfile(fo, '*_dFF.mat'));
    if isempty(dh), fprintf(2,'  no _dFF.mat: %s\n', rel); continue; end
    L = load(fullfile(dh(1).folder, dh(1).name), 'dFF');
    if ~isfield(L, 'dFF'), continue; end
    dFF = L.dFF; [T, N] = size(dFF);

    S = load(fullfile(fo, 'ca_spike_data.mat'), 'roi_spikes');
    rs = S.roi_spikes;
    if numel(rs) ~= N
        fprintf(2,'  ROI count mismatch (%d vs %d): %s\n', numel(rs), N, rel);
        continue;
    end

    fps = detect_session_fps(fo, FALLBACK_FPS);

    % robust per-ROI noise from the frame-to-frame difference (transient-proof)
    sig = median(abs(diff(dFF,1,1)), 1) * 1.4826 / sqrt(2);

    pre  = round(ISO_PRE  * fps);
    post = round(ISO_POST * fps);
    if isempty(tGrid), tGrid = (-pre:post)'/fps; end

    nKept = 0;
    for k = 1:N
        idx = rs(k).spike_idx;
        if isempty(idx), continue; end
        idx = sort(idx(:));

        % isolation: nearest labelled neighbour must clear both windows
        dl = [inf; diff(idx)];        % gap to previous event
        dr = [diff(idx); inf];        % gap to next event
        iso = dl > pre & dr > post;

        % keep events with full context inside the recording
        inb = idx > pre & idx <= T - post;
        cand = idx(iso & inb);

        for e = 1:numel(cand)
            j = cand(e);
            w = dFF(j-pre:j+post, k);

            bi = tGrid >= BASE_WIN(1) & tGrid <= BASE_WIN(2);
            b  = median(w(bi));
            y  = w - b;
            amp = y(pre+1);                       % height above local baseline
            if amp <= 0, continue; end
            snr = amp / sig(k);
            if snr < SNR_MIN, continue; end

            [rise, tau, r2] = fit_event(y, amp, tGrid, fps, FIT_WIN, FIT_FLOOR);
            if isnan(tau), continue; end

            E(end+1) = struct('group',grp, 'session',rel, 'roi',k, ...
                'snr',snr, 'amp',amp, 'rise',rise, 'tau',tau, 'r2',r2, ...
                'wave',{y/amp}); %#ok<AGROW>
            nKept = nKept + 1;
        end
    end
    fprintf('%-52s fps=%g  %4d/%4d events kept\n', ...
        rel(1:min(52,end)), fps, nKept, sum([rs.n_spikes]));
end

assert(~isempty(E), 'No events survived the SNR/isolation filters.');

%% ------------------------------ per-group summary --------------------------
groups = unique({E.group});
fprintf('\n=================== KERNEL PER CELL TYPE ===================\n');
fprintf('%-8s %6s %5s  %-16s %-18s %-16s %6s\n', ...
    'group','nEvt','nSes','rise 10-90% (ms)','tau (s)','g_AR1','fit R2');
K = struct();
for gi = 1:numel(groups)
    g  = groups{gi};
    m  = strcmp({E.group}, g);
    ri = [E(m).rise]*1000;  ta = [E(m).tau];  r2 = [E(m).r2];
    nSes = numel(unique({E(m).session}));

    [riMed, riCI] = med_ci(ri);
    [taMed, taCI] = med_ci(ta);
    gAR = exp(-(1/30)/taMed);
    gCI = exp(-(1/30)./taCI);

    fprintf('%-8s %6d %5d  %5.0f [%3.0f %3.0f]  %5.2f [%4.2f %4.2f]  %5.3f [%5.3f %5.3f] %6.2f\n', ...
        g, sum(m), nSes, riMed, riCI(1), riCI(2), taMed, taCI(1), taCI(2), ...
        gAR, gCI(1), gCI(2), median(r2));

    K.(g) = struct('n',sum(m), 'nSessions',nSes, ...
        'rise_ms',riMed, 'rise_ci',riCI, 'tau_s',taMed, 'tau_ci',taCI, ...
        'g_AR1_at30fps',gAR, 'g_ci',gCI, 'fit_r2',median(r2), ...
        'wave_median', median(cat(2, E(m).wave), 2), ...
        'rise_frames_at30fps', riMed/1000*30);
end
fprintf('\nGUI currently hardcodes g = 0.93  ->  tau = %.2f s (calcium_spike_gui.m:761)\n', ...
    -(1/30)/log(0.93));
fprintf(['AR(1) assumes an INSTANTANEOUS rise. Any group whose rise exceeds ~1 frame\n' ...
         '(33 ms at 30 fps) violates that and needs a rise term (AR(2)) or a\n' ...
         'matched filter instead.\n']);

%% ------------------------------ per-session spread -------------------------
fprintf('\n============ PER-SESSION tau (is one kernel per group enough?) ============\n');
fprintf('%-52s %6s %8s %8s\n','session','nEvt','tau_s','rise_ms');
sessions = unique({E.session});
for si = 1:numel(sessions)
    m = strcmp({E.session}, sessions{si});
    if sum(m) < 5, continue; end
    fprintf('%-52s %6d %8.2f %8.0f\n', sessions{si}(1:min(52,end)), ...
        sum(m), median([E(m).tau]), median([E(m).rise])*1000);
end

%% ------------------------------ figure -------------------------------------
f = figure('Color','w','Position',[80 80 1180 460]);
subplot(1,3,1); hold on; grid on; box on;
cols = lines(numel(groups));
for gi = 1:numel(groups)
    plot(tGrid, K.(groups{gi}).wave_median, 'LineWidth', 1.8, 'Color', cols(gi,:));
end
xline(0,'k:'); yline(0,'k:');
xlabel('time from labelled peak (s)'); ylabel('dF/F (peak-normalised)');
title('median event waveform'); legend(groups,'Location','northeast','Box','off');
xlim([-0.5 2.5]);

subplot(1,3,2); hold on; grid on; box on;
for gi = 1:numel(groups)
    m = strcmp({E.group}, groups{gi});
    histogram([E(m).tau], 0:0.05:2.5, 'Normalization','probability', ...
        'DisplayStyle','stairs', 'EdgeColor', cols(gi,:), 'LineWidth', 1.5);
end
xline(-(1/30)/log(0.93), 'k--', 'hardcoded g=0.93');
xlabel('decay \tau (s)'); ylabel('fraction of events'); title('per-event \tau');

subplot(1,3,3); hold on; grid on; box on;
for gi = 1:numel(groups)
    m = strcmp({E.group}, groups{gi});
    histogram([E(m).rise]*1000, 0:16.7:500, 'Normalization','probability', ...
        'DisplayStyle','stairs', 'EdgeColor', cols(gi,:), 'LineWidth', 1.5);
end
xline(1000/30, 'k--', '1 frame');
xlabel('rise 10-90% (ms)'); ylabel('fraction of events'); title('per-event rise');
sgtitle(sprintf('Calcium kernel from %d supervised events (SNR>=%g, isolated)', ...
    numel(E), SNR_MIN));

saveas(f, fullfile(outDir, 'spike_kernel_summary.png'));
save(fullfile(outDir, 'spike_kernel_fit.mat'), 'K', 'E', 'tGrid', ...
     'SNR_MIN', 'ISO_PRE', 'ISO_POST', 'FIT_WIN', 'FIT_FLOOR', '-v7.3');
fprintf('\nSaved: %s\n', outDir);
end

%% =============================== LOCAL FUNCTIONS ===============================
function [rise, tau, r2] = fit_event(y, amp, tGrid, fps, fitWin, floorFrac)
% Per-event rise (10-90%) and single-exponential decay tau.
% Measured on the individual event, so peak-click jitter cannot blur it the
% way it blurs an event-triggered average.
    rise = NaN; tau = NaN; r2 = NaN;
    ip = find(tGrid >= 0, 1);          % labelled peak sample

    % ---- rise: walk back from the peak to the 10% and 90% crossings ----
    up = y(1:ip);
    i90 = find(up <= 0.9*amp, 1, 'last');
    i10 = find(up <= 0.1*amp, 1, 'last');
    if ~isempty(i10) && ~isempty(i90) && i90 >= i10
        rise = (i90 - i10) / fps;
    end

    % ---- decay: exponential fit from the peak forward ----
    jEnd = find(tGrid <= fitWin, 1, 'last');
    yd = y(ip:jEnd);
    xd = (0:numel(yd)-1)'/fps;
    ok = yd > floorFrac*amp;
    % only the leading contiguous run: once it drops below the floor, stop
    stop = find(~ok, 1, 'first');
    if ~isempty(stop), ok(stop:end) = false; end
    if sum(ok) < 4, return; end

    % weighted log-linear fit; weights = y, which approximates a nonlinear
    % least-squares fit on the raw trace instead of over-weighting the tail
    xw = xd(ok); yw = yd(ok); w = yw;
    X  = [xw, ones(size(xw))];
    Wm = diag(w);
    b  = (X' * Wm * X) \ (X' * Wm * log(yw));
    if b(1) >= 0, return; end
    tau = -1/b(1);
    if tau <= 0 || tau > 10, tau = NaN; return; end

    yhat = exp(X * b);
    r2   = 1 - sum(w.*(yw-yhat).^2) / sum(w.*(yw-mean(yw)).^2);
end

function [m, ci] = med_ci(x)
% Median with a bootstrap 95% CI (fixed resampling, no RNG dependence issues).
    x = x(:); x = x(isfinite(x));
    m = median(x);
    if numel(x) < 8, ci = [NaN NaN]; return; end
    nB = 2000; n = numel(x);
    rs = RandStream('mt19937ar', 'Seed', 260722);
    bm = zeros(nB,1);
    for b = 1:nB
        bm(b) = median(x(randi(rs, n, n, 1)));
    end
    ci = prctile(bm, [2.5 97.5]);
end
