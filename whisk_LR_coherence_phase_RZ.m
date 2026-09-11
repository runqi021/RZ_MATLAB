% whisk_LR_coherence_phase_RZ.m
%
% Standalone: bilateral whisking COHERENCE SPECTRUM + DELTA PHASE vs frequency.
%
% Reuses the whisker processing standard from
% whisk_LR_integral_xcorr_filtered_RZ.m:
%
%   Whisker:
%     La/Ra (DLC angles) -> fill NaNs -> demean
%     BP 5-50 Hz, zero-phase filtfilt
%     whisk epoch mask from Hilbert envelope of LEFT BP whisker
%     hard envelope threshold = 5 deg
%
% Coherence:
%   Multitaper (Chronux coherencyc). Whisk epochs are tiled into fixed-length
%   non-overlapping windows; every window from every session is stacked as a
%   "trial" and the coherence is trial-averaged with jackknife error bars.
%   -> ONE population coherence spectrum C(f) and phase phi(f).
%
% Phase convention:
%   data1 = L, data2 = R, phi = angle( <conj(J_L) .* J_R> ).
%   phi > 0  ->  R leads L   (L lags R)
%   phi < 0  ->  L leads R
%
% Outputs (figure + .mat) land in <dataRoot>\whisk_LR_coherence_phase\

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW = 400;            % whisker camera frame rate (Hz)

WHISK_BP = [5 50];     % whisker band-pass (Hz)

ENV_THR    = 5;        % hard envelope threshold, deg (whisk-epoch gate)
ENV_SMOOTH = 0.05;     % s
MIN_DUR    = 1;        % s, minimum whisk-epoch duration
MERGE_GAP  = 0.2;      % s, merge whisk epochs closer than this

% --- multitaper coherence settings ---
WIN_SEC = 1.0;         % length of each coherence window (s) -> freq res = TW/WIN_SEC*2
TW      = 3;           % time-bandwidth product (tapers = [TW 2*TW-1])
FPASS   = [1 60];      % frequency band to display (Hz)
ALPHA   = 0.001;        % significance level for jackknife error / confidence

EXCLUDE = "5840027";

doSave  = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

outDir = fullfile(char(dataRoot), 'whisk_LR_coherence_phase');
if doSave && ~isfolder(outDir), mkdir(outDir); end

ord = 4;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2), 'bandpass');

winSamp  = round(WIN_SEC * fpsW);          % samples per coherence window
% require each window long enough to host the tapers
assert(winSamp >= 2*TW, 'WIN_SEC too short for TW: need WIN_SEC*fpsW >= 2*TW');

% Chronux params for the final trial-averaged coherence
pc.Fs       = fpsW;
pc.tapers   = [TW, 2*TW-1];
pc.pad      = 0;
pc.fpass    = FPASS;
pc.err      = [2 ALPHA];   % jackknife error bars
pc.trialave = 1;           % average across stacked windows

% ============================ ACCUMULATORS ============================

L_win = zeros(winSamp, 0);   % stacked L windows  [winSamp x nTrials]
R_win = zeros(winSamp, 0);   % stacked R windows
nsess = 0;
nwin  = 0;

% ============================ SESSION LOOP ============================

ad = dir(char(dataRoot));

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit'))
        continue;
    end
    if any(strcmp(ad(ai).name, EXCLUDE))
        continue;
    end

    animal = ad(ai).name;

    rr = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    [~,o] = sort({rr.name});
    rr = rr(o);

    for kk = 1:numel(rr)

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            continue;
        end

        % ============================ WHISKER ============================
        try
            M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
        catch ME
            warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

        t = (0:numel(La)-1)' / fpsW;
        N = numel(t);

        La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
        Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');
        La0 = La0 - mean(La0, 'omitnan');
        Ra0 = Ra0 - mean(Ra0, 'omitnan');

        xL = filtfilt(b_wbp, a_wbp, La0);   % L BP 5-50 Hz
        xR = filtfilt(b_wbp, a_wbp, Ra0);   % R BP 5-50 Hz

        % ======================= WHISK EPOCH MASK ========================
        env_raw = abs(hilbert(xL));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));
        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, t, MIN_DUR, MERGE_GAP);

        if isempty(ep), continue; end

        % ===================== TILE WINDOWS IN EPOCHS ====================
        sess_added = 0;
        for q = 1:size(ep,1)
            i0 = find(t >= ep(q,1), 1, 'first');
            i1 = find(t <= ep(q,2), 1, 'last');
            if isempty(i0) || isempty(i1), continue; end

            s = i0;
            while s + winSamp - 1 <= i1
                idx = s:(s+winSamp-1);
                wl = xL(idx); wr = xR(idx);
                if all(isfinite(wl)) && all(isfinite(wr))
                    wl = wl - mean(wl);
                    wr = wr - mean(wr);
                    L_win(:, end+1) = wl; %#ok<SAGROW>
                    R_win(:, end+1) = wr; %#ok<SAGROW>
                    sess_added = sess_added + 1;
                end
                s = s + winSamp;            % non-overlapping windows
            end
        end

        if sess_added > 0
            nsess = nsess + 1;
            nwin  = nwin + sess_added;
            fprintf('  %s n%d: %d whisk windows (%.1fs total)\n', ...
                    animal, kk, sess_added, sess_added*WIN_SEC);
        end
    end
end

fprintf('\n%d sessions, %d whisk windows (%.0f s) contributed\n', ...
        nsess, nwin, nwin*WIN_SEC);

if nwin < 2
    error('Too few whisk windows (%d) to estimate coherence.', nwin);
end

% ===================== POPULATION COHERENCE =====================

[~, C, phi, ~, ~, ~, f, confC, phistd, Cerr] = coherencyc(L_win, R_win, pc);

f      = f(:);
C      = C(:);
phi    = phi(:);
phistd = phistd(:);
Clo    = Cerr(1,:)';
Chi    = Cerr(2,:)';

phi_deg    = rad2deg(phi);
phistd_deg = rad2deg(phistd);

sig = C >= confC;                       % significant frequencies

[Cpk, ipk] = max(C);
fprintf('Peak coherence C=%.3f at %.2f Hz, phase=%.1f deg (R-leads-L if +)\n', ...
        Cpk, f(ipk), phi_deg(ipk));

% ============================== FIGURE ==============================

col  = [0.10 0.30 0.85];
colS = [0.6 0.6 0.6];

fig = figure('Color','w','Position',[180 160 620 720]);

% ---- Panel A: coherence magnitude ----
axA = subplot(2,1,1); hold(axA,'on'); grid(axA,'on');
fill(axA, [f; flipud(f)], [Chi; flipud(Clo)], col, ...
     'FaceAlpha',0.20, 'EdgeColor','none');
plot(axA, f, C, '-', 'Color', col, 'LineWidth', 2);
yline(axA, confC, 'k--', sprintf('conf %.2f (\\alpha=%.2g)', confC, ALPHA), ...
      'LabelHorizontalAlignment','left', 'Interpreter','tex');
plot(axA, f(ipk), Cpk, 'o', 'MarkerEdgeColor','k', 'MarkerFaceColor',col, 'MarkerSize',6);
xlim(axA, FPASS); ylim(axA, [0 1]);
xlabel(axA, 'frequency (Hz)'); ylabel(axA, 'coherence  |C_{LR}|');
title(axA, sprintf('Bilateral whisking coherence\nN=%d windows (%.0f s), %d sessions, TW=%d, win=%.1fs', ...
      nwin, nwin*WIN_SEC, nsess, TW, WIN_SEC), 'Interpreter','tex');
box(axA,'off');

% ---- Panel B: delta phase ----
axB = subplot(2,1,2); hold(axB,'on'); grid(axB,'on');
yline(axB, 0, 'k-');

% faint full-band phase + error
fill(axB, [f; flipud(f)], ...
     [phi_deg+1.96*phistd_deg; flipud(phi_deg-1.96*phistd_deg)], ...
     colS, 'FaceAlpha',0.20, 'EdgeColor','none');
plot(axB, f, phi_deg, '-', 'Color', colS, 'LineWidth', 1);

% bold phase where coherence is significant
ph_sig = phi_deg; ph_sig(~sig) = NaN;
plot(axB, f, ph_sig, '-', 'Color', col, 'LineWidth', 2.2);

xlim(axB, FPASS); ylim(axB, [-180 180]);
yticks(axB, -180:90:180);
xlabel(axB, 'frequency (Hz)');
ylabel(axB, '\Delta phase L,R  (deg)');
title(axB, '\Delta phase  ( + : R leads L,   - : L leads R );  bold = significant', ...
      'Interpreter','tex');
box(axB,'off');

% ============================== SAVE ==============================

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LR_coherence_phase.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LR_coherence_phase.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'whisk_LR_coherence_phase_data.mat'), ...
         'f','C','phi','phi_deg','phistd','phistd_deg','Clo','Chi', ...
         'confC','sig','Cpk','ipk','nsess','nwin', ...
         'WIN_SEC','TW','FPASS','ALPHA','WHISK_BP','ENV_THR','ENV_SMOOTH', ...
         'MIN_DUR','MERGE_GAP','fpsW','EXCLUDE');
    fprintf('Saved coherence figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:));
    d = diff([false; a; false]);
    s = find(d == 1);
    e = find(d == -1) - 1;
    ep = [t(s), t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap
            ep2(end,2) = ep(i,2);
        else
            ep2(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end
    ep = ep2;
    ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
