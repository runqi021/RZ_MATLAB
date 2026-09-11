% whisk_breath_phase_reset_RZ.m
%
% Phase resetting of whisking by breathing -- replicates Moore, Deschenes,
% Kleinfeld 2013 Nature SI Fig S4/S5 on the 260615 thermal+Basler data.
%
% Inspiration is treated as a perturbation to the ongoing whisk oscillator.
% For each whisk cycle that contains an inspiration:
%   t_reset  = inspiration time  -  protraction onset of that whisk
%   T0       = duration of the PRECEDING whisk  (expected period)
%   T1       = duration of the perturbed whisk
%   phi_reset  = 2*pi * t_reset / T0           (breath phase within whisk cycle)
%   dphi_whisk = 2*pi * (T1 - T0) / T1         (+ delayed, - advanced)
%
% Conventions (Moore):
%   whisker BP 3-25 Hz, 3-pole Butterworth filtfilt -> Hilbert
%   whisk cycle accepted if peak-to-peak > 5 deg AND duration < 250 ms
%   PHASE 0 = PROTRACTION ONSET = crossing of 10% of the cycle min-to-max
%   breath BP 1-15 Hz; inspiration event detected as in the working pipeline
%   inspiration shifted earlier by INSP_LEAD_MS (drive reaches face first)
%   basal vs sniff by instantaneous breath rate (mouse: basal<5, sniff>6 Hz)
%
% Figures:
%   FIG 1 (S4): (a) PDF of phi_reset, basal;  (b) PRC dphi_whisk vs phi_reset, basal
%   FIG 2 (S5): (a) whisk period shift (T1-T0) vs t_reset, basal+sniff
%               (b) breath period shift (T1-T0) vs t_reset (whisk perturbs breath)

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [3 25];   % Hz (Moore)
BREATH_BP = [1 15];   % Hz (Moore)

% whisk-epoch gate (only analyse active whisking)
ENV_THR    = 5;       % deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

% whisk cycle detection / acceptance (Moore)
WHISK_MINDIST_S = 0.030;   % min trough spacing (<= 1/maxfreq-ish)
WHISK_PROM_DEG  = 1.0;     % trough prominence (deg) on BP angle
AMP_MIN_DEG     = 5.0;     % accept whisk if peak-to-peak > this
DUR_MAX_S       = 0.250;   % accept whisk if cycle duration < this
ONSET_FRAC      = 0.10;    % protraction onset = 10% of min-to-max

% inspiration detection
MIN_INSP_DIST_S = 0.05;
PROM_MAD        = 0.5;
INSP_LEAD_MS    = 30;      % shift inspiration earlier by this (drive lead)

% breathing-state thresholds (mouse)
BASAL_HZ = 5;              % instantaneous breath rate < this = basal
SNIFF_HZ = 6;              % > this = sniff; between = unclassified

% period-stability guard for T0/T1 (drop pathological cycles)
T_MIN_S = 0.040;
T_MAX_S = 0.300;

N_PHI_BINS = 18;           % phi_reset bins for the PRC mean curve

EXCLUDE = "5840027";
doSave  = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_breath_phase_reset');
if doSave && ~isfolder(outDir), mkdir(outDir); end

ord = 3;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

% ---- accumulators ----
% S4 (basal only): phi_reset (rad), dphi_whisk (rad)
PHI_b = []; DPHI_b = [];
% S5a whisk-reset-by-breath: t_reset (s), T1-T0 (s), state (1=basal 2=sniff)
WR_treset = []; WR_dT = []; WR_state = [];
% S5b breath-reset-by-whisk: t_reset (s), T1-T0 (s), state
BR_treset = []; BR_dT = []; BR_state = [];

nSess = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    % ---------------------------- WHISKER ----------------------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)' / fpsW; N = numel(tW);

    La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
    Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');
    La0 = La0 - mean(La0,'omitnan'); Ra0 = Ra0 - mean(Ra0,'omitnan');
    xL = filtfilt(b_wbp, a_wbp, La0);
    xR = filtfilt(b_wbp, a_wbp, Ra0);

    % whisk-epoch mask (from L envelope)
    env = movmean(abs(hilbert(xL)), max(1, round(ENV_SMOOTH*fpsW)));
    ep = bool_to_epochs(env > ENV_THR, tW, MIN_DUR, MERGE_GAP);
    m = false(N,1);
    for q = 1:size(ep,1), m = m | (tW>=ep(q,1) & tW<=ep(q,2)); end
    if nnz(m) < 50, continue; end

    % ---------------------------- BREATH -----------------------------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun)), dataRoot);
        if ~isfile(Pn.breath), continue; end
        Bs = load(Pn.breath);
        br = Bs.breath(:); fb = double(Bs.fps); tBr = (0:numel(br)-1)'/fb;
        brw = interp1(tBr, br, tW, 'linear', NaN);
        breath_finite = isfinite(brw);
        if nnz(breath_finite) < 50, continue; end
        brw = fillmissing(brw,'linear','EndValues','nearest');
        brw_bp = filtfilt(b_br, a_br, brw);
        prom = PROM_MAD * mad(brw_bp,1);
        if ~isfinite(prom) || prom<=0, prom = 0.5*std(brw_bp,'omitnan'); end
        [~, ilocs] = findpeaks(-brw_bp, 'MinPeakDistance', round(MIN_INSP_DIST_S*fpsW), ...
                                        'MinPeakProminence', prom);
    catch ME
        warning('breath failed %s n%d: %s', animal, kRun, ME.message); continue;
    end
    if numel(ilocs) < 3, continue; end

    % inspiration times with drive-lead correction, restricted to recorded breath
    insp_t = tW(ilocs) - INSP_LEAD_MS/1000;
    insp_t = insp_t(insp_t > tW(1) & insp_t < tBr(end));
    if numel(insp_t) < 3, continue; end

    % per-inspiration instantaneous breath rate (1/mean adjacent interval)
    ibi = diff(insp_t);
    insp_rate = nan(size(insp_t));
    insp_rate(2:end-1) = 2 ./ (insp_t(3:end) - insp_t(1:end-2));
    insp_rate(1) = 1/ibi(1); insp_rate(end) = 1/ibi(end);

    % ------------------- WHISK CYCLES (per side) ---------------------
    added = 0;
    for side = 1:2
        if side==1, xS = xL; else, xS = xR; end
        [on_t, dur] = detect_whisk_cycles(xS, tW, fpsW, ...
            WHISK_MINDIST_S, WHISK_PROM_DEG, AMP_MIN_DEG, DUR_MAX_S, ONSET_FRAC, m);
        if numel(on_t) < 3, continue; end

        % ----- S4 / S5a : breathing perturbs whisking -----
        for k = 2:numel(on_t)-1
            T0 = on_t(k)   - on_t(k-1);     % preceding whisk
            T1 = on_t(k+1) - on_t(k);       % perturbed whisk
            if T0<T_MIN_S||T0>T_MAX_S||T1<T_MIN_S||T1>T_MAX_S, continue; end
            % inspiration(s) inside the perturbed whisk [on_t(k), on_t(k+1)]
            ins = insp_t(insp_t>=on_t(k) & insp_t<on_t(k+1));
            if isempty(ins), continue; end
            t_reset = ins(1) - on_t(k);
            rate = insp_rate(find(insp_t==ins(1),1));
            st = state_of(rate, BASAL_HZ, SNIFF_HZ);
            if st==0, continue; end

            dT = T1 - T0;
            WR_treset(end+1,1) = t_reset; %#ok<SAGROW>
            WR_dT(end+1,1)     = dT;      %#ok<SAGROW>
            WR_state(end+1,1)  = st;      %#ok<SAGROW>

            if st==1   % basal -> S4 PRC
                PHI_b(end+1,1)  = 2*pi*t_reset/T0;       %#ok<SAGROW>
                DPHI_b(end+1,1) = 2*pi*(T1-T0)/T1;       %#ok<SAGROW>
            end
            added = added + 1;
        end

        % ----- S5b : whisking perturbs breathing -----
        for j = 2:numel(insp_t)-1
            T0b = insp_t(j)   - insp_t(j-1);
            T1b = insp_t(j+1) - insp_t(j);
            wks = on_t(on_t>=insp_t(j) & on_t<insp_t(j+1));
            if isempty(wks), continue; end
            t_reset_b = wks(1) - insp_t(j);
            st = state_of(insp_rate(j), BASAL_HZ, SNIFF_HZ);
            if st==0, continue; end
            BR_treset(end+1,1) = t_reset_b; %#ok<SAGROW>
            BR_dT(end+1,1)     = T1b - T0b; %#ok<SAGROW>
            BR_state(end+1,1)  = st;        %#ok<SAGROW>
        end
    end

    if added > 0
        nSess = nSess + 1;
        fprintf('  %s n%d: %d reset events (basal so far %d)\n', animal, kRun, added, numel(PHI_b));
    end
end

fprintf('\n%d sessions. whisk-reset events: %d (basal %d, sniff %d). breath-reset events: %d\n', ...
    nSess, numel(WR_dT), nnz(WR_state==1), nnz(WR_state==2), numel(BR_dT));

if isempty(PHI_b), warning('No basal reset events -- S4 panels will be empty.'); end

% ============================== FIGURE 1 (S4) ==============================

fig1 = figure('Color','w','Position',[120 220 980 420]);

axA = subplot(1,2,1); hold(axA,'on'); grid(axA,'on');
if ~isempty(PHI_b)
    edges = linspace(0,2*pi,25);
    histogram(axA, mod(PHI_b,2*pi), edges, 'Normalization','pdf', ...
              'FaceColor',[0.3 0.3 0.3], 'EdgeColor','w');
end
xlim(axA,[0 2*pi]); xticks(axA,[0 pi 2*pi]); xticklabels(axA,{'0','\pi','2\pi'});
xlabel(axA, '\phi_{reset}  (breath phase within whisk cycle)');
ylabel(axA, 'probability density');
title(axA, sprintf('S4b: \\phi_{reset} PDF (basal), n=%d', numel(PHI_b)), 'Interpreter','tex');

axB = subplot(1,2,2); hold(axB,'on'); grid(axB,'on');
if ~isempty(PHI_b)
    phi = mod(PHI_b,2*pi);
    plot(axB, phi, DPHI_b, '.', 'Color',[0.7 0.7 0.7], 'MarkerSize',4);
    % binned circular-x mean
    be = linspace(0,2*pi,N_PHI_BINS+1); bc = 0.5*(be(1:end-1)+be(2:end));
    mu = nan(N_PHI_BINS,1); se = nan(N_PHI_BINS,1);
    for ib=1:N_PHI_BINS
        v = phi>=be(ib) & phi<be(ib+1);
        if nnz(v)>=5, mu(ib)=mean(DPHI_b(v)); se(ib)=std(DPHI_b(v))/sqrt(nnz(v)); end
    end
    errorbar(axB, bc, mu, se, '-o', 'Color',[0.85 0.1 0.1], 'LineWidth',1.6, ...
             'MarkerFaceColor',[0.85 0.1 0.1], 'MarkerSize',4, 'CapSize',0);
end
yline(axB,0,'k--');
xlim(axB,[0 2*pi]); xticks(axB,[0 pi 2*pi]); xticklabels(axB,{'0','\pi','2\pi'});
ylim(axB,[-pi pi]); yticks(axB,[-pi 0 pi]); yticklabels(axB,{'-\pi','0','\pi'});
xlabel(axB, '\phi_{reset}   (protraction \rightarrow retraction)');
ylabel(axB, '\Delta\phi_{whisk}   (- advanced / + delayed)');
title(axB, 'S4c: phase response curve (basal)', 'Interpreter','tex');

% ============================== FIGURE 2 (S5) ==============================

fig2 = figure('Color','w','Position',[140 160 980 420]);
cB = [0.15 0.15 0.15]; cS = [0.85 0.10 0.10];

axC = subplot(1,2,1); hold(axC,'on'); grid(axC,'on');
scatter_state(axC, WR_treset*1000, WR_dT*1000, WR_state, cB, cS);
yline(axC,0,'k--');
xlabel(axC, 't_{reset} (ms): breath time within whisk cycle');
ylabel(axC, 'whisk period shift  T1-T0 (ms)');
title(axC, sprintf('S5a: whisking reset by breathing (basal n=%d, sniff n=%d)', ...
      nnz(WR_state==1), nnz(WR_state==2)), 'Interpreter','tex');
legend(axC, {'basal','sniff'}, 'Location','best');

axD = subplot(1,2,2); hold(axD,'on'); grid(axD,'on');
scatter_state(axD, BR_treset*1000, BR_dT*1000, BR_state, cB, cS);
yline(axD,0,'k--');
xlabel(axD, 't_{reset} (ms): whisk time within breath cycle');
ylabel(axD, 'breath period shift  T1-T0 (ms)');
title(axD, sprintf('S5b: breathing reset by whisking (basal n=%d, sniff n=%d)', ...
      nnz(BR_state==1), nnz(BR_state==2)), 'Interpreter','tex');

% ============================== SAVE ==============================

if doSave
    exportgraphics(fig1, fullfile(outDir,'phase_reset_S4.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig1, fullfile(outDir,'phase_reset_S4.pdf'),'ContentType','vector','BackgroundColor','white');
    exportgraphics(fig2, fullfile(outDir,'phase_reset_S5.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig2, fullfile(outDir,'phase_reset_S5.pdf'),'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'phase_reset_data.mat'), ...
         'PHI_b','DPHI_b','WR_treset','WR_dT','WR_state','BR_treset','BR_dT','BR_state', ...
         'WHISK_BP','BREATH_BP','AMP_MIN_DEG','DUR_MAX_S','ONSET_FRAC','INSP_LEAD_MS', ...
         'BASAL_HZ','SNIFF_HZ','nSess');
    fprintf('Saved S4/S5 figures + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function [on_t, dur] = detect_whisk_cycles(xBP, t, fs, minDistS, promDeg, ampMin, durMax, onsetFrac, mask)
% Whisk cycles from BP angle. Cycle start = trough (protraction start).
% Accept if peak-to-peak > ampMin AND duration < durMax. Protraction onset =
% crossing of onsetFrac of the cycle min-to-max. Onset must be within mask.
    xBP = xBP(:); t = t(:); mask = logical(mask(:));
    [~, tr] = findpeaks(-xBP, 'MinPeakDistance', max(1,round(minDistS*fs)), ...
                              'MinPeakProminence', promDeg);
    on_t = []; dur = [];
    for i = 1:numel(tr)-1
        a = tr(i); b = tr(i+1);
        if b<=a, continue; end
        seg = xBP(a:b);
        [pk, rel] = max(seg);
        mn = xBP(a);
        amp = pk - mn;
        cdur = t(b) - t(a);
        if amp < ampMin || cdur >= durMax, continue; end
        thr = mn + onsetFrac*amp;
        oi = a + find(seg(1:rel) >= thr, 1, 'first') - 1;   % first upward crossing
        if isempty(oi) || ~mask(oi), continue; end
        on_t(end+1,1) = t(oi);  %#ok<AGROW>
    end
    on_t = sort(on_t);
    if numel(on_t) >= 2, dur = diff(on_t); else, dur = []; end
end

function st = state_of(rate, basalHz, sniffHz)
    if rate < basalHz,      st = 1;   % basal
    elseif rate > sniffHz,  st = 2;   % sniff
    else,                   st = 0;   % unclassified
    end
end

function scatter_state(ax, x, y, st, cB, cS)
    v1 = st==1; v2 = st==2;
    if any(v1), scatter(ax, x(v1), y(v1), 10, cB, 'filled', 'MarkerFaceAlpha',0.35); end
    if any(v2), scatter(ax, x(v2), y(v2), 10, cS, 'filled', 'MarkerFaceAlpha',0.45); end
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:)); d = diff([false; a; false]);
    s = find(d==1); e = find(d==-1)-1; ep = [t(s), t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1)-ep2(end,2) <= mergeGap, ep2(end,2)=ep(i,2);
        else, ep2(end+1,:)=ep(i,:); end %#ok<AGROW>
    end
    ep = ep2; ep = ep(ep(:,2)-ep(:,1) >= minDur, :);
end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv')); S = {}; key = {};
    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]); S = S(o);
    end
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
