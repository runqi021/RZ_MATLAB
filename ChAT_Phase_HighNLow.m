% ChAT_Phase_HighNLow.m
% -----------------------------------------------------------------------
%  Piecewise-linear phase coherence polar for the ChAT folder.
%
%  Two reference points per cycle are used:
%    HIGH = breath PEAK    (insp_onset_idx, from breath_peak_data.mat)
%    LOW  = insp START    (insp_start_idx, from breath_insp_start_data.mat,
%                          = foot of the rising flank)
%
%  A continuous phase reference signal phi(t) is built so that
%    phi(peak) = 0 (mod 2pi)   and   phi(foot) = pi (mod 2pi)
%  with linear ramps between consecutive events.  cos(phi(t)) is then
%  fed into Chronux coherencyc against each ROI's spike train.  This
%  makes the polar phase angle anatomically meaningful regardless of
%  inspiration/expiration duty cycle:  th=0   -> spike at peak,
%                                       th=pi  -> spike at foot.
%
%  Plot style mirrors chat_breath_coherence_polar_260526.m:
%    - verified ChAT cells (chat_list) drawn in color, with magnitude
%      and phase 95% CI bars
%    - all other ROIs as hollow black background points
%    - dashed circle = confC at alpha_sig
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE =================================
rootPath = 'D:\Ventral_surface_summary\ChAT';

% verified ChAT cells (from chat_breath_coherence_polar_260526.m)
chat_list = { 'roi5_7x_x-1200y200z-30_3000f_23lp', 1; ...
              'roi3_8x_x-1070y730z0_3000f_15lp',   1; ...
              'roi1_4x_x-900y700z-15_6000f_13lp',  1 };
chat_colors = [0.85 0.10 0.10;    % red
               0.10 0.45 0.85;    % blue
               0.10 0.65 0.20];   % green

nDrop           = 30;
fallback_fps    = 30;
minSpikes       = 2;
TW              = 5;
alpha_sig       = 0.001;

f_breath_search = [0.2 4];   % Hz
fwhm_factor     = 0.6;
min_bw          = 0.05;
fmin            = 0.05;
fmax            = 15;

doSave          = false;     % set true to write png/pdf + .mat
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

%% ---- walk every recording with both peak + insp-start data ----
PP = init_coll(); labels = {}; confC = NaN;

allMat = dir(fullfile(rootPath, '**', 'ca_spike_data.mat'));
fprintf('Found %d recording(s) under %s\n', numel(allMat), rootPath);

for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    [~, recName] = fileparts(folderPath);
    try
        bp_file = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        ip_file = dir(fullfile(folderPath, '*breath_insp_start_data.mat'));
        if isempty(bp_file)
            fprintf('  skip (no peak data): %s\n', recName); continue;
        end
        if isempty(ip_file)
            fprintf('  skip (no insp_start data): %s\n', recName); continue;
        end

        fps = detect_session_fps(folderPath, fallback_fps);
        BP  = load(fullfile(bp_file(1).folder, bp_file(1).name));
        IP  = load(fullfile(ip_file(1).folder, ip_file(1).name));
        CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));

        % --- drop leading nDrop frames; shift event indices ---
        breath = double(BP.breath(:));
        breath = breath(nDrop+1:end);
        T_b    = numel(breath);

        peak_idx = double(BP.insp_onset_idx(:)) - nDrop;
        peak_idx(peak_idx < 1 | peak_idx > T_b) = [];

        foot_idx = double(IP.insp_start_idx(:)) - nDrop;
        foot_idx(foot_idx < 1 | foot_idx > T_b) = [];

        % --- truncate to common length with spike train ---
        nCa = numel(CA.roi_spikes(1).spike_train);
        T   = min(T_b, nCa);
        breath = breath(1:T);
        peak_idx(peak_idx > T) = [];
        foot_idx(foot_idx > T) = [];

        if numel(peak_idx) < 2 || numel(foot_idx) < 1
            fprintf('  skip (too few events): %s  (P=%d, F=%d)\n', recName, numel(peak_idx), numel(foot_idx));
            continue;
        end

        % --- piecewise-linear phase reference: peak=0, foot=pi ---
        phi = piecewise_phase(peak_idx, foot_idx, T);
        if all(isnan(phi))
            fprintf('  skip (phase build failed): %s\n', recName); continue;
        end
        ref = cos(phi);
        ref(isnan(ref)) = 0;

        % --- breath PSD -> coherence band ---
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        bw_ws = detrend(breath); bw_ws = bw_ws - mean(bw_ws);
        [Sb,fb] = mtspectrumc(bw_ws, pB); Sb=Sb(:); fb=fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        if ~any(m), fprintf('  skip (no PSD in band): %s\n', recName); continue; end
        [~,rl]=max(Sb(m)); ip=find(m,1)+rl-1; f_pk=fb(ip);
        h=Sb(ip)/2; lo=ip; while lo>1&&Sb(lo)>h, lo=lo-1; end
        hi=ip;            while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

        pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
        pc.fpass=band; pc.err=[2,alpha_sig];
        ref_m = ref - mean(ref);

        % --- coherence per ROI ---
        nROI = numel(CA.roi_spikes);
        nInc = 0;
        for rid = 1:nROI
            st = double(CA.roi_spikes(rid).spike_train(:));
            st = st(1:min(T,numel(st)));
            if numel(st) < T, st(end+1:T) = 0; end
            if sum(st) < minSpikes, continue; end
            nInc = nInc + 1;

            cidx = 0;     % verified-ChAT lookup
            for c = 1:size(chat_list,1)
                if contains(recName, chat_list{c,1}) && rid == chat_list{c,2}
                    cidx = c; break;
                end
            end

            [PP, confC] = add_coh(PP, ref_m, st - mean(st), pc, band, cidx, confC);
            labels{end+1} = sprintf('%s#%d', recName, rid); %#ok<SAGROW>
        end
        fprintf('[%d] %-50s band [%.2f %.2f] Hz  P=%d F=%d  %d/%d ROI\n', ...
                kk, recName, band(1), band(2), numel(peak_idx), numel(foot_idx), nInc, nROI);

    catch ME
        warning('  ERROR %s: %s', recName, ME.message);
    end
end

if isempty(PP.r), error('No ROIs collected.'); end
nROIs = numel(PP.r);
fprintf('\nTotal: %d ROIs, confC=%.3f, %d >= confC.\n', nROIs, confC, sum(PP.r >= confC));

%% ============================ POLAR =================================
fig = figure('Color','w','Name','ChAT phase polar (peak=0, foot=pi)', ...
             'Units','centimeters','Position',[2 2 15 14]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
ax = polaraxes(fig,'Position',[0.10 0.10 0.78 0.78]); hold(ax,'on');

thc = linspace(0,2*pi,360);
polarplot(ax, thc, repmat(confC,1,360), 'k--','LineWidth',1);

% background ROIs first (hollow black)
for k = 1:numel(PP.r)
    if PP.colorIdx(k) > 0 || isnan(PP.th(k)) || isnan(PP.r(k)), continue; end
    polarplot(ax, PP.th(k), PP.r(k), 'o', ...
              'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
              'MarkerSize', 5, 'LineWidth', 0.5);
end
% verified ChAT cells on top, with CI
for k = 1:numel(PP.r)
    ci = PP.colorIdx(k);
    if ci == 0 || isnan(PP.th(k)) || isnan(PP.r(k)), continue; end
    col = chat_colors(ci,:);
    polarplot(ax, [PP.th(k) PP.th(k)], [PP.rlo(k) PP.rhi(k)], '-', 'Color',col, 'LineWidth',1.5);
    if ~isnan(PP.dphi(k))
        arc = linspace(PP.th(k)-PP.dphi(k), PP.th(k)+PP.dphi(k), 30);
        polarplot(ax, arc, PP.r(k)*ones(size(arc)), '-', 'Color',col, 'LineWidth',1.5);
    end
    polarplot(ax, PP.th(k), PP.r(k), 'o', ...
              'MarkerFaceColor',col, 'MarkerEdgeColor','k', 'MarkerSize', 9);
end

ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
ax.RAxisLocation=180; ax.FontSize=8;
title(ax, sprintf('ChAT phase polar  (peak=0, foot=\\pi)   N=%d ROI   confC=%.2f   \\alpha=%.3f', ...
                  nROIs, confC, alpha_sig), 'Interpreter','tex');

% direction-of-cycle text guides
text(ax, 0,    1.05, 'peak (high)',  'HorizontalAlignment','center','FontSize',9);
text(ax, pi,   1.05, 'foot (low)',   'HorizontalAlignment','center','FontSize',9);

if doSave
    outDir = fullfile(rootPath, 'phase_HighNLow');
    if ~isfolder(outDir), mkdir(outDir); end
    exportgraphics(fig, fullfile(outDir,'ChAT_Phase_HighNLow_polar.png'), 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'ChAT_Phase_HighNLow_polar.pdf'), 'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'ChAT_Phase_HighNLow_data.mat'), ...
         'PP','labels','confC','chat_list','chat_colors','TW','alpha_sig');
    fprintf('Saved figure + .mat to %s\n', outDir);
end

%% ========================= LOCAL FUNCTIONS ==========================
function phi = piecewise_phase(peak_idx, foot_idx, T)
% Continuous phase reference, peaks at 0 (mod 2pi), feet at pi (mod 2pi),
% linear ramps between consecutive events. NaN outside the event range.
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
% strict alternation: drop adjacent same-type duplicates (keep the first)
keep = true(size(events));
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
% phase per event: peak=0, foot=pi, then advance by pi each event
phases = nan(size(events));
phi_cur = (1 - types(1)) * pi;     % type=1 (peak) -> 0; type=0 (foot) -> pi
for i = 1:numel(events)
    phases(i) = phi_cur;
    phi_cur   = phi_cur + pi;
end
% linear interpolation between consecutive events
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end

function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[],'colorIdx',[]);
end

function [C, confC] = add_coh(C, x, y, pc, band, cidx, confC)
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f = f(:); mb = f>=band(1) & f<=band(2); if ~any(mb), mb = true(size(f)); end
    C.th(end+1,1)   = angle(mean(exp(1i*(-phi(mb)))));
    C.r(end+1,1)    = mean(Cxy(mb));
    C.rlo(end+1,1)  = max(0, mean(Cerr(1,mb)));
    C.rhi(end+1,1)  = min(1, mean(Cerr(2,mb)));
    C.dphi(end+1,1) = 1.96*mean(phistd(mb));
    C.colorIdx(end+1,1) = cidx;
end
