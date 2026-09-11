% breath_phase_cosine_check_260528.m
% -----------------------------------------------------------------------
%  Visualize the piecewise-linear phase reference used by ChAT_Phase_HighNLow.
%  For every recording under <rootPath> that has both peak + insp-start
%  data, plot 3 panels stacked on a shared time axis:
%
%    R1  breath waveform (z-scored) + ▲ peak + ▼ insp-start markers
%    R2  phi(t)  -- piecewise-linear phase (the "triangle wave")
%        peaks land at 0/2pi/4pi..., feet at pi/3pi/5pi...
%    R3  cos(phi(t)) overlaid with normalized breath (the smooth
%        sinusoid-shaped wave that gets fed to Chronux coherencyc)
%
%  Window length is configurable so the structure is actually visible.
%  Nothing saved by default; toggle doSave.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =================================
rootPath  = 'D:\Ventral_surface_summary\ChAT';
windowSec = 30;            % seconds shown per figure (use Inf for full)
nDrop     = 30;            % match coherence script
fallback_fps = 30;
doSave    = false;
% =====================================================================

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

%% ---- walk ----
allMat = dir(fullfile(rootPath, '**', 'ca_spike_data.mat'));
fprintf('Found %d recording(s) under %s\n', numel(allMat), rootPath);

for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    [~, recName] = fileparts(folderPath);

    bp = dir(fullfile(folderPath, '*breath_peak_data.mat'));
    ip = dir(fullfile(folderPath, '*breath_insp_start_data.mat'));
    if isempty(bp) || isempty(ip)
        fprintf('  skip (missing peak or insp-start): %s\n', recName); continue;
    end
    try
        fps = detect_session_fps(folderPath, fallback_fps);
        BP  = load(fullfile(bp(1).folder, bp(1).name));
        IP  = load(fullfile(ip(1).folder, ip(1).name));

        breath = double(BP.breath(:));
        breath = breath(nDrop+1:end);
        T = numel(breath);

        peak_idx = double(BP.insp_onset_idx(:)) - nDrop;
        peak_idx(peak_idx<1 | peak_idx>T) = [];
        foot_idx = double(IP.insp_start_idx(:)) - nDrop;
        foot_idx(foot_idx<1 | foot_idx>T) = [];

        if numel(peak_idx) < 2 || numel(foot_idx) < 1
            fprintf('  skip (too few events): %s\n', recName); continue;
        end

        phi = piecewise_phase_unwrap(peak_idx, foot_idx, T);
        ref = cos(phi);
        t   = (0:T-1)' / fps;

        % time window
        t_end = min(t(end), windowSec);
        m     = t <= t_end;

        % normalize breath to [-1, 1] for overlay
        b = breath(m);
        b_norm = 2 * (b - min(b)) / (max(b) - min(b) + eps) - 1;

        %% --- figure ---
        fig = figure('Color','w','Name',recName, ...
                     'Units','centimeters','Position',[1 1 26 16]);
        tl  = tiledlayout(fig, 3, 1, 'TileSpacing','compact','Padding','compact');
        title(tl, sprintf('%s  |  fps=%.2f Hz  |  P=%d F=%d  |  first %.0f s', ...
              recName, fps, numel(peak_idx), numel(foot_idx), t_end), ...
              'Interpreter','none','FontWeight','bold');

        % --- R1: breath + peak/foot markers ---
        ax1 = nexttile(tl); hold(ax1,'on');
        plot(ax1, t(m), breath(m), 'k-', 'LineWidth', 0.8);
        pk = peak_idx(peak_idx <= sum(m));
        ft = foot_idx(foot_idx <= sum(m));
        plot(ax1, t(pk), breath(pk), '^', 'MarkerFaceColor',[0.95 0.5 0], ...
             'MarkerEdgeColor','k','MarkerSize',7);
        plot(ax1, t(ft), breath(ft), 'v', 'MarkerFaceColor',[0 0.55 1], ...
             'MarkerEdgeColor','k','MarkerSize',7);
        ylabel(ax1,'breath (z)');
        title(ax1, 'breath waveform + peak (▲) + insp start (▼)');
        grid(ax1,'on'); box(ax1,'on'); xlim(ax1,[0 t_end]);
        set(ax1,'XTickLabel',[]);

        % --- R2: phi(t) WRAPPED to [0, 2pi) -- continuous, with a vertical
        %     connector at each wrap so the line "touches" 2pi and restarts at 0.
        %     At every wrap we linear-interpolate the time at which the UNWRAPPED
        %     phase = 2pi, then insert two points at that time: (t_cross, 2pi)
        %     and (t_cross, 0). MATLAB then draws a clean vertical between them.
        phi_w   = mod(phi, 2*pi);
        d       = diff(phi_w);
        wraps   = find(d < -pi);           % drops > pi = wrap (just crossed 2pi)
        phi_plot = phi_w; t_plot = t;
        for j = numel(wraps):-1:1
            ix = wraps(j);
            t1 = t(ix); t2 = t(ix+1);
            p1 = phi_w(ix);
            p2_u = phi_w(ix+1) + 2*pi;     % unwrapped value just past the wrap
            alpha   = (2*pi - p1) / (p2_u - p1);
            t_cross = t1 + alpha*(t2 - t1);
            phi_plot = [phi_plot(1:ix); 2*pi; 0; phi_plot(ix+1:end)]; %#ok<AGROW>
            t_plot   = [t_plot(1:ix);   t_cross; t_cross; t_plot(ix+1:end)]; %#ok<AGROW>
        end
        m_plot = t_plot <= t_end;

        ax2 = nexttile(tl); hold(ax2,'on');
        plot(ax2, t_plot(m_plot), phi_plot(m_plot), 'Color',[0.4 0.4 0.4], 'LineWidth', 0.9);
        yline(ax2, 0,    ':', 'Color',[0.7 0.7 0.7]);
        yline(ax2, pi,   ':', 'Color',[0.7 0.7 0.7]);
        yline(ax2, 2*pi, ':', 'Color',[0.7 0.7 0.7]);
        plot(ax2, t(pk), mod(phi(pk),2*pi), '^', 'MarkerFaceColor',[0.95 0.5 0], 'MarkerEdgeColor','k','MarkerSize',6);
        plot(ax2, t(ft), mod(phi(ft),2*pi), 'v', 'MarkerFaceColor',[0 0.55 1], 'MarkerEdgeColor','k','MarkerSize',6);
        ylim(ax2,[0 2*pi]);
        set(ax2,'YTick',[0 pi/2 pi 3*pi/2 2*pi], ...
                'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'});
        ylabel(ax2,'\phi(t) (rad)');
        title(ax2, '\phi(t) wrapped to [0, 2\pi)  --  peak \rightarrow 0, foot \rightarrow \pi');
        grid(ax2,'on'); box(ax2,'on'); xlim(ax2,[0 t_end]);
        set(ax2,'XTickLabel',[]);

        % --- R3: cos(phi) overlaid with normalized breath ---
        ax3 = nexttile(tl); hold(ax3,'on');
        plot(ax3, t(m), b_norm, '-', 'Color',[0.55 0.55 0.55], 'LineWidth', 0.8);
        plot(ax3, t(m), ref(m), 'r-', 'LineWidth', 1.1);
        plot(ax3, t(pk), ref(pk), '^', 'MarkerFaceColor',[0.95 0.5 0], 'MarkerEdgeColor','k','MarkerSize',6);
        plot(ax3, t(ft), ref(ft), 'v', 'MarkerFaceColor',[0 0.55 1], 'MarkerEdgeColor','k','MarkerSize',6);
        yline(ax3, 0, 'k:');
        ylim(ax3,[-1.15 1.15]); ylabel(ax3,'amplitude');
        xlabel(ax3, 'time (s)');
        legend(ax3, {'breath (norm)','cos(\phi)'}, 'Location','northeast', 'FontSize',8);
        title(ax3, 'cos(\phi(t)) -- what is fed to Chronux  (overlaid with normalized breath)');
        grid(ax3,'on'); box(ax3,'on'); xlim(ax3,[0 t_end]);

        if doSave
            outBase = fullfile(folderPath, sprintf('breath_phase_check_%.0fs', t_end));
            exportgraphics(fig, [outBase '.png'], 'Resolution',150, 'BackgroundColor','white');
            exportgraphics(fig, [outBase '.pdf'], 'ContentType','vector','BackgroundColor','white');
            fprintf('  [%d] saved: %s.png/.pdf\n', kk, outBase);
        else
            fprintf('  [%d] %s (P=%d, F=%d)\n', kk, recName, numel(peak_idx), numel(foot_idx));
        end
    catch ME
        warning('  ERROR %s: %s', recName, ME.message);
    end
end
fprintf('Done.\n');

%% ========================= LOCAL FUNCTIONS ==========================
function phi = piecewise_phase_unwrap(peak_idx, foot_idx, T)
% Continuous (unwrapped) phase, peaks at 0/2pi/4pi..., feet at pi/3pi/5pi...
% Linear ramps between events. NaN outside the event range.
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
% strict alternation: drop adjacent same-type duplicates
keep = true(size(events));
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
% assign monotonically increasing phase: each event advances by pi
phases  = nan(size(events));
phi_cur = (1 - types(1)) * pi;     % type=1 (peak) -> 0; type=0 (foot) -> pi
for i = 1:numel(events)
    phases(i) = phi_cur;
    phi_cur   = phi_cur + pi;
end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
