function whisk_LRdiff_bp_vs_raw_diag_RZ()
% whisk_LRdiff_bp_vs_raw_diag_RZ  Per-session diagnostic, ALL sessions.
% Same style as whisk_detect_diag_RZ; ABSOLUTE 5-deg-on-L whisk-epoch mask.
%
% Three stacked, x-linked panels per session:
%   (1) L/R BP whisker angle + L envelope (+/-) + ABSOLUTE threshold + shaded
%       whisk epochs.
%   (2) L - R of the BP 5-50 Hz signal.
%   (3) |L| + |R| of the BP 5-50 Hz signal.
%
% One figure per session (loops every *_whisk_n*.csv).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW      = 400;
BP        = [5 50];     % whisk bandpass (Hz)

% ABSOLUTE 5-deg-on-L whisk-epoch mask
ENV_THR    = 5;         % deg
ENV_SMOOTH = 0.05;      % s
MIN_DUR    = 1;         % s
MERGE_GAP  = 0.2;       % s

ASYM_WIN_S = 0.2;       % s, smoothing window for the moving asymmetry index
                        % (per-whisk ~0.1 s, per-bout ~0.5 s; just smoothing)

ZOOM       = [];        % [t0 t1] s; [] = full trace
EXCLUDE    = "5840027"; % skip this animal's runs ("" = all)
SAVE_FIGS  = false;     % also save a PNG per session
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LRdiff_bp_vs_raw_diag');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
winA = max(1, round(ASYM_WIN_S*fpsW));

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);
fprintf('%d sessions found\n', numel(sess));

nDone = 0;
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;

    % RAW demeaned (filled, NOT band-passed)
    LaRaw = fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest');
    RaRaw = fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest');

    % BP 5-50 Hz
    xL = filtfilt(bw,aw, LaRaw);
    xR = filtfilt(bw,aw, RaRaw);

    % ABSOLUTE 5-deg-on-L mask
    env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);

    dBP  = xL - xR;                 % L - R (differential, signed)
    sAbs = abs(xL) + abs(xR);       % |L| + |R| (total; the bounded denominator)

    % MOVING normalized asymmetry index: windowed |L-R| / windowed (|L|+|R|).
    % Bounded [0,1] because |L-R| <= |L|+|R| at every sample. Window = smoothing.
    AI = movmean(abs(dBP), winA) ./ max(movmean(sAbs, winA), eps);
    AI(~isfinite(AI)) = 0;

    fprintf('  %s n%d: dur=%.0fs | %d epochs (%.0f%% of time) | AI(win=%.0fms) mean=%.3f\n', ...
        animal, kRun, t(end), size(ep,1), 100*sum(ep(:,2)-ep(:,1))/max(t(end),eps), 1000*ASYM_WIN_S, mean(AI));

    zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

    fig = figure('Color','w','Position',[60 60 1240 920]);

    % -------- panel 1: BP L/R + envelope + epochs --------
    ax1 = subplot(4,1,1); hold(ax1,'on'); grid(ax1,'on');
    yl1 = [min([xL;xR;-env]) max([xL;xR;env])];
    shade_epochs(ax1, ep, yl1);
    hL=plot(ax1, t, xL, '-', 'Color',[0 0.5 0]);
    hR=plot(ax1, t, xR, '-', 'Color',[0 0.4 0.85]);
    he=plot(ax1, t, env, 'k-', 'LineWidth',1.4); plot(ax1, t, -env, 'k-', 'LineWidth',1.4);
    yline(ax1, ENV_THR, 'r--'); yline(ax1, -ENV_THR, 'r--');
    ylim(ax1, yl1);
    ylabel(ax1,'whisk angle (deg, BP)');
    title(ax1, sprintf('%s n%d  (green = whisk epoch >%g deg on L;  red-- = thr)', animal,kRun,ENV_THR), ...
          'Interpreter','none');

    % -------- panel 2: L-R of BP signal --------
    ax2 = subplot(4,1,2); hold(ax2,'on'); grid(ax2,'on');
    yl2 = nzlim(dBP);
    shade_epochs(ax2, ep, yl2);
    plot(ax2, t, dBP, '-', 'Color',[0.55 0.15 0.55], 'LineWidth',1.0);
    yline(ax2, 0, 'k-'); ylim(ax2, yl2);
    ylabel(ax2,'L - R  (deg, BP 5-50)');
    title(ax2, 'L - R of band-passed (5-50 Hz) signal', 'Interpreter','none');

    % -------- panel 3: |L| + |R| of BP signal --------
    ax3 = subplot(4,1,3); hold(ax3,'on'); grid(ax3,'on');
    yl3 = nzlim(sAbs);
    shade_epochs(ax3, ep, yl3);
    plot(ax3, t, sAbs, '-', 'Color',[0.10 0.55 0.45], 'LineWidth',1.0);
    yline(ax3, 0, 'k-'); ylim(ax3, yl3);
    ylabel(ax3,'|L| + |R|  (deg, BP 5-50)');
    title(ax3, 'sum of |L| and |R|  (band-passed 5-50 Hz)', 'Interpreter','none');

    % -------- panel 4: MOVING normalized asymmetry index --------
    ax4 = subplot(4,1,4); hold(ax4,'on'); grid(ax4,'on');
    shade_epochs(ax4, ep, [0 1]);
    plot(ax4, t, AI, '-', 'Color',[0.85 0.35 0.10], 'LineWidth',1.0);
    ylim(ax4, [0 1]);
    xlabel(ax4,'time (s)'); ylabel(ax4,'movmean|L-R| / movmean(|L|+|R|)');
    title(ax4, sprintf('moving asymmetry index (win=%.0f ms = smoothing only):  0 = symmetric, 1 = antagonistic', ...
          1000*ASYM_WIN_S), 'Interpreter','none');

    linkaxes([ax1 ax2 ax3 ax4],'x'); xlim(ax1, zoomwin);

    nDone = nDone + 1;
    if SAVE_FIGS
        exportgraphics(fig, fullfile(outDir, sprintf('LRdiff_bp_vs_raw_%s_n%d.png',animal,kRun)), ...
                       'Resolution',150,'BackgroundColor','white');
    end
end
fprintf('Plotted %d sessions.\n', nDone);
end

% ================= helpers =================
function shade_epochs(ax, ep, yl)
    for q=1:size(ep,1)
        patch(ax, ep(q,[1 2 2 1]), yl([1 1 2 2]), [0.3 0.75 0.3], 'FaceAlpha',0.12,'EdgeColor','none');
    end
end

function yl = nzlim(x)
    yl = [min(x) max(x)]; if diff(yl)==0, yl = yl + [-1 1]; end
end

function ep = detect_abs(env, t, thrDeg, minDur, mergeGap)
% ABSOLUTE-threshold epoch detection (env > thrDeg), merge + min-dur.
    a = env(:) > thrDeg;
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

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv')); S = {}; key = {};
    for i=1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o]=sortrows([an(:) rn(:)]); S=S(o);
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
