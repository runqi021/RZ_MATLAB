function whisk_LRdiff_cumint_diag_RZ()
% whisk_LRdiff_cumint_diag_RZ  Per-session diagnostic, ALL sessions.
% Same style as whisk_detect_diag_RZ; ABSOLUTE 5-deg-on-L whisk-epoch mask.
%
% Two stacked, x-linked panels per session:
%   (1) L/R BP whisker angle + L envelope (+/-) + ABSOLUTE threshold + shaded
%       whisk epochs.
%   (2) CUMULATIVE integral of (L - R):  cumtrapz(t, xL-xR)  (deg*s) -- a
%       running integral from the start of the recording (the accumulated
%       bilateral imbalance), same epoch shading, zero line.
%
% RESET_PER_EPOCH = true restarts the cumulative integral at each whisk epoch.
%
% One figure per session (loops every *_whisk_n*.csv).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW      = 400;
BP        = [5 50];     % whisk bandpass (Hz)
USE_BP    = true;       % true: BP signal; false: raw demeaned angle

% ABSOLUTE 5-deg-on-L whisk-epoch mask
ENV_THR    = 5;         % deg
ENV_SMOOTH = 0.2;       % s
MIN_DUR    = 0.2;       % s
MERGE_GAP  = 0.2;       % s

RESET_PER_EPOCH = false; % true: restart cumulative integral at each epoch

PROT_PROM  = 1;         % protraction trough prominence (x std of BP)
SHOW_PEAKS = false;     % overlay protraction-onset markers
ZOOM       = [];        % [t0 t1] s; [] = full trace

EXCLUDE    = "";        % "" = all 17 sessions
SAVE_FIGS  = false;     % also save a PNG per session
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LRdiff_cumint_diag');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
CUMABS_BP = [1 100];                                  % Hz, BP for cumulative |L-R|
[bca,aca] = butter(3, CUMABS_BP/(fpsW/2),'bandpass');

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

    LaRaw = fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest');
    RaRaw = fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest');
    xLbp = filtfilt(bw,aw, LaRaw);
    xRbp = filtfilt(bw,aw, RaRaw);
    if USE_BP, xL = xLbp; xR = xRbp; else, xL = LaRaw; xR = RaRaw; end

    % ABSOLUTE 5-deg-on-L mask (from BP envelope)
    env = movmean(abs(hilbert(xLbp)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);

    d = xL - xR;

    % CUMULATIVE SUM of (L-R) and of |L-R| -- plain cumsum, no fps factor (deg)
    if RESET_PER_EPOCH
        cumI = zeros(size(d)); cumA = zeros(size(d));
        for q=1:size(ep,1)
            idx = find(t>=ep(q,1) & t<=ep(q,2));
            if numel(idx)>1
                cumI(idx) = cumsum(d(idx));
                cumA(idx) = cumsum(abs(d(idx)));
            end
        end
    else
        cumI = cumsum(d);
        cumA = cumsum(abs(d));
    end

    % normalize cumulative |L-R| by band-passing 1-100 Hz (removes the ramp)
    cumA_bp = filtfilt(bca, aca, cumA);

    [~,iL] = findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(0.02*fpsW));
    [~,iR] = findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(0.02*fpsW));

    fprintf('  %s n%d: dur=%.0fs | %d epochs (%.0f%% of time) | cumsum(L-R) end=%.3g deg | cumsum|L-R| end=%.4g deg\n', ...
        animal, kRun, t(end), size(ep,1), 100*sum(ep(:,2)-ep(:,1))/max(t(end),eps), cumI(end), cumA(end));

    zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

    fig = figure('Color','w','Position',[60 80 1240 820]);

    % -------- panel 1: BP L/R + envelope + epochs --------
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    yl1 = [min([xL;xR;-env]) max([xL;xR;env])];
    shade_epochs(ax1, ep, yl1);
    plot(ax1, t, xL, '-', 'Color',[0 0.5 0]);
    plot(ax1, t, xR, '-', 'Color',[0 0.4 0.85]);
    plot(ax1, t, env, 'k-', 'LineWidth',1.4); plot(ax1, t, -env, 'k-', 'LineWidth',1.4);
    yline(ax1, ENV_THR, 'r--'); yline(ax1, -ENV_THR, 'r--');
    if SHOW_PEAKS
        plot(ax1, t(iL), xL(iL), 'v', 'Color',[0 0.5 0],   'MarkerFaceColor',[0 0.5 0],   'MarkerSize',5);
        plot(ax1, t(iR), xR(iR), 'v', 'Color',[0 0.4 0.85],'MarkerFaceColor',[0 0.4 0.85],'MarkerSize',5);
    end
    ylim(ax1, yl1);
    ylabel(ax1,'whisk angle (deg, BP)');
    ttl = sprintf('%s n%d  (green = whisk epoch >%g deg on L;  red-- = thr)', animal,kRun,ENV_THR);
    if SHOW_PEAKS, ttl=[ttl '   v = protraction onset']; end
    title(ax1, ttl, 'Interpreter','none');

    if RESET_PER_EPOCH, rs = ' (reset per epoch)'; else, rs = ' (running from t=0)'; end

    % -------- panel 2: cumulative sum of (L-R) --------
    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    yl2 = [min(cumI) max(cumI)]; if diff(yl2)==0, yl2 = yl2 + [-1 1]; end
    shade_epochs(ax2, ep, yl2);
    plot(ax2, t, cumI, '-', 'Color',[0.55 0.15 0.55], 'LineWidth',1.2);
    yline(ax2, 0, 'k-');
    ylim(ax2, yl2);
    ylabel(ax2,'cumsum(L-R)   (deg)');
    title(ax2, ['cumulative sum of signed (L-R):  cumsum(L-R)' rs], 'Interpreter','tex');

    % -------- panel 3: cumulative |L-R|, band-passed 1-100 Hz --------
    ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
    yl3 = [min(cumA_bp) max(cumA_bp)]; if diff(yl3)==0, yl3 = yl3 + [-1 1]; end
    shade_epochs(ax3, ep, yl3);
    plot(ax3, t, cumA_bp, '-', 'Color',[0.85 0.35 0.10], 'LineWidth',1.0);
    yline(ax3, 0, 'k-');
    ylim(ax3, yl3);
    xlabel(ax3,'time (s)'); ylabel(ax3,'BP[1-100] cumsum(|L-R|)');
    title(ax3, sprintf('cumulative |L-R| band-passed %g-%g Hz%s', CUMABS_BP(1), CUMABS_BP(2), rs), ...
          'Interpreter','none');

    linkaxes([ax1 ax2 ax3],'x'); xlim(ax1, zoomwin);

    nDone = nDone + 1;
    if SAVE_FIGS
        exportgraphics(fig, fullfile(outDir, sprintf('LRdiff_cumint_%s_n%d.png',animal,kRun)), ...
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
