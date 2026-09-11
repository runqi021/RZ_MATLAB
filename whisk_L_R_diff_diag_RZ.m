function whisk_L_R_diff_diag_RZ()
% whisk_L_R_diff_diag_RZ  Per-session diagnostic, ALL sessions.
% Same style as whisk_detect_diag_RZ; ABSOLUTE 5-deg-on-L whisk-epoch mask.
%
% Two stacked, x-linked panels per session:
%   (1) L and R BP whisker angle overlaid + L envelope (+/-) + threshold
%       + shaded whisk epochs.
%   (2) diff(L) and diff(R) overlaid (sample-to-sample difference of the BP
%       angle, deg/sample) + shaded whisk epochs.
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

ZOOM       = [];        % [t0 t1] s; [] = full trace
EXCLUDE    = "";        % "" = all 17 sessions; "5840027" = drop that animal
SAVE_FIGS  = false;     % also save a PNG per session
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_L_R_diff_diag');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

colL = [0 0.5 0];        % L green
colR = [0 0.4 0.85];     % R blue

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

    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest'));

    % ABSOLUTE 5-deg-on-L mask
    env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);

    dL = [0; diff(xL)];      % diff(L), aligned to t
    dR = [0; diff(xR)];      % diff(R)

    fprintf('  %s n%d: dur=%.0fs | %d epochs (%.0f%% of time)\n', ...
        animal, kRun, t(end), size(ep,1), 100*sum(ep(:,2)-ep(:,1))/max(t(end),eps));

    zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

    fig = figure('Color','w','Position',[60 70 1240 880]);

    % -------- panel 1: L and R overlaid --------
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    yl1 = [min([xL;xR;-env]) max([xL;xR;env])];
    shade_epochs(ax1, ep, yl1);
    hL = plot(ax1, t, xL, '-', 'Color',colL);
    hR = plot(ax1, t, xR, '-', 'Color',colR);
    plot(ax1, t, env, 'k-', 'LineWidth',1.4); plot(ax1, t, -env, 'k-', 'LineWidth',1.4);
    yline(ax1, ENV_THR, 'r--'); yline(ax1, -ENV_THR, 'r--');
    ylim(ax1, yl1);
    ylabel(ax1,'whisk angle (deg, BP)');
    legend(ax1, [hL hR], {'L','R'}, 'Location','northeastoutside');
    title(ax1, sprintf('%s n%d  -  L & R BP angle  (green = whisk epoch >%g deg on L)', animal,kRun,ENV_THR), ...
          'Interpreter','none');

    % -------- panel 2: diff(L) and diff(R) overlaid --------
    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    yl2 = nzlim([dL;dR]);
    shade_epochs(ax2, ep, yl2);
    gL = plot(ax2, t, dL, '-', 'Color',colL);
    gR = plot(ax2, t, dR, '-', 'Color',colR);
    yline(ax2, 0, 'k-'); ylim(ax2, yl2);
    ylabel(ax2,'diff (deg/sample)');
    legend(ax2, [gL gR], {'diff(L)','diff(R)'}, 'Location','northeastoutside');
    title(ax2, 'diff(L) & diff(R)  (sample-to-sample difference of BP angle)', 'Interpreter','none');

    % -------- panel 3: L and diff(L) overlaid (twin y-axis) --------
    ax3 = subplot(3,1,3); colD = [0.85 0.35 0.10];
    yyaxis(ax3,'left'); hold(ax3,'on'); grid(ax3,'on');
    ylL = nzlim(xL); shade_epochs(ax3, ep, ylL);
    pL = plot(ax3, t, xL, '-', 'Color',colL);
    ylim(ax3, ylL); ylabel(ax3,'L (deg, BP)'); ax3.YAxis(1).Color = colL;
    yyaxis(ax3,'right');
    pD = plot(ax3, t, dL, '-', 'Color',colD);
    yline(ax3, 0, 'k-'); ylabel(ax3,'diff(L) (deg/sample)'); ax3.YAxis(2).Color = colD;
    xlabel(ax3,'time (s)');
    legend(ax3, [pL pD], {'L','diff(L)'}, 'Location','northeastoutside');
    title(ax3, 'L and diff(L) overlay  (diff = derivative, leads L by ~90 deg)', 'Interpreter','none');

    % NOTE: do NOT linkaxes a yyaxis panel (it collapses the x-limits). Set
    % x-limits explicitly; link only the two single-axis panels.
    linkaxes([ax1 ax2],'x');
    xlim(ax1, zoomwin); xlim(ax2, zoomwin);
    yyaxis(ax3,'left');  xlim(ax3, zoomwin);
    yyaxis(ax3,'right'); xlim(ax3, zoomwin);

    nDone = nDone + 1;
    if SAVE_FIGS
        exportgraphics(fig, fullfile(outDir, sprintf('L_R_diff_%s_n%d.png',animal,kRun)), ...
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
