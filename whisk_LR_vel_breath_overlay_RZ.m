function whisk_LR_vel_breath_overlay_RZ()
% whisk_LR_vel_breath_overlay_RZ  Per-session, ALL sessions.
% Overlay whisker angular VELOCITY of both sides with the breathing trace.
%
%   diff(L)/dt = [0; diff(xL)] * fps   (deg/s)   -- left y-axis
%   diff(R)/dt = [0; diff(xR)] * fps   (deg/s)   -- left y-axis
%   breathing  (canonical thermal nostril, inhale-up)  -- right y-axis
%
% xL/xR = BP 5-50 Hz whisker angle. Whisk epochs (ABSOLUTE 5-deg-on-L) shaded.
% One figure per session; sessions without a breath file are skipped.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW      = 400;
BP        = [5 50];     % whisk bandpass (Hz)
BREATH_BP = [1 20];     % breath bandpass (Hz); [] = use canonical breath as-is

% ABSOLUTE 5-deg-on-L whisk-epoch mask
ENV_THR    = 5;         % deg
ENV_SMOOTH = 0.05;      % s
MIN_DUR    = 1;         % s
MERGE_GAP  = 0.2;       % s

ZOOM       = [];        % [t0 t1] s; [] = full trace
EXCLUDE    = "";        % "" = all sessions (those with breath)
SAVE_FIGS  = false;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_vel_breath_overlay');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
useBrBP = ~isempty(BREATH_BP);
if useBrBP, [bb,ab] = butter(3, BREATH_BP/(fpsW/2),'bandpass'); end

colL = [0 0.5 0]; colR = [0 0.4 0.85]; colB = [0.75 0.1 0.6];

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);
fprintf('%d sessions found\n', numel(sess));

nDone = 0;
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    % ---------------- whisker ----------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW; N = numel(t);
    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest'));

    env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);

    vL = [0; diff(xL)] * fpsW;   % diff(L)/dt  (deg/s)
    vR = [0; diff(xR)] * fpsW;   % diff(R)/dt

    % ---------------- breath ----------------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun)), dataRoot);
        if ~isfile(Pn.breath), fprintf('  %s n%d: no breath file, skipped\n', animal, kRun); continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps); tBr = (0:numel(br)-1)'/fb;
        brw = interp1(tBr, br, t, 'linear', NaN);
        if nnz(isfinite(brw)) < 50, fprintf('  %s n%d: too few breath samples, skipped\n', animal, kRun); continue; end
        brw = fillmissing(brw,'linear','EndValues','nearest');
        if useBrBP, brw = filtfilt(bb,ab, brw); end
    catch ME
        warning('breath load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    fprintf('  %s n%d: dur=%.0fs | %d epochs | vL std=%.0f deg/s\n', ...
        animal, kRun, t(end), size(ep,1), std(vL));

    zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

    fig = figure('Color','w','Position',[60 230 1280 440]);
    ax = axes(fig);

    yyaxis(ax,'left'); hold(ax,'on'); grid(ax,'on');
    ylv = nzlim([vL;vR]);
    for q=1:size(ep,1)
        patch(ax, ep(q,[1 2 2 1]), ylv([1 1 2 2]), [0.3 0.75 0.3], 'FaceAlpha',0.10,'EdgeColor','none');
    end
    pL = plot(ax, t, vL, '-', 'Color',colL, 'LineWidth',0.5);
    pR = plot(ax, t, vR, '-', 'Color',colR, 'LineWidth',0.5);
    ylim(ax, ylv); ylabel(ax,'whisker velocity (deg/s)'); ax.YAxis(1).Color = [0.2 0.2 0.2];

    yyaxis(ax,'right');
    pB = plot(ax, t, brw, '-', 'Color',colB, 'LineWidth',1.6);
    ylabel(ax,'breath (a.u.)'); ax.YAxis(2).Color = colB;

    xlabel(ax,'time (s)');
    legend(ax, [pL pR pB], {'diff(L)/dt','diff(R)/dt','breath'}, 'Location','northeastoutside');
    title(ax, sprintf('%s n%d  -  whisker velocity (L/R) + breathing  (green = whisk epoch >%g deg)', ...
          animal,kRun,ENV_THR), 'Interpreter','none');

    yyaxis(ax,'left');  xlim(ax, zoomwin);
    yyaxis(ax,'right'); xlim(ax, zoomwin);

    nDone = nDone + 1;
    if SAVE_FIGS
        exportgraphics(fig, fullfile(outDir, sprintf('vel_breath_%s_n%d.png',animal,kRun)), ...
                       'Resolution',150,'BackgroundColor','white');
    end
end
fprintf('Plotted %d sessions.\n', nDone);
end

% ================= helpers =================
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
