% whisk_LR_breath_alltraces_RZ  (script)
% Full-recording three-trace overlay for EVERY session (no epoch gating):
%   left  y-axis : L/R fast-whisk BP angle (deg)
%   right y-axis : band-passed breathing (inhale up)
% ONE FIGURE PER SESSION.
% Use this when whisk-epoch detection rejects a session — it draws the raw
% traces regardless of whisking amplitude.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BP        = [6 30];      % fast-whisking bandpass (Hz)
BR_BP     = [2 15];      % breath bandpass (Hz)
TLIM      = [];          % [] = full recording; [t0 t1] (s) to crop every panel
fpsW      = 400;         % whisk camera fps
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% ---- enumerate sessions (digit-named animal dirs, cam1_* runs with a whisk csv) ----
sess = {};
ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess), 'no sessions with whisk csv in %s', whiskDir);
nS = numel(sess);
fprintf('found %d sessions\n', nS);

% ---- one figure per session ----
cL = [0 0.5 0]; cR = [0 0.4 0.85]; cB = [0.10 0.10 0.10];

for e = 1:nS
    animal = sess{e}{1}; kk = sess{e}{2};
    % whisk angles (lik<0.6 -> linear interp, same as the other scripts)
    M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
    t  = (0:numel(xL)-1)'/fpsW;

    % breathing (optional; panel still drawn if missing)
    tB = []; brf = [];
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);
        if isfile(Pn.breath)
            Bs = load(Pn.breath); br = Bs.breath(:); fpsB = double(Bs.fps);
            tB = (0:numel(br)-1)'/fpsB;
            [bb,ab] = butter(2, BR_BP/(fpsB/2), 'bandpass');
            brf = filtfilt(bb,ab, fillmissing(br,'linear'));
        end
    catch
    end

    figure('Color','w','Position',[120 200 1100 440]);
    ax = gca; hold(ax,'on');
    yyaxis(ax,'left');
    p1 = plot(ax, t, xL, '-', 'Color',cL, 'LineWidth',0.9);
    p2 = plot(ax, t, xR, '-', 'Color',cR, 'LineWidth',0.9);
    ylabel(ax,'whisk angle (deg, BP)'); set(ax,'YColor','k');
    yyaxis(ax,'right'); p3 = [];
    if ~isempty(brf)
        p3 = plot(ax, tB, brf, '-', 'Color',cB, 'LineWidth',1.1);
    end
    ylabel(ax,'breathing (BP, inhale up)'); set(ax,'YColor',cB);
    grid(ax,'on'); box(ax,'on'); xlabel(ax,'time (s)');
    if ~isempty(TLIM), xlim(ax,TLIM); else, xlim(ax,[0 t(end)]); end
    hh = [p1 p2]; ll = {'L whisk','R whisk'};
    if ~isempty(p3), hh=[hh p3]; ll=[ll {'breath'}]; end
    legend(hh, ll, 'Orientation','horizontal', 'Location','northoutside');
    title(ax, sprintf('%s n%d  —  L/R whisk (BP %g-%g Hz) + breath (BP %g-%g Hz)', ...
        animal, kk, BP(1),BP(2),BR_BP(1),BR_BP(2)), 'Interpreter','none');
end

% ================= helpers =================
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
