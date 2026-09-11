% whisk_LR_breath_window_RZ  (script)
% Three-trace overlay of ONE session over a FIXED time window (no epoch gating):
%   left  y-axis : L/R fast-whisk BP angle (deg)
%   right y-axis : band-passed breathing (inhale up)
% Just pick the animal/run and the [t0 t1] window — draws whatever is there.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
animal    = "5916296";   % session animal id
kRun      = 2;           % run index (n<k>)
TWIN      = [39 40];     % [t0 t1] window to plot (s)
BP        = [6 30];      % fast-whisking bandpass (Hz)
BR_BP     = [2 15];      % breath bandpass (Hz)
fpsW      = 400;         % whisk camera fps
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths
[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

% ---- whisk angles (lik<0.6 -> linear interp, same as other scripts) ----
M  = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), kRun)), 0.6);
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));   % LEFT  (x mirrored)
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));    % RIGHT
xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
t  = (0:numel(xL)-1)'/fpsW;

% ---- breathing (optional) ----
tB = []; brf = [];
try
    Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', char(animal), kRun)), dataRoot);
    if isfile(Pn.breath)
        Bs = load(Pn.breath); br = Bs.breath(:); fpsB = double(Bs.fps);
        tB = (0:numel(br)-1)'/fpsB;
        [bb,ab] = butter(2, BR_BP/(fpsB/2), 'bandpass');
        brf = filtfilt(bb,ab, fillmissing(br,'linear'));
    end
catch
end

% ---- plot (whisk left axis, breath right axis) ----
cL = [0 0.5 0]; cR = [0 0.4 0.85]; cB = [0.10 0.10 0.10];
figure('Color','w','Position',[120 200 1100 440]);
yyaxis left;  hold on;
hL = plot(t, xL, '-', 'Color',cL, 'LineWidth',1.0);
hR = plot(t, xR, '-', 'Color',cR, 'LineWidth',1.0);
ylabel('whisk angle (deg, BP)'); set(gca,'YColor','k');
yyaxis right; hold on; hB = [];
if ~isempty(brf)
    hB = plot(tB, brf, '-', 'Color',cB, 'LineWidth',1.3);
end
ylabel('breathing (BP, inhale up)'); set(gca,'YColor',cB);
xlim(TWIN); grid on; box on; xlabel('time (s)');
hh = [hL hR]; ll = {'L whisk','R whisk'};
if ~isempty(hB), hh=[hh hB]; ll=[ll {'breath'}]; end
legend(hh, ll, 'Orientation','horizontal','Location','northoutside');
title(sprintf('%s n%d  —  L/R whisk + breath  (%.2f-%.2f s)', char(animal), kRun, TWIN(1), TWIN(2)), ...
      'Interpreter','none');

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
