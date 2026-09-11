% whisk_LR_slow_overlay_RZ  (script)
% Quick check of L/R synchrony of the SLOW whisking component (set-point):
% low-pass (<LP_HZ) the whisker angle, demean, overlay L (green) vs R (blue),
% ONE figure per session, with the L-R correlation in the title.

% ============================ USER-EDITABLE ============================
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
LP_HZ    = 1;           % low-pass cutoff for the slow component (Hz)
fpsW     = 400;         % whisk camera fps
TWIN     = [];          % [] = full recording; [t0 t1] (s) to zoom every figure
EXCLUDE  = "";          % animal ids to skip ("" = none; e.g. "5840027")
% ======================================================================

[bl,al] = butter(3, LP_HZ/(fpsW/2), 'low');
cL=[0 0.55 0]; cR=[0 0.4 0.85];

sess = list_sessions(whiskDir);
assert(~isempty(sess),'no *_whisk_n*.csv in %s', whiskDir);
nfig = 0;
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE)>0 && strcmp(animal,char(EXCLUDE)), continue; end
    M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun)), 0.6);
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;
    sL = filtfilt(bl,al, fillmissing(La,'linear')); sL = sL - mean(sL,'omitnan');   % slow set-point, demeaned
    sR = filtfilt(bl,al, fillmissing(Ra,'linear')); sR = sR - mean(sR,'omitnan');
    g = isfinite(sL)&isfinite(sR); r = corr(sL(g), sR(g));

    figure('Color','w','Position',[60 80 1300 380]); hold on; grid on;
    plot(t, sL, '-', 'Color',cL, 'LineWidth',1.0);
    plot(t, sR, '-', 'Color',cR, 'LineWidth',1.0);
    xlabel('time (s)'); ylabel(sprintf('slow whisk set-point (deg, <%g Hz, demeaned)', LP_HZ));
    legend({'L','R'}, 'Orientation','horizontal','Location','northoutside');
    if ~isempty(TWIN), xlim(TWIN); else, xlim([0 t(end)]); end
    title(sprintf('%s n%d   slow whisking L vs R   r = %.2f', animal, kRun, r), 'Interpreter','none');
    fprintf('%s n%d: L-R slow-component r = %.2f\n', animal, kRun, r);
    nfig = nfig + 1;
end
fprintf('slow L/R overlay figures: %d\n', nfig);

% ================= helpers =================
function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));
    S = {}; key = {};
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
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
