% whisk_LR_scatter_RZ  (script)
% L vs R whisker scatter, ONE panel per session, all sessions in a tiled grid.
% Same signal pipeline as whisk_LR_slow_overlay_RZ (slow set-point: low-pass
% <LP_HZ, demeaned). x = L angle, y = R angle; identity line + L-R corr r in title.
% Animal 5840027 is excluded by default.

% ============================ USER-EDITABLE ============================
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
LP_HZ    = 1;            % low-pass cutoff for the slow component (Hz)
fpsW     = 400;          % whisk camera fps
SLOW     = true;         % true = slow set-point (<LP_HZ); false = full angle, demeaned
TWIN     = [];           % [] = full recording; [t0 t1] (s) to restrict every session
EXCLUDE  = "5840027";    % animal ids to skip ("" = none)
MAXPTS   = 8000;         % max scatter points per panel (random thinning; Inf = all)
% ======================================================================

[bl,al] = butter(3, LP_HZ/(fpsW/2), 'low');

sess = list_sessions(whiskDir);
assert(~isempty(sess),'no *_whisk_n*.csv in %s', whiskDir);

% keep only sessions we will actually plot
keep = true(1,numel(sess));
for e = 1:numel(sess)
    if strlength(EXCLUDE)>0 && strcmp(sess{e}{1},char(EXCLUDE)), keep(e)=false; end
end
sess = sess(keep);
assert(~isempty(sess),'all sessions excluded');

nP = numel(sess);
nc = ceil(sqrt(nP)); nr = ceil(nP/nc);
figure('Color','w','Position',[40 40 320*nc+120 300*nr+80]);
tl = tiledlayout(nr, nc, 'TileSpacing','compact', 'Padding','compact');

allL = []; allR = [];   % pooled across sessions for the second figure

for e = 1:nP
    animal = sess{e}{1}; kRun = sess{e}{2};
    M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun)), 0.6);
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;

    if SLOW
        sL = filtfilt(bl,al, fillmissing(La,'linear'));
        sR = filtfilt(bl,al, fillmissing(Ra,'linear'));
    else
        sL = fillmissing(La,'linear'); sR = fillmissing(Ra,'linear');
    end
    sL = sL - mean(sL,'omitnan');
    sR = sR - mean(sR,'omitnan');

    win = true(size(t)); if ~isempty(TWIN), win = t>=TWIN(1) & t<=TWIN(2); end
    g  = isfinite(sL)&isfinite(sR)&win;
    gi = find(g);
    xL = sL(gi); yR = sR(gi); tg = t(gi);
    r  = corr(xL, yR);

    % thin for plotting only (corr uses all points)
    idx = 1:numel(xL);
    if isfinite(MAXPTS) && numel(xL) > MAXPTS
        idx = round(linspace(1, numel(xL), MAXPTS));
    end

    nexttile; hold on; grid on; axis square;
    scatter(xL(idx), yR(idx), 4, tg(idx), 'filled', ...
            'MarkerFaceAlpha',0.25);                       % colored by time (frame)
    lim = max([abs(xL); abs(yR); eps]);
    plot([-lim lim], [-lim lim], 'k--', 'LineWidth',0.8);  % identity
    xlim([-lim lim]); ylim([-lim lim]);
    xlabel('L (deg)'); ylabel('R (deg)');
    title(sprintf('%s n%d   r = %.2f', animal, kRun, r), 'Interpreter','none');
    fprintf('%s n%d: L-R %s r = %.2f  (n=%d)\n', animal, kRun, ...
            ternary(SLOW,'slow','full'), r, numel(xL));

    allL = [allL; xL]; allR = [allR; yR]; %#ok<AGROW>  pool for figure 2
end
cb = colorbar(nexttile(nP)); cb.Label.String = 'time (s)';
ttl = ternary(SLOW, sprintf('slow set-point (<%g Hz)', LP_HZ), 'full angle');
title(tl, sprintf('L vs R whisker scatter — %s  (excl %s)', ttl, ...
      ternary(strlength(EXCLUDE)>0, char(EXCLUDE), 'none')), 'Interpreter','none');
fprintf('scatter panels: %d\n', nP);

% ---- Figure 2: all sessions pooled, plain black dots, one panel ----
rAll = corr(allL, allR);
idxA = 1:numel(allL);
if isfinite(MAXPTS) && numel(allL) > MAXPTS
    idxA = round(linspace(1, numel(allL), MAXPTS));
end
figure('Color','w','Position',[120 120 560 560]); hold on; grid on; axis square;
scatter(allL(idxA), allR(idxA), 4, 'k', 'filled', 'MarkerFaceAlpha',0.15);
lim = max([abs(allL); abs(allR); eps]);
plot([-lim lim], [-lim lim], 'k--', 'LineWidth',0.8);   % identity
xlim([-lim lim]); ylim([-lim lim]);
xlabel('L (deg)'); ylabel('R (deg)');
title(sprintf('L vs R whisker — all sessions pooled — %s   r = %.2f  (n=%d)', ...
      ttl, rAll, numel(allL)), 'Interpreter','none');
fprintf('pooled L-R r = %.2f  (n=%d)\n', rAll, numel(allL));

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
function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
