% whisk_LR_angle_filterbank_RZ.m
% Per-session L/R whisker angle filter-bank overview.
% For every session, ONE figure with L (blue) + R (red) overlaid in 4 rows:
%   1) raw whisk angle (deg)
%   2) low-pass < LP_HZ          (setpoint / slow envelope of the angle)
%   3) band-pass BP_HZ           (fast whisking band)
%   4) Hilbert amplitude of HB_HZ band-pass (whisking-power envelope)
%
% Per-session EXCLUDE: list animal/run pairs to drop from plotting. The list
% starts EMPTY (off) — fill EXCLUDE_SESS below once you've eyeballed the figures.

% ========================= USER-EDITABLE ==========================
dataRoot   = "D:\260615_thermalNbasler";
whiskDir   = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW       = 400;            % whisk camera frame rate (Hz)
LIK_THR    = 0.6;            % DLC likelihood gate

LP_HZ      = 5;              % row 2: low-pass cutoff (Hz)
BP_HZ      = [5 100];        % row 3: band-pass (Hz)
HB_HZ      = [5 50];         % row 4: Hilbert-amplitude band-pass (Hz)

% Per-session exclude. Each entry is {animal, runIndex}; runIndex matches n<k>.
% Leave EMPTY to plot every session. To drop e.g. animal 5840031 run 2:
%   EXCLUDE_SESS = {{"5840031", 2}};
EXCLUDE_SESS = {};           % <-- currently OFF (plot all sessions)

% Animal-level exclude (whole animal, any run). Leave "" / [] for none.
EXCLUDE_ANIMAL = "";

SAVE_FIG   = false;          % true = also save each figure as PNG next to script
% ==================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

colL = [0.20 0.40 0.80];     % left  = blue
colR = [0.85 0.10 0.10];     % right = red

% --- session discovery (animal folders of digits, paired to whisk DLC csv) ---
ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE_ANIMAL)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*'));
    [~,o] = sort({rr.name}); rr = rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), ...
                sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess), 'no sessions with whisk csv in %s', whiskDir);

% --- design filters once (sample rate is fixed) ---
[bLP,aLP] = butter(4, LP_HZ/(fpsW/2), 'low');
[bBP,aBP] = butter(3, BP_HZ/(fpsW/2), 'bandpass');
[bHB,aHB] = butter(3, HB_HZ/(fpsW/2), 'bandpass');

nPlotted = 0;
for e = 1:numel(sess)
    animal = sess{e}{1}; kk = sess{e}{2};

    % honor per-session exclude
    if is_excluded(animal, kk, EXCLUDE_SESS)
        fprintf('skip (excluded): %s n%d\n', animal, kk);
        continue;
    end

    % --- load whisk DLC ---
    try
        M = dlc_gate_interp(pick_csv(whiskDir, ...
            sprintf('%s_whisk_n%d', animal, kk)), LIK_THR);
    catch ME
        warning('whisk DLC load failed %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    % --- L/R angle (same convention as whisk_LR_sync_combined_10panels_RZ.m) ---
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)' / fpsW;

    La = fillmissing(La, 'linear', 'EndValues', 'nearest');
    Ra = fillmissing(Ra, 'linear', 'EndValues', 'nearest');

    % demean before filtering (LP/BP/Hilbert operate on the AC component)
    Lc = La - mean(La,'omitnan');
    Rc = Ra - mean(Ra,'omitnan');

    % --- filter bank ---
    Llp = filtfilt(bLP,aLP, Lc);          Rlp = filtfilt(bLP,aLP, Rc);
    Lbp = filtfilt(bBP,aBP, Lc);          Rbp = filtfilt(bBP,aBP, Rc);
    Lhb = abs(hilbert(filtfilt(bHB,aHB, Lc)));
    Rhb = abs(hilbert(filtfilt(bHB,aHB, Rc)));

    % --- figure ---
    f = figure('Color','w', 'Position',[40 40 1400 900], ...
        'Name', sprintf('%s n%d  whisk L/R filter bank', animal, kk));
    tlo = tiledlayout(f, 4, 1, 'TileSpacing','compact', 'Padding','compact');
    title(tlo, sprintf('%s  n%d   L/R whisk angle filter bank   (%.1f s @ %g fps)', ...
        animal, kk, t(end), fpsW), 'FontSize', 11, 'FontWeight','bold', ...
        'Interpreter','none');

    ax1 = nexttile(tlo);
    plot(ax1, t, La, 'Color',colL, 'LineWidth',0.6); hold(ax1,'on');
    plot(ax1, t, Ra, 'Color',colR, 'LineWidth',0.6);
    ylabel(ax1,'angle (deg)'); title(ax1,'raw angle');
    legend(ax1, {'L','R'}, 'Location','northeast','Box','off');

    ax2 = nexttile(tlo);
    plot(ax2, t, Llp, 'Color',colL, 'LineWidth',1.0); hold(ax2,'on');
    plot(ax2, t, Rlp, 'Color',colR, 'LineWidth',1.0);
    ylabel(ax2,'deg'); title(ax2, sprintf('low-pass < %g Hz', LP_HZ));

    ax3 = nexttile(tlo);
    plot(ax3, t, Lbp, 'Color',colL, 'LineWidth',0.5); hold(ax3,'on');
    plot(ax3, t, Rbp, 'Color',colR, 'LineWidth',0.5);
    ylabel(ax3,'deg'); title(ax3, sprintf('band-pass %g–%g Hz', BP_HZ(1), BP_HZ(2)));

    ax4 = nexttile(tlo);
    plot(ax4, t, Lhb, 'Color',colL, 'LineWidth',1.0); hold(ax4,'on');
    plot(ax4, t, Rhb, 'Color',colR, 'LineWidth',1.0);
    ylabel(ax4,'amplitude (deg)'); xlabel(ax4,'time (s)');
    title(ax4, sprintf('Hilbert amplitude of %g–%g Hz band', HB_HZ(1), HB_HZ(2)));

    linkaxes([ax1 ax2 ax3 ax4], 'x');
    xlim(ax1, [t(1) t(end)]);
    arrayfun(@(ax) set(ax,'Box','off'), [ax1 ax2 ax3 ax4]);

    nPlotted = nPlotted + 1;

    if SAVE_FIG
        outPng = fullfile(repoRoot, ...
            sprintf('%s_n%d_whisk_LR_filterbank.png', animal, kk));
        exportgraphics(f, outPng, 'Resolution', 150);
        fprintf('saved %s\n', outPng);
    end
end

fprintf('\nplotted %d / %d sessions\n', nPlotted, numel(sess));

% ===================== LOCAL FUNCTIONS =====================

function tf = is_excluded(animal, kk, excl)
    tf = false;
    for q = 1:numel(excl)
        ent = excl{q};
        if strcmp(char(ent{1}), char(animal)) && double(ent{2}) == kk
            tf = true; return;
        end
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv for %s in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix] = max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
