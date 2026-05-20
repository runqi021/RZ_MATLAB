% whisker_L_3w_batch.m
%  Per-whisker analysis with INDEPENDENT per-whisker run selection.
%  For each selected whisker:
%    - X-Y trajectory grid (one panel per run, demeaned)
%    - PSD_X + PSD_Y grid (one panel per run)
%    - PSD_X + PSD_Y average (mean ± SEM across that whisker's runs)
%
%  Set per-whisker run lists at the top (run_idx_vL1, run_idx_vR1, etc.)
%  Set `whiskers` to pick which whiskers to actually process.
%
%  Grid size is auto-sized per whisker based on its own nPanels.
%  xlim/ylim are shared across all plotted whiskers for comparability.
%
%  Outputs (saved into folderPath), per whisker:
%    whisker_<vW>_<N>runs_traces_XY.pdf
%    whisker_<vW>_<N>runs_spectra.pdf
%    whisker_<vW>_<N>runs_spectra_avg.pdf
%    whisker_<vW>_<N>runs.mat
%
%  Dependencies: Chronux (mtspectrumc).

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
% ======================================================================

%% ========================= USER PARAMETERS ===========================
folderPath  = "C:\260407_KA_electro_dorsal_whisking";
fps         = 90;

% Per-whisker run lists (declare ONE per whisker you want to process).
% Any run number from 1 to 30 is valid.
run_idx_vL1 = [1 2 3 4 5 6 7 8 27 28 29 30];
run_idx_vL3 = 1:30;
run_idx_vR1 = [1 2 3 4 5 6 7 8 25 26 27 28];
run_idx_vR2 = 1:30;

% Whiskers to actually process (subset of declared). Order matters for limits.
whiskers    = {'vL1', 'vR1'};

% Side definitions: maps whisker name → CSV pattern + color
sides = struct( ...
    'label',   {'L', 'R'}, ...
    'pattern', {'*DLC_Resnet50_whisker_dorsalApr*.csv', ...
                '*DLC_Resnet50_whisker_dorsal_RApr*.csv'}, ...
    'points',  {{'vL1', 'vL3'}, {'vR1', 'vR2'}}, ...
    'cols',    {[0.00 0.45 0.85;   % vL1 blue
                 0.10 0.65 0.20], ... % vL3 green
                [0.85 0.10 0.10;   % vR1 red
                 0.95 0.55 0.10]}); % vR2 orange

pmin        = 0.5;
fillMeth    = 'linear';
bpBand      = [];

TW          = 6;
fpass       = [5 44];

traceXlim   = [20 35];
maxGridCols = 6;         % grid auto-sized: cols = min(maxGridCols, nPanels)

doSave      = true;
% ======================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

% Build whisker → run_idx map from the declared variables above
run_idx_map = struct();
for ww = {'vL1','vL3','vR1','vR2'}
    varName = ['run_idx_' ww{1}];
    if exist(varName, 'var')
        run_idx_map.(ww{1}) = eval(varName);
    end
end

% Build whisker → side lookup
whiskerToSide = struct();
for s = 1:numel(sides)
    for i = 1:numel(sides(s).points)
        whiskerToSide.(sides(s).points{i}) = s;
    end
end

if ~isempty(bpBand)
    upperCut = min(bpBand(2), 0.999 * fps / 2);
    bpSpec   = [bpBand(1), upperCut];
    fprintf('Bandpass: %.2f - %.2f Hz\n', bpSpec(1), bpSpec(2));
else
    bpSpec = [];
end

fprintf('Processing %d whisker(s): %s\n', numel(whiskers), strjoin(whiskers, ', '));

%% =================== PASS 1: LOAD EACH WHISKER =======================
wData = repmat(struct(), numel(whiskers), 1);

for wi = 1:numel(whiskers)
    wname = whiskers{wi};
    assert(isfield(run_idx_map, wname), 'No run_idx_%s declared', wname);
    assert(isfield(whiskerToSide, wname), 'Unknown whisker %s', wname);

    runs = run_idx_map.(wname);
    sidx = whiskerToSide.(wname);
    sideLabel = sides(sidx).label;
    pattern   = sides(sidx).pattern;
    points    = sides(sidx).points;
    colRow    = find(strcmp(points, wname), 1);
    col       = sides(sidx).cols(colRow, :);

    fprintf('\n=== %s (side %s, %d runs) ===\n', wname, sideLabel, numel(runs));

    % Find CSVs + match to runs
    Files = dir(fullfile(folderPath, pattern));
    assert(~isempty(Files), 'No CSVs for side %s', sideLabel);
    fileRuns = zeros(numel(Files), 1);
    for i = 1:numel(Files)
        tok = regexp(Files(i).name, 'run(\d+)', 'tokens', 'once');
        if isempty(tok), fileRuns(i) = NaN; else, fileRuns(i) = str2double(tok{1}); end
    end

    nPanelsW = numel(runs);
    X_all   = cell(nPanelsW, 1);
    Y_all   = cell(nPanelsW, 1);
    t_all   = cell(nPanelsW, 1);
    run_lbl = strings(nPanelsW, 1);

    for k = 1:nPanelsW
        idx = find(fileRuns == runs(k), 1);
        assert(~isempty(idx), 'No %s CSV for run %03d', sideLabel, runs(k));
        csv = fullfile(Files(idx).folder, Files(idx).name);

        bp_line = local_bodypart_line(csv);
        bp_parts = strsplit(strtrim(bp_line), ',');
        unique_bps = {};
        for c = 2:numel(bp_parts)
            name = strtrim(bp_parts{c});
            if isempty(unique_bps) || ~strcmp(unique_bps{end}, name)
                unique_bps{end+1} = name; %#ok<AGROW>
            end
        end
        pid = find(strcmp(unique_bps, wname), 1);
        assert(~isempty(pid), 'Bodypart %s not in %s', wname, csv);

        data = readmatrix(csv, 'NumHeaderLines', 3);
        T    = size(data, 1);
        tvec = data(:, 1) / fps;
        cx = 2 + 3*(pid - 1);
        cy = 3 + 3*(pid - 1);
        cp = 4 + 3*(pid - 1);
        x = data(:, cx);  y = data(:, cy);  p = data(:, cp);
        bad = (p < pmin) | isnan(p);
        x(bad) = NaN;  y(bad) = NaN;
        x = fillmissing(x, fillMeth, 'EndValues', 'nearest');
        y = fillmissing(y, fillMeth, 'EndValues', 'nearest');
        if ~isempty(bpSpec)
            x = bandpass(x, bpSpec, fps);
            y = bandpass(y, bpSpec, fps);
        end

        X_all{k} = x - mean(x, 'omitnan');
        Y_all{k} = y - mean(y, 'omitnan');
        t_all{k} = tvec;

        [~, stem] = fileparts(csv);
        run_lbl(k) = local_label(stem);
        fprintf('  [%2d/%2d] %s  (%d samples)\n', k, nPanelsW, run_lbl(k), T);
    end

    % PSD
    params_mt.Fs     = fps;
    params_mt.tapers = [TW, 2*TW - 1];
    params_mt.pad    = 0;
    params_mt.fpass  = fpass;
    params_mt.err    = 0;

    Lmin = min(cellfun(@numel, X_all));
    f_psd = [];
    S_x = nan(1, nPanelsW);
    S_y = nan(1, nPanelsW);
    for k = 1:nPanelsW
        [sx, fk] = mtspectrumc(X_all{k}(1:Lmin), params_mt);
        [sy, ~ ] = mtspectrumc(Y_all{k}(1:Lmin), params_mt);
        if isempty(f_psd)
            f_psd = fk(:);
            S_x = nan(numel(fk), nPanelsW);
            S_y = nan(numel(fk), nPanelsW);
        end
        S_x(:, k) = sx(:);
        S_y(:, k) = sy(:);
    end

    wData(wi).name      = wname;
    wData(wi).side      = sideLabel;
    wData(wi).col       = col;
    wData(wi).runs      = runs;
    wData(wi).run_lbl   = run_lbl;
    wData(wi).X_all     = X_all;
    wData(wi).Y_all     = Y_all;
    wData(wi).t_all     = t_all;
    wData(wi).f_psd     = f_psd;
    wData(wi).S_x       = S_x;
    wData(wi).S_y       = S_y;
    wData(wi).params_mt = params_mt;
end

%% ============ PASS 1.5: GLOBAL BOUNDS ACROSS ALL WHISKERS ============
xSliceMin = +inf; xSliceMax = -inf;
ySliceMin = +inf; ySliceMax = -inf;
xFullMin  = +inf; xFullMax  = -inf;
yFullMin  = +inf; yFullMax  = -inf;
specMax    = 0;
specAvgMax = 0;
specAvgMin = +inf;

for wi = 1:numel(wData)
    X_all = wData(wi).X_all;
    Y_all = wData(wi).Y_all;
    S_x   = wData(wi).S_x;
    S_y   = wData(wi).S_y;
    nPw   = numel(wData(wi).runs);

    iS = max(1, round(traceXlim(1) * fps) + 1);
    iE = min(cellfun(@numel, X_all));
    iE = min(iE, round(traceXlim(2) * fps));

    xSliceMin = min(xSliceMin, min(cellfun(@(v) min(v(iS:iE)), X_all)));
    xSliceMax = max(xSliceMax, max(cellfun(@(v) max(v(iS:iE)), X_all)));
    ySliceMin = min(ySliceMin, min(cellfun(@(v) min(v(iS:iE)), Y_all)));
    ySliceMax = max(ySliceMax, max(cellfun(@(v) max(v(iS:iE)), Y_all)));

    xFullMin  = min(xFullMin, min(cellfun(@min, X_all)));
    xFullMax  = max(xFullMax, max(cellfun(@max, X_all)));
    yFullMin  = min(yFullMin, min(cellfun(@min, Y_all)));
    yFullMax  = max(yFullMax, max(cellfun(@max, Y_all)));

    specMax = max(specMax, max(arrayfun(@(kk) max(S_x(:,kk) + S_y(:,kk)), 1:nPw)));

    S_sum = S_x + S_y;
    m_w  = mean(S_sum, 2, 'omitnan');
    se_w = std(S_sum, 0, 2, 'omitnan') ./ sqrt(nPw);
    specAvgMax = max(specAvgMax, max(m_w + se_w));
    specAvgMin = min(specAvgMin, min(max(m_w - se_w, eps)));
end

fprintf('\nShared limits (across %d whiskers):\n', numel(wData));
fprintf('  X-Y box:      x [%.2f, %.2f], y [%.2f, %.2f]\n', xFullMin, xFullMax, yFullMin, yFullMax);
fprintf('  Spec ymax:    %.4g\n', specMax);
fprintf('  SpecAvg ylim: [%.4g, %.4g]\n', specAvgMin, specAvgMax);

%% =================== PASS 2: FIGURES + SAVE ==========================
for wi = 1:numel(wData)
    wname = wData(wi).name;
    col   = wData(wi).col;
    nPw   = numel(wData(wi).runs);
    run_lbl = wData(wi).run_lbl;

    % Auto grid per whisker
    nCols_w = min(maxGridCols, nPw);
    nRows_w = ceil(nPw / nCols_w);
    figW    = 5 * nCols_w;
    figH    = 4 * nRows_w;

    X_all = wData(wi).X_all;
    Y_all = wData(wi).Y_all;
    t_all = wData(wi).t_all;
    f_psd = wData(wi).f_psd;
    S_x   = wData(wi).S_x;
    S_y   = wData(wi).S_y;

    % --- X vs Y trajectory ---
    fT_XY = figure('Units', 'centimeters', 'Position', [1 1 figW figH], 'Color', 'w', ...
        'Name', sprintf('%s X-Y trajectory (%d runs)', wname, nPw));
    tl = tiledlayout(fT_XY, nRows_w, nCols_w, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s — X-Y raw trajectory (demeaned)', wname), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', col);
    for k = 1:nPw
        ax = nexttile(tl);
        plot(ax, X_all{k}, Y_all{k}, 'Color', col, 'LineWidth', 0.3);
        axis(ax, 'equal');
        xlim(ax, [xFullMin, xFullMax]);  ylim(ax, [yFullMin, yFullMax]);
        title(ax, run_lbl(k), 'FontSize', 7);
        box(ax, 'on');
        if k > nPw - nCols_w, xlabel(ax, 'X (px)'); end
        if mod(k - 1, nCols_w) == 0, ylabel(ax, 'Y (px)'); end
    end

    % --- Per-run X+Y PSD grid ---
    fS_sum = figure('Units', 'centimeters', 'Position', [1 1 figW figH], 'Color', 'w', ...
        'Name', sprintf('%s X+Y PSDs (%d runs)', wname, nPw));
    tl = tiledlayout(fS_sum, nRows_w, nCols_w, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s — PSD_X + PSD_Y (lin-lin, TW=%g)', wname, TW), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', col);
    for k = 1:nPw
        ax = nexttile(tl);
        plot(ax, f_psd, S_x(:,k) + S_y(:,k), 'Color', col, 'LineWidth', 0.7);
        xlim(ax, fpass);  ylim(ax, [0, specMax]);
        title(ax, run_lbl(k), 'FontSize', 7);
        box(ax, 'on');
        if k > nPw - nCols_w, xlabel(ax, 'Hz'); end
        if mod(k - 1, nCols_w) == 0, ylabel(ax, 'px^2/Hz'); end
    end

    % --- Average spectrum (lin-y) ---
    fS_avg = figure('Units', 'centimeters', 'Position', [1 1 12 9], 'Color', 'w', ...
        'Name', sprintf('%s mean PSD (%d runs)', wname, nPw));
    axA = axes(fS_avg);  hold(axA, 'on');
    S_sum = S_x + S_y;
    m_w   = mean(S_sum, 2, 'omitnan');
    se_w  = std(S_sum, 0, 2, 'omitnan') ./ sqrt(nPw);
    fill(axA, [f_psd; flipud(f_psd)], [m_w + se_w; flipud(m_w - se_w)], ...
        col, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot(axA, f_psd, m_w, 'Color', col, 'LineWidth', 1.6);
    xlim(axA, fpass);  ylim(axA, [0, specAvgMax]);
    xlabel(axA, 'Hz');  ylabel(axA, 'px^2/Hz');
    title(axA, sprintf('%s — mean \\pm SEM across %d runs (PSD_X + PSD_Y)', wname, nPw), ...
        'FontSize', 9, 'Color', col);
    box(axA, 'on');

    % --- Average spectrum on LOG Y: log(mean ± sem) ---  [DISABLED]
    % fS_avg_log = figure('Units', 'centimeters', 'Position', [1 1 12 9], 'Color', 'w', ...
    %     'Name', sprintf('%s mean PSD log-y (%d runs)', wname, nPw));
    % axAL = axes(fS_avg_log);  hold(axAL, 'on');
    % up_log = m_w + se_w;
    % lo_log = max(m_w - se_w, eps);
    % fill(axAL, [f_psd; flipud(f_psd)], [up_log; flipud(lo_log)], ...
    %     col, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    % plot(axAL, f_psd, m_w, 'Color', col, 'LineWidth', 1.6);
    % set(axAL, 'YScale', 'log');
    % xlim(axAL, fpass);  ylim(axAL, [max(specAvgMin, eps), specAvgMax]);
    % xlabel(axAL, 'Hz');  ylabel(axAL, 'px^2/Hz (log)');
    % title(axAL, sprintf('%s — log(mean \\pm sem), log-y (%d runs)', wname, nPw), ...
    %     'FontSize', 9, 'Color', col);
    % box(axAL, 'on');

    if doSave
        pdfTXY  = fullfile(folderPath, sprintf('whisker_%s_%druns_traces_XY.pdf',   wname, nPw));
        pdfSS   = fullfile(folderPath, sprintf('whisker_%s_%druns_spectra.pdf',     wname, nPw));
        pdfSAvg = fullfile(folderPath, sprintf('whisker_%s_%druns_spectra_avg.pdf', wname, nPw));
        print(fT_XY,  '-dpdf', '-bestfit', char(pdfTXY));   fprintf('Saved: %s\n', pdfTXY);
        print(fS_sum, '-dpdf', '-bestfit', char(pdfSS));    fprintf('Saved: %s\n', pdfSS);
        print(fS_avg, '-dpdf', '-bestfit', char(pdfSAvg));  fprintf('Saved: %s\n', pdfSAvg);

        matPath = fullfile(folderPath, sprintf('whisker_%s_%druns.mat', wname, nPw));
        results = struct( ...
            'folderPath', char(folderPath), ...
            'whisker',    wname, ...
            'side',       wData(wi).side, ...
            'fps',        fps, ...
            'runs',       wData(wi).runs, ...
            'run_lbl',    {cellstr(run_lbl)}, ...
            'X_all',      {X_all}, ...
            'Y_all',      {Y_all}, ...
            't_all',      {t_all}, ...
            'f_psd',      f_psd, ...
            'S_x',        S_x, ...
            'S_y',        S_y, ...
            'params_mt',  wData(wi).params_mt);
        save(matPath, '-struct', 'results', '-v7.3');
        fprintf('Saved: %s\n', matPath);
    end
end

fprintf('\nDone (%d whiskers, %d PDFs total).\n', numel(wData), 3*numel(wData));

%% ==================== LOCAL FUNCTIONS ================================
function bp_line = local_bodypart_line(csv)
    fid = fopen(csv, 'r');
    fgetl(fid);
    bp_line = fgetl(fid);
    fclose(fid);
end

function lbl = local_label(stem)
    rn = regexp(stem, 'run(\d+)', 'tokens', 'once');
    ts = regexp(stem, '_(\d{8})_(\d{6})_', 'tokens', 'once');
    if ~isempty(rn) && ~isempty(ts)
        hh = ts{2}(1:2);  mm = ts{2}(3:4);
        lbl = sprintf('run%s [%s:%s]', rn{1}, hh, mm);
    elseif ~isempty(rn)
        lbl = sprintf('run%s', rn{1});
    else
        lbl = stem;
    end
    lbl = string(lbl);
end
