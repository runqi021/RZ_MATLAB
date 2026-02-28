clear; clc;

%% ------------------------- USER SETTINGS --------------------------------
dataDir = 'E:\RZ\Whisker\260111_whisking_KA_dorsal_20nL_right\260111_whisker_KA_dorsal_R_20nL-RZ-2026-01-12\videos';

fps = 100; fs = fps;

% ---- choose N whiskers (DLC point IDs; each point = x,y,p triplet) ----
nwhisk_to_track = 3;                 % <-- YOU decide
whisk_point_ids = [1 2 3 4 5];       % uses first N entries

% DLC likelihood handling
pmin        = 0.5;
fill_method = 'linear';              % 'linear' or 'pchip'

% Bandpass filter
hp_Hz = 1;
lp_Hz = 49;                          % MUST be < fs/2
filt_order = 3;

% Multitaper params (Chronux)
fpass = [0.1 45];                    % with fs=100, Nyq=50
TW    = 3;                           % time-bandwidth
K     = 2*TW - 1;

% ---- plotting windows ----
t_full_xlim = [30 60];               % like your style (full trace panel)
t_zoom = [10 11];                    % zoom panel

% ---- layout ----
nPerCol = 10;                        % tiles per column (like you)

% ---- band power for whisking ----
bandHz = [8 45];

% ---- time axis for band power (assume each video is 5 min apart) ----
startTime_min = 30;                  % <-- default start time (minutes)
dt_min        = 5;                   % <-- 5 minutes between files (assumption)

% ---- PSD choice ----
use_diff_for_psd = true;             % true = mtspectrumc(diff(s_pos_f)), false = mtspectrumc(s_pos_f)

%% ------------------------- FIND FILES -----------------------------------
files = dir(fullfile(dataDir, "*.csv"));
[~, idx] = sort({files.name});
files = files(idx);

fprintf('Found %d CSV files\n', numel(files));
if isempty(files), error('No CSV files found.'); end

%% ------------------------- DESIGN FILTER --------------------------------
nyq = fs/2;
if ~(hp_Hz > 0 && lp_Hz > hp_Hz && lp_Hz < nyq)
    error('Bad bandpass: need 0<hp<lp<Nyq. hp=%g lp=%g Nyq=%g', hp_Hz, lp_Hz, nyq);
end
[bBP,aBP] = butter(filt_order, [hp_Hz lp_Hz]/nyq, 'bandpass');

%% ------------------------- MULTITAPER PARAMS ----------------------------
params_spec.Fs     = fs;
params_spec.tapers = [TW, K];
params_spec.pad    = 0;
params_spec.fpass  = fpass;
params_spec.err    = [2 0.05];       % Type 2 => jackknife CI, p=0.05 (95% CI)

%% ------------------------- SELECT WHISKERS ------------------------------
whisk_point_ids = whisk_point_ids(:)';
if numel(whisk_point_ids) < nwhisk_to_track
    error('whisk_point_ids must contain at least nwhisk_to_track entries.');
end
whisk_point_ids = whisk_point_ids(1:nwhisk_to_track);

% colors for overlay
col = lines(max(nwhisk_to_track,3));

%% ------------------------- PREALLOC HOLDERS -----------------------------
Nfile = numel(files);
ok = false(Nfile,1);

t_all   = cell(Nfile,1);
p1_all  = cell(Nfile,1);             % plot DLC prob from whisker #1

pos_raw_all = cell(Nfile,1);         % each: [T x nWhisk]
pos_f_all   = cell(Nfile,1);         % each: [T x nWhisk]

fk_all      = cell(Nfile,1);
Sk_pos_all  = cell(Nfile,1);         % each: [F x nWhisk]

% ---- band power per file + CI FOR ALL WHISKERS ----
bandPow_all = nan(Nfile, nwhisk_to_track);
bandPow_lo  = nan(Nfile, nwhisk_to_track);
bandPow_hi  = nan(Nfile, nwhisk_to_track);

fileLabel   = strings(Nfile,1);

%% ------------------------- MAIN LOOP (WHOLE FOLDER) ---------------------
for i = 1:Nfile
    fn = fullfile(files(i).folder, files(i).name);
    fprintf('[%d/%d] %s\n', i, Nfile, files(i).name);

    data = readmatrix(fn, 'NumHeaderLines', 3);
    if size(data,2) < 4
        warning('Skipping %s: not enough columns', files(i).name);
        continue;
    end

    % validate point availability in this file
    nTriplets = floor((size(data,2) - 1)/3);
    if any(whisk_point_ids > nTriplets)
        warning('Skipping %s: requested points exceed available (%d).', files(i).name, nTriplets);
        continue;
    end

    frame = data(:,1);
    t = frame / fps;

    T = numel(t);
    pos_raw = nan(T, nwhisk_to_track);
    pos_f   = nan(T, nwhisk_to_track);

    p_for_plot = nan(T,1);

    Sk_pos = [];
    fk = [];

    % Store PSD + CI per whisker for band integration
    Spos_w = [];   % [F x nWhisk]
    Slo_w  = [];   % [F x nWhisk]
    Shi_w  = [];   % [F x nWhisk]

    for w = 1:nwhisk_to_track
        pid = whisk_point_ids(w);
        cx = 2 + 3*(pid-1);
        cy = 3 + 3*(pid-1);
        cp = 4 + 3*(pid-1);

        x = data(:,cx);
        y = data(:,cy);
        p = data(:,cp);

        bad = (p < pmin) | isnan(p);
        x(bad) = NaN; y(bad) = NaN;

        x = fillmissing(x, fill_method);
        y = fillmissing(y, fill_method);
        p = fillmissing(p, 'nearest');

        if w == 1
            p_for_plot = p;
        end

        % -------- POSITION hypotenuse relative to mean --------
        x0 = x - mean(x,'omitnan');
        y0 = y - mean(y,'omitnan');
        s_pos = hypot(x0, y0);

        % detrend + filter
        s_pos   = detrend(s_pos);
        s_pos_f = filtfilt(bBP,aBP,s_pos);

        pos_raw(:,w) = s_pos;
        pos_f(:,w)   = s_pos_f;

        % PSD input choice
        if use_diff_for_psd
            psd_input = diff(s_pos_f(:));
        else
            psd_input = s_pos_f(:);
        end

        % PSD + jackknife CI (Serr)
        [Spos, fk, Serr] = mtspectrumc(psd_input, params_spec);
        Slo = Serr(1,:).';
        Shi = Serr(2,:).';

        if isempty(Sk_pos)
            Sk_pos  = nan(numel(Spos), nwhisk_to_track);
            Spos_w  = nan(numel(Spos), nwhisk_to_track);
            Slo_w   = nan(numel(Spos), nwhisk_to_track);
            Shi_w   = nan(numel(Spos), nwhisk_to_track);
        end

        Sk_pos(:,w) = Spos(:);

        % Save for band integration
        Spos_w(:,w) = Spos(:);
        Slo_w(:,w)  = Slo;
        Shi_w(:,w)  = Shi;
    end

    % ---- band power + CI (ALL whiskers) ----
    idxBand = (fk >= bandHz(1)) & (fk <= bandHz(2));
    for w = 1:nwhisk_to_track
        bandPow_all(i,w) = trapz(fk(idxBand), Spos_w(idxBand,w));
        bandPow_lo(i,w)  = trapz(fk(idxBand), Slo_w(idxBand,w));
        bandPow_hi(i,w)  = trapz(fk(idxBand), Shi_w(idxBand,w));
    end

    % store
    t_all{i}        = t;
    p1_all{i}       = p_for_plot;

    pos_raw_all{i}  = pos_raw;
    pos_f_all{i}    = pos_f;

    fk_all{i}       = fk;
    Sk_pos_all{i}   = Sk_pos;

    ok(i) = true;
    fileLabel(i) = string(files(i).name);
end

validIdx = find(ok);
if isempty(validIdx), error('No valid CSVs processed.'); end

% compress arrays to valid only
t_all       = t_all(validIdx);
p1_all      = p1_all(validIdx);
pos_raw_all = pos_raw_all(validIdx);
pos_f_all   = pos_f_all(validIdx);
fk_all      = fk_all(validIdx);
Sk_pos_all  = Sk_pos_all(validIdx);

bandPow_all = bandPow_all(validIdx,:);
bandPow_lo  = bandPow_lo(validIdx,:);
bandPow_hi  = bandPow_hi(validIdx,:);

fileLabel   = fileLabel(validIdx);

N = numel(validIdx);

%% ------------------------- ONE BIG SUMMARY: TIME (POS ONLY) --------------
nColBlocks = ceil(N / nPerCol);
nRows = nPerCol;
nCols = nColBlocks * 2;  % [full | zoom]

figT = figure('Color','w');
tlT  = tiledlayout(figT, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

ax_full = gobjects(N,1);
ax_zoom = gobjects(N,1);

for k = 1:N
    colBlock = floor((k-1)/nPerCol);
    row      = mod((k-1), nPerCol) + 1;

    col_full = colBlock*2 + 1;
    col_zoom = colBlock*2 + 2;

    tile_full = (row-1)*nCols + col_full;
    tile_zoom = (row-1)*nCols + col_zoom;

    t = t_all{k};
    Praw = pos_raw_all{k};
    Pf   = pos_f_all{k};

    % FULL
    ax_full(k) = nexttile(tlT, tile_full); hold(ax_full(k),'on');
    for w = 1:nwhisk_to_track
        plot(ax_full(k), t, Praw(:,w), 'Color', [col(w,:) 0.15], 'LineWidth', 0.1);
        plot(ax_full(k), t, Pf(:,w),   'Color', [col(w,:) 0.85], 'LineWidth', 0.2);
    end

    yyaxis(ax_full(k),'left');
    xlim(ax_full(k), t_full_xlim);
    ylim(ax_full(k), [-10 25]);
    if row == ceil(nRows/2), ylabel(ax_full(k), 'pos |(x,y)-mean| (px)'); end

    yyaxis(ax_full(k),'right');
    plot(ax_full(k), t, p1_all{k}, 'k', 'LineWidth', 0.1);
    ylim(ax_full(k), [0 1]);
    if row == ceil(nRows/2), ylabel(ax_full(k), 'DLC prob (whisk1)'); end

    grid(ax_full(k),'on');
    ttl = strrep(regexprep(files(validIdx(k)).name,'DLC.*$',''),'_','\_');
    title(ax_full(k), ttl);

    if row < nRows, set(ax_full(k),'XTickLabel',[]); else, xlabel(ax_full(k),'Time (s)'); end

    % ZOOM
    ax_zoom(k) = nexttile(tlT, tile_zoom); hold(ax_zoom(k),'on');
    for w = 1:nwhisk_to_track
        plot(ax_zoom(k), t, Praw(:,w), 'Color', [col(w,:) 0.15], 'LineWidth', 0.6);
        plot(ax_zoom(k), t, Pf(:,w),   'Color', [col(w,:) 0.85], 'LineWidth', 0.8);
    end
    xlim(ax_zoom(k), t_zoom);
    ylim(ax_zoom(k), [-5 5]);
    grid(ax_zoom(k),'on');
    title(ax_zoom(k), sprintf('pos %g–%g s', t_zoom(1), t_zoom(2)));

    if row < nRows, set(ax_zoom(k),'XTickLabel',[]); else, xlabel(ax_zoom(k),'Time (s)'); end
end

linkaxes(ax_full,'x');
linkaxes(ax_zoom,'x');

%% ------------------------- ONE BIG SUMMARY: PSD (POS ONLY) ---------------
nColsPSD = ceil(N / nPerCol);
nRowsPSD = nPerCol;

figP = figure('Color','w');
tlP  = tiledlayout(figP, nRowsPSD, nColsPSD, 'TileSpacing','compact', 'Padding','compact');
axP = gobjects(N,1);

for k = 1:N
    colIdx = floor((k-1)/nPerCol) + 1;
    rowIdx = mod((k-1), nPerCol) + 1;
    tile   = (rowIdx-1)*nColsPSD + colIdx;

    axP(k) = nexttile(tlP, tile); hold(axP(k),'on');

    fk   = fk_all{k};
    Spos = Sk_pos_all{k};   % [F x nWhisk]

    for w = 1:nwhisk_to_track
        plot(axP(k), fk, Spos(:,w), 'Color', col(w,:), 'LineWidth', 1.0);
    end

    xlim(axP(k), fpass);
    ylim(axP(k), [0 0.35]);  % adjust if needed
    grid(axP(k),'on');
    title(axP(k), strrep(regexprep(files(validIdx(k)).name,'DLC.*$',''),'_','\_'));

    if rowIdx < nRowsPSD
        set(axP(k),'XTickLabel',[]);
    else
        xlabel(axP(k),'Frequency (Hz)');
    end

    if colIdx == 1 && rowIdx == ceil(nRowsPSD/2)
        if use_diff_for_psd
            ylabel(axP(k),'Power (PSD of diff(pos))');
        else
            ylabel(axP(k),'Power (PSD of pos)');
        end
    end
end
linkaxes(axP,'x');

%% ------------------------- BAND POWER vs TIME (+ CI), OVERLAY N ----------

t_minutes = startTime_min + (0:N-1)' * dt_min;

figBP = figure('Color','w'); hold on;

for w = 1:nwhisk_to_track
    y  = bandPow_all(:,w);
    yL = bandPow_lo(:,w);
    yH = bandPow_hi(:,w);

    % shaded CI band for this whisker
    fill([t_minutes; flipud(t_minutes)], [yL; flipud(yH)], ...
        col(w,:), 'EdgeColor','none', 'FaceAlpha', 0.18);

    % mean curve
    plot(t_minutes, y, '-o', 'Color', col(w,:), 'LineWidth', 1.5, 'MarkerSize', 5);
end

grid on;
xlabel('Time (min)');

if use_diff_for_psd
    ylabel(sprintf('Band power %g–%g Hz (∫ PSD(diff(pos)) df)', bandHz(1), bandHz(2)));
    title(sprintf('Whisking band power vs time (jackknife CI) — PSD of diff(pos); Δt=%g min; start=%g min', dt_min, startTime_min));
else
    ylabel(sprintf('Band power %g–%g Hz (∫ PSD(pos) df)', bandHz(1), bandHz(2)));
    title(sprintf('Whisking band power vs time (jackknife CI) — PSD of pos; Δt=%g min; start=%g min', dt_min, startTime_min));
end

leg = strings(nwhisk_to_track,1);
for w = 1:nwhisk_to_track
    leg(w) = sprintf('whisk%d (point %d)', w, whisk_point_ids(w));
end
legend(leg, 'Location','best', 'Interpreter','none');

%% ------------------------- DONE -----------------------------------------
fprintf('Done. Valid files: %d\n', N);
