clear; clc;

%% ------------------------- USER SETTINGS --------------------------------
dataDir = 'E:\RZ\Whisker\251224_whisking_KA\pooled';
fps = 100;     fs = fps;

% DLC likelihood handling
pmin        = 0.5;
fill_method = 'linear';     % 'linear' or 'pchip'

hp_Hz = 0.1;
lp_Hz = 49;                % MUST be < fs/2
filt_order = 3;

% multitaper params (Chronux)
fpass = [1 45];             % note: with fs=100, Nyq=50
TW    = 3;                  % time-bandwidth
% params_spec.tapers = [TW, 2*TW-1];

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

%% ------------------------- PREALLOC HOLDERS -----------------------------
N = numel(files);
ok = false(N,1);

t_all  = cell(N,1);
x_all  = cell(N,1);   % raw chosen signal
xf_all = cell(N,1);   % filtered chosen signal

fk_all = cell(N,1);
Sk_all = cell(N,1);
Slo_all = cell(N,1);
Shi_all = cell(N,1);

%% ------------------------- MAIN LOOP ------------------------------------
for i = 1:N
    fn = fullfile(files(i).folder, files(i).name);
    fprintf('[%d/%d] %s\n', i, N, files(i).name);

    data = readmatrix(fn, 'NumHeaderLines', 3);
    if size(data,2) < 4
        warning('Skipping %s: not enough columns', files(i).name);
        continue;
    end

    frame = data(:,1);
    t = frame / fps;

    % dot1
    x1 = data(:,2); y1 = data(:,3); p1 = data(:,4);
    % dot2
    %x2 = data(:,5); y2 = data(:,6); p2 = data(:,7);

    % likelihood mask
    bad1 = (p1 < pmin) | isnan(p1);
    %bad2 = (p2 < pmin) | isnan(p2);
    x1(bad1)=NaN; y1(bad1)=NaN;
    %x2(bad2)=NaN; y2(bad2)=NaN;

    % fill
    x1 = fillmissing(x1, fill_method);
    y1 = fillmissing(y1, fill_method);
    %x2 = fillmissing(x2, fill_method);
    %y2 = fillmissing(y2, fill_method);


    s1 = hypot(x1,y1);
    %s2 = hypot(x2,y2);
    sig = s1; 

    sig = detrend(sig);

    sig_f = sig;
    sig_f = filtfilt(bBP,aBP,sig);

    % ----- multitaper PSD (Chronux) -----
    params_spec.Fs     = fs;
    params_spec.tapers = [TW, 2*TW - 1];
    params_spec.pad    = 0;
    params_spec.fpass  = fpass;
    params_spec.err    = [2 0.05];

    [Sk, fk, Sconf] = mtspectrumc(diff(sig_f(:)), params_spec);

    % store
    t_all{i}   = t;
    p_all{i}   = p1;
    x_all{i}   = sig;
    xf_all{i}  = sig_f;

    fk_all{i}  = fk;
    Sk_all{i}  = Sk;
    Slo_all{i} = Sconf(1,:);
    Shi_all{i} = Sconf(2,:);

    ok(i) = true;
end

validIdx = find(ok);
if isempty(validIdx), error('No valid CSVs processed.'); end

%%
N = numel(validIdx);

nPerCol = 10;                     % traces per column
nColBlocks = ceil(N / nPerCol);   % how many left→right blocks
nRows = nPerCol;
nCols = nColBlocks * 2;           % each block = [full | zoom]

fig1 = figure('Color','w');
set(fig1,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);

tl1 = tiledlayout(fig1, nRows, nCols, ...
    'TileSpacing','compact', 'Padding','compact');

ax_full = gobjects(N,1);
ax_zoom = gobjects(N,1);

t_zoom = [10 11];

for k = 1:N
    i = validIdx(k);

    % ----- compute grid position -----
    colBlock = floor((k-1) / nPerCol);   % 0-based
    row      = mod(k-1, nPerCol) + 1;    % 1..10

    col_full = colBlock*2 + 1;
    col_zoom = colBlock*2 + 2;

    tile_full = (row-1)*nCols + col_full;
    tile_zoom = (row-1)*nCols + col_zoom;

    % -------- LEFT: full trace --------
    ax_full(k) = nexttile(tl1, tile_full); hold on;

    y_raw = x_all{i};
    y     = xf_all{i};

    plot(t_all{i}, y_raw, 'Color', [0 0 1 0.4], 'LineWidth', 0.1);
    plot(t_all{i}, y,     'Color', [1 0 0 0.4], 'LineWidth', 0.1);
    ylim([-10 25]);

    if row == ceil(nRows/2)
        ylabel('ΔWhisker pos (px)');
    end

    yyaxis right;
    plot(t_all{i}, p_all{i}, 'k', 'LineWidth', 0.1);
    xlim([30 60]);
    ylim([0 1]);

    if row == ceil(nRows/2)
        ylabel('DLC prob');
    end

    grid on;

    title(strrep(regexprep(files(i).name,'DLC.*$',''),'_','\_'));

    if row < nRows
        set(gca,'XTickLabel',[]);
    else
        xlabel('Time (s)');
    end

    % -------- RIGHT: zoom --------
    ax_zoom(k) = nexttile(tl1, tile_zoom); hold on;

    plot(t_all{i}, y_raw, 'Color', [0 0 1 0.7], 'LineWidth', 1);
    plot(t_all{i}, y,     'Color', [1 0 0 0.7], 'LineWidth', 1);
    xlim(t_zoom);
    ylim([-13 13]);
    grid on;

    title(sprintf('%g–%g s', t_zoom(1), t_zoom(2)));

    if row < nRows
        set(gca,'XTickLabel',[]);
    else
        xlabel('Time (s)');
    end
end

linkaxes(ax_full,'x');
linkaxes(ax_zoom,'x');


%% ------------------------- FIG 2: STACKED MT PSD (freq) ------------------
N = numel(validIdx);

nPerCol   = 10;                  % traces per column
nCols     = ceil(N / nPerCol);
nRows     = nPerCol;

fig2 = figure('Color','w');
tl2 = tiledlayout(fig2, nRows, nCols, ...
    'TileSpacing','compact', 'Padding','compact');

ax2 = gobjects(N,1);

for k = 1:N
    i = validIdx(k);

    % ---- grid position (left -> right, 10 per column) ----
    col = floor((k-1)/nPerCol) + 1;   % 1..nCols
    row = mod(k-1, nPerCol) + 1;      % 1..10

    tile = (row-1)*nCols + col;
    ax2(k) = nexttile(tl2, tile); hold on;

    fk  = fk_all{i};
    Sk  = Sk_all{i};
    Slo = Slo_all{i};
    Shi = Shi_all{i};

    plot(fk, Sk, 'k', 'LineWidth', 1.0);
    % plot(fk, Slo, 'k--', 'LineWidth', 0.6);
    % plot(fk, Shi, 'k--', 'LineWidth', 0.6);

    xlim(fpass);
    ylim([0 1]);
    grid on;

    title(strrep(regexprep(files(i).name,'DLC.*$',''),'_','\_'));

    if row < nRows
        set(gca,'XTickLabel',[]);
    else
        xlabel('Frequency (Hz)');
    end

    if col == 1 && row == ceil(nRows/2)
        ylabel('Power');
    end
end

linkaxes(ax2,'x');


% %% Select to plot:
% plotSel = [1 3 4 5 6 7 8 9 10];   % <-- YOU SET THIS (e.g. [1,3,5,7])
% 
% % clamp / validate selection
% plotSel = plotSel(plotSel >= 1 & plotSel <= numel(validIdx));
% if isempty(plotSel)
%     error('plotSel is empty or out of range. validIdx has %d entries.', numel(validIdx));
% end
% 
% plotIdx = validIdx(plotSel);
% nPlot   = numel(plotIdx);
% 
% fig1 = figure('Color','w');
% set(fig1,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
% 
% tl1 = tiledlayout(fig1, nPlot, 2, ...
%     'TileSpacing','compact', 'Padding','compact');
% 
% ax_full = gobjects(nPlot,1);
% ax_zoom = gobjects(nPlot,1);
% 
% t_zoom = [10 15];                 % seconds
% midRow = ceil(nPlot/2);           % row to keep y-labels
% 
% for k = 1:nPlot
%     i = plotIdx(k);
% 
%     % -------- LEFT: full trace --------
%     ax_full(k) = nexttile(tl1, (k-1)*2 + 1); hold(ax_full(k),'on');
% 
%     y_raw = x_all{i};
%     y     = xf_all{i};
% 
%     plot(ax_full(k), t_all{i}, y_raw, 'Color', [0 0 1 0.4], 'LineWidth', 0.1);
%     plot(ax_full(k), t_all{i}, y,     'Color', [1 0 0 0.4], 'LineWidth', 0.1);
% 
%     % LEFT axis
%     yyaxis(ax_full(k),'left');
%     ylim(ax_full(k), [-10 25]);
% 
%     if k == midRow
%         ylabel(ax_full(k), 'ΔWhisker position (pixel)');
%     end
% 
%     yyaxis(ax_full(k),'right');
%     plot(ax_full(k), t_all{i}, p_all{i}, 'Color', 'k', 'LineWidth', 0.1);
%     ylim(ax_full(k), [0 1]);
% 
%     if k == midRow
%         ylabel(ax_full(k), 'probability (DLC)');
%     end
% 
%     grid(ax_full(k),'on');
% 
%     title(ax_full(k), sprintf('%s', strrep(regexprep(strrep(files(i).name,'__','_'),'DLC.*$', ''),'_','\_')));
% 
%     if k < nPlot
%         set(ax_full(k), 'XTickLabel', []);
%     else
%         xlabel(ax_full(k), 'Time (s)');
%     end
% 
%     % -------- RIGHT: zoomed (10–15 s) --------
%     ax_zoom(k) = nexttile(tl1, (k-1)*2 + 2); hold(ax_zoom(k),'on');
% 
%     plot(ax_zoom(k), t_all{i}, y_raw, 'Color', [0 0 1 0.7], 'LineWidth', 1);
%     plot(ax_zoom(k), t_all{i}, y,     'Color', [1 0 0 0.7], 'LineWidth', 1);
% 
%     xlim(ax_zoom(k), t_zoom);
%     ylim(ax_zoom(k), [-13 13]);
%     grid(ax_zoom(k),'on');
% 
%     title(ax_zoom(k), sprintf('%g–%g s', t_zoom(1), t_zoom(2)));
% 
%     if k < nPlot
%         set(ax_zoom(k), 'XTickLabel', []);
%     else
%         xlabel(ax_zoom(k), 'Time (s)');
%     end
% end
% 
% % link axes
% linkaxes(ax_full,'x');
% linkaxes(ax_zoom,'x');
% 
% %
% fig2 = figure('Color','w');
% tl2 = tiledlayout(fig2, numel(plotIdx), 1, 'TileSpacing','compact', 'Padding','compact');
% 
% ax2 = gobjects(numel(plotIdx),1);
% midRow = ceil(numel(plotIdx)/2);   % row to keep y-labels
% 
% for k = 1:numel(plotIdx)
%     i = plotIdx(k);
%     ax2(k) = nexttile(tl2); hold(ax2(k),'on');
% 
%     fk  = fk_all{i};
%     Sk  = Sk_all{i};
%     %Slo = Slo_all{i};
%     %Shi = Shi_all{i};
% 
%     plot(ax2(k), fk, Sk, 'k', 'LineWidth', 1.0);
%     %plot(ax2(k), fk, Slo, 'k--', 'LineWidth', 0.6);
%     %plot(ax2(k), fk, Shi, 'k--', 'LineWidth', 0.6);
% 
%     if k == midRow
%         ylabel(ax2(k), 'Power');
%     end
% 
%     xlim(ax2(k), fpass);
%     ylim(ax2(k), [0 1]);
%     grid(ax2(k),'on');
% 
%     title(ax2(k), sprintf('%s', strrep(regexprep(strrep(files(i).name,'__','_'),'DLC.*$', ''),'_','\_')));
% 
%     if k < numel(plotIdx)
%         set(ax2(k), 'XTickLabel', []);
%     else
%         xlabel(ax2(k), 'Frequency (Hz)');
%     end
% end
% 
% linkaxes(ax2,'x');
% 
% %% select to plot (final)
% % ------------------------- FINAL FIG: TRACE (2 cols) + PSD (1 col) ------
% fig = figure('Color','w');
% set(fig,'defaultAxesColorOrder',[[0,0,0]; [0,0,0]]);
% 
% % 3 columns: [trace | trace | PSD]
% tl = tiledlayout(fig, nPlot, 3, ...
%     'TileSpacing','compact', 'Padding','compact');
% 
% ax_time = gobjects(nPlot,1);
% ax_psd  = gobjects(nPlot,1);
% 
% t_zoom = [10 20];                 % seconds
% midRow = ceil(nPlot/2);           % only label once
% 
% for k = 1:nPlot
%     i = plotIdx(k);
% 
%     y_raw = x_all{i};
%     y     = xf_all{i};
% 
%     % -------- LEFT: zoomed time trace (spans 2 columns) --------
%     ax_time(k) = nexttile(tl, (k-1)*3 + 1, [1 2]);  % span cols 1–2
%     hold(ax_time(k),'on');
% 
%     %plot(ax_time(k), t_all{i}, y_raw, 'Color', [0 0 1 0.7], 'LineWidth', 1);
%     plot(ax_time(k), t_all{i}, y,     'Color', [0 0 0], 'LineWidth', 0.6);
% 
%     xlim(ax_time(k), t_zoom);
%     ylim(ax_time(k), [-13 13]);
%     grid(ax_time(k),'on');
% 
%     if k == midRow
%         ylabel(ax_time(k), 'Δposition (pixel)');
%     end
% 
%     if k < nPlot
%         set(ax_time(k), 'XTickLabel', []);
%     else
%         xlabel(ax_time(k), 'Time (s)');
%     end
% 
%     % -------- RIGHT: multitaper PSD (single column) --------
%     ax_psd(k) = nexttile(tl, (k-1)*3 + 3);
%     hold(ax_psd(k),'on');
% 
%     fk = fk_all{i};
%     Sk = Sk_all{i};
% 
%     plot(ax_psd(k), fk, Sk, 'k', 'LineWidth', 1.0);
%     xlim(ax_psd(k), fpass);
%     ylim(ax_psd(k), [0 1.5]);
%     grid(ax_psd(k),'on');
% 
%     if k == midRow
%         ylabel(ax_psd(k), 'Power');
%     end
% 
%     if k < nPlot
%         set(ax_psd(k), 'XTickLabel', []);
%     else
%         xlabel(ax_psd(k), 'Frequency (Hz)');
%     end
% 
%     % -------- Row title (on trace axis) --------
%     %ttl = strrep(regexprep(strrep(files(i).name,'__','_'),'DLC.*$', ''),'_','\_');
%     %title(ax_time(k), ttl);
%     if k ==1; title(ax_time(k), 'Representative whisker tracking (10s)'); title(ax_psd(k), 'Spectrum'); end;
% end
% 
% % link x-axes
% linkaxes(ax_time,'x');
% linkaxes(ax_psd,'x');
% 
% %% ------------------------- BAND POWER (6–30 Hz) -------------------------
% bandHz = [8 40];
% 
% bandPower = nan(nPlot,1);
% 
% for k = 1:numel(validIdx)
%     i = validIdx(k);
% 
%     fk = fk_all{i};
%     Sk = Sk_all{i};
% 
%     % frequency mask
%     idxBand = fk >= bandHz(1) & fk <= bandHz(2);
% 
%     % integrate power over band (area under PSD)
%     bandPower(k) = trapz(fk(idxBand), Sk(idxBand));
% end
% 
% %% ------------------------- PLOT BAND POWER ------------------------------
% figBP = figure('Color','w');
% bar(bandPower, 'FaceColor', [0.2 0.2 0.2]);
% ylabel('Integrated power (8–40 Hz)');
% xlabel('Selected trace #');
% title('Whisker power (8–40 Hz)');
% grid on;
% 
% xticks(1:nPlot);

