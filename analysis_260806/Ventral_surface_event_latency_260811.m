function Ventral_surface_event_latency_260811()
%% Ventral_surface_event_latency_260811
% TIME-DOMAIN event latency for EVERY ACTIVE CELL in the ventral surface
% archive. NO PHASE ANYWHERE. Same archive scan, same exclusions, same cell
% identity and the same active gate as Ventral_surface_breath_time_summary_260808,
% so the two summaries describe the SAME cell set and can be read side by side.
%
% Output: D:\Ventral_surface_summary\event_latency_260811\
%
% ------------------------------- LAYOUT ----------------------------------
% Identical grammar to breath_time_summary: three columns, row 1 of page 1 is
% the pooled population, cells overflow onto _p02, _p03...  Every slot is a
% PAIR of axes -- but here the pair is the two TRIGGERS, not time and phase:
%
%          ALL              SIG            NON-SIG
%   r1  [onset|peak ]   [onset|peak ]   [onset|peak ]    <- pooled over cells
%   ---------------------------------------------------
%   r2  [cell 1     ]   [cell 2     ]   [cell 3     ]    <- SIG first
%
%   RED   = inspiration ONSET triggered
%   BLUE  = inspiratory PEAK triggered
%
% ------------------------------ THE MEASURE ------------------------------
% Per cell, per trigger:
%   1. WINDOW = one mean IBI, centred on the trigger, and that IBI is THIS
%      CELL'S own -- pooled over every recording the cell appears in. A cell
%      recorded only in a 2.9 s-cycle session gets a 2.9 s window; a 0.5 s cell
%      gets 0.5 s. SECONDS ARE THEREFORE NOT COMPARABLE BETWEEN CELLS of
%      different breath rates. Each panel prints its own IBI.
%   2. every event goes to its NEAREST trigger with a SIGNED latency:
%      < 0 LEADS the trigger, > 0 LAGS it. One event, one trigger.
%   3. PERMUTATION TEST FIRST, then description. Nothing is described until
%      the latencies are shown to be CONCENTRATED:
%           statistic  max(D) - mean(D), where D is the latency histogram on the
%                      frame lattice smoothed by a FIXED-bandwidth Gaussian
%                      (SD = IBI/smoothDiv)
%           null       nPerm circular shifts of the BREATH TRIGGER train
%           p          (1 + #{T_null >= T_obs}) / (1 + nPerm)
%      FIXED bandwidth is essential -- an adaptive one would differ between the
%      observation and every shuffle, so the null would not be comparable.
%
%      WHY NOT THE IQR. An earlier version used the IQR of the latencies as the
%      statistic, on the principle that the tested quantity should be the
%      reported one. That principle cost too much: THE IQR ONLY SEES THE MIDDLE
%      HALF OF THE EVENTS, so it is blind to a sharp peak sitting on a broad
%      background. Measured on ChAT/0521 roi5 r1 (117 events, ~40% of them in a
%      167 ms peak): IQR = 708 ms -> p = 0.77 "not significant", while the peak
%      statistic on the SAME null gives p = 0.0005. The peak was plainly visible
%      in the panel. A statistic that misses what the eye sees is the wrong
%      statistic.
%
%      This is NOT a return to the original correlogram test. That one was
%      broken by its NULL (an FFT correlogram whose ceil() clipping piled the
%      last trigger and last events into a shared bin, inflating lag 0), not by
%      being peak-based. Latencies are still recomputed from scratch against
%      every shifted trigger train, and the flat cell that exposed the old bug
%      -- Vglut2/1124 pFN roi5_1400-1230-0 r4 -- still fails here (p = 0.13
%      onset, 0.05 peak, both far above alpha = 0.001).
%   4. DESCRIPTION = MODE and FWHM of that same smoothed curve, plus the
%      FRACTION OF EVENTS INSIDE THE PEAK. Read off the curve the test scored,
%      so the number describes exactly what was detected.
%
%      MEDIAN and IQR are still written to the CSV but are NOT what the panels
%      report, and for the same reason: they describe the whole in-window
%      distribution. On the ChAT cell above, median/IQR read +400/708 ms while
%      mode/FWHM read +467/167 ms -- only the latter describes the peak that is
%      actually there. The event fraction is printed so a peak carrying 20% of
%      the events can never be mistaken for a tight distribution.
%
% WHY SHIFT THE BREATH. One trigger train serves every cell in a recording, so
% one set of shifts nulls them all. Implemented as its exact equivalent:
% shifting the triggers +v is shifting the EVENTS -v, so the events are shifted
% instead and the trigger midpoints are built once per observation. All nPerm
% shifts are then one matrix operation. Shifts within minShiftIBI cycles of
% no-shift are rejected -- half a breath is the same alignment displaced.
%
% ------------------------- NO LAG COMPENSATION ---------------------------
% caLagSec = 0. The ca_spike_data events sit ON the dF/F peak and so lag the
% true spike by the GCaMP rise, but that correction is NOT applied here: it is
% an assumption, it slides every median by the same amount, and it changes
% whether cells read as leading or following. Left off deliberately. The
% Vglut2/1124 rising-edge +1 frame IS applied -- that is a timing fact about
% that session's trigger, not a model of the indicator.
%
% ----------------------- NO FAMILY-WISE CORRECTION -----------------------
% Every p here is RAW. There is no Bonferroni, no FDR, no correction of any
% kind across the ~450 cells tested. At alpha = 0.001 and 1000 shifts the
% p FLOOR is 1/1001 = 0.000999, so "p < 0.001" means EXACTLY "not one of 1000
% shifted breath trains beat this cell" -- every surviving cell prints the same
% p. Raising nPerm is the only way to separate them.
%
% -------------------------- ACTIVE / TESTED ------------------------------
%   active cell <=> pooled nnz(spike_train>0) > activeMinEv   (the archive's
%                   own criterion, so the cell set matches breath_time_summary)
%   TESTED      <=> pooled in-window events >= minTestEv
% A median and an IQR off six events are noise, so cells between the two gates
% are DRAWN but not tested and not described. They are labelled "n<NN" and
% counted separately -- they are neither significant nor non-significant.
%
% Runqi Zhang / 2026-08-11
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(scriptDir);
% coh_cfg_260727 lives under analysis_260727, NOT here. The old line pointed at a
% folder that does not exist and only worked because the function happened to be
% left on the path by an earlier script -- a fresh MATLAB would have failed.
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath, 'event_latency_260811');

groups    = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};
scan_dirs = {'ChAT', 'Vglut2', 'Vgat', 'Sst', 'Sert'};

cell_link_sources = { ...
    'Sert',   '0721', fullfile(rootPath, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(rootPath, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(rootPath, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0806', fullfile(rootPath, 'Sst',    '0806', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0807', fullfile(rootPath, 'Sst',    '0807', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0810', fullfile(rootPath, 'Vglut2', '0810', 'cell_pooled', 'cell_link.mat')
    };

nDrop        = 30;      % breath frames tossed up front, to align with Ca
fallback_fps = 30;
activeMinEv  = 5;       % ACTIVE  = pooled nnz(spike_train>0) > this
minTestEv    = 20;      % TESTED  = pooled in-window events >= this
caLagSec     = 0;       % NO GCaMP lead compensation.  See header.
fix1124      = true;    % Vglut2/1124 rising-edge trigger: delay breath 1 frame
ampFrac      = 0.20;    % breath-cycle amplitude QC, same as breath_time_summary

nPerm        = 10000;   % circular-shift permutations.  RAISED from 1000 because
                        % the polar figure plots -log10(p) as the RADIUS: with
                        % 1000 shifts the smallest p is 1/1001, so every
                        % significant cell would sit at exactly r = 3 and the
                        % radial axis would carry no information at all. 10000
                        % shifts give a floor of 1e-4 and a usable range 0..4.
                        % alpha is unchanged at 0.001.
alphaPerm    = 0.001;   % SIG <=> p < this.  RAW, no family-wise correction.
minShiftIBI  = 2;       % reject shifts within this many cycles of no-shift
nBins        = 25;      % display bins across the window (one IBI)
smoothDiv    = 25;      % Gaussian SD of the tested/drawn density = IBI / this
rngSeed      = 260811;

rowsPerPage  = 6;
colsPerRow   = 3;
doSave       = true;
% =====================================================================

excludeRecordings = coh_cfg_260727().excludeRecordings;
rng(rngSeed);
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

TRIG    = {'onset','peak'};
COL_T   = [0.85 0.20 0.10; 0.15 0.45 0.85];    % RED onset, BLUE peak
COL_POP = [0.25 0.25 0.25];

fprintf('\n=========== Ventral_surface_event_latency_260811 ===========\n');
fprintf('active = pooled >%d events | tested = >=%d in-window events\n', activeMinEv, minTestEv);
fprintf('SIG = permutation p < %.4g (RAW, no family-wise correction), %d shifts\n', alphaPerm, nPerm);
fprintf('p floor = 1/(nPerm+1) = %.5f  -> "p<%.3g" means ZERO shifts beat the cell\n', ...
        1/(nPerm+1), alphaPerm);
fprintf('window = one mean IBI of THAT CELL | caLagSec = %.3f (no lag comp)\n\n', caLagSec);

%% ---- cell-identity lookup ----------------------------------------------
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
for s = 1:size(cell_link_sources,1)
    pre = sprintf('%s/%s', cell_link_sources{s,1}, cell_link_sources{s,2});
    if ~isfile(cell_link_sources{s,3})
        warning('cell_link missing for %s -- that session stays per-ROI.', pre); continue;
    end
    Lk = load(cell_link_sources{s,3},'link');  Tk = Lk.link.obsT;
    for i = 1:height(Tk)
        kk = sprintf('%s/%s/%d', pre, Tk.rec_name(i), Tk.maskL_label(i));
        if isnan(Tk.cell_id(i)), tossed_set(kk) = true;
        else,                    cell_map(kk) = sprintf('%s#c%d', pre, Tk.cell_id(i));
        end
    end
    fprintf('  cell identity %-14s %4d masks -> %4d cells (%d tossed)\n', pre, height(Tk), ...
            numel(unique(Tk.cell_id(~isnan(Tk.cell_id)))), nnz(isnan(Tk.cell_id)));
end

%% ---- SCAN ---------------------------------------------------------------
REC = struct('fps',{},'T',{},'Trec',{},'trig',{},'gi',{},'folder',{});
OBS = struct('rec',{},'ev',{},'gi',{},'label',{},'cellKey',{},'nEv',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), warning('Folder missing: %s', gdir); continue; end
    allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
    fprintf('\n=== [%s] %d recordings with spikes ===\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;
        recName    = folder_basename(folderPath);
        recDate    = date_from_path(folderPath, gdir);
        gname = sname;
        if is_io_path(folderPath, gdir), gname = 'IO'; end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), continue; end
        if any(strcmp(recName, excludeRecordings))
            fprintf('  skip (excluded): %s\n', recName); continue;
        end
        try
            bpFile = fullfile(folderPath,'breath_peak_pc1.mat');
            ipFile = fullfile(folderPath,'breath_insp_start_pc1.mat');
            if ~isfile(bpFile) || ~isfile(ipFile)
                fprintf('  skip (no breath): %s\n', recName); continue;
            end
            fps = detect_session_fps(folderPath, fallback_fps);
            CA  = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
            nCa = numel(CA.roi_spikes(1).spike_train);

            BP = load(bpFile);  IP = load(ipFile);
            bw = detrend(double(BP.breath(:)));
            bw(1:min(nDrop,numel(bw))) = [];  bw = bw - mean(bw);
            peak_idx = round(BP.insp_onset_idx(:)) - nDrop;   % insp PEAK
            foot_idx = round(IP.insp_start_idx(:)) - nDrop;   % insp ONSET

            if fix1124 && strcmpi(sname,'Vglut2') && strcmp(recDate,'1124')
                peak_idx = peak_idx + 1;  foot_idx = foot_idx + 1;
                bw = [bw(1); bw(1:end-1)];
            end

            T = min([numel(bw), nCa]);
            peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
            foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
            if numel(peak_idx) < 3 || numel(foot_idx) < 3
                fprintf('  skip (too few landmarks): %s\n', recName); continue;
            end
            bw = bw(1:T);

            % ---- amplitude QC on the cycles, same rule as breath_time_summary ----
            bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
            ft  = sort(foot_idx);
            amp = nan(numel(ft)-1,1);  pkf = nan(numel(ft)-1,1);
            for i = 1:numel(ft)-1
                q = peak_idx(peak_idx>ft(i) & peak_idx<ft(i+1));
                if ~isempty(q), amp(i) = bwz(q(1)) - bwz(ft(i)); pkf(i) = q(1); end
            end
            good   = amp > ampFrac*median(amp,'omitnan');
            onTrig = ft(good);
            pkTrig = pkf(good);  pkTrig = pkTrig(~isnan(pkTrig));
            if numel(onTrig) < 3 || numel(pkTrig) < 3
                fprintf('  skip (QC left too few triggers): %s\n', recName); continue;
            end

            REC(end+1) = struct('fps',fps,'T',T,'Trec',T/fps, ...
                                'trig',{{sort(onTrig(:))/fps, sort(pkTrig(:))/fps}}, ...
                                'gi',gi,'folder',folderPath); %#ok<AGROW>
            recIdx = numel(REC);

            lag  = round(caLagSec*fps);
            nInc = 0;
            for rid = 1:numel(CA.roi_spikes)
                st = double(CA.roi_spikes(rid).spike_train(:));
                st = st(1:min(T,numel(st)));
                if numel(st) < T, st(end+1:T,1) = 0; end %#ok<AGROW>
                ckin = sprintf('%s/%s/%s/%d', sname, recDate, recName, rid);
                if isKey(tossed_set, ckin), continue; end
                if lag > 0, st = [st(1+lag:end); zeros(lag,1)]; end %#ok<AGROW>
                lab = sprintf('%s/%s/%s/%d', gname, recDate, recName, rid);
                if isKey(cell_map, ckin), ck = cell_map(ckin); else, ck = ['roi:' lab]; end
                OBS(end+1) = struct('rec',recIdx,'ev',find(st>0),'gi',gi, ...
                                    'label',lab,'cellKey',ck,'nEv',nnz(st>0)); %#ok<AGROW>
                nInc = nInc + 1;
            end
            fprintf('  [%2d] %-52s %-7s fps %.2f  %3d ROI  %3d/%3d cycles pass QC\n', ...
                    kk, recName(1:min(52,end)), gname, fps, nInc, numel(onTrig), numel(ft)-1);
        catch ME
            warning('  ERROR %s: %s', recName, ME.message);
        end
    end
end
if isempty(OBS), error('No observations collected from %s', rootPath); end

%% ---- CELLS --------------------------------------------------------------
[uCell, ~, obsOfCell] = unique({OBS.cellKey}, 'stable');
CELL = struct('key',{},'gi',{},'obs',{},'nEv',{},'nObs',{},'label',{});
for c = 1:numel(uCell)
    m  = find(obsOfCell == c);
    ne = sum([OBS(m).nEv]);
    if ne <= activeMinEv, continue; end
    if numel(m) == 1, lab = OBS(m(1)).label;
    else,             lab = sprintf('%s [+%d rec]', OBS(m(1)).label, numel(m)-1);
    end
    CELL(end+1) = struct('key',uCell{c},'gi',OBS(m(1)).gi,'obs',m(:)', ...
                         'nEv',ne,'nObs',numel(m),'label',lab); %#ok<AGROW>
end
fprintf('\n%d ROI-observations -> %d keys -> %d ACTIVE cells (>%d events)\n', ...
        numel(OBS), numel(uCell), numel(CELL), activeMinEv);

%% ---- per-cell latencies + permutation test ------------------------------
S = struct('ibi',{},'fps',{},'L',{},'mode',{},'fwhm',{},'frac',{},'dens',{},'dctr',{}, ...
           'med',{},'q1',{},'q3',{},'iqr',{},'p',{},'z',{}, ...
           'nEv',{},'nCyc',{},'tested',{},'sig',{});
fprintf('\ntesting %d cells x 2 triggers, %d shifts each...\n', numel(CELL), nPerm);
tStart = tic;
for c = 1:numel(CELL)
    S(c) = cell_latency(CELL(c), OBS, REC, nPerm, minShiftIBI, minTestEv, alphaPerm, smoothDiv); %#ok<AGROW>
    if mod(c,25)==0, fprintf('  %d/%d  (%.0f s)\n', c, numel(CELL), toc(tStart)); end
end
fprintf('  done in %.0f s\n', toc(tStart));

tested = vertcat(S.tested);  sigM = vertcat(S.sig);
for q = 1:2
    fprintf('%-5s: %d tested, %d SIG (p<%.4g)\n', TRIG{q}, nnz(tested(:,q)), nnz(sigM(:,q)), alphaPerm);
end
anySig = any(sigM,2);  anyTest = any(tested,2);
fprintf('cells with EITHER trigger significant: %d of %d tested (%d active)\n', ...
        nnz(anySig), nnz(anyTest), numel(CELL));

%% ---- figures ------------------------------------------------------------
if doSave && ~isfolder(outDir), mkdir(outDir); end
rowsCSV = {};

for gi = 1:numel(groups)
    idx = find([CELL.gi] == gi);
    if isempty(idx), fprintf('\n[%s] no active cells -- skipped\n', groups{gi}); continue; end

    sg  = anySig(idx);  te = anyTest(idx);
    % SIG first (tightest peak IQR first), then tested non-sig, then untested.
    iSig = idx(sg);
    [~,o] = sort(arrayfun(@(k) min_iqr(S(k)), iSig));       iSig = iSig(o);
    iNon = idx(te & ~sg);
    [~,o] = sort(arrayfun(@(k) min([S(k).p Inf]), iNon));   iNon = iNon(o);
    iUnt = idx(~te);
    [~,o] = sort([CELL(iUnt).nEv],'descend');               iUnt = iUnt(o);
    ord  = [iSig, iNon, iUnt];

    medIBI = median(arrayfun(@(k) mean(S(k).ibi,'omitnan'), idx), 'omitnan');
    medFPS = mode(arrayfun(@(k) S(k).fps, idx));
    [pEdges, pCtrs] = lat_axis(medIBI, nBins, medFPS);

    POP = struct('name',{},'H',{},'n',{},'nEv',{});
    POP(1) = struct('name','ALL',    'H',{pool_lat(S(idx),  pEdges)}, 'n',numel(idx), 'nEv',sum(vertcat(S(idx).nEv),'all'));
    POP(2) = struct('name','SIG',    'H',{pool_lat(S(iSig), pEdges)}, 'n',numel(iSig),'nEv',0);
    POP(3) = struct('name','NON-SIG','H',{pool_lat(S([iNon iUnt]), pEdges)},'n',numel(iNon)+numel(iUnt),'nEv',0);

    slotsP1 = colsPerRow*(rowsPerPage-1);
    slotsPn = colsPerRow*rowsPerPage;
    nPage   = 1 + ceil(max(numel(ord)-slotsP1,0)/slotsPn);
    fprintf('\n[%s] %d active (%d sig, %d tested-ns, %d untested), median IBI %.3f s -> %d page(s)\n', ...
            groups{gi}, numel(idx), numel(iSig), numel(iNon), numel(iUnt), medIBI, nPage);

    for pg = 1:nPage
        if pg == 1
            take = ord(1:min(slotsP1,numel(ord)));  nRow = 1 + ceil(numel(take)/colsPerRow);
        else
            a = slotsP1 + (pg-2)*slotsPn + 1;
            take = ord(a:min(a+slotsPn-1,numel(ord)));  nRow = ceil(numel(take)/colsPerRow);
        end
        nRow = max(nRow,1);
        hf = figure('Color','w','Visible','off','Units','centimeters', ...
                    'Position',[1 1 34 3.35*nRow + 2.6]);
        set(hf,'DefaultAxesFontSize',7);
        yTop = 0.915; yBot = 0.050;  rowH = (yTop-yBot)/nRow;
        colW = 0.288; colGap = 0.028; X0 = 0.050;
        subGap = 0.050; subW = (colW-subGap)/2;

        rowsUsed = 0;
        if pg == 1
            for q = 1:3
                x0 = X0 + (q-1)*(colW+colGap);
                y0 = yTop - rowH + rowH*0.30;
                axA = axes('Parent',hf,'Position',[x0,            y0, subW, rowH*0.60]); %#ok<LAXES>
                axB = axes('Parent',hf,'Position',[x0+subW+subGap,y0, subW, rowH*0.60]); %#ok<LAXES>
                if POP(q).n == 0
                    axis(axA,'off'); axis(axB,'off');
                    text(axA,0.5,0.5,sprintf('%s: no cells',POP(q).name), ...
                         'Units','normalized','HorizontalAlignment','center','FontSize',8);
                    continue;
                end
                ttl = sprintf('%s  |  %s POP  |  %d cells', groups{gi}, POP(q).name, POP(q).n);
                draw_pop(axA, pCtrs, POP(q).H{1}, COL_POP, COL_T(1,:), 'ONSET', true, q==1, ttl);
                draw_pop(axB, pCtrs, POP(q).H{2}, COL_POP, COL_T(2,:), 'PEAK',  true, false, '');
            end
            rowsUsed = 1;
        end

        for k = 1:numel(take)
            c  = take(k);
            r  = rowsUsed + ceil(k/colsPerRow);
            q  = mod(k-1,colsPerRow)+1;
            x0 = X0 + (q-1)*(colW+colGap);
            y0 = yTop - r*rowH + rowH*0.30;
            axA = axes('Parent',hf,'Position',[x0,            y0, subW, rowH*0.60]); %#ok<LAXES>
            axB = axes('Parent',hf,'Position',[x0+subW+subGap,y0, subW, rowH*0.60]); %#ok<LAXES>
            isLast = (r == nRow);
            mark = ''; if anySig(c), mark = ' *'; end
            ttl = sprintf('%s%s  n=%d  IBI %.2fs   on %s | pk %s', ...
                          disp_label(CELL(c).label), mark, max(S(c).nEv), mean(S(c).ibi,'omitnan'), ...
                          fmt(S(c),1,minTestEv), fmt(S(c),2,minTestEv));
            draw_cell(axA, S(c), 1, nBins, COL_T(1,:), isLast, q==1, ttl);
            draw_cell(axB, S(c), 2, nBins, COL_T(2,:), isLast, false, '');
            if anySig(c)
                set([axA axB],'LineWidth',1.6);
                set(get(axA,'Title'),'FontWeight','bold');
            end
        end

        sgtitle({sprintf(['%s   |   page %d/%d   |   %d active cells (%d sig, %d tested n.s., %d untested)   ' ...
                          '|   y = epc (events per cycle per bin)'], ...
                         groups{gi}, pg, nPage, numel(idx), numel(iSig), numel(iNon), numel(iUnt)), ...
                 sprintf(['LEFT of each pair = INSPIRATION ONSET triggered (red),   RIGHT = INSPIRATORY PEAK triggered (blue).   ' ...
                          'x = latency, NEGATIVE LEADS the trigger.   window = that cell''s own mean IBI   |   ' ...
                          'curve = KDE, dashed = median']), ...
                 sprintf(['* , bold title, heavy frame = p < %.4g on either trigger.  RAW p, NO family-wise correction.  ' ...
                          '%d shifts -> p floor %.5f.  untested = fewer than %d in-window events.  NO GCaMP lag compensation.'], ...
                         alphaPerm, nPerm, 1/(nPerm+1), minTestEv)}, ...
                'Interpreter','none','FontSize',8);

        if doSave
            base = fullfile(outDir, sprintf('%s_p%02d', groups{gi}, pg));
            exportgraphics(hf,[base '.png'],'Resolution',200,'BackgroundColor','white');
            exportgraphics(hf,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
            fprintf('   saved %s.png/.pdf  (%d rows, %d cells)\n', folder_basename(base), nRow, numel(take));
        end
        close(hf);
    end

    % ---- population-only figure ----
    hp = pop_figure(POP, pCtrs, groups{gi}, numel(idx), numel(iSig), medIBI, ...
                    alphaPerm, nPerm, COL_POP, COL_T, '');
    if doSave
        base = fullfile(outDir, sprintf('%s_population', groups{gi}));
        exportgraphics(hp,[base '.png'],'Resolution',200,'BackgroundColor','white');
        exportgraphics(hp,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
        fprintf('   saved %s.png/.pdf\n', folder_basename(base));
    end
    close(hp);

    for k = 1:numel(ord)
        c = ord(k);
        rowsCSV(end+1,:) = {groups{gi}, CELL(c).key, CELL(c).label, CELL(c).nObs, ...
            mean(S(c).ibi,'omitnan')*1000, S(c).nEv(1), S(c).nEv(2), ...
            S(c).p(1), S(c).p(2), S(c).z(1), S(c).z(2), ...
            S(c).mode(1), S(c).mode(2), S(c).fwhm(1), S(c).fwhm(2), ...
            S(c).frac(1), S(c).frac(2), ...
            S(c).med(1), S(c).med(2), S(c).iqr(1), S(c).iqr(2), ...
            double(S(c).tested(1)), double(S(c).tested(2)), ...
            double(S(c).sig(1)), double(S(c).sig(2)), k}; %#ok<AGROW>
    end
end

%% ---- GRAND population ---------------------------------------------------
% Pooling SECONDS across genotypes mixes 0.5 s and 2.5 s breath cycles, so this
% panel is drawn on the median cycle over ALL cells and is dominated by whichever
% genotypes breathe near it. Read it as a mixed-cycle summary, NOT a latency.
allIBI = arrayfun(@(k) mean(S(k).ibi,'omitnan'), 1:numel(S));
grandIBI = median(allIBI,'omitnan');
grandFPS = mode([S.fps]);
[gEdges, gCtrs] = lat_axis(grandIBI, nBins, grandFPS);
POPG(1) = struct('name','ALL',    'H',{pool_lat(S,          gEdges)},'n',numel(S),'nEv',0);
POPG(2) = struct('name','SIG',    'H',{pool_lat(S(anySig),  gEdges)},'n',nnz(anySig),'nEv',0);
POPG(3) = struct('name','NON-SIG','H',{pool_lat(S(~anySig), gEdges)},'n',nnz(~anySig),'nEv',0);
hp = pop_figure(POPG, gCtrs, 'ALL GENOTYPES', numel(S), nnz(anySig), grandIBI, ...
                alphaPerm, nPerm, COL_POP, COL_T, ...
                'seconds are NOT comparable across genotypes: this mixes 0.5-2.9 s breath cycles');
if doSave
    base = fullfile(outDir,'ALL_GENOTYPES_population');
    exportgraphics(hp,[base '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(hp,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
    fprintf('\n   saved ALL_GENOTYPES_population.png/.pdf\n');
end
close(hp);

%% ---- CSV ----------------------------------------------------------------
if doSave && ~isempty(rowsCSV)
    Tc = cell2table(rowsCSV,'VariableNames', ...
        {'group','cell_key','label','n_recordings','cell_mean_ibi_ms', ...
         'n_events_onset','n_events_peak','p_onset','p_peak','z_onset','z_peak', ...
         'mode_onset_ms','mode_peak_ms','fwhm_onset_ms','fwhm_peak_ms', ...
         'frac_in_peak_onset','frac_in_peak_peak', ...
         'median_onset_ms','median_peak_ms','iqr_onset_ms','iqr_peak_ms', ...
         'tested_onset','tested_peak','sig_onset','sig_peak','panel_order'});
    writetable(Tc, fullfile(outDir,'event_latency_cells.csv'));
    fprintf('Wrote event_latency_cells.csv (%d active cells)\n', height(Tc));
    % Everything the downstream figures need, so nothing is ever recomputed with
    % different parameters than the pages were drawn with.
    prm = struct('nPerm',nPerm,'alphaPerm',alphaPerm,'minTestEv',minTestEv, ...
                 'activeMinEv',activeMinEv,'smoothDiv',smoothDiv,'caLagSec',caLagSec, ...
                 'minShiftIBI',minShiftIBI,'nBins',nBins);
    % OBS and REC go too: downstream figures (cycle dF/F heatmaps) need the
    % recording folder, the trigger times and the per-observation event indices,
    % and re-deriving them would mean duplicating the scan, the exclusions and
    % the cell-identity logic -- three places for them to drift apart.
    save(fullfile(outDir,'event_latency_data.mat'),'S','CELL','OBS','REC','groups','prm','-v7.3');
    fprintf('Wrote event_latency_data.mat\n');
end
fprintf('\nOutput: %s\n', outDir);
fprintf('==============================================================\n\n');
end

% =========================================================================
%                              HELPERS
% =========================================================================

function S = cell_latency(C, OBS, REC, nPerm, minShiftIBI, minTestEv, alphaPerm, smoothDiv)
% Latencies and the permutation test for ONE cell, pooled over its recordings.
%
% THE SHIFT IS DONE ON THE EVENTS, NOT THE TRIGGERS -- they are the same thing.
% Circularly shifting the trigger train by +v gives every event the same latency
% as shifting the events by -v against the unshifted triggers. Doing it on the
% events lets the trigger midpoints be built ONCE per observation, so all nPerm
% shifts collapse into a single discretize on an nEv x nPerm matrix.
S = struct('ibi',[NaN NaN],'fps',NaN,'L',{{[],[]}}, ...
           'mode',[NaN NaN],'fwhm',[NaN NaN],'frac',[NaN NaN], ...
           'dens',{{[],[]}},'dctr',{{[],[]}}, ...
           'med',[NaN NaN],'q1',[NaN NaN],'q3',[NaN NaN],'iqr',[NaN NaN], ...
           'p',[NaN NaN],'z',[NaN NaN],'nEv',[0 0],'nCyc',[0 0], ...
           'tested',[false false],'sig',[false false]);

% The frame lattice this cell's latencies live on. mode() because a cell can be
% re-imaged in sessions at different frame rates; almost every cell is single-fps,
% and the bin grid can only be snapped to one lattice.
S.fps = mode(arrayfun(@(j) REC(OBS(j).rec).fps, C.obs));

for q = 1:2
    % this CELL's own mean IBI, pooled over the recordings it appears in
    dt = []; nCyc = 0;
    for j = C.obs
        ts = REC(OBS(j).rec).trig{q};
        if numel(ts) > 1, dt = [dt; diff(ts)]; nCyc = nCyc + numel(ts) - 1; end %#ok<AGROW>
    end
    if isempty(dt), continue; end
    ibi = mean(dt);  half = ibi/2;
    S.ibi(q) = ibi;  S.nCyc(q) = nCyc;

    % ---- the grid the test lives on, built ONCE before the observation loop --
    % Fine bins on the frame lattice (one frame each -- the finest the data can
    % support), then a Gaussian of FIXED width. Fixed is essential: an adaptive
    % bandwidth would differ between the observed data and every shuffle, so the
    % null would not be comparable to the observation.
    dtm  = 1000/S.fps;                                   % ms per frame
    nf   = max(2, floor(1000*half/dtm));
    ctr  = (-nf:nf)*dtm;
    eg   = [ctr - dtm/2, ctr(end) + dtm/2];
    sg   = max(1000*ibi/smoothDiv, dtm);                 % smoothing SD, scales with cycle
    gx   = -3*sg : dtm : 3*sg;
    g    = exp(-0.5*(gx/sg).^2);  g = g/sum(g);

    % The null counts are ACCUMULATED per observation rather than storing every
    % shuffled latency: at nPerm = 10000 a 500-event cell would otherwise hold a
    % 500 x 10000 matrix. This keeps the footprint at (bins x nPerm).
    Lobs = [];  cnt = zeros(numel(ctr), nPerm);
    for j = C.obs
        r  = REC(OBS(j).rec);
        ts = r.trig{q};  Trec = r.Trec;
        tE = OBS(j).ev(:)/r.fps;
        tE = tE(tE >= 0 & tE <= Trec);
        if isempty(tE), continue; end

        te   = [ts(:).'-Trec, ts(:).', ts(:).'+Trec];
        edg  = [-inf, (te(1:end-1)+te(2:end))/2, inf];

        L0 = tE - te(discretize(tE, edg)).';
        L0(abs(L0) > half) = NaN;
        Lobs = [Lobs; 1000*L0]; %#ok<AGROW>

        % nPerm shifts, all at once. Shifts too close to no-shift are rejected:
        % half a breath is the same alignment displaced, not a new null.
        sh = draw_shifts(nPerm, Trec, minShiftIBI*ibi);
        M  = mod(tE - sh(:).', Trec);            % nEv x nPerm
        Ls = M - te(discretize(M, edg));
        Ls(abs(Ls) > half) = NaN;
        b  = discretize(1000*Ls, eg);
        ok = ~isnan(b);
        if any(ok(:))
            [~, cc] = ind2sub(size(b), find(ok));
            cnt = cnt + accumarray([b(ok), cc], 1, [numel(ctr) nPerm]);
        end
    end

    L = Lobs(~isnan(Lobs));
    S.L{q} = L;  S.nEv(q) = numel(L);
    if numel(L) < minTestEv, continue; end

    hO = conv(histcounts(L, eg), g, 'same');
    TO = max(hO) - mean(hO);
    HS = conv2(cnt, g(:), 'same');
    TN = max(HS,[],1) - mean(HS,1);

    S.p(q)      = (1 + nnz(TN >= TO)) / (1 + nPerm);
    % EFFECT SIZE, in null standard deviations. p SATURATES -- it cannot go
    % below 1/(nPerm+1), and measured on this archive 100% of the onset-
    % significant cells and 97% of the peak-significant ones sit exactly on that
    % floor, so -log10(p) cannot rank them at all. z keeps separating cells long
    % after p has bottomed out, and is monotonic with p so the ordering agrees
    % wherever p can still discriminate.
    S.z(q)      = (TO - mean(TN)) / max(std(TN), eps);
    S.tested(q) = true;
    S.sig(q)    = S.p(q) < alphaPerm;
    S.dens{q}   = hO / max(nCyc,1);      % epc per frame-bin, for drawing
    S.dctr{q}   = ctr;

    % ---- describe the PEAK, not the whole window ---------------------------
    % Computed for EVERY tested cell, not only the significant ones, so the
    % polar figure can place non-significant cells too. MODE and FWHM come off
    % the SAME smoothed curve the test scored, so the number describes exactly
    % what was detected. On ChAT/0521 r1 the onset panel is a 167 ms peak
    % holding 40% of the events: mode/FWHM read +467/167 ms, while a plain
    % median/IQR over the whole window reads +400/708 ms.
    [pk, im] = max(hO);
    base = median(hO);  hm = base + (pk - base)/2;
    lo = im; while lo > 1         && hO(lo-1) >= hm, lo = lo - 1; end
    hi = im; while hi < numel(hO) && hO(hi+1) >= hm, hi = hi + 1; end
    S.mode(q) = ctr(im);
    S.fwhm(q) = ctr(hi) - ctr(lo);
    S.frac(q) = nnz(L >= ctr(lo) & L <= ctr(hi)) / numel(L);

    % CIRCULAR median and quartiles, anchored on the mode. The window is one
    % whole breath cycle, so it WRAPS: a cell peaking near +IBI/2 has half its
    % events at -IBI/2, and a plain linear median would land on the opposite
    % side of the cycle from the peak. Re-expressing every latency relative to
    % the mode (wrapped into +/-IBI/2) removes the wrap, and the quartiles are
    % then mapped back. For a cell peaking mid-window this is identical to the
    % linear median; for one peaking at the edge it is the difference between
    % right and wrong.
    W  = 1000*ibi;
    d  = mod(L - S.mode(q) + W/2, W) - W/2;
    Qd = prctile(d, [25 50 75]);
    S.med(q) = S.mode(q) + Qd(2);
    S.q1(q)  = S.mode(q) + Qd(1);
    S.q3(q)  = S.mode(q) + Qd(3);
    S.iqr(q) = Qd(3) - Qd(1);
end
end

function sh = draw_shifts(nPerm, Trec, shMin)
% Uniform circular shifts, excluding a band of +/- shMin around no-shift (and
% around Trec, which wraps to the same alignment).
if 2*shMin >= Trec, sh = Trec*rand(nPerm,1); return; end
sh = shMin + (Trec - 2*shMin)*rand(nPerm,1);
end

function v = min_iqr(Sc)
% sort key: the TIGHTEST significant peak first (FWHM, not IQR -- the IQR
% describes the whole window and a sharp peak on a background has a wide one)
v = min(Sc.fwhm(Sc.sig));
if isempty(v) || isnan(v), v = Inf; end
end

function [edges, ctrs] = lat_axis(ibi, nBins, fps)
% BINS MUST BE A WHOLE NUMBER OF FRAMES WIDE.
%
% A latency is (event frame - trigger frame)/fps, so it can ONLY take values on
% an integer frame lattice: 0, +/-1/fps, +/-2/fps ... At 30 fps that is 33.3 ms.
% Binning at an arbitrary width (the old IBI/nBins) puts a different NUMBER OF
% LATTICE POINTS in each bin -- measured across this archive, every combination
% of fps and IBI gave bins holding [1 2] or [2 3] or [3 4] possible values. A
% bin that can hold 3 values collects ~50% more events than its neighbour that
% can hold 2, purely from the grid, so the bars alternate high-low and the
% histogram looks jagged. That sawtooth is the BIN GRID, not the neuron.
%
% So the bin width is snapped to k whole frames, with a bin CENTRED on zero
% lag. Every bin then spans exactly k lattice points and the sampling is
% uniform. k is chosen to land as close as possible to the requested IBI/nBins,
% so the resolution is unchanged -- only the alignment is.
if ~isfinite(ibi) || ibi <= 0, ibi = 1; end
if nargin < 3 || ~isfinite(fps) || fps <= 0, fps = 30; end
% k MUST BE ODD. With a bin centred on zero lag, an even k puts the bin edges at
% +/-(k/2)*dt -- an exact multiple of the frame period, i.e. ON a lattice point.
% Values landing exactly on an edge go to one side or the other by floating-point
% luck, which reintroduces the very unevenness this is meant to remove ([1 2 3]
% points per bin was measured at 42 and 47 fps). An odd k puts every edge at a
% half-integer multiple of dt, strictly BETWEEN two frames, so each bin holds
% exactly k lattice points and no value is ambiguous.
dt    = 1000/fps;                              % ms per frame
k     = max(1, round(1000*ibi/nBins/dt));      % whole frames per bin
if mod(k,2) == 0                               % snap to the nearer odd value
    if 1000*ibi/nBins/dt >= k, k = k + 1; else, k = max(1, k - 1); end
end
nb    = max(1, floor(((1000*ibi/2) - k*dt/2)/(k*dt)));   % whole bins only
ctrs  = (-nb:nb)*k*dt;
edges = [ctrs - k*dt/2, ctrs(end) + k*dt/2];
end

function H = pool_lat(Sset, edges)
% Pooled histogram per trigger, in epc: every cell's latencies counted, divided
% by the total breath cycles those cells were exposed to. A cell contributes in
% proportion to its events, which is the point of a population panel.
H = {zeros(1,numel(edges)-1), zeros(1,numel(edges)-1)};
if isempty(Sset), return; end
for q = 1:2
    cnt = zeros(1,numel(edges)-1);  cyc = 0;
    for k = 1:numel(Sset)
        if isempty(Sset(k).L{q}), continue; end
        cnt = cnt + histcounts(Sset(k).L{q}, edges);
        cyc = cyc + Sset(k).nCyc(q);
    end
    H{q} = cnt / max(cyc, 1);
end
end

function draw_cell(ax, Sc, q, nBins, col, isLast, showY, ttl)
hold(ax,'on'); box(ax,'on');
L = Sc.L{q};
if isempty(L) || ~isfinite(Sc.ibi(q))
    axis(ax,'off');
    if ~isempty(ttl), title(ax,ttl,'Interpreter','none','FontSize',6, ...
        'HorizontalAlignment','left','Units','normalized','Position',[0 1.03 0]); end
    return;
end
[edges, ctrs] = lat_axis(Sc.ibi(q), nBins, Sc.fps);
y   = histcounts(L, edges) / max(Sc.nCyc(q),1);      % epc
top = max([y, eps]);
bar(ax, ctrs, y, 1, 'FaceColor',col, 'EdgeColor','none', 'FaceAlpha',0.55);
if Sc.sig(q) && ~isempty(Sc.dens{q})
    % The curve drawn IS the curve the permutation test scored -- same fixed
    % bandwidth, same frame lattice -- rescaled from per-frame to the display
    % bin so it sits on the bars. Dashed line = the MODE, which is the number
    % reported, not the median.
    kb = diff(edges(1:2)) / (1000/Sc.fps);           % display bins per frame bin
    plot(ax, Sc.dctr{q}, Sc.dens{q}*kb, '-', 'Color',col, 'LineWidth',1.3);
    xline(ax, Sc.mode(q), '--', 'Color',col, 'LineWidth',1.1);
    top = max(top, max(Sc.dens{q}*kb));
end
xlim(ax, edges([1 end])); ylim(ax, [0 top*1.12]);
if isLast, xlabel(ax,'latency (ms)','FontSize',6.5); else, set(ax,'XTickLabel',[]); end
if showY, ylabel(ax,'epc','FontSize',6.5); end
if ~isempty(ttl)
    title(ax, ttl, 'Interpreter','none','FontSize',6,'FontWeight','normal', ...
          'HorizontalAlignment','left','Units','normalized','Position',[0 1.03 0]);
end
set(ax,'TickDir','out','FontSize',6);
end

function draw_pop(ax, ctrs, y, colBar, colEdge, name, isLast, showY, ttl)
hold(ax,'on'); box(ax,'on');
top = max([y, eps]);
bar(ax, ctrs, y, 1, 'FaceColor',colBar, 'EdgeColor','none');
plot(ax, ctrs, y, '-', 'Color',colEdge, 'LineWidth',1.1);
xlim(ax, [ctrs(1)-diff(ctrs(1:2))/2, ctrs(end)+diff(ctrs(1:2))/2]);
ylim(ax, [0 top*1.12]);
if isLast, xlabel(ax, sprintf('latency from %s (ms)', name), 'FontSize',6.5); end
if showY, ylabel(ax,'epc','FontSize',6.5); end
if ~isempty(ttl)
    title(ax, ttl, 'Interpreter','none','FontSize',6,'FontWeight','normal', ...
          'HorizontalAlignment','left','Units','normalized','Position',[0 1.03 0]);
end
set(ax,'TickDir','out','FontSize',6);
end

function hp = pop_figure(POP, ctrs, name, nAll, nSig, medIBI, alphaPerm, nPerm, colPop, colT, note)
hp = figure('Color','w','Visible','off','Units','centimeters','Position',[1 1 34 8.4]);
set(hp,'DefaultAxesFontSize',8);
colW = 0.288; colGap = 0.028; X0 = 0.050; subGap = 0.050; subW = (colW-subGap)/2;
for q = 1:3
    x0  = X0 + (q-1)*(colW+colGap);
    axA = axes('Parent',hp,'Position',[x0,             0.19, subW, 0.55]); %#ok<LAXES>
    axB = axes('Parent',hp,'Position',[x0+subW+subGap, 0.19, subW, 0.55]); %#ok<LAXES>
    if POP(q).n == 0
        axis(axA,'off'); axis(axB,'off');
        text(axA,0.5,0.5,sprintf('%s: no cells',POP(q).name), ...
             'Units','normalized','HorizontalAlignment','center','FontSize',9);
        continue;
    end
    ttl = sprintf('%s  |  %s  |  %d cells', name, POP(q).name, POP(q).n);
    draw_pop(axA, ctrs, POP(q).H{1}, colPop, colT(1,:), 'ONSET', true, q==1, ttl);
    draw_pop(axB, ctrs, POP(q).H{2}, colPop, colT(2,:), 'PEAK',  true, false, '');
end
axAll = findall(hp,'Type','axes'); set(axAll,'FontSize',8);
for a = axAll(:)'
    set(get(a,'Title'),'FontSize',8.5);
    set(get(a,'XLabel'),'FontSize',8.5); set(get(a,'YLabel'),'FontSize',8.5);
end
lines = {sprintf('%s POPULATION   |   %d cells (%d sig)   |   epc, window +/- half the median cycle (%.3f s)', ...
                 name, nAll, nSig, medIBI), ...
         'LEFT of each pair = INSPIRATION ONSET triggered,   RIGHT = INSPIRATORY PEAK triggered.   NEGATIVE latency LEADS the trigger.', ...
         sprintf('SIG = permutation p < %.4g, RAW (no family-wise correction), %d shifts, p floor %.5f.  NO GCaMP lag compensation.', ...
                 alphaPerm, nPerm, 1/(nPerm+1))};
if ~isempty(note), lines{end+1} = note; end
sgtitle(lines,'Interpreter','none','FontSize',9);
end

function s = fmt(Sc, q, minTestEv)
if ~Sc.tested(q),  s = sprintf('n<%d', minTestEv);
elseif ~Sc.sig(q), s = sprintf('n.s. p=%.3f', Sc.p(q));
else,              s = sprintf('%+.0f ms (w%.0f, %.0f%%) p=%.4f', ...
                               Sc.mode(q), Sc.fwhm(q), 100*Sc.frac(q), Sc.p(q));
end
end

function d = date_from_path(folderPath, groupRoot)
rel = regexprep(strrep(folderPath, groupRoot, ''), '^[\\/]+','');
parts = regexp(rel,'[\\/]','split');
if isempty(parts), d = ''; else, d = parts{1}; end
end

function tf = is_io_path(folderPath, groupRoot)
rel = regexprep(strrep(folderPath, groupRoot, ''), '^[\\/]+','');
parts = regexp(rel,'[\\/]','split');
tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2},'IO','once'));
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; rebuild the full last segment
% for folders with dots in the name (e.g. "1.7x", "15.5lp").
p = char(p);
while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
[~,n,e] = fileparts(p);  name = [n e];
end

function s = disp_label(x)
x = char(x);  p = regexp(x,'/','split');
if numel(p) >= 4
    rec = p{end-1};
    if numel(rec) > 24, rec = [rec(1:12) '..' rec(end-9:end)]; end
    s = sprintf('%s/%s/%s', p{end-2}, rec, p{end});
else
    s = x;
end
if numel(s) > 42, s = [s(1:20) '..' s(end-19:end)]; end
end
