% perm_hist_sig_260903.m
% -----------------------------------------------------------------------
%  The PEAK-TRIGGERED event histogram of every SIGNIFICANT cell, drawn over the
%  permutation null it was tested against. This is the picture behind the p that
%  polar_recheck_perm_260901.m prints -- nothing new is measured here, the same
%  test is shown one cell at a time.
%
%  THE TRIGGER IS THE INSPIRATORY PEAK. Onset appears only as the red landmark at
%  -dOnPk (median onset-to-peak lag); it is not the trigger.
%
%  THE BAND DRAWN IS THE POINTWISE ONE: 2.5-97.5 percentile of the shuffled
%  histograms bin by bin, the same band the per-cell summary figure draws. It is
%  NOT the criterion the p came from. The test's own criterion is a SIMULTANEOUS
%  envelope, null mean +/- sCrit, sCrit being the (1-alpha) quantile of
%  max|null - nullMu| over the window; that envelope sits HIGHER than the
%  pointwise band, and it was drawn here until 2026-09-03, when it was taken off
%  at the user's request. So do not read "sticks out of the grey" as
%  "significant" -- with ~60 bins about 3 such excursions are expected by chance.
%  sCrit is still computed, and the console table prints maxDev / sCrit per cell.
%
%  THE X AXIS IS THE TESTED WINDOW AND NOTHING ELSE: |t| <= IBI/2, one full
%  breath cycle. The histogram is still BUILT over +/-1 IBI -- the null needs the
%  wider record -- but the wings are not drawn, because no part of them enters
%  the statistic and they are just the neighbouring cycles' peaks.
%
%  sCrit at alpha = 0.001 from 2000 shuffles is the 99.9th percentile of 2000
%  numbers, i.e. it is pinned by the top two draws and is correspondingly noisy.
%  The envelope's HEIGHT is a rough number; whether the observation clears it is
%  the same yes/no the p already gave.
%
%  p IS REPRODUCED, NOT RE-DRAWN. The loop walks all cells in the same order with
%  the same seed and the same randi() calls as the polar script, so the p
%  computed here is bit-identical to the CSV; the script asserts that before it
%  plots anything. Coherence is the only thing dropped -- it draws no random
%  numbers, so removing it cannot move the stream.
%
%  Runqi Zhang / 2026-09-03
% -----------------------------------------------------------------------
clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
sumRoot   = 'D:\Ventral_surface_summary';
bundleDir = fullfile(sumRoot,'per-cell-summary_active_260812','spike_recheck_260901');
outDir    = fullfile(sumRoot,'polar_recheck_260901');

nShuffle    = 2000;      % must match polar_recheck_perm_260901.m
shiftMinCyc = 3;
histBinFr   = 2;
pCrit       = 0.001;     % -log10(p) >= 3
rngSeed     = 260901;

TW              = 4;
minSpikes       = 2;
nBins           = 36;    % phase bins; only used to reproduce the polar script's
                         % "too thin" skip exactly, so the RNG stream matches
f_breath_search = [0.2 4];
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
rng(rngSeed);

CU = ca_recheck_load_curation_260901(bundleDir);
d  = dir(fullfile(bundleDir,'cells','*.mat'));
GC = genotype_colors_260817();

H = struct('stem',{},'group',{},'p',{},'ctrs',{},'obs',{},'nullMu',{}, ...
           'lo',{},'hi',{},'sCrit',{},'testHalf',{},'winSec',{}, ...
           'nTrig',{},'nSpk',{},'dOnPk',{});
tAll = tic;

for k = 1:numel(d)
    L = load(fullfile(d(k).folder, d(k).name),'cellInfo','RECc');
    Ci = L.cellInfo;  REC = L.RECc;

    Q = struct('T',{},'fps',{},'ev',{},'peak',{},'foot',{},'f_pk',{},'bidx',{});
    tossedWhole = false;  nEvC = 0;
    for i = 1:numel(REC)
        r = REC(i);
        g = CU.get(r.folder, r.roi);
        if g.found && g.cell_toss, tossedWhole = true; break; end
        if g.found && g.toss, continue; end
        ev = r.spike_idx(:);
        if g.found, ev = g.spike_idx(:); end

        T = r.T; fps = r.fps;
        ev = ev(ev >= 1 & ev <= T);
        nEvC = nEvC + numel(ev);

        phi   = piecewise_phase_local(r.peak, r.foot, T);
        phiW  = mod(phi, 2*pi);
        valid = ~isnan(phiW);
        if nnz(valid) < 10 || isempty(ev), continue; end
        edgesP = linspace(0, 2*pi, nBins+1);
        bidx = nan(T,1);
        bidx(valid) = min(discretize(phiW(valid), edgesP), nBins);

        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[0.05,min(15,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(r.breath(1:T), pB); Sb=Sb(:); fb=fb(:);
        mm = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        if ~any(mm), continue; end
        [~,rl]=max(Sb(mm)); ipk=find(mm,1)+rl-1; f_pk=fb(ipk);

        Q(end+1) = struct('T',T,'fps',fps,'ev',ev,'peak',r.peak, ...
                          'foot',r.foot,'f_pk',f_pk,'bidx',bidx); %#ok<SAGROW>
    end
    if tossedWhole || isempty(Q), continue; end

    nPool = 0;    % events with a defined phase bin -- polar_recheck_perm's aPool
    for i = 1:numel(Q)
        b = Q(i).bidx(Q(i).ev);
        nPool = nPool + nnz(~isnan(b));
    end
    if nPool < minSpikes, continue; end

    allIBI = [];
    for i = 1:numel(Q), allIBI = [allIBI; diff(Q(i).foot)/Q(i).fps]; end %#ok<AGROW>
    IBI = median(allIBI);
    if ~isfinite(IBI) || IBI <= 0, IBI = 1/max(median([Q.f_pk]),eps); end
    fpsRef  = median([Q.fps]);
    winSecH = IBI;  testHalf = IBI/2;
    binW    = histBinFr/fpsRef;
    Mb      = floor(winSecH/binW);
    ctrsC   = (-Mb:Mb)*binW;
    edgesC  = ((-Mb-0.5):(Mb+0.5))*binW;
    nB      = numel(ctrsC);

    cnt = zeros(1,nB); nTrig = 0; nullH = zeros(nShuffle, nB); okNull = true;
    dOn = [];
    for i = 1:numel(Q)
        T = Q(i).T; fps = Q(i).fps; ev = Q(i).ev;
        wH  = max(1, round(winSecH*fps));
        trg = Q(i).peak(Q(i).peak-wH>=1 & Q(i).peak+wH<=T);
        if isempty(trg), continue; end
        cnt   = cnt + trig_hist_local(ev, trg, wH, edgesC, fps);
        nTrig = nTrig + numel(trg);
        % median onset-to-peak lag, used only to place the red landmark
        ft = Q(i).foot(:);
        for tq = trg(:).'
            pr = ft(ft < tq);
            if ~isempty(pr), dOn(end+1) = (tq - pr(end))/fps; end %#ok<SAGROW>
        end
        if isempty(ev), continue; end
        minShift = max([wH+1, round(shiftMinCyc/max(Q(i).f_pk,eps)*fps), round(fps)]);
        if T - minShift <= minShift, okNull = false; continue; end

        m  = mod(double(ev(:)).' - double(trg(:)), T);
        c  = accumarray(m(:)+1, 1, [T 1]);
        lags   = (-wH:wH).';
        binLag = discretize(lags/fps, edgesC);
        okL    = ~isnan(binLag);
        cRev = c([1, T:-1:2]);
        Hs = zeros(T, nB);
        Ls = lags(okL);  Bs = binLag(okL);
        for q = 1:numel(Ls)
            Hs(:,Bs(q)) = Hs(:,Bs(q)) + circshift(cRev, Ls(q));
        end
        sh = randi([minShift, T-minShift], nShuffle, 1);
        nullH = nullH + Hs(sh+1, :);
    end

    p = NaN; spkH = []; nullMu = []; lo = []; hi = []; sCrit = NaN;
    if nTrig > 0 && okNull
        spkH   = 100*cnt/nTrig;
        nullH  = 100*nullH/nTrig;
        nullMu = mean(nullH,1);
        lo = prctile(nullH,2.5,1);  hi = prctile(nullH,97.5,1);
        tm = abs(ctrsC) <= testHalf;
        if any(tm)
            sObs  = max(abs(spkH(tm) - nullMu(tm)));
            sNull = max(abs(nullH(:,tm) - nullMu(tm)), [], 2);
            p     = (1 + nnz(sNull >= sObs))/(1 + numel(sNull));
            sCrit = prctile(sNull, 100*(1-pCrit));
        end
    end

    H(end+1) = struct('stem',Ci.stem,'group',Ci.group,'p',p,'ctrs',ctrsC, ...
        'obs',spkH,'nullMu',nullMu,'lo',lo,'hi',hi,'sCrit',sCrit, ...
        'testHalf',testHalf,'winSec',winSecH,'nTrig',nTrig,'nSpk',nEvC, ...
        'dOnPk',median(dOn)); %#ok<SAGROW>

    if mod(numel(H),50)==0
        fprintf('  %d cells done (%.1f min)\n', numel(H), toc(tAll)/60);
    end
end
fprintf('\n%d cells recomputed in %.1f min\n', numel(H), toc(tAll)/60);

%% ---- the p must be identical to the polar CSV, else this is not that test ----
csv = readtable(fullfile(outDir,'polar_recheck_perm_percell.csv'),'TextType','string');
[tf, loc] = ismember(string({H.stem})', csv.stem);
assert(all(tf), 'a recomputed cell is missing from the CSV');
dP = abs([H.p]' - csv.p_perm_peak(loc));
fprintf('max |p_here - p_csv| = %.3g over %d cells\n', max(dP), numel(H));
assert(max(dP) < 1e-12, 'the permutation did not reproduce -- do not trust these bands');

%% ===================== SELECT AND ORDER =====================
sub = strings(numel(H),1);
for q = 1:numel(H)
    parts = split(string(H(q).stem),'_');  dt = parts(2);
    switch H(q).group
        case 'Vglut2'
            if dt=="0224",     sub(q) = "intact FAST";
            elseif dt=="0824", sub(q) = "vagotomised";
            else,              sub(q) = "intact rest"; end
        case 'Sert'
            if any(dt==["0828","0829"]), sub(q) = "vagotomised"; else, sub(q) = "intact"; end
        case 'ChAT'
            if dt=="0826", sub(q) = "vagotomised"; else, sub(q) = "intact"; end
        otherwise
            sub(q) = "";
    end
end

sig = [H.p]' <= pCrit;
ORD = ["Vglut2|intact FAST","Vglut2|intact rest","Vglut2|vagotomised", ...
       "Sert|intact","Sert|vagotomised","ChAT|intact","ChAT|vagotomised", ...
       "IO|","Sst|","Vgat|"];
key  = string({H.group})' + "|" + sub;
rank = zeros(numel(H),1);
for q = 1:numel(H)
    ii = find(ORD==key(q),1);
    if isempty(ii), rank(q) = 99; else, rank(q) = ii; end
end

PAGES = { 'Vglut2', "Vglut2"
          'Sert',   "Sert"
          'other',  ["ChAT","Vgat","Sst","IO"] };

for pg = 1:size(PAGES,1)
    sel = find(sig & ismember(string({H.group})', PAGES{pg,2}));
    if isempty(sel), continue; end
    [~,o] = sortrows([rank(sel), (1:numel(sel))']);  sel = sel(o);
    n = numel(sel);
    if     n <= 8,  nCols = 4;
    elseif n <= 18, nCols = 5;
    else,           nCols = 6;
    end
    nRows = ceil(n/nCols);

    figW = min(6.3*nCols + 2.0, 40);
    figH = min(4.6*nRows + 3.6, 24);
    fig = figure('Color','w','Units','centimeters','Position',[1 1 figW figH], ...
                 'Name',sprintf('peak-triggered histogram vs null: %s', PAGES{pg,1}));
    set(fig,'DefaultAxesFontSize',7.5);

    xL = 0.055; xR = 0.995; yB = 0.055; yT = 0.885;
    xw = (xR-xL)/nCols;  yh = (yT-yB)/nRows;
    for q = 1:n
        h = H(sel(q));
        col = GC.(h.group);  if strcmp(h.group,'IO'), col = [0.5 0.5 0.5]; end
        rr = ceil(q/nCols);  cc = q - (rr-1)*nCols;
        ax = axes(fig,'Position',[xL+(cc-1)*xw, yT-rr*yh+0.16*yh, xw*0.80, yh*0.58]); %#ok<LAXES>
        hold(ax,'on');

        % Everything is cropped to the TESTED window, |t| <= IBI/2, which spans
        % exactly one breath cycle. Nothing outside it is drawn, because nothing
        % outside it enters the statistic; the +/-1 IBI the histogram was built
        % over only existed so the edge bins had neighbours.
        tm = abs(h.ctrs) <= h.testHalf;
        ct = h.ctrs(tm);
        fill(ax,[ct fliplr(ct)],[h.hi(tm) fliplr(h.lo(tm))],[0.6 0.6 0.6], ...
             'FaceAlpha',0.30,'EdgeColor','none');
        bar(ax, ct, h.obs(tm), 1, 'FaceColor','k','FaceAlpha',0.85,'EdgeColor','none');

        yMax = max([h.obs(tm), h.hi(tm), 1])*1.12;
        xline(ax, 0, 'Color',[0.10 0.35 0.85],'LineWidth',1);
        if isfinite(h.dOnPk) && h.dOnPk <= h.testHalf
            xline(ax, -h.dOnPk, '-','Color',[0.90 0.15 0.10],'LineWidth',0.9);
        end
        xlim(ax,[-h.testHalf h.testHalf]); ylim(ax,[0 yMax]); box(ax,'on');
        if cc == 1, ylabel(ax,'spk/cyc %'); end
        if rr == nRows, xlabel(ax,'time from insp PEAK (s)'); end
        ttl = strrep(h.stem,[h.group '_'],'');
        title(ax, {sprintf('%s  %s', ttl, sub(sel(q))), ...
                   sprintf('p = %.4f   %d spk / %d cyc', h.p, h.nSpk, h.nTrig)}, ...
              'FontSize',7,'FontWeight','normal','Interpreter','none');
    end

    annotation(fig,'textbox',[0.02 0.945 0.96 0.045],'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontWeight','bold','FontSize',10,'Interpreter','none', ...
        'String',sprintf(['%s -- SIGNIFICANT cells (n = %d), inspiratory PEAK-triggered ' ...
                          'event histogram vs its circular-shift null'], PAGES{pg,1}, n));
    annotation(fig,'textbox',[0.02 0.900 0.96 0.042],'EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',8,'Interpreter','tex', ...
        'String',{sprintf(['grey band = the permutation null, pointwise 2.5-97.5%% of %d circular-shift shuffles'], nShuffle), ...
                  ['x axis = the TESTED window only, +/- IBI/2, i.e. ONE breath cycle     ' ...
                   'blue = peak (trigger),  red = median inspiration onset']});

    stem = fullfile(outDir, sprintf('perm_hist_sig_%s', PAGES{pg,1}));
    exportgraphics(fig, [stem '.png'], 'Resolution',300, 'BackgroundColor','white');
    exportgraphics(fig, [stem '.pdf'], 'ContentType','vector','BackgroundColor','white');
    fi = dir([stem '.png']);
    fprintf('%-8s %2d cells -> %s.png (%.0f kB)\n', PAGES{pg,1}, n, stem, fi.bytes/1024);
end

%% ---- how far outside the envelope each significant cell sits ----
fprintf('\n%-24s %-22s %8s %8s %7s %6s\n', ...
        'cell','group','p','maxDev','sCrit','ratio');
selAll = find(sig);  [~,o] = sortrows([rank(selAll), (1:numel(selAll))']);
selAll = selAll(o);
for q = selAll(:).'
    h = H(q);  tm = abs(h.ctrs) <= h.testHalf;
    md = max(abs(h.obs(tm) - h.nullMu(tm)));
    fprintf('%-24s %-22s %8.4f %8.2f %7.2f %6.2f\n', h.stem, ...
            sprintf('%s %s',h.group,sub(q)), h.p, md, h.sCrit, md/h.sCrit);
end

%% ===================== LOCAL FUNCTIONS =====================
function h = trig_hist_local(evIdx, trigIdx, winH, edgesC, fps)
h = zeros(1, numel(edgesC)-1);
if isempty(evIdx) || isempty(trigIdx), return; end
d = double(evIdx(:)).' - double(trigIdx(:));
d = d(abs(d) <= winH);
if isempty(d), return; end
h = histcounts(d(:)/fps, edgesC);
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
[events, ord] = sort(events); types = types(ord);
keep = true(size(events));
for i = 2:numel(events), if types(i) == types(i-1), keep(i) = false; end, end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases = nan(size(events)); phi_cur = types(1) * pi;
for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
