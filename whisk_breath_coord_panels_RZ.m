% whisk_breath_coord_panels_RZ  (script)
% Protraction-onset timing relative to BREATHING across 5 conditions, pooled
% over all sessions (whisk epochs only), L & R onsets COMBINED.
% Rows (each a PAIR: peri-inspiration PSTH | breath-phase distribution):
%   1 population   2 fast whisk   3 slow whisk   4 fast breathing   5 slow breathing
% Classification per onset:  whisk freq = 1/ITI ; breath = rate of its cycle.
% The inspiration-onset autocorrelogram (blue) is overlaid ONLY on the
% population PSTH (panel 1), on a second y-axis (yyaxis).
% Phase: 0 = insp onset, 0.5 = breath peak, 1 = next onset.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % protraction-onset (trough) prominence (x std of BP whisk)
MINDIST   = 0.02;       % s, min distance between whisk onsets (findpeaks)
THR_FRAC  = 0; MIN_DUR = 0.5; MERGE_GAP = 0.10; ENV_SMOOTH = 0;   % whisk-epoch gating
BR_BP     = [2 15];     % breath bandpass (Hz) + inspiration detection band
BR_PROM   = 0.5;        % inspiration-onset (trough) prominence (x std of BP breath)
SLOW_HZ   = 10;         % whisk freq < this = SLOW whisk
FAST_HZ   = 20;         % whisk freq > this = FAST whisk
BASAL_HZ  = 5;          % breath rate < this = SLOW breathing
SNIFF_HZ  = 8;          % breath rate > this = FAST breathing
WHISK_ITI_MAX = 0.25;   % s, drop whisk ITIs longer than this
CYC_MAX   = 0.5;        % s, drop breath cycles longer than this
PERI_WIN  = 0.50;       % s, +/- pairing window for PSTH/ACG (bin-width source; must be >= every WIN_T*)
% per-panel DISPLAY x-window (+/- s) for the TIME panels, laid out as the 3x3 grid:
%             col1          col2 (slow breath)   col3 (fast breath)
%   row1:   population      slow breathing       fast breathing
%   row2:   slow whisk      sBsW                 fBsW
%   row3:   fast whisk      sBfW                 fBfW
gridWin   = [ 0.50 0.50 0.20 ;
              0.20 0.20 0.20 ;
              0.1 0.1 0.1 ];
BIN_W     = 0.0050;     % s, FIXED PSTH bin width (consistent across panels); 1/fpsW = 1 frame = 2.5 ms
NBINS_PH  = 31;         % breath-phase bins
EXCLUDE   = "5840027";  % animal ids to drop ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
PAIR_WIN = max([PERI_WIN, gridWin(:)']);   % collect pairs over the widest display window

% ---------------- categorical state palette ----------------
hex2rgb = @(h) reshape(sscanf(char(erase(h,"#")),'%2x')/255,1,3);
% Single-condition colors
C.pop   = hex2rgb("#333333");  % population / reference: dark gray
C.bSlow = hex2rgb("#0B2C6B");  % slow breathing: navy
C.bFast = hex2rgb("#F04A24");  % fast breathing / sniffing: orange-red
C.wSlow = hex2rgb("#7B3FB3");  % slow whisking: violet
C.wFast = hex2rgb("#009B7A");  % fast whisking: emerald
% Joint breath x whisk states
C.sBsW  = hex2rgb("#20B7E5");  % slow breath + slow whisk: cyan
C.sBfW  = hex2rgb("#0052CC");  % slow breath + fast whisk: royal blue
C.fBsW  = hex2rgb("#B11264");  % fast breath + slow whisk: burgundy-magenta
C.fBfW  = hex2rgb("#F2B93B");  % fast breath + fast whisk: sand-yellow
% Reference/annotation colors
C.insp  = hex2rgb("#333333");  % inspiration reference / ACG overlay
FALPHA  = 0.70;                % histogram fill alpha

D = struct('pop',[],'fastW',[],'slowW',[],'fastB',[],'slowB',[], 'sBsW',[],'sBfW',[],'fBsW',[],'fBfW',[]);   % PSTH lags
P = struct('pop',[],'fastW',[],'slowW',[],'fastB',[],'slowB',[], 'sBsW',[],'sBfW',[],'fBsW',[],'fBfW',[]);   % phases
O = struct('pop',[],'fastW',[],'slowW',[],'fastB',[],'slowB',[], 'sBsW',[],'sBfW',[],'fBsW',[],'fBfW',[]);   % phase OCCUPANCY (all frames)
dt_in=[]; bcyc=[]; nsess=0;
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    if any(strcmp(ad(ai).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(ai).name,kk)), 0.6);
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = movmean((abs(hilbert(xL))+abs(hilbert(xR)))/2, max(1,round(ENV_SMOOTH*fpsW)));
        ep  = detect(env, t, THR_FRAC, MIN_DUR, MERGE_GAP);
        m   = false(size(t)); for q=1:size(ep,1), m=m|(t>=ep(q,1)&t<=ep(q,2)); end
        [~,oL] = findpeaks(-xL,'MinPeakProminence',PROT_PROM*std(xL),'MinPeakDistance',round(MINDIST*fpsW)); oL=oL(m(oL));
        [~,oR] = findpeaks(-xR,'MinPeakProminence',PROT_PROM*std(xR),'MinPeakDistance',round(MINDIST*fpsW)); oR=oR(m(oR));
        [twL,fwL] = onsetFreq((subsamp(-xL,oL)-1)/fpsW, WHISK_ITI_MAX);   % sub-frame onset times
        [twR,fwR] = onsetFreq((subsamp(-xR,oR)-1)/fpsW, WHISK_ITI_MAX);
        tw=[twL;twR]; fw=[fwL;fwR];
        if isempty(tw), continue; end
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(ai).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
        catch
            continue
        end
        [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        [~,il]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        if numel(il) < 2, continue; end
        tInsp=(subsamp(-brf,il)-1)/fb; cycdur=diff(tInsp);   % sub-frame inspiration onsets
        tPeak=zeros(numel(il)-1,1); for k=1:numel(il)-1, [~,rel]=max(brf(il(k):il(k+1))); tPeak(k)=(subsamp(brf,il(k)+rel-1)-1)/fb; end

        % per-onset breath rate (cycle containing the onset)
        bn=discretize(tw,tInsp); ok=~isnan(bn)&bn<numel(tInsp);
        cyc_tw=nan(size(tw)); cyc_tw(ok)=cycdur(bn(ok)); cyc_tw(cyc_tw>CYC_MAX)=NaN; brate=1./cyc_tw;
        % phase + aligned whisk-freq / breath-rate
        [ph,fwk,brk] = coordPhaseFB(tw, fw, tInsp, tPeak, CYC_MAX);

        % ---- accumulate PSTH (time) per condition ----
        D.pop  =[D.pop;  psth(tw,                 tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.fastW=[D.fastW;psth(tw(fw>FAST_HZ),     tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.slowW=[D.slowW;psth(tw(fw<SLOW_HZ),     tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.fastB=[D.fastB;psth(tw(brate>SNIFF_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.slowB=[D.slowB;psth(tw(brate<BASAL_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        % ---- accumulate phase per condition ----
        P.pop  =[P.pop;  ph]; %#ok<AGROW>
        P.fastW=[P.fastW;ph(fwk>FAST_HZ)]; %#ok<AGROW>
        P.slowW=[P.slowW;ph(fwk<SLOW_HZ)]; %#ok<AGROW>
        P.fastB=[P.fastB;ph(brk>SNIFF_HZ)]; %#ok<AGROW>
        P.slowB=[P.slowB;ph(brk<BASAL_HZ)]; %#ok<AGROW>
        % ---- joint breath-state x whisk-speed (figure 2) ----
        D.sBsW=[D.sBsW; psth(tw(brate<BASAL_HZ & fw<SLOW_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.sBfW=[D.sBfW; psth(tw(brate<BASAL_HZ & fw>FAST_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.fBsW=[D.fBsW; psth(tw(brate>SNIFF_HZ & fw<SLOW_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        D.fBfW=[D.fBfW; psth(tw(brate>SNIFF_HZ & fw>FAST_HZ), tInsp, PAIR_WIN)]; %#ok<AGROW>
        P.sBsW=[P.sBsW; ph(brk<BASAL_HZ & fwk<SLOW_HZ)]; %#ok<AGROW>
        P.sBfW=[P.sBfW; ph(brk<BASAL_HZ & fwk>FAST_HZ)]; %#ok<AGROW>
        P.fBsW=[P.fBsW; ph(brk>SNIFF_HZ & fwk<SLOW_HZ)]; %#ok<AGROW>
        P.fBfW=[P.fBfW; ph(brk>SNIFF_HZ & fwk>FAST_HZ)]; %#ok<AGROW>
        % ---- phase OCCUPANCY: breath phase of ALL whisk-epoch frames, same phase def ----
        tAll = t(m);                                     % whisk-epoch frame times (fpsW grid)
        [twU,iu] = unique(tw);                           % per-frame whisk freq from onset freqs
        if numel(twU) >= 2
            fwAll = interp1(twU, fw(iu), tAll, 'nearest', 'extrap');   % same freq estimator as events
            [phO,fwO,brO] = coordPhaseFB(tAll, fwAll, tInsp, tPeak, CYC_MAX);
            O.pop  =[O.pop;  phO]; %#ok<AGROW>
            O.fastW=[O.fastW;phO(fwO>FAST_HZ)]; %#ok<AGROW>
            O.slowW=[O.slowW;phO(fwO<SLOW_HZ)]; %#ok<AGROW>
            O.fastB=[O.fastB;phO(brO>SNIFF_HZ)]; %#ok<AGROW>
            O.slowB=[O.slowB;phO(brO<BASAL_HZ)]; %#ok<AGROW>
            O.sBsW=[O.sBsW; phO(brO<BASAL_HZ & fwO<SLOW_HZ)]; %#ok<AGROW>
            O.sBfW=[O.sBfW; phO(brO<BASAL_HZ & fwO>FAST_HZ)]; %#ok<AGROW>
            O.fBsW=[O.fBsW; phO(brO>SNIFF_HZ & fwO<SLOW_HZ)]; %#ok<AGROW>
            O.fBfW=[O.fBfW; phO(brO>SNIFF_HZ & fwO>FAST_HZ)]; %#ok<AGROW>
        end
        dt_in=[dt_in; psthACG(tInsp, PAIR_WIN)]; %#ok<AGROW>
        bcyc =[bcyc; cycdur(cycdur<=CYC_MAX)]; %#ok<AGROW>
        nsess = nsess + 1;
    end
end
medCyc = median(bcyc);
fprintf('%d sessions | breath %.0f ms (%.1f Hz)\n', nsess, 1000*medCyc, 1/medCyc);

% ============================ ONE WIDE FIGURE: TIME LEFT, PHASE RIGHT ============================
gridKeys = {
    'pop',   'slowB', 'fastB';
    'slowW', 'sBsW',  'fBsW';
    'fastW', 'sBfW',  'fBfW'
};

gridLabs = {
    'population', ...
    sprintf('slow breathing <%g Hz',BASAL_HZ), ...
    sprintf('fast breathing >%g Hz',SNIFF_HZ);

    sprintf('slow whisk <%g Hz',SLOW_HZ), ...
    sprintf('slow breath <%g + slow whisk <%g Hz',BASAL_HZ,SLOW_HZ), ...
    sprintf('fast breath >%g + slow whisk <%g Hz',SNIFF_HZ,SLOW_HZ);

    sprintf('fast whisk >%g Hz',FAST_HZ), ...
    sprintf('slow breath <%g + fast whisk >%g Hz',BASAL_HZ,FAST_HZ), ...
    sprintf('fast breath >%g + fast whisk >%g Hz',SNIFF_HZ,FAST_HZ)
};

gridCols = {
    C.pop,   C.bSlow, C.bFast;
    C.wSlow, C.sBsW,  C.fBsW;
    C.wFast, C.sBfW,  C.fBfW
};

pedges = linspace(0,1,NBINS_PH+1);

fig = figure('Color','w','Position',[30 80 1900 650]);

annotation(fig,'textbox',[0.055 0.935 0.39 0.045], ...
    'String',sprintf('peri-inspiration whisk onset timing  (%d sessions)',nsess), ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

annotation(fig,'textbox',[0.565 0.935 0.39 0.045], ...
    'String','breath phase of whisk protraction onsets', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontWeight','bold','FontSize',12);

for r = 1:3
    for c = 1:3
        key = gridKeys{r,c};
        lab = gridLabs{r,c};
        col = gridCols{r,c};
        W   = gridWin(r,c);

        % ---------- LEFT BLOCK: peri-inspiration time histogram ----------
        axT = axes('Position',panelpos(1,r,c)); hold(axT,'on'); grid(axT,'on');
        dt = D.(key);
        edges = binedges(W, BIN_W);

        if strcmp(key,'pop')
            yyaxis(axT,'left');
            histogram(axT,dt,edges,'Normalization','probability', ...
                'FaceColor',col,'FaceAlpha',FALPHA,'EdgeColor','none');
            ylabel(axT,'whisk prob');
            set(axT,'YColor','k');

            yyaxis(axT,'right');
            histogram(axT,dt_in,edges,'Normalization','probability', ...
                'DisplayStyle','stairs','EdgeColor',C.insp,'LineWidth',1.4);
            ylabel(axT,'insp ACG');
            set(axT,'YColor',C.insp);
        else
            histogram(axT,dt,edges,'Normalization','probability', ...
                'FaceColor',col,'FaceAlpha',FALPHA,'EdgeColor','none');
            ylabel(axT,'prob');
        end

        xline(axT,0,'k--','LineWidth',1);
        xlim(axT,W*[-1 1]);
        axis(axT,'square');
        title(axT,sprintf('%s  (n=%d)',lab,numel(dt)), ...
            'Color',col,'FontWeight','bold','FontSize',8,'Interpreter','none');

        if r == 3
            xlabel(axT,'time from inspiration onset (s)');
        end

        % ---------- RIGHT BLOCK: breath phase histogram ----------
        axP = axes('Position',panelpos(2,r,c)); hold(axP,'on'); grid(axP,'on');
        ph = P.(key);

        hEvt = histogram(axP,ph,pedges,'Normalization','probability', ...
            'FaceColor',col,'FaceAlpha',FALPHA,'EdgeColor','none');

        % phase occupancy (all frames in this state) as a thin black line
        octr = (pedges(1:end-1)+pedges(2:end))/2;
        occ  = histcounts(O.(key), pedges, 'Normalization','probability');
        hOcc = plot(axP, octr, occ, 'k-', 'LineWidth',1.4);

        [mu,Rv,pr] = circstat(2*pi*ph);
        meanPhase = mod(mu/(2*pi),1);

        xline(axP,0.5,'k:','LineWidth',1);      % breath peak
        xline(axP,meanPhase,'--','Color',0.65*col,'LineWidth',1.6);
        if r==1 && c==1
            legend(axP,[hEvt hOcc],{'whisk onset phase','phase occupancy'}, ...
                'Location','northwest','Box','off','FontSize',7);
        end

        xlim(axP,[0 1]);
        axis(axP,'square');
        ylabel(axP,'prob');

        title(axP,sprintf('%s  phase %.2f (R %.2f, p=%.1g, n=%d)', ...
            lab, meanPhase, Rv, pr, numel(ph)), ...
            'Color',col,'FontWeight','bold','FontSize',8,'Interpreter','none');

        if r == 3
            xlabel(axP,'breath phase (0=insp, 0.5=peak, 1=next)');
        end
    end
end

sgtitle(fig,sprintf('whisk protraction onset vs breathing  (%d sessions, breath %.0f ms)', ...
    nsess,1000*medCyc),'FontSize',14,'FontWeight','bold');

% ================= helpers =================
function e = binedges(W, bw)
    % symmetric edges of fixed width bw, with one bin CENTERED on 0, covering +/-W
    n = max(1, round(W/bw)); e = (-(n+0.5):(n+0.5)) * bw;
end
function [tw, fw] = onsetFreq(ton, itimax)
    % CENTERED local whisk frequency: 1/(mean of the intervals before & after the
    % onset) = the cycle straddling the event (matches the bracketing breath cycle)
    tw=[]; fw=[]; if numel(ton)<3, return; end
    d=diff(ton); dprev=d(1:end-1); dnext=d(2:end);          % adjacent intervals of interior onsets
    per=(dprev+dnext)/2;                                    % centered cycle period
    keep = dprev>0 & dnext>0 & dprev<=itimax & dnext<=itimax;
    tw=ton(2:end-1); tw=tw(keep); fw=1./per(keep); tw=tw(:); fw=fw(:);
end
function [ph, fwk, brk] = coordPhaseFB(tw, fw, tInsp, tPeak, cycmax)
    ph=[]; fwk=[]; brk=[];
    if isempty(tw) || numel(tInsp)<2, return; end
    bn=discretize(tw,tInsp); ok=~isnan(bn)&bn<numel(tInsp);
    bn=bn(ok); tw2=tw(ok); fw2=fw(ok);
    on0=tInsp(bn); on1=tInsp(bn+1); pk=tPeak(bn); cyc=on1-on0;
    phh=nan(size(tw2)); rise=tw2<pk;
    phh(rise)  = 0.5*(tw2(rise)-on0(rise))./(pk(rise)-on0(rise));
    phh(~rise) = 0.5 + 0.5*(tw2(~rise)-pk(~rise))./(on1(~rise)-pk(~rise));
    good=cyc>0 & cyc<=cycmax & pk>on0 & pk<on1 & phh>=0 & phh<1;
    ph=phh(good); fwk=fw2(good); brk=1./cyc(good); ph=ph(:); fwk=fwk(:); brk=brk(:);
end
function d = psth(tw, tInsp, win)
    d=[]; if isempty(tw)||isempty(tInsp), return; end
    tw=tw(:);
    for i=1:numel(tInsp)
        sel = tw>=tInsp(i)-win & tw<=tInsp(i)+win;
        if any(sel), d=[d; tw(sel)-tInsp(i)]; end %#ok<AGROW>
    end
    d=d(:);
end
function d = psthACG(tev, win)
    d=[]; if numel(tev)<2, return; end
    tev=tev(:);
    for i=1:numel(tev)
        dd = tev - tev(i); dd = dd(dd>=-win & dd<=win & abs(dd)>1e-9);
        d=[d; dd]; %#ok<AGROW>
    end
    d=d(:);
end
function [mu, R, pRay] = circstat(phr)
    if isempty(phr), mu=NaN; R=0; pRay=NaN; return; end
    N=numel(phr); C=mean(cos(phr)); S=mean(sin(phr)); R=hypot(C,S); mu=mod(atan2(S,C),2*pi);
    Z=N*R^2; pRay=exp(-Z)*(1+(2*Z-Z^2)/(4*N));
end
function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac*prctile(env,95);
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1;
    ep = [t(s) t(e)];
    if ~isempty(ep)
        m = ep(1,:);
        for i=2:size(ep,1)
            if ep(i,1)-m(end,2) <= mergeGap, m(end,2)=ep(i,2); else, m(end+1,:)=ep(i,:); end %#ok<AGROW>
        end
        ep = m(m(:,2)-m(:,1) >= minDur, :);
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
function p = subsamp(x, idx)
    % parabolic sub-sample refinement of peak locations (removes frame quantization)
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx>1 & idx<numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp;                 % <0 at a true maximum
    delta = 0.5*(ym - yp) ./ den;
    delta(~isfinite(delta) | abs(delta)>0.5) = 0;   % guard flat / ill-conditioned
    p(in) = i + delta;
end
function pos = panelpos(block,row,col)
% block 1 = left time block
% block 2 = right phase block
% row/col are 1-based in the 3x3 grid.

    panelW = 0.105;
    panelH = 0.210;
    colGap = 0.038;
    rowGap = 0.075;

    xStart = [0.055, 0.565];
    yTop   = 0.680;

    x = xStart(block) + (col-1)*(panelW + colGap);
    y = yTop - (row-1)*(panelH + rowGap);

    pos = [x y panelW panelH];
end
