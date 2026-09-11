% whisk_breath_coord_pop_RZ  (script)
% POPULATION whisking-breathing coordination, all protraction onsets (L & R
% combined), pooled over all sessions (whisk epochs only). One figure (1x3):
%   1) DETECTION method on an example window: breath (onsets v / peaks ^) and
%      L/R whisk angle (protraction onsets o).
%   2) peri-inspiration PSTH of whisk onsets (grey) with the INSPIRATION-ONSET
%      autocorrelogram overlaid in BLUE (the breath rhythm reference).
%   3) breath-phase distribution of whisk onsets (0=insp onset,0.5=peak,1=next).
% Detection params match whisk_breath_coord_RZ.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BP        = [6 30];     % whisk bandpass (Hz)
PROT_PROM = 1;          % protraction-onset (trough) prominence (x std of BP whisk)
MINDIST   = 0.02;       % s, min distance between whisk onsets (findpeaks)
THR_FRAC  = 0.2; MIN_DUR = 0.5; MERGE_GAP = 0.10; ENV_SMOOTH = 0;   % whisk-epoch gating
BR_BP     = [2 15];     % breath bandpass (Hz) + inspiration detection band
BR_PROM   = 0.5;        % inspiration-onset (trough) prominence (x std of BP breath)
CYC_MAX   = 1.0;        % s, drop breath cycles longer than this (gaps)
PERI_WIN  = 0.50;       % s, +/- window for the peri-inspiration PSTH / ACG
NBINS_T   = 61;         % peri-inspiration time-histogram bins
NBINS_PH  = 48;         % breath-phase bins
EX_DUR    = 3;          % s, example window length for the detection panel
EXCLUDE   = "5840027";  % animal ids to drop ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
cL=[0 0.55 0]; cR=[0 0.4 0.85]; cB=[0.1 0.1 0.1]; cIN=[0.1 0.3 0.9];

dt_wh=[]; dt_in=[]; ph_wh=[]; bcyc=[]; nsess=0; nonset=0; ex=[];
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
        tw = sort([oL;oR] - 1)/fpsW;                         % all protraction onsets (L+R)
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
        tInsp=(il-1)/fb; cycdur=diff(tInsp);
        tPeak=zeros(numel(il)-1,1); ip=zeros(numel(il)-1,1);
        for k=1:numel(il)-1, [~,rel]=max(brf(il(k):il(k+1))); ip(k)=il(k)+rel-1; tPeak(k)=(ip(k)-1)/fb; end

        dt_wh = [dt_wh; psth(tw, tInsp, PERI_WIN)]; %#ok<AGROW>
        dt_in = [dt_in; psthACG(tInsp, PERI_WIN)]; %#ok<AGROW>
        ph_wh = [ph_wh; coordPhase(tw, tInsp, tPeak, CYC_MAX)]; %#ok<AGROW>
        bcyc  = [bcyc; cycdur(cycdur<=CYC_MAX)]; %#ok<AGROW>
        nonset = nonset + numel(tw); nsess = nsess + 1;
        if isempty(ex)                                       % keep first session as detection example
            ex.t=t; ex.xL=xL; ex.xR=xR; ex.oL=oL; ex.oR=oR;
            ex.tB=(0:numel(brf)-1)'/fb; ex.brf=brf; ex.il=il; ex.ip=ip; ex.name=sprintf('%s n%d',ad(ai).name,kk);
            ex.t0 = max(0, median(tw)-EX_DUR/2);
        end
    end
end
medCyc = median(bcyc); [mu,R,p] = circstat(2*pi*ph_wh);
fprintf('%d sessions, %d whisk onsets | breath %.0f ms (%.1f Hz) | phase mean %.2f, Rvec %.3f, p=%.2g\n', ...
    nsess, nonset, 1000*medCyc, 1/medCyc, mu/(2*pi), R, p);

% ============================ FIGURE (1x3) ============================
figure('Color','w','Position',[60 120 1500 460]);

% ---- 1) detection method on example window ----
ax0 = subplot(1,3,1); hold(ax0,'on'); grid(ax0,'on');
w = [ex.t0 ex.t0+EX_DUR];
yyaxis(ax0,'left');
plot(ax0, ex.t, ex.xL, '-', 'Color',cL, 'LineWidth',0.8);
plot(ax0, ex.t, ex.xR, '-', 'Color',cR, 'LineWidth',0.8);
plot(ax0, (ex.oL-1)/fpsW, ex.xL(ex.oL), 'o', 'Color',cL, 'MarkerSize',4, 'MarkerFaceColor','w');
plot(ax0, (ex.oR-1)/fpsW, ex.xR(ex.oR), 'o', 'Color',cR, 'MarkerSize',4, 'MarkerFaceColor','w');
ylabel(ax0,'whisk angle (deg)'); set(ax0,'YColor','k');
yyaxis(ax0,'right');
plot(ax0, ex.tB, ex.brf, '-', 'Color',cB, 'LineWidth',1.0);
plot(ax0, ex.tB(ex.il), ex.brf(ex.il), 'v', 'Color',[0.85 0.1 0.1],'MarkerFaceColor',[0.85 0.1 0.1],'MarkerSize',5);  % inspiration onset
plot(ax0, ex.tB(ex.ip), ex.brf(ex.ip), '^', 'Color',[0.9 0.4 0],  'MarkerFaceColor',[0.9 0.4 0],  'MarkerSize',5);   % breath peak
ylabel(ax0,'breath (BP)'); set(ax0,'YColor',cB);
xlim(ax0, w); xlabel(ax0,'time (s)');
title(ax0, sprintf('detection (example %s)', ex.name), 'Interpreter','none');

% ---- 2) peri-inspiration PSTH (whisk) + inspiration-onset ACG (blue) ----
W = PERI_WIN; edges = linspace(-W, W, NBINS_T+1);
ax1 = subplot(1,3,2); hold(ax1,'on'); grid(ax1,'on');
histogram(ax1, dt_wh, edges, 'Normalization','probability','FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
histogram(ax1, dt_in, edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor',cIN,'LineWidth',1.6);
xline(ax1,0,'k--','LineWidth',1.2);
xlim(ax1,[-W W]); xlabel(ax1,'time from inspiration onset (s)'); ylabel(ax1,'prob');
legend(ax1, {'whisk onsets','inspiration ACG'}, 'Location','northeast');
title(ax1, sprintf('peri-inspiration PSTH (n=%d, %d sess)  breath %.0f ms', nonset, nsess, 1000*medCyc));
axis(ax1,'square');

% ---- 3) breath-phase distribution of whisk onsets ----
ax2 = subplot(1,3,3); hold(ax2,'on'); grid(ax2,'on'); pedges=linspace(0,1,NBINS_PH+1);
histogram(ax2, ph_wh, pedges, 'Normalization','probability','FaceColor',[0.45 0.45 0.45],'EdgeColor','none');
xline(ax2,0.5,'k:','LineWidth',1); xline(ax2,mu/(2*pi),'r--','LineWidth',1.4);
xlim(ax2,[0 1]); xlabel(ax2,'breath phase (0=insp onset, 0.5=peak, 1=next onset)'); ylabel(ax2,'prob');
title(ax2, sprintf('whisk-onset breath phase   mean %.2f (Rvec %.2f, p=%.1g, n=%d)', mu/(2*pi), R, p, numel(ph_wh)));
axis(ax2,'square');

% ================= helpers =================
function ph = coordPhase(tw, tInsp, tPeak, cycmax)
    ph=[]; if isempty(tw) || numel(tInsp)<2, return; end
    bn = discretize(tw, tInsp); ok = ~isnan(bn) & bn < numel(tInsp);
    bn=bn(ok); tw=tw(ok);
    on0=tInsp(bn); on1=tInsp(bn+1); pk=tPeak(bn); cyc=on1-on0;
    phh=nan(size(tw)); rise=tw<pk;
    phh(rise)  = 0.5*(tw(rise)-on0(rise))./(pk(rise)-on0(rise));
    phh(~rise) = 0.5 + 0.5*(tw(~rise)-pk(~rise))./(on1(~rise)-pk(~rise));
    good = cyc>0 & cyc<=cycmax & pk>on0 & pk<on1 & phh>=0 & phh<1;
    ph = phh(good); ph=ph(:);
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
    % inspiration-onset autocorrelogram: other onsets within +/-win (excl. self)
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
