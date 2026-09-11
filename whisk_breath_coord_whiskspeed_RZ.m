% whisk_breath_coord_whiskspeed_RZ  (script)
% Protraction-onset timing relative to BREATHING, split by WHISKING SPEED,
% pooled over all sessions (whisk epochs only), L & R onsets COMBINED.
% Each protraction onset is classified by its local whisk frequency (1/ITI):
%   SLOW : whisk freq < SLOW_HZ      FAST : whisk freq > FAST_HZ
% TWO figures (slow, fast), each 1x2:
%   LEFT  : peri-inspiration PSTH (onset - inspiration, +/-PERI_WIN)
%   RIGHT : breath-phase distribution (0=insp onset, 0.5=peak, 1=next onset)
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
SLOW_HZ   = 10;         % whisk freq (1/ITI) < this = SLOW
FAST_HZ   = 24;         % whisk freq (1/ITI) > this = FAST  (between = dropped)
WHISK_ITI_MAX = 0.25;   % s, drop whisk ITIs longer than this (cross-bout gaps)
CYC_MAX   = 1.0;        % s, drop breath cycles longer than this (gaps)
PERI_WIN  = 0.25;       % s, +/- window for the peri-inspiration PSTH
NBINS_T   = 51;         % peri-inspiration time-histogram bins
NBINS_PH  = 48;         % breath-phase bins
EXCLUDE   = "5840027";  % animal ids to drop ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

dt_slo=[]; dt_fas=[]; ph_slo=[]; ph_fas=[]; bcyc=[]; nsess=0;
ad = dir(char(dataRoot));
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    if any(strcmp(ad(ai).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        % ---- whisk protraction onsets (troughs), L & R, in epochs ----
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
        % per-onset local whisk frequency = 1/ITI (interval to next onset, same side)
        [twL,fwL] = onsetFreq((oL-1)/fpsW, WHISK_ITI_MAX);
        [twR,fwR] = onsetFreq((oR-1)/fpsW, WHISK_ITI_MAX);
        tw = [twL;twR]; fw = [fwL;fwR];
        if isempty(tw), continue; end
        % ---- breath inspiration onsets + peaks ----
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
        tPeak=zeros(numel(il)-1,1); for k=1:numel(il)-1, [~,rel]=max(brf(il(k):il(k+1))); tPeak(k)=(il(k)+rel-2)/fb; end
        bcyc=[bcyc; cycdur(cycdur<=CYC_MAX)]; %#ok<AGROW>

        % PSTH: onsets (by whisk speed) around all inspirations
        dt_slo=[dt_slo; psth(tw(fw<SLOW_HZ), tInsp, PERI_WIN)]; %#ok<AGROW>
        dt_fas=[dt_fas; psth(tw(fw>FAST_HZ), tInsp, PERI_WIN)]; %#ok<AGROW>
        % breath phase of each onset, classified by whisk speed
        [ph, fr] = coordPhaseF(tw, fw, tInsp, tPeak, CYC_MAX);
        ph_slo=[ph_slo; ph(fr<SLOW_HZ)]; ph_fas=[ph_fas; ph(fr>FAST_HZ)]; %#ok<AGROW>
        nsess = nsess + 1;
    end
end
fprintf('%d sessions | slow: %d pairs, %d phase | fast: %d pairs, %d phase\n', ...
    nsess, numel(dt_slo), numel(ph_slo), numel(dt_fas), numel(ph_fas));

plotState(dt_slo, ph_slo, bcyc, PERI_WIN, NBINS_T, NBINS_PH, [0.2 0.4 0.8], sprintf('SLOW whisk <%g Hz', SLOW_HZ));
plotState(dt_fas, ph_fas, bcyc, PERI_WIN, NBINS_T, NBINS_PH, [0.85 0.1 0.1], sprintf('FAST whisk >%g Hz', FAST_HZ));

% ================= plotting =================
function plotState(dt, ph, bcyc, PERI_WIN, NBINS_T, NBINS_PH, col, name)
    if isempty(ph) && isempty(dt), warning('%s: no data', name); return; end
    medCyc = median(bcyc); [mu,R,p] = circstat(2*pi*ph);
    figure('Color','w','Position',[120 160 1080 470],'Name',name,'NumberTitle','off');
    W=PERI_WIN; edges=linspace(-W,W,NBINS_T+1);
    ax1=subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
    histogram(ax1, dt, edges, 'Normalization','probability','FaceColor',col,'EdgeColor','none');
    xline(ax1,0,'k--','LineWidth',1.2);
    for jj=1:floor(W/medCyc), xline(ax1,jj*medCyc,'Color',[.6 .6 .6],'LineStyle',':'); xline(ax1,-jj*medCyc,'Color',[.6 .6 .6],'LineStyle',':'); end
    xlim(ax1,[-W W]); xlabel(ax1,'whisk protraction onset - inspiration (s)'); ylabel(ax1,'prob');
    title(ax1, sprintf('%s  peri-insp PSTH (n=%d)  breath %.0f ms (%.1f Hz)', name, numel(dt), 1000*medCyc, 1/medCyc));
    axis(ax1,'square');
    ax2=subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); pedges=linspace(0,1,NBINS_PH+1);
    histogram(ax2, ph, pedges, 'Normalization','probability','FaceColor',col,'EdgeColor','none');
    xline(ax2,0.5,'k:','LineWidth',1); xline(ax2,mu/(2*pi),'--','Color',col,'LineWidth',1.4);
    xlim(ax2,[0 1]); xlabel(ax2,'breath phase (0=insp onset, 0.5=peak, 1=next onset)'); ylabel(ax2,'prob');
    title(ax2, sprintf('%s  phase mean %.2f (Rvec %.2f, p=%.1g, n=%d)', name, mu/(2*pi), R, p, numel(ph)));
    axis(ax2,'square');
end

% ================= helpers =================
function [tw, fw] = onsetFreq(ton, itimax)
    % per-onset local whisk frequency = 1/(interval to next onset)
    tw = []; fw = [];
    if numel(ton)<2, return; end
    d = diff(ton); keep = d>0 & d<=itimax;
    tw = ton(1:end-1); tw = tw(keep); fw = 1./d(keep); tw=tw(:); fw=fw(:);
end
function [ph, fr] = coordPhaseF(tw, fw, tInsp, tPeak, cycmax)
    % piecewise breath phase (0=insp onset, 0.5=peak, 1=next onset), carry per-onset freq
    ph=[]; fr=[];
    if isempty(tw) || numel(tInsp)<2, return; end
    bn = discretize(tw, tInsp); ok = ~isnan(bn) & bn < numel(tInsp);
    bn=bn(ok); tw2=tw(ok); fw2=fw(ok);
    on0=tInsp(bn); on1=tInsp(bn+1); pk=tPeak(bn); cyc=on1-on0;
    phh=nan(size(tw2)); rise=tw2<pk;
    phh(rise)  = 0.5*(tw2(rise)-on0(rise))./(pk(rise)-on0(rise));
    phh(~rise) = 0.5 + 0.5*(tw2(~rise)-pk(~rise))./(on1(~rise)-pk(~rise));
    good = cyc>0 & cyc<=cycmax & pk>on0 & pk<on1 & phh>=0 & phh<1;
    ph = phh(good); ph = ph(:);
    fr = fw2(good); fr = fr(:);
end
function d = psth(tw, tInsp, win)
    d = []; if isempty(tw)||isempty(tInsp), return; end
    tw = tw(:);
    for i = 1:numel(tInsp)
        sel = tw>=tInsp(i)-win & tw<=tInsp(i)+win;
        if any(sel), d = [d; tw(sel)-tInsp(i)]; end %#ok<AGROW>
    end
    d = d(:);
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
