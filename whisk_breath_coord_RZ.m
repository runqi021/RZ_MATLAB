% whisk_breath_coord_RZ  (script)
% Coordination of whisking with breathing, pooled over all sessions (whisk
% epochs only). For every whisker PROTRACTION ONSET (trough of BP angle, L & R
% pooled) we measure its timing relative to the breathing cycle defined by
% INSPIRATION ONSETS (troughs of the BP nose-thermal breath signal).
%   LEFT  : peri-inspiration PSTH (whisk onset - inspiration, all pairs within +/-PERI_WIN)
%   RIGHT : breath-PHASE distribution (0 = insp onset, 0.5 = breath peak, 1 = next
%           onset; piecewise-linear), linear axis, with circular mean + Rayleigh.
% Parameters match whisk_LR_sync_RZ (whisk side) + breath_freq_whisk_psd_RZ (breath).

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
PERI_WIN  = 0.2;        % s, +/- window for the peri-inspiration PSTH (play with this)
NBINS_T   = 41;         % peri-inspiration time-histogram bins (over +/-PERI_WIN)
NBINS_PH  = 48;         % breath-phase bins
EXCLUDE   = "5840027";  % animal ids to drop from the analysis ("0027")
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
cL=[0 0.55 0]; cR=[0 0.4 0.85];   % L green / R blue (match video trace)

dt_L=[]; dt_R=[]; ph_L=[]; ph_R=[]; bcyc=[]; nsess=0; nonsetL=0; nonsetR=0;
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
        tL=(oL-1)/fpsW; tR=(oR-1)/fpsW;                     % L / R protraction-onset times (s)
        if isempty(tL) && isempty(tR), continue; end
        % ---- breath inspiration onsets (troughs of BP breath) ----
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(ai).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
        catch
            continue
        end
        [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        [~,il]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        tInsp=(il-1)/fb;                                     % inspiration-onset times (s)
        if numel(tInsp) < 2, continue; end
        bcyc = [bcyc; diff(tInsp)]; %#ok<AGROW>
        tPeak = zeros(numel(il)-1,1);                        % breath PEAK time within each cycle
        for k = 1:numel(il)-1, [~,rel]=max(brf(il(k):il(k+1))); tPeak(k)=(il(k)+rel-2)/fb; end
        pL = coordPhase(tL, tInsp, tPeak, CYC_MAX); pR = coordPhase(tR, tInsp, tPeak, CYC_MAX);   % breath phase
        dL = psth(tL, tInsp, PERI_WIN);      dR = psth(tR, tInsp, PERI_WIN);        % peri-insp times (all pairs)
        dt_L=[dt_L;dL]; ph_L=[ph_L;pL]; dt_R=[dt_R;dR]; ph_R=[ph_R;pR]; %#ok<AGROW>
        nonsetL=nonsetL+numel(tL); nonsetR=nonsetR+numel(tR); nsess=nsess+1;
    end
end
assert(~isempty(ph_L)||~isempty(ph_R),'no whisk onsets paired to breaths');
medCyc = median(bcyc);
[muL,RL,pRayL] = circstat(2*pi*ph_L);
[muR,RR,pRayR] = circstat(2*pi*ph_R);
fprintf('%d sessions | breath cycle %.0f ms (%.1f Hz)\n', nsess, 1000*medCyc, 1/medCyc);
fprintf('L: n=%d  phase mean %.2f (%.0f deg)  Rvec=%.3f  Rayleigh p=%.2g\n', numel(ph_L), muL/(2*pi), rad2deg(muL), RL, pRayL);
fprintf('R: n=%d  phase mean %.2f (%.0f deg)  Rvec=%.3f  Rayleigh p=%.2g\n', numel(ph_R), muR/(2*pi), rad2deg(muR), RR, pRayR);

% ============================ FIGURE (1x2) ============================
figure('Color','w','Position',[120 160 1080 470]);

% ---- LEFT: peri-inspiration PSTH (all onset-inspiration pairs within +/-PERI_WIN) ----
W = PERI_WIN; edges = linspace(-W, W, NBINS_T+1);
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
hh(ax1, dt_L, edges, cL); hh(ax1, dt_R, edges, cR);
xline(ax1,0,'k--','LineWidth',1.2);                        % inspiration onset
for jj = 1:floor(W/medCyc)                                 % grey marks at +/- breath-period multiples
    xline(ax1, jj*medCyc,'Color',[.6 .6 .6],'LineStyle',':'); xline(ax1,-jj*medCyc,'Color',[.6 .6 .6],'LineStyle',':');
end
xlim(ax1,[-W W]); xlabel(ax1,'whisk protraction onset - inspiration (s)'); ylabel(ax1,'prob');
legend(ax1,{'L','R'},'Location','northeast');
title(ax1, sprintf('peri-inspiration PSTH (+/-%.2f s)  L n=%d / R n=%d, %d sess  breath %.0f ms', W, nonsetL, nonsetR, nsess, 1000*medCyc));
axis(ax1,'square');

% ---- RIGHT: breath-phase distribution (LINEAR), L vs R + mean phase ----
ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
pedges = linspace(0, 1, NBINS_PH+1);
hh(ax2, ph_L, pedges, cL); hh(ax2, ph_R, pedges, cR);
xline(ax2, 0.5, 'k:', 'LineWidth',1);                       % breath peak
xline(ax2, muL/(2*pi), '--', 'Color',cL, 'LineWidth',1.4);
xline(ax2, muR/(2*pi), '--', 'Color',cR, 'LineWidth',1.4);
xlim(ax2,[0 1]); xlabel(ax2,'breath phase of whisk onset (0=insp onset, 0.5=peak, 1=next onset)'); ylabel(ax2,'prob');
legend(ax2,{'L','R'},'Location','northeast');
title(ax2, sprintf('breath phase   L %.2f (Rvec %.2f, p=%.1g) / R %.2f (Rvec %.2f, p=%.1g)', ...
    muL/(2*pi), RL, pRayL, muR/(2*pi), RR, pRayR));
axis(ax2,'square');

% ================= helpers =================
function ph = coordPhase(tw, tInsp, tPeak, cycmax)
    % piecewise breath phase: 0 = inspiration onset, 0.5 = breath peak, 1 = next onset
    ph = [];
    if isempty(tw) || numel(tInsp)<2, return; end
    bn = discretize(tw, tInsp); ok = ~isnan(bn) & bn < numel(tInsp);
    bn = bn(ok); tw = tw(ok);
    on0 = tInsp(bn); on1 = tInsp(bn+1); pk = tPeak(bn); cyc = on1 - on0;
    ph = nan(size(tw)); rise = tw < pk;                                        % onset->peak vs peak->next onset
    ph(rise)  = 0.5 * (tw(rise)  - on0(rise))  ./ (pk(rise)  - on0(rise));
    ph(~rise) = 0.5 + 0.5 * (tw(~rise) - pk(~rise)) ./ (on1(~rise) - pk(~rise));
    good = cyc>0 & cyc<=cycmax & pk>on0 & pk<on1 & ph>=0 & ph<1;
    ph = ph(good); ph = ph(:);
end
function d = psth(tw, tInsp, win)
    % all (onset - inspiration) lags for pairs within +/-win (onsets may count >once)
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
    Z=N*R^2; pRay=exp(-Z)*(1+(2*Z-Z^2)/(4*N));            % Rayleigh test (approx)
end
function hh(ax, x, edges, col)
    histogram(ax, x, edges, 'Normalization','probability','DisplayStyle','stairs','EdgeColor',col,'LineWidth',1.5);
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
