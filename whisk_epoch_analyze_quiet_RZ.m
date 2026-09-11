% whisk_epoch_analyze_quiet_RZ  (run as a script)
% QUIET (non-whisking) epochs: envelope BELOW THR_FRAC*95th-pct for >= MIN_DUR.
% Control for whisk_epoch_analyze_auto_RZ.m (same engine; only detection inverted).
% In-epoch samples split by PER-BREATH rate (1/ITI from inspiration onsets):
% BASAL (f < BASAL_HZ, BLACK) vs SNIFFING (f > SNIFF_HZ, RED); between = dropped.
%
% SINGLE 3x3 figure:
%   row1 BASAL scatters (black) : xL vs xR | breath vs xL | breath vs xR
%   row2 SNIFF scatters (red)   : xL vs xR | breath vs xL | breath vs xR
%   row3 triggered averages     : R-peak-triggered L | insp-triggered L | insp-triggered R
% Scatters: linear fit + r,p; shared lims (WLIM/BLIM crop percentiles); all axes square.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
THR_FRAC  = 0.50;       % QUIET = envelope BELOW THR_FRAC * 95th-pct(envelope)
MIN_DUR   = 1.0;        % s, minimum quiet-epoch duration
MERGE_GAP = 0.20;       % s, merge epochs closer than this
BASAL_HZ  = 5;          % Hz, per-breath 1/ITI <  this = basal,  BLACK
SNIFF_HZ  = 7;          % Hz, per-breath 1/ITI >  this = sniff,  RED  (between = dropped)
BP        = [6 30];     % whisk bandpass (Hz)
BR_BP     = [2 15];     % breathing bandpass (Hz) — also breath onset detection band
TRIG_WIN  = 0.6;        % s, +/- window for triggered averages
WLIM      = [0.005 99.95];   % whisker scatter CROP percentiles [lo hi]; [] = [0.1 99.9]
BLIM      = [0.005 99.95];   % breath  scatter CROP percentiles [lo hi]; [] = [0.1 99.9]
TLIM      = [];         % s, x-limit for ALL row3 triggered averages [lo hi]; [] = full +/-TRIG_WIN
fpsW      = 400;
MAXPTS    = 8000;       % subsample for scatter display
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

[XL,XR,BR,FB,T,X] = pool_epochs(@(env,t) detect_quiet(env,t, THR_FRAC, MIN_DUR, MERGE_GAP), ...
    dataRoot, whiskDir, noseDir, BP, BR_BP, BASAL_HZ, SNIFF_HZ, TRIG_WIN, fpsW);

base  = ~isnan(FB) & FB <  BASAL_HZ;
sniff = ~isnan(FB) & FB >  SNIFF_HZ;
durB = nnz(base)/fpsW; durS = nnz(sniff)/fpsW;
fprintf('pooled %d samples, %d QUIET epochs, %d sessions\n', numel(XL), X.nep, X.nsess);
fprintf('  basal %.1f s (med %.1f Hz, <%g) | sniff %.1f s (med %.1f Hz, >%g)\n', ...
    durB, median(FB(base),'omitnan'), BASAL_HZ, durS, median(FB(sniff),'omitnan'), SNIFF_HZ);

ttl = sprintf('QUIET epochs: basal (black, <%g Hz, %.1f s) vs sniff (red, >%g Hz, %.1f s) — %d ep, %d sess', ...
        BASAL_HZ, durB, SNIFF_HZ, durS, X.nep, X.nsess);
render_overlay(XL,XR,BR, base, sniff, MAXPTS, T, WLIM, BLIM, TLIM, ttl);

% =================================================================
function [XL,XR,BR,FB,T,X] = pool_epochs(detectFn, dataRoot, whiskDir, noseDir, ...
        BP, BR_BP, BASAL_HZ, SNIFF_HZ, TRIG_WIN, fpsW)
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess),'no sessions with whisk csv');

nW = round(TRIG_WIN*fpsW); W = 2*nW+1;
XL=[];XR=[];BR=[];FB=[]; nep=0; nsess=0;
T.lags_t = (-nW:nW)/fpsW;
[T.sXL,T.sXL2,T.sXR,T.sXR2] = deal(zeros(2,W));
[T.sLR,T.sLR2] = deal(zeros(2,W)); T.nLR = zeros(2,1);
T.nT = zeros(2,1);

for e = 1:numel(sess)
    animal = sess{e}{1}; kk = sess{e}{2};
    M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kk)), 0.6);    % lik<0.6 -> linear interp
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;
    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
    env  = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;
    phiR = angle(hilbert(xR));                       % for R-peak triggering
    brw = nan(size(t)); phiB = nan(size(t)); fbr = nan(size(t));
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kk)), dataRoot);
        if isfile(Pn.breath)
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
            [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
            brw = interp1(tBr, brf, t, 'linear', NaN);
            valid = ~isnan(brw);
            phiB(valid) = angle(hilbert(brw(valid)));
            [~,iloc] = findpeaks(-brf, 'MinPeakProminence', 0.5*std(brf), 'MinPeakDistance', round(fb/BR_BP(2)));
            tInsp = (iloc-1)/fb; bfreq = 1./diff(tInsp);
            bn = discretize(t, tInsp); okb = ~isnan(bn) & bn<=numel(bfreq);
            fbr(okb) = bfreq(bn(okb));        % each sample tagged with its breath's 1/ITI
        end
    catch
    end
    ep = detectFn(env, t);
    m = false(size(t));
    for q=1:size(ep,1), m = m | (t>=ep(q,1) & t<=ep(q,2)); end
    XL=[XL;xL(m)]; XR=[XR;xR(m)]; BR=[BR;brw(m)]; FB=[FB;fbr(m)]; %#ok<AGROW>
    nep=nep+size(ep,1); nsess=nsess+1;
    % inspiration-triggered whisk averages (t=0 = insp peak = phiB upward zero-cross)
    pk = find(phiB(1:end-1)<0 & phiB(2:end)>=0);
    for p = pk'
        if p-nW<1 || p+nW>numel(t) || isnan(fbr(p)), continue; end
        if ~all(m(p-nW:p+nW)), continue; end
        if     fbr(p) <  BASAL_HZ, c = 1;
        elseif fbr(p) >  SNIFF_HZ, c = 2;
        else,  continue; end
        idx = p-nW:p+nW; l=xL(idx)'; r=xR(idx)';
        T.sXL(c,:)=T.sXL(c,:)+l; T.sXL2(c,:)=T.sXL2(c,:)+l.^2;
        T.sXR(c,:)=T.sXR(c,:)+r; T.sXR2(c,:)=T.sXR2(c,:)+r.^2; T.nT(c)=T.nT(c)+1;
    end
    % R-peak-triggered whisk L (t=0 = R protraction peak = phiR upward zero-cross)
    pkR = find(phiR(1:end-1)<0 & phiR(2:end)>=0);
    for p = pkR'
        if p-nW<1 || p+nW>numel(t) || isnan(fbr(p)), continue; end
        if ~all(m(p-nW:p+nW)), continue; end
        if     fbr(p) <  BASAL_HZ, c = 1;
        elseif fbr(p) >  SNIFF_HZ, c = 2;
        else,  continue; end
        l = xL(p-nW:p+nW)';
        T.sLR(c,:)=T.sLR(c,:)+l; T.sLR2(c,:)=T.sLR2(c,:)+l.^2; T.nLR(c)=T.nLR(c)+1;
    end
end
assert(~isempty(XL),'no samples in quiet epochs');
X.nep = nep; X.nsess = nsess;
end

% ================= rendering =================
function render_overlay(XL,XR,BR, mb, ms, MAXPTS, T, WLIM, BLIM, TLIM, ttl)
    mb = mb(:); ms = ms(:);
    if isempty(WLIM), WLIM=[0.1 99.9]; end
    if isempty(BLIM), BLIM=[0.1 99.9]; end
    vw=[XL(mb|ms);XR(mb|ms)]; vw=vw(~isnan(vw)); wlim=max(abs(prctile(vw,WLIM)))*[-1 1];   % symmetric crop
    vb=BR(mb|ms); vb=vb(~isnan(vb)); blim=max(abs(prctile(vb,BLIM)))*[-1 1];
    figure('Color','w','Position',[40 40 1180 980]);
    amp1(1, XL, XR, mb, [0 0 0], 'xL','xR','basal  xL vs xR',         wlim, wlim);
    amp1(2, BR, XL, mb, [0 0 0], 'breath','xL','basal  breath vs xL', blim, wlim);
    amp1(3, BR, XR, mb, [0 0 0], 'breath','xR','basal  breath vs xR', blim, wlim);
    amp1(4, XL, XR, ms, [0.85 0.1 0.1], 'xL','xR','sniff  xL vs xR',         wlim, wlim);
    amp1(5, BR, XL, ms, [0.85 0.1 0.1], 'breath','xL','sniff  breath vs xL', blim, wlim);
    amp1(6, BR, XR, ms, [0.85 0.1 0.1], 'breath','xR','sniff  breath vs xR', blim, wlim);
    [mLR,sLR]=msem(T.sLR,T.sLR2,T.nLR); trigplot(subplot(3,3,7), T.lags_t, mLR,sLR, T.nLR, 'R-triggered L', 'xL', 'time from R peak (s)', TLIM);
    [mL,sL]=msem(T.sXL,T.sXL2,T.nT);    trigplot(subplot(3,3,8), T.lags_t, mL,sL, T.nT, 'insp-triggered L', 'xL', 'time from insp peak (s)', TLIM);
    [mR,sR]=msem(T.sXR,T.sXR2,T.nT);    trigplot(subplot(3,3,9), T.lags_t, mR,sR, T.nT, 'insp-triggered R', 'xR', 'time from insp peak (s)', TLIM);
    sgtitle(ttl);
    function amp1(p, Xx, Yy, mask, col, xlab, ylab, ttl1, xl, yl)
        [r,pp]=corrp(Xx(mask),Yy(mask));
        ax=subplot(3,3,p); scat1(ax,Xx,Yy,mask,[col 0.15],MAXPTS); fitline(ax,Xx(mask),Yy(mask),xl,col);
        xlim(ax,xl); ylim(ax,yl); axis(ax,'square'); xlabel(ax,xlab); ylabel(ax,ylab);
        title(ax, sprintf('%s\nr=%.2f  p=%.1e', ttl1, r, pp));
    end
end
function trigplot(ax, x, m, s, n, ttl, ylab, xlab, tlim)
    hold(ax,'on'); grid(ax,'on');
    shade(ax,x,m(1,:),s(1,:),[0 0 0]); shade(ax,x,m(2,:),s(2,:),[0.85 0.1 0.1]);
    plot(ax,x,m(1,:),'k-','LineWidth',1.4); plot(ax,x,m(2,:),'-','Color',[0.85 0.1 0.1],'LineWidth',1.4);
    xline(ax,0,'b:');
    if isempty(tlim), mg=0.02*(x(end)-x(1)); xlim(ax,[x(1)-mg x(end)+mg]); else, xlim(ax,tlim); end
    yy=[m+s; m-s]; yy=yy(isfinite(yy)); if ~isempty(yy), pad=0.05*range(yy)+eps; ylim(ax,[min(yy)-pad max(yy)+pad]); end
    axis(ax,'square'); xlabel(ax,xlab); ylabel(ax,ylab);
    title(ax, sprintf('%s  (basal n=%d, sniff n=%d)', ttl, n(1), n(2)));
end
function shade(ax,x,m,s,col)
    if any(isnan(m)) || any(isnan(s)), return; end
    patch(ax,[x fliplr(x)],[m+s fliplr(m-s)],col,'FaceAlpha',0.15,'EdgeColor','none');
end

% ================= helpers =================
function ep = detect_quiet(env, t, thrFrac, minDur, mergeGap)
    a = env(:) < thrFrac*prctile(env,95);          % BELOW threshold = quiet
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
function [mu,se] = msem(s, s2, n)
    mu = s ./ n; se = sqrt(max(s2./n - mu.^2, 0)) ./ sqrt(n);   % SEM = STD/sqrt(n)
    mu(n<1,:) = NaN; se(n<1,:) = NaN;
end
function [r,p] = corrp(x,y), g=~isnan(x)&~isnan(y); if nnz(g)<3, r=NaN; p=NaN; else, [r,p]=corr(x(g),y(g)); end, end
function fitline(ax, x, y, xl, col)
    g = ~isnan(x) & ~isnan(y); if nnz(g) < 3, return; end
    pf = polyfit(x(g), y(g), 1); xs = linspace(xl(1), xl(2), 2);
    plot(ax, xs, polyval(pf, xs), '-', 'Color', col, 'LineWidth', 1.6);
end
function scat1(ax, x, y, m, col, MAXPTS)
    i = find(m & ~isnan(x) & ~isnan(y)); i = i(subidx(numel(i),MAXPTS));
    plot(ax, x(i),y(i),'.','MarkerSize',1,'Color',col); hold(ax,'on'); grid(ax,'on');
end
function idx = subidx(N, n), if N<=n, idx=1:N; else, idx=randperm(N,n); end, end
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
