% whisk_epoch_analyze_auto_RZ  (run as a script)
% Auto-detected WHISKING epochs (envelope ABOVE THR_FRAC*95th-pct for >= MIN_DUR).
% In-epoch samples split by PER-BREATH rate (1/ITI from thermal inspiration troughs):
% BASAL (f < BASAL_HZ, BLACK) vs SNIFFING (f > SNIFF_HZ, RED); between = dropped.
%
% Two 3×3 figures:
%   Figure 1 (angle):    row1 BASAL  : xL vs xR | breath vs xL | breath vs xR
%                        row2 SNIFF
%                        row3 triggered averages: R-peak->L | insp->L | insp->R
%   Figure 2 (velocity): same 3×3 layout but angular velocity vL/vR (deg/s)

% ============================ USER-EDITABLE ============================
dataRoot   = "D:\260615_thermalNbasler";
whiskDir   = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
THR_FRAC   = 0.2;            % epoch threshold = THR_FRAC * 95th-pct(envelope); 0 = whole recording
MIN_DUR    = 0;          % s, minimum epoch duration
MERGE_GAP  = 0.30;         % s, merge epochs closer than this
ENV_SMOOTH = 0;            % frames, movmean smoothing on envelope; 0 = no smoothing
BASAL_HZ   = 5;            % Hz, per-breath rate <  this = basal  (BLACK)
SNIFF_HZ   = 8;            % Hz, per-breath rate >  this = sniff  (RED)
SLOW_HZ    = 10;           % Hz, whisk frequency threshold for slow-whisking condition
BP         = [5 30];       % whisk bandpass (Hz)
BR_BP      = [2 15];       % breathing bandpass (Hz)
TRIG_WIN   = 0.6;          % s, +/- window for triggered averages
PROT_PROM  = 1;            % protraction peak prominence (x std of BP whisk)
MINDIST    = 0.025;        % s, min distance between protraction peaks
BR_PROM    = 0.5;          % inspiration trough prominence (x std of BP breath)
WLIM       = [0.005 99.95]; % whisker scatter crop percentiles
BLIM       = [0.05 99.5]; % breath  scatter crop percentiles
VLIM       = [0.005 99.95]; % velocity scatter crop percentiles
TLIM       = [-.1 .1];      % s, x-limits for row-3 triggered averages
fpsW       = 400;
MAXPTS     = 800000;
EXCLUDE    = "5840027";
SAVE_FIGS  = true;
OUTDIR     = fullfile(pwd, 'whisk_epoch_exports');
% ======================================================================

if SAVE_FIGS && ~exist(OUTDIR,'dir')
    mkdir(OUTDIR);
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

[XL,XR,VL,VR,BR,FB,FW,T,Tvel,X] = pool_epochs( ...
    @(env,t) detect(env,t, THR_FRAC, MIN_DUR, MERGE_GAP), ...
    dataRoot, whiskDir, noseDir, BP, BR_BP, ...
    BASAL_HZ, SNIFF_HZ, TRIG_WIN, PROT_PROM, MINDIST, BR_PROM, ENV_SMOOTH, fpsW, EXCLUDE);

base  = ~isnan(FB) & FB <  BASAL_HZ;
sniff = ~isnan(FB) & FB >  SNIFF_HZ;
durB = nnz(base)/fpsW; durS = nnz(sniff)/fpsW;
fprintf('pooled %d samples, %d epochs, %d sessions\n', numel(XL), X.nep, X.nsess);
fprintf('  basal %.1f s (med %.1f Hz, <%g) | sniff %.1f s (med %.1f Hz, >%g)\n', ...
    durB, median(FB(base),'omitnan'), BASAL_HZ, durS, median(FB(sniff),'omitnan'), SNIFF_HZ);

ttl = sprintf('WHISK epochs: basal (black, <%g Hz, %.1f s) vs sniff (red, >%g Hz, %.1f s) — %d ep, %d sess', ...
    BASAL_HZ, durB, SNIFF_HZ, durS, X.nep, X.nsess);

% Figure 1: whisker angle
fig1 = render_overlay(XL, XR, BR, base, sniff, MAXPTS, T, WLIM, BLIM, TLIM, ['Angle — ' ttl]);

% Figure 2: angular velocity
fig2 = render_overlay_vel(VL, VR, BR, base, sniff, MAXPTS, T, Tvel, VLIM, BLIM, TLIM, ['Velocity — ' ttl]);

% Figure 3: L vs R scatter — all whisking vs slow-whisk + fast-sniff
fig3 = render_LR_scatter_two_conditions(XL, XR, FB, FW, SNIFF_HZ, SLOW_HZ, MAXPTS, WLIM);

if SAVE_FIGS
    exportVectorFig(fig1, fullfile(OUTDIR, 'whisk_epoch_angle_overlay'));
    exportVectorFig(fig2, fullfile(OUTDIR, 'whisk_epoch_angular_velocity_overlay'));
    exportVectorFig(fig3, fullfile(OUTDIR, 'whisk_LR_scatter_all_vs_slowWhisk_fastSniff'));
end

% =================================================================
% pool_epochs: collect scatter samples + triggered-average accumulators.
%   XL/XR      = epoch-masked angle samples
%   VL/VR      = epoch-masked angular velocity samples (deg/s)
%   BR         = breath signal interpolated to whisk frame rate
%   FB         = per-sample breath frequency tag (1/ITI)
%   T          = angle triggered-average struct
%   Tvel       = velocity triggered-average struct (same triggers as T)
%   X          = {nep, nsess}
% =================================================================
function [XL,XR,VL,VR,BR,FB,FW,T,Tvel,X] = pool_epochs(detectFn, dataRoot, whiskDir, noseDir, ...
        BP, BR_BP, BASAL_HZ, SNIFF_HZ, TRIG_WIN, PROT_PROM, MINDIST, BR_PROM, ENV_SMOOTH, fpsW, EXCLUDE)

[bw,aw] = butter(3, BP/(fpsW/2), 'bandpass');

ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess), 'no sessions with whisk csv');

nW = round(TRIG_WIN*fpsW); W = 2*nW+1;
XL=[]; XR=[]; VL=[]; VR=[]; BR=[]; FB=[]; FW=[]; nep=0; nsess=0;

T.lags_t = (-nW:nW)/fpsW;
[T.sXL, T.sXL2, T.sXR, T.sXR2] = deal(zeros(2,W));
[T.sLR,  T.sLR2]  = deal(zeros(2,W));
[T.sRR,  T.sRR2]  = deal(zeros(2,W));
T.nT  = zeros(2,1);
T.nLR = zeros(2,1);

[Tvel.sXL, Tvel.sXL2, Tvel.sXR, Tvel.sXR2] = deal(zeros(2,W));
[Tvel.sLR,  Tvel.sLR2]  = deal(zeros(2,W));
[Tvel.sRR,  Tvel.sRR2]  = deal(zeros(2,W));

for e = 1:numel(sess)
    animal = sess{e}{1}; kk = sess{e}{2};

    % --- whisk DLC ---
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
    catch ME
        warning('whisk load failed: %s n%d: %s', animal, kk, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),   M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)' / fpsW;

    xL = filtfilt(bw,aw, fillmissing(La - mean(La,'omitnan'), 'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra - mean(Ra,'omitnan'), 'linear'));

    vL = gradient(xL) * fpsW;   % deg/s
    vR = gradient(xR) * fpsW;

    % --- local whisk frequency from L protraction peaks ---
    [~,iLpeak] = findpeaks(xL, ...
        'MinPeakProminence', PROT_PROM*std(xL), ...
        'MinPeakDistance',   round(MINDIST*fpsW));
    fwFrame = nan(size(t));
    if numel(iLpeak) >= 3
        tLpeak = (subsamp(xL, iLpeak) - 1) / fpsW;
        fLpeak = nan(size(tLpeak));
        for k = 2:numel(tLpeak)-1
            dt = tLpeak(k+1) - tLpeak(k-1);
            if dt > 0 && dt <= 0.5
                fLpeak(k) = 2/dt;
            end
        end
        ok = isfinite(fLpeak);
        if nnz(ok) >= 2
            fwFrame = interp1(tLpeak(ok), fLpeak(ok), t, 'linear', NaN);
        end
    end

    env = movmean((abs(hilbert(xL)) + abs(hilbert(xR))) / 2, max(1, round(ENV_SMOOTH*fpsW)));

    % --- breath ---
    brw = nan(size(t)); fbr = nan(size(t)); pkInsp = [];
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);
        if isfile(Pn.breath)
            Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
            tBr = (0:numel(br)-1)' / fb;
            [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
            brf = filtfilt(b2,a2, fillmissing(br,'linear'));

            brw = interp1(tBr, brf, t, 'linear', NaN);

            [~,iloc] = findpeaks(-brf, ...
                'MinPeakProminence', BR_PROM*std(brf), ...
                'MinPeakDistance',   round(fb/BR_BP(2)));

            if numel(iloc) >= 2
                tInsp  = (subsamp(-brf, iloc) - 1) / fb;
                tInsp  = tInsp(:);
                cycdur = diff(tInsp);
                bfreq  = 1 ./ cycdur;
                bn  = discretize(t, tInsp); okb = ~isnan(bn) & bn <= numel(bfreq);
                fbr(okb) = bfreq(bn(okb));
                pkInsp = round(tInsp * fpsW) + 1;
                pkInsp = pkInsp(pkInsp >= 1 & pkInsp <= numel(t));
            end
        end
    catch
    end

    % --- epoch mask ---
    ep = detectFn(env, t);
    m = false(size(t));
    for q = 1:size(ep,1), m = m | (t >= ep(q,1) & t <= ep(q,2)); end

    % pooled scatter samples
    XL = [XL; xL(m)]; XR = [XR; xR(m)]; %#ok<AGROW>
    VL = [VL; vL(m)]; VR = [VR; vR(m)];
    BR = [BR; brw(m)]; FB = [FB; fbr(m)];
    FW = [FW; fwFrame(m)];
    nep = nep + size(ep,1); nsess = nsess + 1;

    % --- inspiration-triggered averages ---
    for p = pkInsp'
        if p-nW<1 || p+nW>numel(t) || isnan(fbr(p)), continue; end
        if ~all(m(p-nW:p+nW)), continue; end
        if     fbr(p) <  BASAL_HZ, c = 1;
        elseif fbr(p) >  SNIFF_HZ, c = 2;
        else,  continue; end
        idx = p-nW:p+nW;
        l=xL(idx)'; r=xR(idx)'; vl=vL(idx)'; vr=vR(idx)';
        T.sXL(c,:) = T.sXL(c,:) + l;   T.sXL2(c,:) = T.sXL2(c,:) + l.^2;
        T.sXR(c,:) = T.sXR(c,:) + r;   T.sXR2(c,:) = T.sXR2(c,:) + r.^2;
        T.nT(c) = T.nT(c) + 1;
        Tvel.sXL(c,:) = Tvel.sXL(c,:) + vl;  Tvel.sXL2(c,:) = Tvel.sXL2(c,:) + vl.^2;
        Tvel.sXR(c,:) = Tvel.sXR(c,:) + vr;  Tvel.sXR2(c,:) = Tvel.sXR2(c,:) + vr.^2;
    end

    % --- R-protraction-peak-triggered averages ---
    [~,iRpeak] = findpeaks(xR, ...
        'MinPeakProminence', PROT_PROM*std(xR), ...
        'MinPeakDistance',   round(MINDIST*fpsW));
    iRpeak = iRpeak(m(iRpeak));

    for p = iRpeak'
        if p-nW<1 || p+nW>numel(t) || isnan(fbr(p)), continue; end
        if ~all(m(p-nW:p+nW)), continue; end
        if     fbr(p) <  BASAL_HZ, c = 1;
        elseif fbr(p) >  SNIFF_HZ, c = 2;
        else,  continue; end
        idx = p-nW:p+nW;
        l=xL(idx)'; vl=vL(idx)';
        r=xR(idx)'; vr=vR(idx)';
        T.sLR(c,:)    = T.sLR(c,:)    + l;   T.sLR2(c,:)    = T.sLR2(c,:)    + l.^2;
        T.sRR(c,:)    = T.sRR(c,:)    + r;   T.sRR2(c,:)    = T.sRR2(c,:)    + r.^2;
        T.nLR(c) = T.nLR(c) + 1;
        Tvel.sLR(c,:) = Tvel.sLR(c,:) + vl;  Tvel.sLR2(c,:) = Tvel.sLR2(c,:) + vl.^2;
        Tvel.sRR(c,:) = Tvel.sRR(c,:) + vr;  Tvel.sRR2(c,:) = Tvel.sRR2(c,:) + vr.^2;
    end
end

assert(~isempty(XL), 'no samples in detected epochs');
X.nep = nep; X.nsess = nsess;
end

% ================= rendering =================
function fig = render_LR_scatter_two_conditions(XL, XR, FB, FW, SNIFF_HZ, SLOW_HZ, MAXPTS, WLIM)
    validAll = isfinite(XL) & isfinite(XR);
    slowFast = validAll & isfinite(FB) & isfinite(FW) & FB > SNIFF_HZ & FW < SLOW_HZ;

    vw = [XL(validAll); XR(validAll)];
    vw = vw(isfinite(vw));
    if isempty(vw)
        wlim = [-1 1];
    else
        wlim = max(abs(prctile(vw, WLIM))) * [-1 1];
    end

    fig = figure('Color','w','Position',[80 80 850 420], ...
        'Renderer','painters','InvertHardcopy','off');
    set(fig,'GraphicsSmoothing','off');

    ax1 = subplot(1,2,1);
    scatterPanel(ax1, XL, XR, validAll, [0.20 0.20 0.20], MAXPTS, ...
        'all whisking samples', wlim);

    ax2 = subplot(1,2,2);
    scatterPanel(ax2, XL, XR, slowFast, [0.94 0.29 0.14], MAXPTS, ...
        sprintf('slow whisk <%.0f Hz + sniff >%.0f Hz', SLOW_HZ, SNIFF_HZ), wlim);

    sgtitle('L vs R whisker angle scatter');
end

function scatterPanel(ax, XL, XR, mask, col, MAXPTS, ttl, wlim)
    hold(ax,'on'); grid(ax,'on');

    idx = find(mask & isfinite(XL) & isfinite(XR));
    nTotal = numel(idx);

    if nTotal > MAXPTS
        idx = idx(randperm(nTotal, MAXPTS));
    end

    lightCol = 0.75*[1 1 1] + 0.25*col;
    plot(ax, XL(idx), XR(idx), '.', 'MarkerSize', 2, 'Color', lightCol);

    [r,p] = corrp(XL(mask), XR(mask));
    fitline(ax, XL(mask), XR(mask), wlim, col);

    xlim(ax,wlim); ylim(ax,wlim);
    axis(ax,'square'); box(ax,'off');
    xlabel(ax,'xL angle'); ylabel(ax,'xR angle');
    title(ax, sprintf('%s\nn=%d  r=%.2f  p=%.1e', ttl, nTotal, r, p));
end

function fig = render_overlay(XL, XR, BR, mb, ms, MAXPTS, T, WLIM, BLIM, TLIM, ttl)
    mb = mb(:); ms = ms(:);
    if isempty(WLIM), WLIM=[0.1 99.9]; end
    if isempty(BLIM), BLIM=[0.1 99.9]; end
    vw=[XL(mb|ms);XR(mb|ms)]; vw=vw(~isnan(vw)); wlim=max(abs(prctile(vw,WLIM)))*[-1 1];
    vb=BR(mb|ms); vb=vb(~isnan(vb)); blim=max(abs(prctile(vb,BLIM)))*[-1 1];
    fig = figure('Color','w','Position',[40 40 1180 980], ...
        'Renderer','painters','InvertHardcopy','off');
    set(fig,'GraphicsSmoothing','off');
    amp1(1, XL, XR, mb, [0 0 0],       'xL','xR','basal  xL vs xR',         wlim, wlim);
    amp1(2, BR, XL, mb, [0 0 0],       'breath','xL','basal  breath vs xL', blim, wlim);
    amp1(3, BR, XR, mb, [0 0 0],       'breath','xR','basal  breath vs xR', blim, wlim);
    amp1(4, XL, XR, ms, [0.85 0.1 0.1],'xL','xR','sniff  xL vs xR',         wlim, wlim);
    amp1(5, BR, XL, ms, [0.85 0.1 0.1],'breath','xL','sniff  breath vs xL', blim, wlim);
    amp1(6, BR, XR, ms, [0.85 0.1 0.1],'breath','xR','sniff  breath vs xR', blim, wlim);
    [mLR,sLR]=msem(T.sLR, T.sLR2, T.nLR);
    [mRR,sRR]=msem(T.sRR, T.sRR2, T.nLR);
    trigplot_selfR(subplot(3,3,7), T.lags_t, mLR, sLR, mRR, sRR, T.nLR, ...
        'R-triggered L + R-self', 'angle', 'time from R peak (s)', TLIM);
    [mL,sL]=msem(T.sXL, T.sXL2, T.nT);
    trigplot(subplot(3,3,8), T.lags_t, mL,  sL,  T.nT,  'insp-triggered L','xL', 'time from insp peak (s)', TLIM);
    [mR,sR]=msem(T.sXR, T.sXR2, T.nT);
    trigplot(subplot(3,3,9), T.lags_t, mR,  sR,  T.nT,  'insp-triggered R','xR', 'time from insp peak (s)', TLIM);
    sgtitle(ttl);
    function amp1(p, Xx, Yy, mask, col, xlab, ylab, ttl1, xl, yl)
        [r,pp]=corrp(Xx(mask),Yy(mask));
        ax=subplot(3,3,p); scat1(ax,Xx,Yy,mask,col,MAXPTS); fitline(ax,Xx(mask),Yy(mask),xl,col);
        xlim(ax,xl); ylim(ax,yl); axis(ax,'square'); xlabel(ax,xlab); ylabel(ax,ylab);
        title(ax, sprintf('%s\nr=%.2f  p=%.1e', ttl1, r, pp));
    end
end

function fig = render_overlay_vel(VL, VR, BR, mb, ms, MAXPTS, T, Tvel, VLIM, BLIM, TLIM, ttl)
    mb = mb(:); ms = ms(:);
    if isempty(VLIM), VLIM=[0.005 99.95]; end
    if isempty(BLIM), BLIM=[0.1 99.9]; end
    vv=[VL(mb|ms);VR(mb|ms)]; vv=vv(~isnan(vv)); vlim=max(abs(prctile(vv,VLIM)))*[-1 1];
    vb=BR(mb|ms); vb=vb(~isnan(vb)); blim=max(abs(prctile(vb,BLIM)))*[-1 1];
    fig = figure('Color','w','Position',[60 60 1180 980], ...
        'Renderer','painters','InvertHardcopy','off');
    set(fig,'GraphicsSmoothing','off');
    amp1(1, VL, VR, mb, [0 0 0],       'vL (deg/s)','vR (deg/s)','basal  vL vs vR',         vlim, vlim);
    amp1(2, BR, VL, mb, [0 0 0],       'breath','vL (deg/s)','basal  breath vs vL',         blim, vlim);
    amp1(3, BR, VR, mb, [0 0 0],       'breath','vR (deg/s)','basal  breath vs vR',         blim, vlim);
    amp1(4, VL, VR, ms, [0.85 0.1 0.1],'vL (deg/s)','vR (deg/s)','sniff  vL vs vR',         vlim, vlim);
    amp1(5, BR, VL, ms, [0.85 0.1 0.1],'breath','vL (deg/s)','sniff  breath vs vL',         blim, vlim);
    amp1(6, BR, VR, ms, [0.85 0.1 0.1],'breath','vR (deg/s)','sniff  breath vs vR',         blim, vlim);
    [mLR,sLR]=msem(Tvel.sLR, Tvel.sLR2, T.nLR);
    [mRR,sRR]=msem(Tvel.sRR, Tvel.sRR2, T.nLR);
    trigplot_selfR(subplot(3,3,7), T.lags_t, mLR, sLR, mRR, sRR, T.nLR, ...
        'R-triggered vL + vR-self', 'angular velocity', 'time from R peak (s)', TLIM);
    [mVL,sVL]=msem(Tvel.sXL, Tvel.sXL2, T.nT);
    trigplot(subplot(3,3,8), T.lags_t, mVL, sVL, T.nT,  'insp-triggered vL','vL (deg/s)', 'time from insp peak (s)', TLIM);
    [mVR,sVR]=msem(Tvel.sXR, Tvel.sXR2, T.nT);
    trigplot(subplot(3,3,9), T.lags_t, mVR, sVR, T.nT,  'insp-triggered vR','vR (deg/s)', 'time from insp peak (s)', TLIM);
    sgtitle(ttl);
    function amp1(p, Xx, Yy, mask, col, xlab, ylab, ttl1, xl, yl)
        [r,pp]=corrp(Xx(mask),Yy(mask));
        ax=subplot(3,3,p); scat1(ax,Xx,Yy,mask,col,MAXPTS); fitline(ax,Xx(mask),Yy(mask),xl,col);
        xlim(ax,xl); ylim(ax,yl); axis(ax,'square'); xlabel(ax,xlab); ylabel(ax,ylab);
        title(ax, sprintf('%s\nr=%.2f  p=%.1e', ttl1, r, pp));
    end
end

function trigplot_selfR(ax, x, mL, sL, mR, sR, n, ttl, ylab, xlab, tlim)
    hold(ax,'on'); grid(ax,'on');
    shade(ax,x,mL(1,:),sL(1,:),[0 0 0]);
    shade(ax,x,mL(2,:),sL(2,:),[0.85 0.1 0.1]);
    plot(ax,x,mL(1,:),'k-','LineWidth',1.4);
    plot(ax,x,mL(2,:),'-','Color',[0.85 0.1 0.1],'LineWidth',1.4);
    plot(ax,x,mR(1,:),'k--','LineWidth',1.0);
    plot(ax,x,mR(2,:),'--','Color',[0.85 0.1 0.1],'LineWidth',1.0);
    xline(ax,0,'b:');
    if isempty(tlim), mg=0.02*(x(end)-x(1)); xlim(ax,[x(1)-mg x(end)+mg]); else, xlim(ax,tlim); end
    yy=[mL+sL; mL-sL; mR+sR; mR-sR]; yy=yy(isfinite(yy));
    if ~isempty(yy), pad=0.05*range(yy)+eps; ylim(ax,[min(yy)-pad max(yy)+pad]); end
    axis(ax,'square'); xlabel(ax,xlab); ylabel(ax,ylab);
    title(ax, sprintf('%s  (basal n=%d, sniff n=%d)', ttl, n(1), n(2)));
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

function shade(ax, x, m, s, col)
    if any(isnan(m)) || any(isnan(s)), return; end
    lightCol = 0.85*[1 1 1] + 0.15*col;
    patch(ax, [x fliplr(x)], [m+s fliplr(m-s)], lightCol, 'EdgeColor','none');
end

% ================= helpers =================
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

function p = subsamp(x, idx)
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx > 1 & idx < numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp; delta = 0.5*(ym-yp)./den;
    delta(~isfinite(delta) | abs(delta)>0.5) = 0;
    p(in) = i + delta;
end

function [mu,se] = msem(s, s2, n)
    mu = s ./ n; se = sqrt(max(s2./n - mu.^2, 0)) ./ sqrt(n);
    mu(n<1,:) = NaN; se(n<1,:) = NaN;
end

function [r,p] = corrp(x,y)
    g = ~isnan(x) & ~isnan(y);
    if nnz(g) < 3, r=NaN; p=NaN; else, [r,p]=corr(x(g),y(g)); end
end

function fitline(ax, x, y, xl, col)
    g = ~isnan(x) & ~isnan(y); if nnz(g) < 3, return; end
    pf = polyfit(x(g), y(g), 1); xs = linspace(xl(1), xl(2), 2);
    plot(ax, xs, polyval(pf, xs), '-', 'Color', col, 'LineWidth', 1.6);
end

function scat1(ax, x, y, m, col, MAXPTS)
    i = find(m & ~isnan(x) & ~isnan(y)); i = i(subidx(numel(i),MAXPTS));
    lightCol = 0.75*[1 1 1] + 0.25*col;
    plot(ax, x(i), y(i), '.', 'MarkerSize', 2, 'Color', lightCol);
    hold(ax,'on'); grid(ax,'on');
end

function idx = subidx(N, n)
    if N<=n, idx=1:N; else, idx=randperm(N,n); end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens');
    if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end

function exportVectorFig(fig, baseName)
    drawnow;
    set(fig, 'Renderer',          'painters');
    set(fig, 'InvertHardcopy',    'off');
    set(fig, 'GraphicsSmoothing', 'off');

    oldUnits   = fig.Units;
    fig.Units  = 'inches';
    pos        = fig.Position;

    fig.PaperUnits        = 'inches';
    fig.PaperSize         = [pos(3) pos(4)];
    fig.PaperPosition     = [0 0 pos(3) pos(4)];
    fig.PaperPositionMode = 'manual';

    pdfFile = baseName + ".pdf";
    epsFile = baseName + ".eps";

    print(fig, pdfFile, '-dpdf',  '-painters');
    print(fig, epsFile, '-depsc', '-painters');

    fig.Units = oldUnits;

    fprintf('Exported vector figure:\n  %s\n  %s\n', pdfFile, epsFile);
end
