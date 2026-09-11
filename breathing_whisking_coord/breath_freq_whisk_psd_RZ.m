function breath_freq_whisk_psd_RZ()
% breath_freq_whisk_psd_RZ  Two-panel summary (pooled over all sessions):
%   TOP    : histogram of instantaneous breathing frequency (1/ITI) with the
%            classification basal < BASAL_HZ (black) | sniff > SNIFF_HZ (red) |
%            unclassified in between (grey).
%   BOTTOM : whisking spectral power (pwelch of BP whisk, L+R) during BASAL-
%            respiration runs (black) vs SNIFFING runs (red), restricted to
%            WHISKING EPOCHS (smoothed mean L/R envelope > ENV_THR deg).
% Thresholds: basal <5 Hz, sniff >7 Hz. Top histogram = all breaths (whole
% recording); bottom PSD = whisk epochs only.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BASAL_HZ  = 5;          % < this = basal (black)
SNIFF_HZ  = 7;          % > this = sniff (red); between = unclassified (grey)
BP        = [5 60];     % whisk bandpass (Hz)
BR_BP     = [2 20];     % breath bandpass (Hz) + onset detection
BR_PROM   = 0.5;        % inspiration trough prominence (x std of BP breath)
ENV_THR   = 5;          % deg, hard threshold on smoothed mean L/R envelope (PSD only)
ENV_SMOOTH= 0.05;       % s, envelope smoothing (movmean) for whisk-epoch detection
MIN_DUR   = 0.5;        % s, min whisk-epoch duration
MERGE_GAP = 0.00;       % s, merge whisk epochs closer than this
EDGES     = 0:0.25:30;  % breath-frequency histogram bins (Hz)
PROT_PROM = 1;          % whisk protraction trough prominence (x std of BP whisk)
ITI_MAX   = 0.2;       % s, drop whisk ITIs longer than this (cross-bout gaps)
WHZ_EDGES = 0:1:30;     % whisk frequency (1/ITI) histogram bins (Hz)
PSD_NFFT  = 512; PSD_WIN = 102; PSD_OVERLAP = 101; PSD_FPASS = [0 40];   % pwelch: seg=102 (~0.25s), step 1 samp
MINRUN    = 128;        % min samples in a class run for a PSD estimate
fpsW      = 400;
EXCLUDE   = "5840027";  % animal ids to drop from the analysis ("0027")
EXCLUDE_SESS = {{"5916297",3}, {"5916297",4}, {"5916300",4}};  % {animal,run} drops
% ======================================================================

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
nf = PSD_NFFT/2+1; fpsd = (0:PSD_NFFT/2)/PSD_NFFT*fpsW;
Pw=zeros(2,nf); Pw2=zeros(2,nf); nPw=zeros(2,1);   % 1=basal 2=sniff
freq = []; nsess = 0; whzB = []; whzS = [];   % centered sub-frame whisk freq by breath class (basal/sniff)
ampB = []; ampS = [];                          % per-cycle whisk amplitude (deg), aligned with whzB/whzS

ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE)), continue; end       % skip excluded animals
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if is_excluded_sess(ad(a).name, kk, EXCLUDE_SESS), continue; end
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk)))), continue; end
        % whisk (BP)
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(a).name,kk)), 0.6);    % lik<0.6 -> linear interp
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;     % whisk envelope (epoch detection)
        % breath -> per-breath 1/ITI on whisk grid
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(a).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
        catch
            continue
        end
        [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        [~,iloc]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        tInsp=(iloc-1)/fb; bfreq=1./diff(tInsp);
        freq = [freq; bfreq(:)]; nsess = nsess+1; %#ok<AGROW>   % for the histogram
        fbr = nan(size(t)); bin = discretize(t, tInsp); ok = ~isnan(bin) & bin<=numel(bfreq);
        fbr(ok) = bfreq(bin(ok));                      % per-sample breath rate
        ep = detect(env, t, ENV_THR, ENV_SMOOTH, MIN_DUR, MERGE_GAP, fpsW);   % whisk epochs
        whmask = false(size(t)); for q=1:size(ep,1), whmask = whmask | (t>=ep(q,1)&t<=ep(q,2)); end
        % whisking PSD + protraction-onset ITI over the SAME contiguous basal/sniff
        % runs (>= MINRUN samples), WITHIN whisk epochs only
        for c = 1:2
            if c==1, cls = whmask & ~isnan(fbr) & fbr <  BASAL_HZ;
            else,    cls = whmask & ~isnan(fbr) & fbr >  SNIFF_HZ; end
            dd = diff([0; cls(:); 0]); rs = find(dd==1); re = find(dd==-1)-1;
            for j = 1:numel(rs)
                seg = rs(j):re(j); if numel(seg) < MINRUN, continue; end
                win = hamming(min(numel(seg), PSD_WIN)); nov = min(PSD_OVERLAP, numel(win)-1);
                pw = (pwelch(xL(seg),win,nov,PSD_NFFT,fpsW) + pwelch(xR(seg),win,nov,PSD_NFFT,fpsW))'/2;
                Pw(c,:)=Pw(c,:)+pw; Pw2(c,:)=Pw2(c,:)+pw.^2; nPw(c)=nPw(c)+1;
                % Row 3: centered local whisk frequency from sub-frame protraction-onset timing.
                for s = {xL, xR}
                    xsg = s{1}(seg);
                    [~,iloc] = findpeaks(-xsg,'MinPeakProminence',PROT_PROM*std(s{1}), ...
                        'MinPeakDistance',round(0.02*fpsW));
                    if numel(iloc) < 3, continue; end

                    % sub-frame protraction-onset timing inside this run
                    pLoc = subsamp(-xsg, iloc);
                    tOn  = (seg(1)-1)/fpsW + (pLoc(:)-1)/fpsW;

                    % centered local frequency and cycle amplitude
                    fw  = nan(numel(tOn),1);
                    amp = nan(numel(tOn),1);

                    for k = 2:numel(tOn)-1
                        dprev = tOn(k)   - tOn(k-1);
                        dnext = tOn(k+1) - tOn(k);

                        if dprev > 0 && dnext > 0 && dprev <= ITI_MAX && dnext <= ITI_MAX
                            per = (dprev + dnext)/2;
                            fw(k) = 1/per;
                        end

                        i1 = max(1, round(pLoc(k)));
                        i2 = min(numel(xsg), round(pLoc(k+1)));
                        if i2 > i1
                            amp(k) = max(xsg(i1:i2)) - min(xsg(i1:i2));
                        end
                    end

                    keep = isfinite(fw) & isfinite(amp);
                    if c==1, whzB=[whzB; fw(keep)]; ampB=[ampB; amp(keep)]; %#ok<AGROW>
                    else,    whzS=[whzS; fw(keep)]; ampS=[ampS; amp(keep)]; end %#ok<AGROW>
                end
            end
        end
    end
end
freq = freq(isfinite(freq) & freq<=EDGES(end));
fprintf('%d breaths | basal<%g %.0f%% | mid %.0f%% | sniff>%g %.0f%%\n', numel(freq), ...
    BASAL_HZ, 100*mean(freq<BASAL_HZ), 100*mean(freq>=BASAL_HZ&freq<=SNIFF_HZ), SNIFF_HZ, 100*mean(freq>SNIFF_HZ));
fprintf('whisk PSD runs: basal=%d, sniff=%d\n', nPw(1), nPw(2));
fprintf('whisk freq centered/subframe: basal n=%d (med %.1f Hz) | sniff n=%d (med %.1f Hz)\n', ...
    numel(whzB), median(whzB,'omitnan'), numel(whzS), median(whzS,'omitnan'));

figure('Color','w','Position',[100 30 620 1100]);
% ---- TOP: breath-frequency histogram, classified ----
ax1=subplot(4,1,1); hold(ax1,'on'); grid(ax1,'on');
histogram(ax1, freq(freq>=BASAL_HZ&freq<=SNIFF_HZ), EDGES, 'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
histogram(ax1, freq(freq<BASAL_HZ),                 EDGES, 'FaceColor',[0 0 0],      'EdgeColor','none');
histogram(ax1, freq(freq>SNIFF_HZ),                 EDGES, 'FaceColor',[0.85 0.1 0.1],'EdgeColor','none');
%xline(ax1,BASAL_HZ,'k--'); xline(ax1,SNIFF_HZ,'r--'); xlim(ax1,[EDGES(1) EDGES(end)]);
xlabel(ax1,'instantaneous breathing frequency (Hz)'); ylabel(ax1,'number of breaths');
legend(ax1, {sprintf('unclassified (%.0f%%)',100*mean(freq>=BASAL_HZ&freq<=SNIFF_HZ)), ...
    sprintf('basal <%g Hz (%.0f%%)',BASAL_HZ,100*mean(freq<BASAL_HZ)), ...
    sprintf('sniff >%g Hz (%.0f%%)',SNIFF_HZ,100*mean(freq>SNIFF_HZ))}, 'Location','northeast');
title(ax1, sprintf('breathing frequency  (%d breaths, %d sessions)', numel(freq), nsess));

% ---- MIDDLE: whisking PSD during basal vs sniffing ----
ax2=subplot(4,1,2); hold(ax2,'on'); grid(ax2,'on'); cols={[0 0 0],[0.85 0.1 0.1]}; nm={'basal','sniff'};
h=gobjects(0); lab={};
for c=1:2
    if nPw(c)<1, continue; end
    mu=Pw(c,:)/nPw(c); sd=sqrt(max(Pw2(c,:)/nPw(c)-mu.^2,0));
    below=fpsd<BP(1); mu(below)=NaN; sd(below)=NaN;   % signal band-passed BP(1)-BP(2); PSD undefined below
    fk=fpsd(~below); muk=mu(~below); sdk=sd(~below);
    patch(ax2,[fk fliplr(fk)],[muk+sdk fliplr(max(muk-sdk,0))],cols{c},'FaceAlpha',0.15,'EdgeColor','none');
    hh=plot(ax2,fpsd,mu,'-','Color',cols{c},'LineWidth',1.6);
    in=fpsd>=max(PSD_FPASS(1),BP(1))&fpsd<=PSD_FPASS(2); ff=fpsd(in); mm=mu(in); [~,ip]=max(mm);
    h(end+1)=hh; lab{end+1}=sprintf('%s n=%d (pk %.1f Hz)',nm{c},nPw(c),ff(ip)); %#ok<AGROW>
end
xlim(ax2,[0 30]); xlabel(ax2,'frequency (Hz)'); ylabel(ax2,'whisking power');
if ~isempty(h), legend(ax2,h,lab,'Location','northeast'); end
title(ax2, sprintf('whisking power (whisk epochs, thr %g deg): basal (black) vs sniffing (red)  [BP %g-%g Hz]',ENV_THR,BP(1),BP(2)));

% ---- THIRD: centered whisk frequency from sub-frame protraction-onset timing ----
ax3=subplot(4,1,3); hold(ax3,'on'); grid(ax3,'on');
histogram(ax3, whzB, WHZ_EDGES, 'Normalization','probability','FaceColor',[0 0 0],      'EdgeColor','none','FaceAlpha',0.6);
histogram(ax3, whzS, WHZ_EDGES, 'Normalization','probability','FaceColor',[0.85 0.1 0.1],'EdgeColor','none','FaceAlpha',0.6);
xlim(ax3,[WHZ_EDGES(1) WHZ_EDGES(end)]); xlabel(ax3,'centered whisk frequency from sub-frame onset timing (Hz)'); ylabel(ax3,'prob');
legend(ax3, {sprintf('basal n=%d (med %.1f Hz)',numel(whzB),median(whzB,'omitnan')), ...
    sprintf('sniff n=%d (med %.1f Hz)',numel(whzS),median(whzS,'omitnan'))}, 'Location','northeast');
title(ax3,'whisking frequency from centered sub-frame protraction-onset timing');

% ---- BOTTOM: mean whisk amplitude vs whisk-frequency bin ----
ax4=subplot(4,1,4); hold(ax4,'on'); grid(ax4,'on');
ctr = WHZ_EDGES(1:end-1) + diff(WHZ_EDGES)/2;
[mB,eB] = binstat(whzB, ampB, WHZ_EDGES);
[mS,eS] = binstat(whzS, ampS, WHZ_EDGES);
hB=errorbar(ax4, ctr, mB, eB, '-', 'Color',[0 0 0],      'LineWidth',1.4, 'CapSize',0);
hS=errorbar(ax4, ctr, mS, eS, '-', 'Color',[0.85 0.1 0.1],'LineWidth',1.4, 'CapSize',0);
xlim(ax4,[WHZ_EDGES(1) WHZ_EDGES(end)]); xlabel(ax4,'centered whisk frequency from sub-frame onset timing (Hz)'); ylabel(ax4,'mean whisk amplitude (deg)');
legend(ax4,[hB hS],{'basal','sniff'},'Location','northeast');
title(ax4, 'mean whisk-cycle amplitude (peak-to-peak) per frequency bin  \pm SEM');
end

function [m,e] = binstat(x, y, edges)
% mean +/- SEM of y in each bin of x (defined by edges)
    b = discretize(x(:), edges); n = numel(edges)-1; m = nan(n,1); e = nan(n,1);
    for i = 1:n
        yi = y(b==i); if isempty(yi), continue; end
        m(i) = mean(yi); e(i) = std(yi)/sqrt(numel(yi));
    end
end

% ================= helpers =================
function ep = detect(env, t, envThr, envSmooth, minDur, mergeGap, fps)
    a = movmean(env(:), max(1, round(envSmooth*fps))) > envThr;   % smoothed hard threshold (deg)
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
function tf = is_excluded_sess(animal, kRun, excl)
    tf = false;
    for q = 1:numel(excl)
        if strcmp(char(excl{q}{1}), char(animal)) && double(excl{q}{2}) == kRun
            tf = true; return;
        end
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
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx > 1 & idx < numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp; delta = 0.5*(ym - yp) ./ den;
    delta(~isfinite(delta) | abs(delta) > 0.5) = 0;
    p(in) = i + delta;
end
