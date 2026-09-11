function dlc_whisker_bilateral(csvPath)
% dlc_whisker_bilateral  Whisker trajectory + bilateral (L vs R) coordination from
% a DLC csv with dots vL0,vL1 (left whisker) and vR0,vR1 (right whisker).
%
% Whisker angle per side = orientation of the base->tip vector (vL0->vL1, vR0->vR1).
% Plots dot trajectories, the two angle traces, their spectra, and the L-R
% cross-correlation + coherence (bilateral whisking synchrony).
%
%   dlc_whisker_bilateral('C:\...\xxxDLC_Resnet50_...csv')

[folder, stem] = fileparts(csvPath);
PMIN = 0.5;            % DLC likelihood gate
BP   = [4 30];        % whisking band (Hz)

% ---- fps from the run's timestamps.csv (camera clock), else 400 ----
fps = 400;
ts = fullfile(folder, 'timestamps.csv');
if isfile(ts)
    A = readmatrix(ts); cam = A(:,2); fps = (numel(cam)-1)/((cam(end)-cam(1))/1e9);
end

% ---- read DLC csv (3 header rows) ----
M = readmatrix(csvPath, 'NumHeaderLines', 3);
get = @(c) M(:, c);
xyL0=[get(2) get(3)]; pL0=get(4);
xyL1=[get(5) get(6)]; pL1=get(7);
xyR0=[get(8) get(9)]; pR0=get(10);
xyR1=[get(11) get(12)]; pR1=get(13);
T = size(M,1); t = (0:T-1)'/fps;

% gate low-likelihood points -> NaN -> interpolate
gate=@(xy,p) fillmissing(xy + 0./[p>=PMIN p>=PMIN], 'linear');
xyL0=gate(xyL0,pL0); xyL1=gate(xyL1,pL1); xyR0=gate(xyR0,pR0); xyR1=gate(xyR1,pR1);

% whisking signal = TIP dot position projected on its own sweep axis.
% (Robust: a linear projection of coordinates can't flip, unlike the base->tip
% angle, which jumps +-180 deg whenever a dot's likelihood dips.)
thL = proj_sweep(xyL1); thR = proj_sweep(xyR1);   % left/right whisker tip sweep (px)

% ---- bilateral coupling ----
bpf=@(x) filtfilt(butter(2,BP/(fps/2),'bandpass'),1,double(x));
zL=zscore(bpf(thL)); zR=zscore(bpf(thR));
[c,lags]=xcorr(zL,zR,round(0.06*fps),'coeff'); [rho,im]=max(abs(c)); lag=lags(im);
[Cxy,fc]=mscohere(thL,thR,hann(min(8192,T)),[],[],fps);
[~,ic]=max(Cxy.*(fc>BP(1)&fc<BP(2)));

% ---- figure ----
f=figure('Color','w','Position',[60 60 1200 760],'Visible','off');
tl=tiledlayout(f,3,3,'TileSpacing','compact','Padding','compact');

nexttile([1 1]); hold on; ss=1:10:T;   % dot trajectories (image coords, y down)
plot(xyL0(ss,1),xyL0(ss,2),'.','Color',[0 .6 0]); plot(xyL1(ss,1),xyL1(ss,2),'.','Color',[0 1 0]);
plot(xyR0(ss,1),xyR0(ss,2),'.','Color',[.7 0 0]); plot(xyR1(ss,1),xyR1(ss,2),'.','Color',[1 .5 0]);
set(gca,'YDir','reverse'); axis equal tight; legend({'vL0','vL1','vR0','vR1'},'Location','best');
title('dot trajectories');

% busiest 4 s window for display
win=round(4*fps); v=movvar(zL,win)+movvar(zR,win); [~,jc]=max(v); i0=max(1,jc-round(win/2)); i1=min(T,i0+win);
nexttile([1 2]); plot(t(i0:i1),thL(i0:i1)-median(thL(i0:i1)),'g'); hold on;
plot(t(i0:i1),thR(i0:i1)-median(thR(i0:i1)),'r');
legend({'left tip','right tip'}); xlabel('s'); ylabel('sweep (px)');
title(sprintf('whisker tip sweep (busiest 4 s @ %.0fs)',t(i0)));

[PL,fr]=pwelch(bpf(thL),hann(min(8192,T)),[],[],fps);
[PR,~]=pwelch(bpf(thR),hann(min(8192,T)),[],[],fps);
nexttile; plot(fr,10*log10(PL),'g'); hold on; plot(fr,10*log10(PR),'r'); xlim([0 40]);
[~,ipL]=max(PL.*(fr>BP(1)&fr<BP(2))); xlabel('Hz'); ylabel('dB');
title(sprintf('PSD (L peak %.1f Hz)',fr(ipL)));

nexttile; plot(lags*1e3/fps,c,'k'); hold on; xline(lag*1e3/fps,'b');
xlabel('lag ms (L-R)'); ylabel('corr'); title(sprintf('L-R xcorr %.2f @ %.0f ms',rho,lag*1e3/fps));

nexttile; plot(fc,Cxy,'b'); xlim([0 40]); xlabel('Hz'); ylabel('coherence');
title(sprintf('L-R coherence %.2f @ %.1f Hz',Cxy(ic),fc(ic)));

ax=nexttile; axis(ax,'off');
text(ax,0,1,{sprintf('fps %.2f, T %d (%.1fs)',fps,T,t(end)); '';
  sprintf('L-R tip-sweep corr %.2f',rho); sprintf('  lag %.0f ms',lag*1e3/fps);
  sprintf('coherence %.2f @ %.1f Hz',Cxy(ic),fc(ic)); '';
  sprintf('mean likelihood L %.2f R %.2f',mean([pL0;pL1]),mean([pR0;pR1]))}, ...
  'VerticalAlignment','top','FontName','Consolas','FontSize',10);

title(tl, stem, 'Interpreter','none');
outPng=fullfile(folder,[stem '_bilateral.png']);
exportgraphics(f,outPng,'Resolution',150); close(f);
save(fullfile(folder,[stem '_bilateral.mat']),'t','thL','thR','fps','rho','lag','Cxy','fc');
fprintf('L-R tip-sweep corr %.2f @ %.0f ms; coherence %.2f @ %.1f Hz\nsaved %s\n', rho, lag*1e3/fps, Cxy(ic), fc(ic), outPng);
end

function s = proj_sweep(xy)
% project a [T x 2] dot trajectory onto its largest-variance (sweep) axis -> [T x 1]
C = xy - mean(xy,1);
[V,~] = eig(cov(C));
s = C * V(:,end);
end
