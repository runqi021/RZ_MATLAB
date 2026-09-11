function whisk_epoch_analyze_RZ()
% whisk_epoch_analyze_RZ  Behavior of the SELECTED whisking epochs.
% Pools samples inside the saved epochs (across sessions, or one if ONLY_FOLDER
% set) and plots, for whisk L/R (BP) + breath (BR_BP, resampled to whisk grid):
%   row1  amplitude scatters : xL vs xR | breath vs xL | breath vs xR
%   row2  Hilbert phase-difference distributions (polar + PLV): L-R | br-L | br-R
%   row3  phase scatters     : phiL vs phiR | phiB vs phiL | phiB vs phiR

% ============================ USER-EDITABLE ============================
dataRoot   = "D:\260615_thermalNbasler";
whiskDir   = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
ONLY_FOLDER = "";       % "" = pool ALL sessions with epochs; or a run folder for one
BP        = [5 30];     % whisk bandpass (Hz)
BR_BP     = [1 15];     % breathing bandpass (Hz)
fpsW      = 400;
MAXPTS    = 8000;       % subsample for scatter display
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

epf = dir(fullfile(char(dataRoot),'*','cam1_*','*_whisk_epochs.mat'));
assert(~isempty(epf),'no *_whisk_epochs.mat found');

XL=[];XR=[];BR=[];PL=[];PR=[];PB=[]; nep=0; nsess=0;
for e = 1:numel(epf)
    rf = epf(e).folder;
    if ~isempty(char(ONLY_FOLDER)) && ~strcmp(rf,char(ONLY_FOLDER)), continue; end
    S = load(fullfile(rf, epf(e).name));
    if isempty(S.epochs), continue; end
    animal = regexp(epf(e).name,'_(\d+)_n','tokens','once'); animal = animal{1};
    kk = str2double(regexp(epf(e).name,'_n(\d+)_','tokens','once'));
    % whisk
    wcsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk));
    M = dlc_gate_interp(wcsv, 0.6);    % lik<0.6 -> linear interp
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;
    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
    phiL = angle(hilbert(xL)); phiR = angle(hilbert(xR));
    % breath, resampled to whisk grid
    brw = nan(size(t)); phiB = nan(size(t));
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kk)), dataRoot);
        if isfile(Pn.breath)
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
            [b2,a2]=butter(2,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
            brw = interp1(tBr, brf, t, 'linear', NaN);
            phiB = angle(hilbert(fillmissing(brw,'linear')));
        end
    catch
    end
    % accumulate samples inside selected epochs
    m = false(size(t));
    for q=1:size(S.epochs,1), m = m | (t>=S.epochs(q,1) & t<=S.epochs(q,2)); end
    XL=[XL;xL(m)]; XR=[XR;xR(m)]; BR=[BR;brw(m)]; PL=[PL;phiL(m)]; PR=[PR;phiR(m)]; PB=[PB;phiB(m)]; %#ok<AGROW>
    nep=nep+size(S.epochs,1); nsess=nsess+1;
end
assert(~isempty(XL),'no samples in selected epochs');
fprintf('pooled %d samples from %d epochs across %d sessions\n', numel(XL), nep, nsess);

okB = ~isnan(BR) & ~isnan(PB);          % samples with breath
i1 = subidx(numel(XL), MAXPTS);         % display subsample (all samples)
ib = find(okB); ib = ib(subidx(numel(ib), MAXPTS));   % subsample of breath samples

figure('Color','w','Position',[40 40 1280 880]);
% ---- row1: amplitude scatters ----
subplot(3,3,1); scat(XL(i1),XR(i1)); xlabel('xL'); ylabel('xR'); title(sprintf('xL vs xR  (r=%.2f)',cc(XL,XR)));
subplot(3,3,2); scat(BR(ib),XL(ib)); xlabel('breath'); ylabel('xL'); title(sprintf('breath vs xL  (r=%.2f)',cc(BR(okB),XL(okB))));
subplot(3,3,3); scat(BR(ib),XR(ib)); xlabel('breath'); ylabel('xR'); title(sprintf('breath vs xR  (r=%.2f)',cc(BR(okB),XR(okB))));
% ---- row2: phase-difference distributions (polar + PLV) ----
linhist(subplot(3,3,4), PL-PR,           '\phi_L - \phi_R');
linhist(subplot(3,3,5), PL(okB)-PB(okB), '\phi_L - \phi_{br}  (0=insp peak)');
linhist(subplot(3,3,6), PR(okB)-PB(okB), '\phi_R - \phi_{br}  (0=insp peak)');
% ---- row3: phase scatters ----
subplot(3,3,7); scat(PL(i1),PR(i1)); phaxes; xlabel('\phi_L'); ylabel('\phi_R'); title('\phi_L vs \phi_R');
subplot(3,3,8); scat(PB(ib),PL(ib)); phaxes; xlabel('\phi_{br}'); ylabel('\phi_L'); title('\phi_{br} vs \phi_L');
subplot(3,3,9); scat(PB(ib),PR(ib)); phaxes; xlabel('\phi_{br}'); ylabel('\phi_R'); title('\phi_{br} vs \phi_R');
sgtitle(sprintf('Selected whisking epochs: %d samples, %d epochs, %d sessions', numel(XL), nep, nsess));

% ================= local helpers =================
    function phaxes(), xlim([-pi pi]); ylim([-pi pi]); xticks([-pi 0 pi]); yticks([-pi 0 pi]);
        xticklabels({'-\pi','0','\pi'}); yticklabels({'-\pi','0','\pi'}); axis square; grid on; end
end

function scat(x,y), plot(x,y,'.','MarkerSize',3,'Color',[0.2 0.4 0.8]); grid on; end
function r = cc(x,y), g=~isnan(x)&~isnan(y); r=corr(x(g),y(g)); end
function linhist(ax, dphi, ttl)
    dphi = angle(exp(1i*dphi(~isnan(dphi))));
    histogram(ax, dphi, linspace(-pi,pi,37), 'Normalization','probability', ...
        'FaceColor',[0.3 0.5 0.9],'EdgeColor','none'); hold(ax,'on'); grid(ax,'on');
    plv = abs(mean(exp(1i*dphi))); mu = angle(mean(exp(1i*dphi)));
    xline(ax, mu, 'r-', 'LineWidth',1.5);
    xlim(ax,[-pi pi]); xticks(ax,[-pi 0 pi]); xticklabels(ax,{'-\pi','0','\pi'});
    xlabel(ax,'\Delta\phi'); ylabel(ax,'prob');
    title(ax, sprintf('%s  PLV=%.2f (\\mu=%+.0f\\circ)', ttl, plv, rad2deg(mu)));
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
