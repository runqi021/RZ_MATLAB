function whisker_crop_track(cropMat)
% whisker_crop_track  Track the bright whisker marker in each cropped ROI from
% whisker_crop_extract.py and produce a 1-D whisking trace per side, plus the
% bilateral (L vs R) cross-correlation / coherence.
%
% Method (per side): isolate the MOVING marker as max(frame - meanImage, 0) so
% static bright structures (eye, fur) cancel; take its intensity-weighted
% centroid each frame -> (cx,cy); project onto the principal axis of the centroid
% cloud -> whisk(t) (the sweep). Confidence = per-frame marker energy.
%
%   whisker_crop_track('C:\...\cam1_..._run001_whiskcrop.mat')

S = load(cropMat);
[outDir, stem] = fileparts(cropMat);
names = cellstr(S.names);
fps = double(S.fps);
t = double(S.t_s(:));

fWhisk = [4 25];        % plausible whisking/sniff-coupled band (Hz) for PSD shading
res = struct();
f = figure('Color','w','Position',[40 40 1280 820],'Visible','off');
tl = tiledlayout(f, numel(names)+2, 3, 'TileSpacing','compact','Padding','compact');

BRIGHT_PCT = 96;        % absolute brightness percentile for the marker pixels
for s = 1:numel(names)
    nm = names{s};
    mov = single(S.(['mov_' nm]));          % [H W T]
    [H,W,T] = size(mov);
    X = reshape(mov, H*W, T);                % [P T]
    mu = mean(X, 2);
    stdv = std(X,0,2);
    [Xc,Yc] = meshgrid(1:W, 1:H);
    px = Xc(:); py = Yc(:);

    % moving marker only = max(frame - mean, 0). Its intensity-weighted centroid
    % is a convex combination of pixel coords -> ALWAYS inside the ROI, so it can
    % never flip; projecting it on the sweep axis gives a bounded whisking trace
    % (noisy at rest = not whisking; clean oscillation during bouts).
    D = max(X - mu, 0);                       % [P T]
    wsum = sum(D,1) + eps;                    % [1 T] marker motion energy
    cx = (px' * D) ./ wsum;
    cy = (py' * D) ./ wsum;
    C = [cx - mean(cx); cy - mean(cy)];       % [2 T]
    [V,~] = eig(cov(C'));
    pc1 = V(:,end);                           % sweep axis (largest variance)
    whisk = (pc1' * C)';                      % [T 1] sweep position (px)
    whisk = whisk - median(whisk);

    res.(nm).cx = cx(:); res.(nm).cy = cy(:);
    res.(nm).whisk = whisk(:); res.(nm).energy = wsum(:);
    res.(nm).meanImg = reshape(mu,H,W);
    res.(nm).stdImg = reshape(stdv,H,W);

    % --- row of panels for this side ---
    nexttile; imagesc(res.(nm).stdImg); axis image off; colormap(gca,parula);
    hold on; plot(cx(1:20:end), cy(1:20:end), '.r','MarkerSize',2);
    title(sprintf('%s: std + sweep-centroid track', nm),'Interpreter','none');

    nexttile;
    % auto-pick the 4 s window with the most whisking (highest angle variance)
    win = round(4*fps); step=round(0.5*fps);
    v = movvar(whisk, win); [~,ic]=max(v); i0=max(1,ic-round(win/2)); i1=min(T,i0+win);
    plot(t(i0:i1), whisk(i0:i1), 'k'); xlabel('s'); ylabel('sweep (px)');
    title(sprintf('%s whisking (busiest 4 s @ %.0fs)', nm, t(i0)),'Interpreter','none');

    nexttile;
    [P,fr] = pwelch(detrend(whisk), hann(min(8192,T)), [], [], fps);
    plot(fr, 10*log10(P+eps), 'k'); xlim([0 40]); hold on;
    yl=ylim; patch([fWhisk fliplr(fWhisk)],[yl(1) yl(1) yl(2) yl(2)],[1 .9 .6],'EdgeColor','none','FaceAlpha',.3);
    plot(fr,10*log10(P+eps),'k'); xlabel('Hz'); ylabel('dB');
    [~,ip]=max(P.*(fr>fWhisk(1)&fr<fWhisk(2)));
    title(sprintf('%s PSD (peak %.1f Hz)', nm, fr(ip)),'Interpreter','none');
end

% --- bilateral L vs R (use first two sides) ---
a = double(res.(names{1}).whisk); b = double(res.(names{2}).whisk);
n = min(numel(a),numel(b)); a=zscore(a(1:n)); b=zscore(b(1:n));
bpf = @(x) filtfilt(butter(2,fWhisk/(fps/2),'bandpass'), 1, x);
af=bpf(a); bf=bpf(b);
[c,lags]=xcorr(af,bf,round(0.06*fps),'coeff');   % +/-60 ms (< half a whisk cycle)
[rho,im]=max(c); lag=lags(im);
[Cxy,fc]=mscohere(a,b,hann(8192),4096,8192,fps);

nexttile([1 2]);
tt=(0:n-1)/fps; plot(tt,af,'b'); hold on; plot(tt,bf,'r');
v2 = movvar(af,round(4*fps))+movvar(bf,round(4*fps)); [~,jc]=max(v2);
w0=max(1,jc-round(2*fps)); xlim([tt(w0) tt(w0)+4]);
legend(names(1:2)); title('Bilateral whisking (band-passed, busiest 4 s)'); xlabel('s');
nexttile;
plot(lags*1e3/fps,c,'k'); xline(lag*1e3/fps,'r'); xlabel('lag ms (L-R)'); ylabel('corr');
title(sprintf('L-R xcorr peak %.2f @ %.0f ms', rho, lag*1e3/fps));

nexttile([1 3]);
plot(fc,Cxy,'b'); xlim([0 40]); xlabel('Hz'); ylabel('coherence');
[~,ic]=max(Cxy.*(fc>fWhisk(1)&fc<fWhisk(2)));
title(sprintf('L-R coherence (peak %.2f @ %.1f Hz)', Cxy(ic), fc(ic)));

title(tl, sprintf('Whisker tracking  |  %s', stem),'Interpreter','none');
outPng = fullfile(outDir,[stem '_track.png']);
exportgraphics(f,outPng,'Resolution',150); close(f);

whisk_result = res; whisk_result_meta = struct('fps',fps,'lag_LR_ms',lag*1e3/fps,'rho_LR',rho); %#ok<NASGU>
save(fullfile(outDir,[stem '_track.mat']),'whisk_result','whisk_result_meta','t','names');
fprintf('L-R whisking: peak corr %.2f at %.0f ms; coherence peak %.2f @ %.1f Hz\n', rho, lag*1e3/fps, Cxy(ic), fc(ic));
fprintf('saved %s\n', outPng);
end
