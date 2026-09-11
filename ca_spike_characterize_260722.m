function ca_spike_characterize_260722()
%CA_SPIKE_CHARACTERIZE_260722  Label-free survey of every ROI's signal regime.
%
% Deep-inspection pass for the automatic detector. Uses NO hand labels: every
% quantity here is computed from the trace itself, so uncurated ROIs count too
% (868 ROIs rather than the 362 that carry clicks).
%
% Per ROI it measures:
%   sn         robust noise, from the frame-to-frame difference
%   skew       skewness of dF/F. Real calcium is one-sided-positive; symmetric
%              means noise or motion, negative means something is wrong.
%   asym       (excursions above +k*sn) / (excursions below -k*sn). The key
%              label-free quality number: calcium cannot go down, so the
%              downward count is an FDR estimate for the upward count.
%   tauAC      decay time from the autocorrelation, i.e. g without needing a
%              kernel fit or any label
%   fPeak      dominant frequency in 0.05-8 Hz and the fraction of band power
%   rhythm     concentration of power at fPeak -> separates sparse transients
%              from continuous rhythmic modulation (the regime that broke the
%              AR(1) detector)
%   rate       events/min at a fixed 3*sn prominence, just for scale
%
% Writes <ROOT_DIR>\_spike_characterize_260722\ (mat + summary figure).

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = 'D:\Ventral_surface_summary';
OUT_SUB    = '_spike_characterize_260722';
K_EXC      = 3;         % sigma for the up/down excursion count
F_BAND     = [0.05 8];  % Hz: band searched for a rhythm
MIN_DIST_S = 0.2;
FALLBACK_FPS = 30;
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
warning('off','signal:findpeaks:largeMinPeakHeight');
outDir = fullfile(ROOT_DIR, OUT_SUB);
if ~isfolder(outDir), mkdir(outDir); end

d = dir(fullfile(ROOT_DIR, '**', 'ca_spike_data.mat'));   % locates the sessions
assert(~isempty(d), 'No sessions under %s', ROOT_DIR);

Q = struct('session',{},'group',{},'roi',{},'curated',{},'nLab',{}, ...
           'sn',{},'skew',{},'asym',{},'tauAC',{},'fPeak',{},'rhythm',{}, ...
           'rate',{},'amp90',{},'T',{},'fps',{});

for i = 1:numel(d)
    fo  = d(i).folder;
    rel = fo(numel(ROOT_DIR)+2:end);
    grp = strtok(rel, filesep);
    if strcmp(grp,'Vglut2_test'), grp = 'Vglut2'; end

    dh = dir(fullfile(fo,'*_dFF.mat'));
    if isempty(dh), continue; end
    L = load(fullfile(dh(1).folder,dh(1).name),'dFF');
    if ~isfield(L,'dFF'), continue; end
    dFF = L.dFF; [T,N] = size(dFF);
    fps = detect_session_fps(fo, FALLBACK_FPS);

    rs = [];
    Sp = load(fullfile(fo,'ca_spike_data.mat'),'roi_spikes');
    if numel(Sp.roi_spikes) == N, rs = Sp.roi_spikes; end

    for k = 1:N
        x = dFF(:,k);
        if ~all(isfinite(x)) || std(x) <= 0, continue; end
        sn = median(abs(diff(x)))*1.4826/sqrt(2);
        if sn <= 0, sn = std(x); end
        xz = (x - median(x)) / sn;

        % one-sidedness: up vs down excursion counts at the same threshold
        md = max(1, round(MIN_DIST_S*fps));
        nUp = numel(pk_idx( xz, K_EXC, md));
        nDn = numel(pk_idx(-xz, K_EXC, md));
        asym = (nUp + 1) / (nDn + 1);      % +1 keeps it finite when nDn = 0

        % autocorrelation decay -> tau, no kernel fit, no labels
        tauAC = ac_tau(x, fps);

        % rhythmicity
        [fPk, rhy] = rhythm_score(x, fps, F_BAND);

        nLab = 0;
        if ~isempty(rs) && ~isempty(rs(k).spike_idx), nLab = numel(rs(k).spike_idx); end

        Q(end+1) = struct('session',rel,'group',grp,'roi',k, ...
            'curated', nLab>0, 'nLab',nLab, 'sn',sn, 'skew',skewness(x), ...
            'asym',asym, 'tauAC',tauAC, 'fPeak',fPk, 'rhythm',rhy, ...
            'rate', nUp/(T/fps)*60, 'amp90', prctile(x,90)/sn, ...
            'T',T, 'fps',fps); %#ok<AGROW>
    end
    fprintf('%-52s %3d ROIs\n', rel(1:min(52,end)), N);
end

fprintf('\n%d ROIs characterised (%d curated, %d not)\n', ...
    numel(Q), sum([Q.curated]), sum(~[Q.curated]));

%% ------------------------------ summary ------------------------------------
asym = [Q.asym]; rhy = [Q.rhythm]; rate = [Q.rate];
tau  = [Q.tauAC]; skw = [Q.skew];

fprintf('\n--- one-sidedness (up/down excursions at %g sigma) ---\n', K_EXC);
fprintf('percentiles 10/25/50/75/90: %s\n', mat2str(round(prctile(asym,[10 25 50 75 90]),2)));
fprintf('ROIs with asym < 1.5 (no calcium-like asymmetry -> likely silent/noise): %d (%.0f%%)\n', ...
    sum(asym<1.5), 100*mean(asym<1.5));

fprintf('\n--- rhythmicity ---\n');
fprintf('percentiles 10/25/50/75/90: %s\n', mat2str(round(prctile(rhy,[10 25 50 75 90]),3)));
fprintf('ROIs with rhythm > 0.25 : %d (%.0f%%)\n', sum(rhy>0.25), 100*mean(rhy>0.25));

fprintf('\n--- autocorrelation tau (s) by group ---\n');
gs = unique({Q.group});
for gi = 1:numel(gs)
    m = strcmp({Q.group}, gs{gi}) & isfinite(tau);
    fprintf('  %-8s n=%4d  tau 25/50/75: %s\n', gs{gi}, sum(m), ...
        mat2str(round(prctile(tau(m),[25 50 75]),2)));
end

fprintf('\n--- activity regime (rate at %g sigma) ---\n', K_EXC);
edges = [0 2 10 20 40 80 1e9];
for b = 1:numel(edges)-1
    m = rate>=edges(b) & rate<edges(b+1);
    if ~any(m), continue; end
    fprintf('  %4g-%4g /min : %4d ROIs (%4.1f%%)  median rhythm %.2f  median asym %5.1f\n', ...
        edges(b), edges(b+1), sum(m), 100*mean(m), median(rhy(m)), median(asym(m)));
end

%% ------------------------------ figure -------------------------------------
f = figure('Color','w','Position',[50 50 1500 460]);
subplot(1,3,1);
histogram(log10(max(asym,0.1)), 40); grid on; box on;
xlabel('log_{10} up/down excursion ratio'); ylabel('ROIs');
title(sprintf('one-sidedness (median %.1f)', median(asym)));
xline(0,'r--','symmetric');

subplot(1,3,2);
scatter(rate, rhy, 12, log10(max(asym,0.1)), 'filled'); grid on; box on;
set(gca,'XScale','log'); colormap(gca,parula); cb=colorbar; cb.Label.String='log_{10} asym';
xlabel('event rate /min (3\sigma)'); ylabel('rhythmicity');
title('regime map');

subplot(1,3,3);
histogram(tau(isfinite(tau)), 0:0.05:2.5); grid on; box on;
xlabel('autocorrelation \tau (s)'); ylabel('ROIs');
title(sprintf('label-free \\tau (median %.2f s)', median(tau,'omitnan')));

sgtitle(sprintf('Label-free ROI survey: %d ROIs, %d sessions', numel(Q), numel(d)));
exportgraphics(f, fullfile(outDir,'characterize_summary.png'), 'Resolution',150);
save(fullfile(outDir,'characterize.mat'), 'Q', 'K_EXC', 'F_BAND', '-v7.3');
fprintf('\nSaved: %s\n', outDir);
end

%% =============================== LOCAL FUNCTIONS ===============================
function loc = pk_idx(z, thr, minDistFr)
    [~, loc] = findpeaks(z, 'MinPeakHeight', thr, 'MinPeakDistance', minDistFr);
    loc = loc(:);
end

function tau = ac_tau(x, fps)
% AR(1) decay time from the autocovariance, label-free and kernel-free.
%
% Must NOT use lag 0: for calcium + white noise, c(0) = var(c) + var(noise)
% while c(L) = var(c)*g^L for L >= 1. Lag 0 is the only lag the noise touches,
% so including it makes the ACF collapse within one frame and reports tau ~ 1
% frame regardless of the true decay (which is what the first version did).
% Taking the ratio across lags >= 1 cancels the noise entirely.
    x = x - movmedian(x, round(20*fps));       % kill slow drift
    x = x - mean(x);
    L = min(max(3, round(0.6*fps)), floor(numel(x)/4));   % lags 1..L
    c = xcorr(x, L, 'biased');
    c = c(L+1:end);                             % lags 0..L
    c1 = c(2:end);                              % lags 1..L, noise-free
    ok = c1 > 0;
    if sum(ok) < 3, tau = NaN; return; end
    % log c(L) = log(var) + L*log(g)  ->  slope gives g
    lags = (1:numel(c1))';
    p = polyfit(lags(ok), log(c1(ok)), 1);
    g = exp(p(1));
    if ~isfinite(g) || g <= 0 || g >= 1, tau = NaN; return; end
    tau = -1/(fps*log(g));
    if tau <= 0 || tau > 5, tau = NaN; end
end

function [fPk, score] = rhythm_score(x, fps, band)
% Dominant in-band frequency and how concentrated the power is around it.
% score = power within +-25% of the peak / total in-band power. Sparse
% transients spread power broadly (low score); a continuous oscillation
% concentrates it (high score).
    x = x - mean(x);
    n = numel(x);
    w = hann(n);
    P = abs(fft(x.*w)).^2;
    nf = floor(n/2)+1;
    P = P(1:nf); fr = (0:nf-1)'*(fps/n);
    m = fr >= band(1) & fr <= band(2);
    if ~any(m), fPk = NaN; score = NaN; return; end
    Pb = P; Pb(~m) = 0;
    [~, j] = max(Pb); fPk = fr(j);
    nearBand = fr >= 0.75*fPk & fr <= 1.25*fPk & m;
    score = sum(P(nearBand)) / max(sum(P(m)), eps);
end
