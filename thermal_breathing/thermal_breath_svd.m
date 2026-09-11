function thermal_breath_svd()
% thermal_breath_svd  Breathing readout from a FLIR thermal (.ats) sequence.
%
% Consumes the compact *_thermal.mat produced by thermal_ats_to_mat.py
% (radiometric temperature [T x H x W] in deg C, with the true fps recovered
% from the camera's per-frame timestamps -- nothing hardcoded). For each file:
%   1. SVD of mean-subtracted (per-pixel) temperature -> spatial eigenimages x PCs.
%   2. PC1 (sign-stabilized) = the breathing trace; auto-flag which PC is most
%      breath-like by breath-band power.
%   3. Per-pixel breath-band POWER MAP -- the QC image showing WHERE the thermal
%      breathing signal lives (nostrils / airflow), independent of the SVD.
%
% Mirrors the SVD / PSD / sign conventions of breath_svd_pc1.m, which is a
% SEPARATE script and is NOT modified by this one.
%
% Per file it writes:  *_thermal_svd.mat  +  *_thermal_svd.png  (next to input).

%% ----------------------------- USER PARAMETERS -----------------------------
MAT_PATH = '';            % '' = glob ROOT_DIR for *_thermal.mat; else a single file
ROOT_DIR = 'D:\260611_thermal_breathing';
nSVD     = 10;            % SVD components to compute
fBreath  = [0.3 12];      % breath band (Hz): PSD shading, breath-PC pick, power map
SVD_MAX_BYTES = 1.5e9;    % auto spatial-bin the SVD matrix below this (double)
%% ---------------------------------------------------------------------------

repoRoot = fileparts(fileparts(mfilename('fullpath')));  % script lives in thermal_breathing/
addpath(repoRoot);
chronuxDir = fullfile(repoRoot, 'chronux_2_12');
if isfolder(chronuxDir), addpath(genpath(chronuxDir)); end

if isempty(MAT_PATH)
    d = dir(fullfile(ROOT_DIR, '**', '*_thermal.mat'));
    files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'uni', 0);
else
    files = {MAT_PATH};
end
assert(~isempty(files), 'No *_thermal.mat found.');
fprintf('thermal_breath_svd: %d file(s) to process.\n', numel(files));

summary = {};
for fi = 1:numel(files)
    matPath = files{fi};
    [outDir, stem] = fileparts(matPath);
    fprintf('\n===== [%d/%d] %s =====\n', fi, numel(files), stem);

    S = load(matPath);
    fps = double(S.fps);
    H = double(S.H); W = double(S.W);
    stack = single(S.stack);                 % [T x H x W] deg C
    T = size(stack, 1);
    if isfield(S, 't'), t = double(S.t(:)); else, t = (0:T-1)'/fps; end
    fprintf('  %d frames, %dx%d px, fps=%.3f Hz (nominal %.0f, %d dropped of %d)\n', ...
        T, H, W, fps, getfielddef(S, 'fps_nominal', fps), ...
        getfielddef(S, 'n_dropped', 0), getfielddef(S, 'n_expected', T));

    % [P x T] temperature, mean-subtracted per pixel (remove static thermal map)
    X = reshape(permute(stack, [2 3 1]), H*W, T);   % [P x T]
    Xc = X - mean(X, 2);

    % --- per-pixel breath-band POWER MAP (QC: where breathing lives) ---
    powMap = bandpower_map(Xc, fps, fBreath);       % [P x 1]
    powImg = reshape(powMap, H, W);
    meanImg = reshape(mean(X, 2), H, W);            % static mean temperature

    % --- SVD of mean-subtracted temperature (spatial-bin if huge) ---
    sbin = 1;
    while (H*W / sbin^2) * T * 8 > SVD_MAX_BYTES, sbin = sbin + 1; end
    cube = reshape(Xc, H, W, T);
    if sbin > 1, cube = binCube(cube, sbin); end
    [Hb, Wb, ~] = size(cube);
    Yd = double(reshape(cube, Hb*Wb, T));
    Yd = Yd - mean(Yd, 2);
    k  = min(nSVD, min(size(Yd)) - 1);
    [U, Sv, V] = svds(Yd, k);
    sv = diag(Sv); varExp = sv.^2 / sum(sv.^2);

    % --- which PC is breathing? sign-stabilized PC1 as the trace ---
    [pcBreath, ~] = pick_breath_pc(V, sv, fps, fBreath);
    pc = pcBreath;                                  % use the breath-band winner
    sgn = sign(sum(U(:,pc).^3) + eps);
    breathTrace = sgn * sv(pc) * V(:,pc);
    eigImg = sgn * reshape(U(:,pc), Hb, Wb);
    [fa, Pa] = mt_psd(breathTrace, fps);
    inb = fa >= fBreath(1) & fa <= fBreath(2);
    [~, jp] = max(Pa .* inb); pcPeak = fa(jp);

    % --- figure: power map | eigenimage | PSD | trace ---
    f = figure('Name', ['thermal breath: ' stem], 'Color', 'w', ...
               'Position', [60 60 1200 760]);
    subplot(2,3,1);
    imagesc(meanImg); axis image off; colormap(gca, gray); colorbar;
    title('mean temperature (\circC)');
    subplot(2,3,2);
    imagesc(powImg); axis image off; colormap(gca, hot); colorbar;
    title(sprintf('breath-band power  [%.1f-%.1f Hz]', fBreath(1), fBreath(2)));
    subplot(2,3,3);
    imagesc(eigImg); axis image off; colormap(gca, parula); colorbar;
    title(sprintf('PC%d eigenimage (%.0f%% var)', pc, 100*varExp(pc)));
    subplot(2,3,4);
    plot(fa, Pa, 'k', 'LineWidth', 1.3); grid on; hold on;
    area(fa(inb), Pa(inb), 'FaceColor', [1 .8 .8], 'EdgeColor', 'none');
    plot(fa, Pa, 'k', 'LineWidth', 1.3);
    xline(pcPeak, 'r--', sprintf('%.2f Hz', pcPeak));
    xlim([0 min(20, fps/2)]); xlabel('frequency (Hz)'); ylabel('power');
    title(sprintf('PC%d PSD', pc));
    subplot(2,3,[5 6]);
    plot(t, breathTrace, 'LineWidth', 0.7); grid on; xlim([t(1) t(end)]);
    xlabel('time (s)'); ylabel(sprintf('PC%d', pc));
    title(sprintf('PC%d temporal trace = breathing readout (%.2f Hz)', pc, pcPeak));
    sgtitle(sprintf('%s    fps=%.2f Hz', strrep(stem, '_', '\_'), fps));

    % --- save ---
    outMat = fullfile(outDir, [stem '_svd.mat']);
    outPng = fullfile(outDir, [stem '_svd.png']);
    save(outMat, 'breathTrace', 't', 'fps', 'pc', 'pcBreath', 'pcPeak', ...
         'eigImg', 'powImg', 'meanImg', 'sv', 'varExp', 'U', 'V', ...
         'fBreath', 'sbin', '-v7.3');
    try
        exportgraphics(f, outPng, 'Resolution', 150);
    catch
        saveas(f, outPng);
    end

    summary(end+1, :) = {stem, fps, pcPeak, pc, 100*varExp(pc)}; %#ok<AGROW>
    fprintf('  breathing PC%d, peak %.2f Hz, %.0f%% var. Saved *_svd.mat/.png\n', ...
        pc, pcPeak, 100*varExp(pc));
end

%% summary table
fprintf('\n=============== SUMMARY ===============\n');
fprintf('%-45s %7s %8s %5s %7s\n', 'file', 'fps', 'breathHz', 'PC', 'var%');
for i = 1:size(summary, 1)
    fprintf('%-45s %7.2f %8.2f %5d %6.0f%%\n', summary{i,1}(1:min(45,end)), ...
        summary{i,2}, summary{i,3}, summary{i,4}, summary{i,5});
end
end

%% =============================== LOCAL FUNCTIONS ===============================
function v = getfielddef(S, name, def)
if isfield(S, name), v = double(S.(name)); else, v = def; end
end

function P = bandpower_map(X, Fs, band)
% Per-row (pixel) power in [band] via one FFT of the whole [P x T] matrix.
n = size(X, 2);
w = hann(n)';                          % taper each pixel's time series
Xw = X .* w;
F = fft(Xw, [], 2);
nf = floor(n/2) + 1;
Pf = (abs(F(:, 1:nf)).^2) / (Fs * sum(w.^2));   % one-sided PSD, per pixel
f = (0:nf-1) * (Fs / n);
m = f >= band(1) & f <= band(2);
P = sum(Pf(:, m), 2);                  % integrated band power per pixel
end

function cb = binCube(c, b)
[H, W, T] = size(c);
H2 = floor(H/b)*b; W2 = floor(W/b)*b;
c = c(1:H2, 1:W2, :);
c = reshape(c, b, H2/b, b, W2/b, T);
cb = squeeze(mean(mean(c, 1), 3));
end

function [f, P] = mt_psd(x, Fs)
x = x(:) - mean(x(:));
if exist('mtspectrumc', 'file') == 2
    pr.Fs = Fs; pr.tapers = [3 5]; pr.pad = 0;
    pr.fpass = [0.05, min(20, Fs/2)];
    [P, f] = mtspectrumc(x, pr); f = f(:); P = P(:);
else
    n = numel(x); nf = floor(n/2) + 1;
    Pf = abs(fft(x .* hann(n))).^2;
    P = Pf(1:nf); f = (0:nf-1)' * (Fs / n);
end
end

function [pc, fpk] = pick_breath_pc(V, sv, Fs, band)
K = size(V, 2); inbPow = zeros(K, 1); pkf = zeros(K, 1);
for kk = 1:K
    [f, P] = mt_psd(sv(kk) * V(:, kk), Fs);
    m = f >= band(1) & f <= band(2);
    if any(m)
        inbPow(kk) = max(P(m));
        [~, j] = max(P .* m); pkf(kk) = f(j);
    end
end
[~, pc] = max(inbPow); fpk = pkf(pc);
end
