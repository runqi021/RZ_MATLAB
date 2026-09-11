function breath_svd_pc1()
% breath_svd_pc1  DLC-free breathing readout from the 2P-triggered Basler video.
%
% For each breath-video folder under ROOT_DIR (or a hand-listed FOLDERS set):
%   1. fps  = detect_session_fps(folder)   -- the breath cam is 2P-frame-
%             triggered, so its rate = the calcium imaging frame rate stored in
%             the ScanImage metadata (_meta.mat / TIFF / _dFF.mat).
%   2. crop -- you drag an ROI box on a preview (OpenCV window). One per video.
%   3. SVD  -- mean-subtracted intensity -> spatial eigenimages x temporal PCs.
%   4. PC1  -- shown + saved as the breathing trace (signed, sign-stabilized).
%
% Reuses: detect_session_fps.m, orofacial_crop_extract.py (the cropper/cuber).
% Per folder it writes:  breath_crop.mat (cube), breath_pc1.mat, breath_pc1.png
%
% RUN THIS IN YOUR OWN MATLAB (not headless) so the crop window can appear.

%% ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = 'D:\260728_vglut2_soma-g8s\phys';
FOLDERS    = {};          % {} = auto-find every breath-video folder under ROOT_DIR
                          % else list specific folders to do a subset, e.g.:
                          % FOLDERS = {'C:\...\cell1\roi3_...00001'};
AVI_PATTERNS = {'Basler_*.avi', 'cam*.avi'};  % breath-cam AVI name patterns:
                          %   'Basler_*.avi' = legacy Pylon-Viewer name
                          %   'cam*.avi'     = new dual-cam GUI (cam1_YYYYMMDD_HHMMSS_runNNN.avi)
                          % one video per recording folder; first match is used.
PYTHON_EXE = 'C:\Program Files\Python314\python.exe';

RE_CROP    = true;        % re-run the interactive crop each time (false = reuse cube)
SHOW_PC    = 1;           % which PC to show as the breathing trace (1 = PC1)
nSVD       = 10;          % SVD components to compute
fBreath    = [0.3 8];     % breath band: PSD shading + which-PC-is-breathing check (Hz)
fallbackFps= 30;          % only if a folder has no calcium metadata
SVD_MAX_BYTES = 1.5e9;    % auto spatial-bin the SVD matrix below this (double)
%% ---------------------------------------------------------------------------

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
chronuxDir = fullfile(repoRoot, 'chronux_2_12');
if isfolder(chronuxDir), addpath(genpath(chronuxDir)); end
extractPy = fullfile(repoRoot, 'orofacial_crop_extract.py');
assert(isfile(extractPy), 'Missing extractor: %s', extractPy);

% Build the folder list
if isempty(FOLDERS)
    avis = [];
    for pi = 1:numel(AVI_PATTERNS)
        d = dir(fullfile(ROOT_DIR, '**', AVI_PATTERNS{pi}));
        if ~isempty(d)
            if isempty(avis), avis = d; else, avis = [avis; d]; end %#ok<AGROW>
        end
    end
    assert(~isempty(avis), 'No AVIs matching {%s} under %s', ...
        strjoin(AVI_PATTERNS, ', '), ROOT_DIR);
    FOLDERS = unique({avis.folder});
end
fprintf('breath_svd_pc1: %d folder(s) to process.\n', numel(FOLDERS));

summary = {};
for fi = 1:numel(FOLDERS)
    folder = FOLDERS{fi};
    fprintf('\n===== [%d/%d] %s =====\n', fi, numel(FOLDERS), folder);

    av = find_avis(folder, AVI_PATTERNS);
    if isempty(av)
        fprintf(2, '  no breath .avi (%s) here -- skip\n', strjoin(AVI_PATTERNS, ', '));
        continue;
    end
    if numel(av) > 1
        fprintf('  %d AVIs here; using "%s"\n', numel(av), av(1).name);
    end
    aviPath = fullfile(folder, av(1).name);
    outMat  = fullfile(folder, 'breath_crop.mat');

    % --- (1) fps from the calcium metadata in the SAME folder ---
    fps = detect_session_fps(folder, fallbackFps);
    fprintf('  imaging fps (2P-triggered breath cam) = %g Hz\n', fps);

    % --- (2) interactive crop -> cube (via the Python extractor) ---
    if RE_CROP && isfile(outMat), delete(outMat); end
    if ~isfile(outMat)
        cmd = sprintf('"%s" "%s" --avi "%s" --ts "%s" --out "%s" --fps %g', ...
            PYTHON_EXE, extractPy, aviPath, ...
            fullfile(folder, 'timestamps.csv'), outMat, fps);
        fprintf('  Crop window will pop up -- drag a box over the snout, press Enter.\n');
        st = system(cmd, '-echo');
        if st ~= 0 || ~isfile(outMat)
            fprintf(2, '  crop/extract failed -- skip\n'); continue;
        end
    end

    % --- (3) load cube + SVD of mean-subtracted intensity ---
    S = load(outMat);
    mov = S.mov; [H, W, T] = size(mov); t = double(S.t_s(:));
    roi = S.roi_xywh(:)';                                % crop box [x y w h]
    X = single(reshape(mov, H*W, T))';
    X = X - mean(X, 1);                                  % mean-sub intensity
    cube = reshape(X', H, W, T);

    sbin = 1;
    while (H*W / sbin^2) * T * 8 > SVD_MAX_BYTES, sbin = sbin + 1; end
    if sbin > 1, cubeB = binCube(cube, sbin); else, cubeB = cube; end
    [Hb, Wb, ~] = size(cubeB);
    Yd = double(reshape(cubeB, Hb*Wb, T));
    Yd = Yd - mean(Yd, 2);
    k  = min(nSVD, min(size(Yd)) - 1);
    [U, Sv, V] = svds(Yd, k);
    sv = diag(Sv); varExp = sv.^2 / sum(sv.^2);

    % --- (4) PC1 (or SHOW_PC) as the breathing trace + diagnostics ---
    pc  = min(SHOW_PC, k);
    sgn = sign(sum(U(:,pc).^3) + eps);                  % stabilize sign
    breathTrace = sgn * sv(pc) * V(:,pc);
    eigImg = sgn * reshape(U(:,pc), Hb, Wb);
    [fa, Pa] = mt_psd(breathTrace, fps);
    inb = fa >= fBreath(1) & fa <= fBreath(2);
    [~, jp] = max(Pa .* inb); pcPeak = fa(jp);          % breath-band peak of this PC
    [pcBreath, fpk] = pick_breath_pc(V, sv, fps, fBreath);  % which PC is most breath-like

    % --- figure: eigenimage + PC trace + PSD ---
    f = figure('Name', ['PC1 breath: ' av(1).name], 'Color', 'w', ...
               'Position', [80 80 1120 720]);
    subplot(2,2,1);
    imagesc(eigImg); axis image off; colormap(gca, parula); colorbar;
    title(sprintf('PC%d eigenimage  (%.0f%% var)', pc, 100*varExp(pc)));
    subplot(2,2,2);
    plot(fa, Pa, 'k', 'LineWidth', 1.3); grid on; hold on;
    xline(pcPeak, 'r--', sprintf('%.2f Hz', pcPeak));
    xlim([0 min(15, fps/2)]); xlabel('frequency (Hz)'); ylabel('power');
    title(sprintf('PC%d PSD', pc));
    subplot(2,1,2);
    plot(t, breathTrace, 'LineWidth', 0.7); grid on; xlim([t(1) t(end)]);
    xlabel('time (s)'); ylabel(sprintf('PC%d', pc));
    title(sprintf('PC%d temporal trace = breathing readout  (%.2f Hz)', pc, pcPeak));
    sgtitle(sprintf('%s    fps=%g Hz', strrep(av(1).name, '_', '\_'), fps));

    if pcBreath ~= pc
        fprintf(2, ['  NOTE: PC%d carries the most breath-band power (peak %.2f Hz), ' ...
            'not PC%d. Breathing may live in PC%d here.\n'], pcBreath, fpk, pc, pcBreath);
    end

    % --- cached diff map for the peak GUI: top - bottom frames by breath signal ---
    % (full-res inspiration - baseline; oriented so inspiration is the high side)
    btf = breathTrace; zz = (btf-mean(btf))/std(btf); if mean(zz.^3) < 0, btf = -btf; end
    nfD = min(100, floor(T/2));
    [~, ordD] = sort(btf, 'descend'); topD = ordD(1:nfD); botD = ordD(end-nfD+1:end);
    camImg  = mean(single(mov), 3);
    diffImg = mean(single(mov(:,:,topD)), 3) - mean(single(mov(:,:,botD)), 3);

    % --- save ---
    save(fullfile(folder, 'breath_pc1.mat'), 'breathTrace', 't', 'fps', ...
         'pc', 'pcBreath', 'pcPeak', 'eigImg', 'diffImg', 'camImg', 'sv', ...
         'varExp', 'U', 'V', 'roi', '-v7.3');
    try
        saveas(f, fullfile(folder, 'breath_pc1.png'));
    catch
    end
    summary(end+1, :) = {av(1).name, fps, pcPeak, pcBreath, 100*varExp(pc)}; %#ok<AGROW>
    fprintf('  PC%d peak %.2f Hz, %.0f%% var. Saved breath_pc1.mat/.png\n', ...
        pc, pcPeak, 100*varExp(pc));
end

%% summary table
fprintf('\n=============== SUMMARY ===============\n');
fprintf('%-45s %5s %8s %7s %7s\n', 'video', 'fps', 'PC1_Hz', 'brPC', 'PC1var%');
for i = 1:size(summary, 1)
    fprintf('%-45s %5g %8.2f %7d %6.0f%%\n', summary{i,1}(1:min(45,end)), ...
        summary{i,2}, summary{i,3}, summary{i,4}, summary{i,5});
end
end

%% =============================== LOCAL FUNCTIONS ===============================
function av = find_avis(dirPath, patterns)
% Breath-cam AVIs in dirPath matching any name pattern (legacy 'Basler_*.avi'
% or new 'cam*.avi'). Deduped by name, sorted so selection is deterministic.
    av = [];
    for pi = 1:numel(patterns)
        d = dir(fullfile(dirPath, patterns{pi}));
        if ~isempty(d)
            if isempty(av), av = d; else, av = [av; d]; end %#ok<AGROW>
        end
    end
    if ~isempty(av)
        [~, ix] = unique({av.name}, 'stable'); av = av(ix);   % drop dup pattern hits
        [~, ord] = sort({av.name}); av = av(ord);
    end
end

function cb = binCube(c, b)
% Block-mean spatial binning of an [H,W,T] cube by integer factor b.
    [H, W, T] = size(c);
    H2 = floor(H/b)*b; W2 = floor(W/b)*b;
    c = c(1:H2, 1:W2, :);
    c = reshape(c, b, H2/b, b, W2/b, T);
    cb = squeeze(mean(mean(c, 1), 3));
end

function [f, P] = mt_psd(x, Fs)
% Temporal PSD: Chronux multitaper if available, else windowed FFT.
    x = x(:) - mean(x(:));
    if exist('mtspectrumc', 'file') == 2
        pr.Fs = Fs; pr.tapers = [3 5]; pr.pad = 0;
        pr.fpass = [0.05, min(15, Fs/2)];
        [P, f] = mtspectrumc(x, pr); f = f(:); P = P(:);
    else
        n = numel(x); nf = floor(n/2) + 1;
        Pf = abs(fft(x .* hann(n))).^2;
        P = Pf(1:nf); f = (0:nf-1)' * (Fs / n);
    end
end

function [pc, fpk] = pick_breath_pc(V, sv, Fs, band)
% Among the computed PCs, return the one with the most power in the breath band.
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
