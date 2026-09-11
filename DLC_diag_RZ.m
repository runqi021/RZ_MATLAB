% DLC_diag_RZ.m
% DLC statistics check: pooled keypoint-LIKELIHOOD histogram for a folder of DLC
% analyzed csvs. Run it once on the thermal (nose) folder and once on the whisk
% folder to eyeball tracking quality.
%
%   - percentage histogram of all likelihoods (pooled over bodyparts + videos)
%   - total frames + time (frames/FPS), %% below each threshold
%   - per-bodypart median likelihood + %% below 0.6 (console + bar)
%
% Right-shifted mass = good tracking; a bump near 0 = frames to refine/relabel.

vidDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";  % DLC csv folder
FPS    = 400;             % for the time label (frames / FPS)
THRESH = [0.6 0.9];      % reference lines + reported fractions
BINS   = 50;

files = dir(fullfile(char(vidDir), '*DLC*.csv'));
assert(~isempty(files), 'no DLC csvs in %s', vidDir);

pooled = [];                 % all likelihoods
bpAll  = containers.Map('KeyType','char','ValueType','any');   % per-bodypart pools
totalFrames = 0;
for f = 1:numel(files)
    csv = fullfile(files(f).folder, files(f).name);
    % bodypart names from header row 2: bodyparts,bp,bp,bp,bp2,...
    fid = fopen(csv); fgetl(fid); hdr = fgetl(fid); fclose(fid);
    parts = strsplit(hdr, ','); bpnames = parts(2:3:end);
    M = readmatrix(csv, 'NumHeaderLines', 3);
    likCols = 4:3:size(M,2);                 % likelihood columns
    totalFrames = totalFrames + size(M,1);
    for b = 1:numel(likCols)
        v = M(:, likCols(b)); v = v(~isnan(v));
        pooled = [pooled; v]; %#ok<AGROW>
        nm = bpnames{b};
        if isKey(bpAll, nm), bpAll(nm) = [bpAll(nm); v]; else, bpAll(nm) = v; end
    end
end
n = numel(pooled); med = median(pooled);

% ---- percentage histogram ----
edges = linspace(0, 1, BINS+1); ctr = (edges(1:end-1)+edges(2:end))/2;
pc = histcounts(pooled, edges, 'Normalization','probability') * 100;
figure('Color','w','Position',[80 80 900 520]);
bar(ctr, pc, 1, 'FaceColor',[0.30 0.50 0.90], 'EdgeColor','w'); hold on; grid on;
xline(med,'k-','LineWidth',1.5); text(med, max(pc)*0.97, sprintf(' median %.2f', med), 'VerticalAlignment','top');
for th = THRESH
    frac = 100*mean(pooled < th);
    xline(th,'r--');
    text(th, max(pc)*0.88, sprintf(' <%.2g: %.1f%%', th, frac), 'Color','r','VerticalAlignment','top');
end
xlim([0 1]); xlabel('likelihood'); ylabel('% of keypoint detections');
[~,stem] = fileparts(char(vidDir));
title(sprintf('%s — %d videos, %d frames, %.1f s   (%d detections, %d bodyparts)', ...
    stem, numel(files), totalFrames, totalFrames/FPS, n, bpAll.Count), 'Interpreter','none');

% ---- console + per-bodypart summary ----
fprintf('\nPOOLED: %d detections | median %.3f | %%<0.6 %.1f | %%<0.9 %.1f\n', ...
    n, med, 100*mean(pooled<0.6), 100*mean(pooled<0.9));
keys = bpAll.keys; medbp = zeros(1,numel(keys)); lowbp = zeros(1,numel(keys));
fprintf('per-bodypart  median  %%<0.6:\n');
for k = 1:numel(keys)
    v = bpAll(keys{k}); medbp(k) = median(v); lowbp(k) = 100*mean(v<0.6);
    fprintf('  %-8s %.2f   %.1f%%\n', keys{k}, medbp(k), lowbp(k));
end
figure('Color','w','Position',[120 120 700 420]);
subplot(1,2,1); bar(categorical(keys), medbp); ylim([0 1]); grid on;
ylabel('median likelihood'); title('per-bodypart median');
subplot(1,2,2); bar(categorical(keys), lowbp); grid on;
ylabel('% < 0.6'); title('per-bodypart % low-confidence');
sgtitle(sprintf('DLC diag — %s', stem), 'Interpreter','none');
