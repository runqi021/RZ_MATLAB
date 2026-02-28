lambda = diag(SVD_result.S).^2;

cut     = 120;     % <-- your cut
tailFrac = 0.80;   % <-- your setting
z       = 3;

[kStar, fitInfo] = pick_k_from_tail_line_cut(lambda, cut, tailFrac, z);

%% -------- Plot 1: FULL (1:nAll), fit computed ONLY on 1:cut ----------
figure('Color','w'); hold on;
semilogy((1:1200), lambda(1:1200), 'k-', 'LineWidth', 1.2);

%xline(cut, 'k--', 'cut', 'LabelVerticalAlignment','bottom');

xlabel('Mode number');
ylabel('Eigenvalue, \lambda^2');
grid on; box on;
axis square;
hold off;

%% -------- Plot 2: CUT ZOOM (1:cut) ----------
idxCut = (1:cut)';

figure('Color','w'); hold on;
semilogy(idxCut, lambda(1:cut), 'k-', 'LineWidth', 1.2);

baselineCut = 10.^fitInfo.yhat;
thrCurveCut = 10.^fitInfo.thr;

semilogy(idxCut, baselineCut, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
semilogy(idxCut, thrCurveCut, ':',  'Color',[0.4 0.4 0.4], 'LineWidth',1.5);

plot(kStar, lambda(kStar), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7);

xlabel('Mode number');
ylabel('Eigenvalue (\lambda_k = s_k^2)');
title(sprintf('ZOOM 1:%d (tailFrac=%.2f, z=%g), k^*=%d', cut, tailFrac, z, kStar));
grid on; box on;
legend({'Eigenvalues','Tail linear fit','Fit + z\sigma','k^*'}, ...
       'Location','northeast');
hold off;

%%
function [kStar, fitInfo] = pick_k_from_tail_line_cut(lambda, cut, tailFrac, z)
lambda = lambda(:);
nAll = numel(lambda);

cut = min(cut, nAll);
lambda_use = lambda(1:cut);
n = numel(lambda_use);
idx = (1:n)';

lambda_use(lambda_use<=0) = eps;

nTail = max(20, round(tailFrac*n));
tailIdx = (n-nTail+1):n;

y = log10(lambda_use);
p = polyfit(idx(tailIdx), y(tailIdx), 1);
yhat = polyval(p, idx);

resTail = y(tailIdx) - yhat(tailIdx);
sigma = 1.4826 * median(abs(resTail - median(resTail)));
if ~isfinite(sigma) || sigma == 0
    sigma = std(resTail) + eps;
end

thr = yhat + z*sigma;

above = find(y > thr);
if isempty(above), kStar = 1; else, kStar = above(end); end

fitInfo.p = p;
fitInfo.sigma = sigma;
fitInfo.tailIdx = tailIdx;
fitInfo.cut = cut;
fitInfo.tailFrac = tailFrac;
fitInfo.z = z;
fitInfo.lambda_use = lambda_use;
fitInfo.y = y;
fitInfo.yhat = yhat;
fitInfo.thr = thr;
end
