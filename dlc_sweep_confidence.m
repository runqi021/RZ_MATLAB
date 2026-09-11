function dlc_sweep_confidence(csvPath)
% dlc_sweep_confidence  Plot whisker tip SWEEP with DLC CONFIDENCE overlaid on time,
% for left (vL1) and right (vR1) whisker tips from a DLC csv.
%
%   dlc_sweep_confidence('C:\...\xxxDLC_Resnet50_...csv')

if nargin < 1
    csvPath = "D:\260615_thermalNbasler\260615_thermalNbasler-RZ-2026-06-15\videos\5916296_whisk_n4DLC_Resnet50_260615_thermalNbaslerJun15shuffle1_snapshot_best-20.csv";
end
[folder, stem] = fileparts(csvPath);

M = readmatrix(csvPath, 'NumHeaderLines', 3);     % cols: frame, vL0(xyl), vL1, vR0, vR1
fps = 400;
ts = fullfile(folder, 'timestamps.csv');
if isfile(ts), A = readmatrix(ts); fps = (size(A,1)-1)/((A(end,2)-A(1,2))/1e9); end
t = (0:size(M,1)-1)'/fps;

L1 = [M(:,5) M(:,6)];  pL = M(:,7);     % left tip x,y + likelihood
R1 = [M(:,11) M(:,12)]; pR = M(:,13);   % right tip
swL = projsweep(L1);  swR = projsweep(R1);

f = figure('Color','w','Position',[80 80 1100 600]);
subplot(2,1,1);
yyaxis left;  plot(t, swL, '-'); ylabel('left sweep (px)');
yyaxis right; plot(t, pL, 'k-'); ylabel('confidence'); ylim([0 1]);
title('LEFT tip (vL1)'); xlabel('s'); grid on;

subplot(2,1,2);
yyaxis left;  plot(t, swR, '-'); ylabel('right sweep (px)');
yyaxis right; plot(t, pR, 'k-'); ylabel('confidence'); ylim([0 1]);
title('RIGHT tip (vR1)'); xlabel('s'); grid on;

exportgraphics(f, fullfile(folder, [stem '_sweepconf.png']), 'Resolution', 150);
fprintf('saved %s\n', fullfile(folder, [stem '_sweepconf.png']));
end

function s = projsweep(xy)
C = xy - mean(xy,1); [V,~] = eig(cov(C)); s = C*V(:,end);   % project on sweep axis
end
