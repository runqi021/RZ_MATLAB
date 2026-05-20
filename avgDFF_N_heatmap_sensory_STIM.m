clear; close all; clc;

currentPath = "C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys\IO\sst_dao-z45_2x_256x512_1000f_sound_5k_00001";
sam_hits = dir(fullfile(currentPath, '*cpSAM_output.mat'));
SAM     = fullfile(currentPath, sam_hits(1).name);
load(SAM);

%%
n_pulse = 10;
T_base = 10;
T_pulse = .25;
T_cycle = 10;
BaselineWinSec  = 20;

fs = detect_session_fps(currentPath, 30);
toss = 4; % Sec
F(1:fs*toss, :) = [];
F(end-98:end, :) = [];

[T, N] = size(F);
t = (0:T-1) / fs;

dFF = helper.dFF_RZ(F, 'BaselineWinSec', BaselineWinSec).dFF;

t_pulse = linspace(T_base,T_base+T_cycle-1, n_pulse);

% %% stackDFF and stim overlay
% stackDFF(dFF);
% 
% for i = 1:numel(t_pulse)
%     x1 = t_pulse(i);
%     x2 = x1 + T_pulse;
%     yl = ylim;
%     h = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1 1 0], ...
%     'EdgeColor', 'none', 'FaceAlpha', 0.7);
%     uistack(h,'bottom');
% end
% 
% outFig = fullfile(currentPath, 'stackdFF_Stim.png')
% exportgraphics(gcf, outFig, 'ContentType','vector');
% outFig = fullfile(currentPath, 'stackDFF_Stim.fig')
% savefig(gcf, outFig);
% 
% close all;
% %fprintf('Saved:\n  %s\n', outFig);

%%
figure('Position',[100 100 720 560], 'Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% full trace average
dFF_avg = mean(dFF, 2);
dFF_std = std(dFF,0,2);

ax1 = nexttile(1);
hold(ax1,'on');
plot(ax1, t, dFF_avg, 'Color','k'); hold on;
y1 = dFF_avg - dFF_std; y2 = dFF_avg + dFF_std;
y1 = y1'; y2 = y2';
X = [t, fliplr(t)];
Y = [y1, fliplr(y2)];
fill(ax1, X, Y, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
xlim(ax1, [0 25]);
ylim(ax1, [-0.3, 0.6]);
for i = 1:numel(t_pulse)
    x1 = t_pulse(i);
    x2 = x1 + T_pulse;
    yl = ylim;
    h = patch(ax1, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.4 0.4 0.4], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7);
    uistack(h,'top');
end
ylabel(ax1, "Population average dF/F")

ax2 = nexttile(2, [2 1]);
hold(ax2,'on');

% heatmap
prcLim = [0.1 99.9];
cmap = bluewhitered;
% percentile-based color limits
cLow  = prctile(dFF(:), prcLim(1));
cHigh = prctile(dFF(:), prcLim(2));
cAbs = max(abs(cLow), abs(cHigh));

hHeat = imagesc(ax2, t, 1:size(dFF,2), dFF');   % rows = ROIs
hold on;

set(ax2, 'YDir', 'reverse');

xlim(ax2, [0 25]);
ylim(ax2, [1 N]);
xlabel('Time (s)');
ylabel('ROI');
colormap(cmap);
caxis([-cAbs cAbs]);

cb = colorbar;
ylabel(cb, 'dF/F');

%%
% --- hide non-heatmap stuff ---
ax1.Visible = 'off';
cb.Visible  = 'off';
set(ax2, 'XColor','none','YColor','none');  % hide axes ticks/box
set(ax2, 'Color','none');              % transparent background
% axis(ax2,'tight');  % removed: overrides xlim to full data range

outTif = fullfile(currentPath,'AVGdFF_heatmap.tif');
exportgraphics(ax2, outTif, 'Resolution', 600, 'BackgroundColor','none');

% --- restore ---
ax1.Visible = 'on';
cb.Visible  = 'on';
set(ax2, 'XColor','k','YColor','k');   % or whatever you had

set(hHeat,'Visible','off');
outPDF = fullfile(currentPath,'AVGdFF_overlay.pdf');
exportgraphics(gcf, outPDF, 'ContentType','vector', 'BackgroundColor','none');
set(hHeat,'Visible','on');
fprintf('Saved:\n  %s\n', outPDF);

% outPdf = fullfile(currentPath, 'AVGdFF_all_raster.pdf');
% exportgraphics(gcf, outPdf, 'ContentType','vector','BackgroundColor','none');


% %% triggered average
% pulse_onsets = t_pulse * fs;
% 
% win = 1*fs;
% nEv = numel(pulse_onsets);
% 
% fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);
% 
% dff_seg_all = cell(N,1);
% 
% for i = 1:N
%     roi_trace = dFF(:, i);          % RAW ΔF/F
% 
%     seg = zeros(nEv, win+1);
%     for e = 1:nEv
%         c = pulse_onsets(e);
%         seg(e,:) = roi_trace(c : c+win);
%     end
% 
%     dff_seg_all{i} = seg;
% end
% 
% t = (0:win)/fs;
% ct=0;
% figure('Color', 'white');
% for i = 1:N
%     seg = dff_seg_all{i};
%     avg = mean(seg,1);
%     plot(t, avg, 'color', [0.8, 0.8, 0.8]); hold on;
%     ct=ct+1
%     seg_avg(i,:) = avg;
% end
% 
% pop_avg = mean(seg_avg,1); hold on;
% pop_std = std(seg_avg, 0, 1); hold on;
% plot(t, pop_avg, 'Color', 'k');
% y1 = pop_avg - pop_std; y2 = pop_avg + pop_std;
% 
% X = [t, fliplr(t)];
% Y = [y1, fliplr(y2)];
% 
% fill(X, Y, [0.7 0.7 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% xlabel('Time (s)');
% ylabel('ΔF/F');
% title('Event-Triggered Average with Standard Deviation');
% axis square
% 
% outFig = fullfile(currentPath, 'triggered_avg.png')
% exportgraphics(gcf, outFig, 'ContentType','vector');
% outFig = fullfile(currentPath, 'triggered_avg.fig')
% savefig(gcf, outFig)
% 
% close all;
% %fprintf("%d ROI average overlay \n", ct);


