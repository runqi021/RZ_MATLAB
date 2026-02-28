%% path to your DLC csv
fn = "C:\Users\zhang\Desktop\251224_whisking_KA\bright\whisker_KA_bright-RZ-2025-12-28\videos\Basler_acA1300-200um__22607866__20251224_180941807DLC_Resnet101_whisker_KA_brightDec28shuffle1_snapshot_best-10.csv";

% DLC CSV has 3 header lines: scorer / bodyparts / coords
data = readmatrix(fn, 'NumHeaderLines', 3);

% columns:
% 1: frame index
% 2–4: dot1 (x, y, likelihood)
% 5–7: dot2 (x, y, likelihood)
frame = data(:,1);

x_dot1 = data(:,2);
y_dot1 = data(:,3);
p_dot1 = data(:,4);

x_dot2 = data(:,5);
y_dot2 = data(:,6);
p_dot2 = data(:,7);

%%
s1 = sqrt(x_dot1.^2 + y_dot1.^2);
s2 = sqrt(x_dot2.^2 + y_dot2.^2);
%%
fps = 100;       % <<< change to your actual FPS
t = frame / fps; % time in seconds

%%
figure('Color','w'); hold on;

%plot(t, x_dot1, 'LineWidth', 1, 'Color', 'k');
%plot(t, x_dot2, 'LineWidth', 1, 'Color', 'k');

%plot(t, y_dot1, 'Color', 'b');
%plot(t, y_dot2, 'Color', 'b');

plot(t, s1, 'LineWidth', 1, 'Color', 'r');
plot(t, s2, 'LineWidth', 1, 'Color', 'r');

set(gca, 'YDir', 'reverse');  % so “up” in the image is up in the plot
xlabel('Time (s)');
ylabel('Vertical position (pixels)');
%legend({'dot1','dot2'});
title('Whisker motion (DLC)');
grid on;

%%
figure('Color','w'); hold on; plot(t, p_dot1); plot(t, p_dot2);