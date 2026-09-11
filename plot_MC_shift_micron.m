%% plot_MC_shift_micron.m
% Plot the rigid motion-correction |shift| (microns) as a time series.
% Reads a NoRMCorre *_MC_output.mat (mcOut.shifts) and the matching
% *_meta.mat for pixelSize_um and fps.

clear; clc;

%% ---------------- user params ----------------
mcFile   = "c:\Users\Admin\Desktop\260521_ChAT_soma_G8s\roi7_8x_x-1150y100z-35_3000f_25lp_00001\roi7_8x_x-1150y100z-35_3000f_25lp_00001_ch1_preproc_MC_output.mat";
metaFile = "c:\Users\Admin\Desktop\260521_ChAT_soma_G8s\roi7_8x_x-1150y100z-35_3000f_25lp_00001\roi7_8x_x-1150y100z-35_3000f_25lp_00001_ch1_meta.mat";
% ----------------------------------------------

%% Load
S  = load(mcFile);   mcOut = S.mcOut;
Me = load(metaFile);
px  = Me.pixelSize_um;   % microns per pixel
fps = Me.fps;            % Hz

%% Extract rigid shifts -> [y x] pixels per frame
sh  = mcOut.shifts;                 % struct array, one per frame
nF  = numel(sh);
yx  = zeros(nF, 2);
for k = 1:nF
    s = squeeze(sh(k).shifts);      % [y x] in pixels
    yx(k,:) = s(:)';
end

%% Magnitude in microns
shift_um = sqrt(sum(yx.^2, 2)) * px;     % |shift| (microns)
t = (0:nF-1)' / fps;                     % seconds

%% Plot
figure('Color','w','Position',[100 100 1100 380]);
plot(t, shift_um, 'k', 'LineWidth', 0.8);
xlabel('Time (s)');
ylabel('|shift| (\mum)');
title(sprintf('Rigid MC |shift|  (px size %.4f \\mum, %.0f fps)', px, fps));
box off; xlim([t(1) t(end)]); ylim([0 max(shift_um)*1.05]);

fprintf('|shift| (um): mean %.3f, median %.3f, max %.3f over %d frames\n', ...
    mean(shift_um), median(shift_um), max(shift_um), nF);
