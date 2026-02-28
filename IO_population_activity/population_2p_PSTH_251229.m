clear; clc;

SAM = "C:\Users\zhang\Desktop\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001_ch1\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001_ch1_ch1_minusDark_MC_cpSAM_output.mat";
SAMload = load(SAM);

%%
F_raw = SAMload.F;
F = F_raw;
F(1:30, :)=[];
%%
dFFout = dFF_RZ(F);
dFF = dFFout.dFF;

%%
stackDFF(dFF);
hmapDFF(dFF);
%hmapDFF(dFF, 30, [0 99.5], bluewhitered);
%%
g = 0.93;
[F_oasis_deconv, dFF_oasis_deconv, spikes_oasis, baseline_oasis] = oasis_deconv_and_dff_AR1(F, g);
%% event detection
% assume: dFF is T x N, fs = 30 Hz
[T,N] = size(dFF_oasis_deconv);

t = (1:1:T)/30;
% --- peak detection (per ROI) ---
pks  = cell(N,1);
locs = cell(N,1);

for i = 1:N
    
    [pks{i}, locs{i}] = findpeaks(dFF_oasis_deconv(:, i), ...
        'MinPeakProminence', 0.35, 'MinPeakHeight', 0.3);
end

% --- binary spike/event trains S(t, ROI) ---
S = zeros(T, N);   % 0/1 events

for i = 1:N
    if ~isempty(locs{i})
        S(locs{i}, i) = 1;
    end
end

figure('Name','Peak detection (good ROIs only)','Color','w');
tiledlayout('flow');  % automatically fills panels

for i = 1:N
    
    nexttile;
    plot(t, dFF_oasis_deconv(:, i), 'k'); hold on;

    if ~isempty(locs{i})
        plot(t(locs{i}), dFF_oasis_deconv(locs{i}, i), 'ro', 'MarkerFaceColor','r', 'MarkerSize',3);
    end

    title(sprintf('ROI %d', i));   % keep original ROI index in title
    xlabel('Time (s)');
    ylabel('\DeltaF/F');
end

% Save peak-detection figure
%outFig = fullfile(tdir, sprintf('%s_peakDetection.pdf', tname));
%exportgraphics(gcf, outFig, 'ContentType','vector');
%fprintf('Saved peak-detection figure:\n  %s\n', outFig);

%%
fs=30;

[T, N] = size(dFF);

% ---- find ROIs that have at least one event ----
hasEvent = any(S ~= 0, 1);     % 1 x N logical
idxActive = find(hasEvent);    % indices of ROIs with events

fprintf('Active ROIs: %s\n', mat2str(idxActive));


dFF_active        = dFF(:, idxActive);
dFF_denoise_active= dFF_oasis_deconv(:, idxActive);
S_active          = S(:, idxActive);

% now all are T x N_active, same columns = same ROIs
stackDFF(dFF_active, dFF_denoise_active, fs, 'krb', 0.5, 0.6, S_active, 0.2);
hmapDFF(dFF_active);
%hmapDFF(dFF_active, 30, [0 99.5], bluewhitered);
%% Raster (vertical bar style) for good ROIs only
Ngood    = numel(idxActive);

figRaster = figure('Name','ROI Event Raster (Good ROIs, Bars)','Color','w'); 
hold on;

barHeight = 0.8;   % length of each vertical bar

for k = 1:Ngood
    r = idxActive(k);               % actual ROI index
    tt = t(S(:, r) == 1);          % event times for that ROI
    
    % draw a vertical line for each detected event
    for i = 1:numel(tt)
        plot([tt(i) tt(i)], [k - barHeight/2, k + barHeight/2], ...
            'k-', 'LineWidth', 2);
    end
end

xlabel('Time (s)');
ylabel('ROI#');
title('Event Detection Raster');

ylim([0 Ngood+1]);
set(gca,'YDir','normal');

% y-tick labels = original ROI IDs
%set(gca,'YTick', 1:Ngood, 'YTickLabel', string(idxActive));

grid off;

% ---- SAVE RASTER FIGURE ----
%outRaster = fullfile(tdir, sprintf('%s_raster.pdf', tname));
%exportgraphics(figRaster, outRaster, 'ContentType','vector');

%fprintf('Saved raster plot:\n  %s\n', outRaster);

%%
binSec = 0.1;              
bin    = round(binSec * fs);

[T, N] = size(S);

t = (0:T-1)'/fs;

pop_raw = sum(S, 2);             % instantaneous population spike count

nBins = floor(T / bin);
PSTH = zeros(nBins, 1);

for k = 1:nBins
    idx = ( (k-1)*bin + 1 ) : (k*bin);
    PSTH(k) = sum(pop_raw(idx));
end

tPSTH = ((0:nBins-1)' + 0.5) * binSec;

PSTH_rate = PSTH / binSec;

figure('Color','w');
bar(tPSTH, PSTH_rate, 'histc');
xlabel('Time (s)');
ylabel('Population rate (spikes/sec)');
title(sprintf('Population PSTH firing rate (bin = %.2f s)', binSec));
grid on;

%%
g = gausswin(1);
g = g/sum(g);

PSTH_sm = conv(PSTH_rate, g, 'same');

figure('Color','w');
plot(tPSTH, PSTH_sm, 'LineWidth',1.2, 'Color', 'k');
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title('PSTH');
grid on;
