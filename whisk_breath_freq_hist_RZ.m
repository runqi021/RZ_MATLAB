function whisk_breath_freq_hist_RZ()
% whisk_breath_freq_hist_RZ  Histogram of instantaneous breathing frequency
% (1 / inter-inspiration interval), pooled over all sessions, with the
% basal / sniff classification:
%   < BASAL_HZ (3 Hz)  -> basal respiration (black)
%   > SNIFF_HZ (5 Hz)  -> sniffs           (red)
%   in between          -> unclassified     (grey)
% Inspiration onset = trough of the (inhale-up) breath signal.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BASAL_HZ  = 5;          % < this = basal respiration
SNIFF_HZ  = 7;          % > this = sniff
BR_BP     = [2 15];     % breath bandpass before onset detection (Hz)
BR_PROM   = 0.5;        % inspiration trough prominence (x std of bandpassed breath)
EDGES     = 0:0.25:15;  % frequency histogram bin edges (Hz)
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

ad = dir(char(dataRoot)); sess = {};
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sess{end+1} = {ad(a).name, kk}; %#ok<AGROW>
        end
    end
end
assert(~isempty(sess),'no sessions with whisk csv');

freq = [];
for e = 1:numel(sess)
    animal = sess{e}{1}; kk = sess{e}{2};
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,kk)), dataRoot);
        if ~isfile(Pn.breath), continue; end
        Bs=load(Pn.breath); br=Bs.breath(:); fbB=double(Bs.fps);
    catch
        continue
    end
    [b2,a2] = butter(2, BR_BP/(fbB/2), 'bandpass');
    brf = filtfilt(b2, a2, fillmissing(br,'linear'));
    [~,iloc] = findpeaks(-brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fbB/BR_BP(2)));
    tInsp = (iloc-1)/fbB;
    freq = [freq; 1./diff(tInsp)]; %#ok<AGROW>
end
freq = freq(isfinite(freq) & freq<=EDGES(end));
nB = nnz(freq<BASAL_HZ); nS = nnz(freq>SNIFF_HZ); nU = nnz(freq>=BASAL_HZ & freq<=SNIFF_HZ);
fprintf('%d breaths: basal %.0f%% | unclassified %.0f%% | sniff %.0f%%\n', ...
    numel(freq), 100*nB/numel(freq), 100*nU/numel(freq), 100*nS/numel(freq));

figure('Color','w','Position',[100 100 620 360]); hold on; grid on;
histogram(freq(freq>=BASAL_HZ & freq<=SNIFF_HZ), EDGES, 'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
histogram(freq(freq<BASAL_HZ),                   EDGES, 'FaceColor',[0 0 0],      'EdgeColor','none');
histogram(freq(freq>SNIFF_HZ),                   EDGES, 'FaceColor',[0.85 0.1 0.1],'EdgeColor','none');
xline(BASAL_HZ,'k--','LineWidth',1); xline(SNIFF_HZ,'r--','LineWidth',1);
xlim([EDGES(1) EDGES(end)]);
xlabel('instantaneous breathing frequency (Hz)'); ylabel('number of breaths');
legend({sprintf('unclassified (%.0f%%)',100*nU/numel(freq)), ...
        sprintf('basal <%g Hz (%.0f%%)',BASAL_HZ,100*nB/numel(freq)), ...
        sprintf('sniff >%g Hz (%.0f%%)',SNIFF_HZ,100*nS/numel(freq))}, 'Location','northeast');
title(sprintf('breathing frequency  —  %d breaths, %d sessions', numel(freq), numel(sess)));
end

% ================= helpers =================
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
