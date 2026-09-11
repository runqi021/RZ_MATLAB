function whisk_breath_amp_freq_RZ()
% whisk_breath_amp_freq_RZ  Per-breath AMPLITUDE vs FREQUENCY scatter to check
% that detected sniffs are faster AND lower-amplitude than basal breaths.
%   x = breath frequency (1/ITI, Hz)
%   y = breath amplitude (peak-to-trough of the raw breath within the cycle)
%   colour = class: basal <BASAL_HZ (black), sniff >SNIFF_HZ (red), mid (grey)
% One dot per breath, pooled over all sessions. Big markers = class medians.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BR_BP     = [2 15];     % breath bandpass for onset detection (Hz)
BASAL_HZ  = 5; SNIFF_HZ = 7;
BR_PROM   = 0.5;        % inspiration trough prominence (x std of BP breath)
FMAX      = 15;         % x-axis max (Hz)
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
F = []; A = [];
ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk)))), continue; end
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(a).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
        catch
            continue
        end
        [b2,a2]=butter(2,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        [~,iloc]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        for i=1:numel(iloc)-1
            seg = iloc(i):iloc(i+1);
            F(end+1,1) = fb/(iloc(i+1)-iloc(i));            % 1/ITI (Hz) %#ok<AGROW>
            A(end+1,1) = max(br(seg)) - min(br(seg));       % raw peak-to-trough (deg C) %#ok<AGROW>
        end
    end
end
assert(~isempty(F),'no breaths');
F=F(:); A=A(:); keep=F<=FMAX; F=F(keep); A=A(keep);
mB = F<BASAL_HZ; mS = F>SNIFF_HZ; mU = ~mB & ~mS;
fprintf('basal n=%d (amp %.3f) | mid n=%d | sniff n=%d (amp %.3f)\n', ...
    nnz(mB), median(A(mB)), nnz(mU), nnz(mS), median(A(mS)));

figure('Color','w','Position',[100 100 640 520]); hold on; grid on;
plot(F(mU), A(mU), '.', 'Color',[0.6 0.6 0.6], 'MarkerSize',3);
plot(F(mB), A(mB), '.', 'Color',[0 0 0],       'MarkerSize',3);
plot(F(mS), A(mS), '.', 'Color',[0.85 0.1 0.1],'MarkerSize',3);
plot(median(F(mB)), median(A(mB)), 'ko', 'MarkerFaceColor','k','MarkerSize',10);
plot(median(F(mS)), median(A(mS)), 'o', 'Color',[0.85 0.1 0.1],'MarkerFaceColor',[0.85 0.1 0.1],'MarkerSize',10);
xline(BASAL_HZ,'k--'); xline(SNIFF_HZ,'r--');
set(gca,'YScale','log'); xlim([0 FMAX]);
xlabel('breath frequency  1/ITI (Hz)'); ylabel('breath amplitude  peak-to-trough (\circC)');
legend({'unclassified','basal <5 Hz','sniff >7 Hz','basal median','sniff median'}, 'Location','northeast');
title(sprintf('per-breath amplitude vs frequency  (%d breaths)', numel(F)));
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function nn = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), nn=0; else, nn=str2double(tok{1}{1}); end
end
