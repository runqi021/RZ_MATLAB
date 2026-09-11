function whisk_breath_trig_avg_RZ()
% whisk_breath_trig_avg_RZ  Inspiration-triggered average of the breath waveform,
% split into BASAL (<BASAL_HZ, black) and SNIFFING (>SNIFF_HZ, red) breaths.
% Per-breath rate = 1/ITI; trigger = inspiration onset (trough of BP breath).
% Pooled over all sessions, resampled to a common grid; shows +/-1 s (2 s).

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BR_BP     = [2 15];     % breath bandpass (Hz)
BASAL_HZ  = 5; SNIFF_HZ = 7;
BR_PROM   = 0.5;        % inspiration trough prominence (x std of BP breath)
TWIN      = 1.0;        % s, +/- window (=> 2 s shown)
FS        = 200;        % Hz, common resample grid
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
lags = -TWIN:1/FS:TWIN;  L = numel(lags);
sumv = zeros(2,L); sumsq = zeros(2,L); n = zeros(2,1);   % 1=basal 2=sniff

ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk)))), continue; end
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(a).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
        catch
            continue
        end
        [b2,a2]=butter(2,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        [~,iloc]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        tInsp=(iloc-1)/fb; freq=1./diff(tInsp);
        for i=1:numel(freq)
            if     freq(i)<BASAL_HZ, c=1;
            elseif freq(i)>SNIFF_HZ, c=2;
            else,  continue; end
            w = interp1(tBr, brf, tInsp(i)+lags, 'linear', NaN);
            if any(isnan(w)), continue; end
            sumv(c,:)=sumv(c,:)+w; sumsq(c,:)=sumsq(c,:)+w.^2; n(c)=n(c)+1;
        end
    end
end
assert(any(n>0),'no breaths');
mu = sumv./n; sd = sqrt(max(sumsq./n - mu.^2,0));
sem = sd ./ sqrt(n);                          % SEM (std/sqrt(n)) -- tighter than STD
fprintf('basal n=%d | sniff n=%d breaths\n', n(1), n(2));

figure('Color','w','Position',[100 100 620 440]); hold on; grid on;
band(lags, mu(1,:), sem(1,:), [0 0 0]);        h1=plot(lags, mu(1,:), 'k-', 'LineWidth',1.8);
band(lags, mu(2,:), sem(2,:), [0.85 0.1 0.1]); h2=plot(lags, mu(2,:), '-', 'Color',[0.85 0.1 0.1],'LineWidth',1.8);
xline(0,'b:','LineWidth',1); xlim([-TWIN TWIN]);
xlabel('time from inspiration onset (s)'); ylabel(sprintf('breath (BP %g-%g Hz)',BR_BP(1),BR_BP(2)));
legend([h1 h2], {sprintf('basal <%g Hz (n=%d)',BASAL_HZ,n(1)), sprintf('sniff >%g Hz (n=%d)',SNIFF_HZ,n(2))}, ...
    'Location','northeast');
title('inspiration-triggered breath: basal vs sniffing (mean \pm SEM)');
end

function band(x, m, s, col)
    patch([x fliplr(x)], [m+s fliplr(m-s)], col, 'FaceAlpha',0.15, 'EdgeColor','none');
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
