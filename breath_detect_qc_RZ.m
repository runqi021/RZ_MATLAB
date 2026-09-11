% breath_detect_qc_RZ  (script)
% QC the inspiration-onset (trough) and breath-peak detection used by
% whisk_breath_coord_RZ, for EVERY session (one figure each). Plots the band-
% passed breath trace with onsets (red v) and peaks (blue ^); raw breath faint.
% Detection params match whisk_breath_coord_RZ.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BR_BP     = [2 15];     % breath bandpass (Hz)
BR_PROM   = 0.5;        % inspiration-onset (trough) prominence (x std of BP breath)
TWIN      = [];         % [] = full recording; [t0 t1] (s) to zoom every figure
EXCLUDE   = "5840027";  % animal ids to skip ("0027")
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

ad = dir(char(dataRoot)); nfig = 0;
for ai = 1:numel(ad)
    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    if any(strcmp(ad(ai).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(ai).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(ai).name, kk)))), continue; end
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(ai).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps);
        catch
            continue
        end
        t = (0:numel(br)-1)'/fb;
        [b2,a2] = butter(3, BR_BP/(fb/2), 'bandpass');
        brf = filtfilt(b2,a2, fillmissing(br,'linear'));
        [~,il] = findpeaks(-brf, 'MinPeakProminence', BR_PROM*std(brf), 'MinPeakDistance', round(fb/BR_BP(2)));
        if numel(il) < 2, continue; end
        tInsp = (il-1)/fb;
        ip = zeros(numel(il)-1,1);
        for k = 1:numel(il)-1, [~,rel]=max(brf(il(k):il(k+1))); ip(k)=il(k)+rel-1; end
        tPeak = (ip-1)/fb; cyc = diff(tInsp);
        fprintf('%s n%d: fps=%.2f, %.0fs | %d onsets | breath %.0f ms (%.1f Hz)\n', ...
            ad(ai).name, kk, fb, t(end), numel(il), 1000*median(cyc), 1/median(cyc));

        figure('Color','w','Position',[60 80 1300 420]); hold on; grid on;
        plot(t, normalize(br,'range',[-1 1])*max(abs(brf)), '-', 'Color',[0.82 0.82 0.82]);   % raw (scaled)
        hbf = plot(t, brf, 'k-', 'LineWidth',1);
        ho  = plot(tInsp, brf(il), 'v', 'Color',[0.85 0.1 0.1], 'MarkerFaceColor',[0.85 0.1 0.1], 'MarkerSize',5);
        hp  = plot(tPeak, brf(ip), '^', 'Color',[0 0.4 0.85],  'MarkerFaceColor',[0 0.4 0.85],  'MarkerSize',5);
        xlabel('time (s)'); ylabel('breath (BP)');
        legend([hbf ho hp], {'BP breath','inspiration onset','breath peak'}, 'Location','northeast');
        if ~isempty(TWIN), xlim(TWIN); else, xlim([0 t(end)]); end
        title(sprintf('%s n%d  breath QC   %d onsets, median %.0f ms (%.1f Hz)   [BR %g-%g Hz, prom %.2f]', ...
            ad(ai).name, kk, numel(il), 1000*median(cyc), 1/median(cyc), BR_BP(1), BR_BP(2), BR_PROM), 'Interpreter','none');
        nfig = nfig + 1;
    end
end
fprintf('QC figures: %d\n', nfig);

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
