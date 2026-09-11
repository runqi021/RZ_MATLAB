function whisk_LR_diff_then_int_RZ()
% whisk_LR_diff_then_int_RZ  Per-session, ALL sessions.
% Differentiate THEN integrate the (L - R) trace and plot the result.
%
%   d   = xL - xR                  (BP 5-50 Hz, deg)
%   v   = gradient(d, t)           % differentiate  (deg/s)
%   rec = cumtrapz(t, v)           % integrate back (deg) -> recovers d - d(1)
%
% NOTE: differentiate-then-integrate is the identity up to a constant, so
% `rec` is just d with its initial offset removed. The original d is overlaid
% (faint) to show the round-trip. Whisk epochs (ABSOLUTE 5-deg-on-L) shaded.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW      = 400;
BP        = [5 50];     % whisk bandpass (Hz)

ENV_THR    = 5;         % deg, whisk-epoch gate
ENV_SMOOTH = 0.05;      % s
MIN_DUR    = 1;         % s
MERGE_GAP  = 0.2;       % s

ZOOM       = [];        % [t0 t1] s; [] = full trace
EXCLUDE    = "";        % "" = all 17 sessions
SAVE_FIGS  = false;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_diff_then_int');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);
fprintf('%d sessions found\n', numel(sess));

nDone = 0;
for e = 1:numel(sess)
    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW;

    xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest'));
    xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest'));

    env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);

    d   = xL - xR;                 % the trace
    v   = gradient(d, t);          % differentiate  (deg/s)
    rec = cumtrapz(t, v);          % integrate back (deg) ~= d - d(1)

    fprintf('  %s n%d: dur=%.0fs | %d epochs | max|rec-(d-d1)|=%.2g deg\n', ...
        animal, kRun, t(end), size(ep,1), max(abs(rec-(d-d(1)))));

    zoomwin = ZOOM; if isempty(zoomwin), zoomwin = [t(1) t(end)]; end

    fig = figure('Color','w','Position',[60 230 1280 440]); hold on; grid on;
    yl = nzlim([rec; d-d(1)]);
    for q=1:size(ep,1)
        patch(ep(q,[1 2 2 1]), yl([1 1 2 2]), [0.3 0.75 0.3], 'FaceAlpha',0.12,'EdgeColor','none');
    end
    ho = plot(t, d-d(1), '-', 'Color',[0.6 0.6 0.6], 'LineWidth',1.6);   % original (faint)
    hr = plot(t, rec,    '-', 'Color',[0.55 0.15 0.55], 'LineWidth',0.8); % diff-then-int
    yline(0,'k-'); ylim(yl); xlim(zoomwin);
    xlabel('time (s)'); ylabel('(L-R)  (deg)');
    legend([ho hr], {'(L-R) - offset','differentiate then integrate'}, 'Location','northeastoutside');
    title(sprintf('%s n%d  -  differentiate then integrate of (L-R)  (round-trip = identity)', ...
          animal,kRun), 'Interpreter','none');

    nDone = nDone + 1;
    if SAVE_FIGS
        exportgraphics(fig, fullfile(outDir, sprintf('diff_then_int_%s_n%d.png',animal,kRun)), ...
                       'Resolution',150,'BackgroundColor','white');
    end
end
fprintf('Plotted %d sessions.\n', nDone);
end

% ================= helpers =================
function yl = nzlim(x)
    yl = [min(x) max(x)]; if diff(yl)==0, yl = yl + [-1 1]; end
end

function ep = detect_abs(env, t, thrDeg, minDur, mergeGap)
    a = env(:) > thrDeg;
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1;
    ep = [t(s) t(e)];
    if ~isempty(ep)
        m = ep(1,:);
        for i=2:size(ep,1)
            if ep(i,1)-m(end,2) <= mergeGap, m(end,2)=ep(i,2); else, m(end+1,:)=ep(i,:); end %#ok<AGROW>
        end
        ep = m(m(:,2)-m(:,1) >= minDur, :);
    end
end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv')); S = {}; key = {};
    for i=1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o]=sortrows([an(:) rn(:)]); S=S(o);
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
