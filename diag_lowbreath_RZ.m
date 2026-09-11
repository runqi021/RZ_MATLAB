function diag_lowbreath_RZ()
% diag_lowbreath_RZ  Scan all sessions for epochs where the instantaneous
% breathing rate (1/ITI of inspiration onsets) drops below LOW_HZ, report how
% common they are, and PLOT the clearest (longest) low-rate run found.
% Same breath detection convention as breath_freq_whisk_psd_RZ + sub-frame onsets.

% ============================ USER-EDITABLE ============================
dataRoot  = "D:\260615_thermalNbasler";
noseDir   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
whiskDir  = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
BR_BP     = [2 15];     % breath bandpass (Hz) + onset detection
BR_PROM   = 0.5;        % inspiration trough prominence (x std of BP breath)
LOW_HZ    = 5;          % "slow breathing" threshold (Hz)
MINRUN    = 2;          % min # consecutive cycles below LOW_HZ to call a run
PAD_S     = 1.0;        % s of context to show around the plotted run
% whisk-epoch gating (a breath cycle counts only if its midpoint is in a whisk epoch)
BP        = [5 30];     % whisk bandpass (Hz)
THR_FRAC  = 0.5;        % whisk-epoch threshold = THR_FRAC * 95th-pct(envelope)
MIN_DUR   = 0.5;        % s, min whisk-epoch duration
MERGE_GAP = 0.10;       % s, merge whisk epochs closer than this
fpsW      = 400;
EXCLUDE   = "5840027";
% ======================================================================
[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));

allFreq = [];                 % every cycle frequency (Hz), pooled
R = {};                       % all low-rate runs (browsable)
ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    if any(strcmp(ad(a).name, EXCLUDE)), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*')); [~,o]=sort({rr.name}); rr=rr(o);
    for kk = 1:numel(rr)
        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk)))), continue; end
        % whisk epochs (gate)
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',ad(a).name,kk)), 0.6);
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        tW = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bw,aw, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bw,aw, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;
        ep  = detect(env, tW, THR_FRAC, MIN_DUR, MERGE_GAP);   % whisk epochs [t0 t1]
        if isempty(ep), continue; end
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',ad(a).name,kk)), dataRoot);
            if ~isfile(Pn.breath), continue; end
            Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps);
        catch
            continue
        end
        [b2,a2]=butter(3,BR_BP/(fb/2),'bandpass'); brf=filtfilt(b2,a2,fillmissing(br,'linear'));
        tBr = (0:numel(brf)-1)'/fb;
        [~,il]=findpeaks(-brf,'MinPeakProminence',BR_PROM*std(brf),'MinPeakDistance',round(fb/BR_BP(2)));
        if numel(il) < 3, continue; end
        tInsp = (subsamp(-brf,il)-1)/fb;                  % sub-frame onsets
        f = 1./diff(tInsp);                               % per-cycle rate (Hz), length numel(il)-1
        tMid = (tInsp(1:end-1)+tInsp(2:end))/2;           % cycle midpoint time
        inWh = any(tMid(:) >= ep(:,1)' & tMid(:) <= ep(:,2)', 2);   % cycle midpoint inside a whisk epoch
        allFreq = [allFreq; f(inWh)]; %#ok<AGROW>          % pool ONLY whisk-epoch cycles

        % runs of consecutive whisk-epoch cycles with f < LOW_HZ -> collect each
        low = (f(:) < LOW_HZ) & inWh(:);
        dd = diff([0; low(:); 0]); rs = find(dd==1); re = find(dd==-1)-1;
        for j = 1:numel(rs)
            L = re(j)-rs(j)+1; if L < MINRUN, continue; end
            i0 = rs(j); i1 = re(j)+1;                      % onset indices bounding the run
            t0 = tInsp(i0)-PAD_S; t1 = tInsp(i1)+PAD_S;
            selB = tBr>=t0 & tBr<=t1; selW = tW>=t0 & tW<=t1; onIn = tInsp>=t0 & tInsp<=t1;
            R{end+1} = struct('animal',ad(a).name,'kk',kk,'len',L,'minf',min(f(i0:i1-1)), ...
                't0',t0,'t1',t1,'run0',tInsp(i0),'run1',tInsp(i1), ...
                'tB',tBr(selB),'brf',brf(selB),'tW',tW(selW),'xL',xL(selW),'xR',xR(selW), ...
                'onsets',tInsp(onIn),'fcyc',f(i0:i1-1),'ftimes',(tInsp(i0:i1-1)+tInsp(i0+1:i1))/2); %#ok<AGROW>
        end
    end
end

% ---------------- report ----------------
allFreq = allFreq(isfinite(allFreq));
fprintf('pooled breath cycles: %d\n', numel(allFreq));
fprintf('rate min %.2f Hz | 1st pct %.2f | 5th pct %.2f | median %.2f Hz\n', ...
    min(allFreq), prctile(allFreq,1), prctile(allFreq,5), median(allFreq));
fprintf('fraction < %g Hz: %.2f%%   (%d cycles)\n', LOW_HZ, 100*mean(allFreq<LOW_HZ), sum(allFreq<LOW_HZ));

if isempty(R)
    fprintf('NO run of >=%d consecutive whisk-epoch cycles below %g Hz found.\n', MINRUN, LOW_HZ);
    return
end
[~,ord] = sort(cellfun(@(x) x.len, R), 'descend'); R = R(ord);   % longest runs first
nR = numel(R);
fprintf('found %d low-rate runs (>=%d cycles, in whisk epochs).  arrows / buttons to browse.\n', nR, MINRUN);

% ---------------- interactive browser ----------------
idx = 1;
fig = figure('Color','w','Position',[120 180 1150 460]);
uicontrol(fig,'Style','pushbutton','String','<< Prev','Units','normalized', ...
    'Position',[0.32 0.02 0.12 0.06],'FontSize',10,'Callback',@(s,e)step(-1));
uicontrol(fig,'Style','pushbutton','String','Next >>','Units','normalized', ...
    'Position',[0.56 0.02 0.12 0.06],'FontSize',10,'Callback',@(s,e)step(1));
set(fig,'KeyPressFcn',@onkey);
draw();

    function step(d), idx = mod(idx-1+d, nR)+1; draw(); end
    function onkey(~,ev)
        switch ev.Key
            case {'rightarrow','d'}, step(1);
            case {'leftarrow','a'},  step(-1);
        end
    end
    function draw()
        delete(findall(fig,'Type','axes'));        % keep the uicontrols, redraw axes
        r = R{idx};
        ax = axes('Parent',fig,'Position',[0.08 0.20 0.86 0.66]); hold(ax,'on'); grid(ax,'on');
        % --- breath (left axis) ---
        yyaxis(ax,'left');
        plot(ax, r.tB, r.brf, '-', 'Color',[0.15 0.15 0.15], 'LineWidth',1.1);
        plot(ax, r.onsets, interp1(r.tB,r.brf,r.onsets), 'v', ...
            'MarkerFaceColor',[0 0.4 0.9],'MarkerEdgeColor','none','MarkerSize',7);
        ylabel(ax,'breath (bandpassed)'); set(ax,'YColor',[0.15 0.15 0.15]);
        yl = ylim(ax);
        patch(ax,[r.run0 r.run1 r.run1 r.run0],[yl(1) yl(1) yl(2) yl(2)], ...
            [0.2 0.5 1],'FaceAlpha',0.10,'EdgeColor','none');
        for q = 1:numel(r.fcyc)
            text(ax, r.ftimes(q), yl(2)*0.92, sprintf('%.1f',r.fcyc(q)), ...
                'HorizontalAlignment','center','Color',[0 0.3 0.8],'FontSize',8);
        end
        % --- whisker L & R (right axis) ---
        yyaxis(ax,'right');
        plot(ax, r.tW, r.xL, '-', 'Color',[0 0.62 0.45], 'LineWidth',0.9);   % L = green
        plot(ax, r.tW, r.xR, '-', 'Color',[0.70 0.35 0.10], 'LineWidth',0.9); % R = orange
        ylabel(ax,'whisker angle (deg, BP)'); set(ax,'YColor',[0.3 0.3 0.3]);
        xlim(ax,[r.t0 r.t1]); xlabel(ax,'time (s)');
        title(ax, sprintf('run %d/%d   %s n%d : %d cycles < %g Hz (slowest %.1f Hz)   [v=insp onset, green=whiskL, orange=whiskR]', ...
            idx, nR, r.animal, r.kk, r.len, LOW_HZ, r.minf), 'Interpreter','none');
    end
end

% ================= helpers =================
function p = subsamp(x, idx)
    x = x(:); idx = double(idx(:)); p = idx;
    in = idx>1 & idx<numel(x); i = idx(in);
    ym = x(i-1); y0 = x(i); yp = x(i+1);
    den = ym - 2*y0 + yp;
    delta = 0.5*(ym - yp) ./ den;
    delta(~isfinite(delta) | abs(delta)>0.5) = 0;
    p(in) = i + delta;
end
function ep = detect(env, t, thrFrac, minDur, mergeGap)
    a = env(:) > thrFrac*prctile(env,95);
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
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
