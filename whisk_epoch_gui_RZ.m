function whisk_epoch_gui_RZ()
% whisk_epoch_gui_RZ  Select a few whisking epochs per session (to focus on later).
% Auto-detects candidate bouts from the envelope; you keep/delete/adjust/add to
% end with the few you want. No good/bad marking.
%
% Top: whisk envelope (bandpass BP -> Hilbert amplitude, L/R combined) + a
% threshold line (visual guide only). Bottom: L/R whisk (BP) + breath (BR_BP),
% filtered. Add epochs by clicking 2 points; navigate/adjust/delete; jump
% sessions (auto-saves first; reloads a session's saved epochs to resume).
%
% Keys:  a add (2 clicks)   e adjust (2 clicks)   d delete   <-/-> prev/next epoch
%        f full view   s save   [ / ]  prev/next SESSION
%
% RUN IN MATLAB (interactive).

% ============================ USER-EDITABLE ============================
startFolder = "D:\260615_thermalNbasler\5916297\cam1_20260615_195807_run001";  % open this first
dataRoot    = "D:\260615_thermalNbasler";
whiskDir    = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir     = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";
BP        = [5 30];     % whisk bandpass (Hz)
BR_BP     = [1 15];     % breathing bandpass (Hz)
THR_FRAC  = 0.25;       % envelope guide line = THR_FRAC * 95th-pct(envelope)
PAD       = 1;          % s, zoom padding
fpsW      = 400;
% ======================================================================

addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));
dstr = string(regexp(char(dataRoot),'\d+','match','once'));
[bb,aa] = butter(3, BP/(fpsW/2),'bandpass');

% ---- session list: every cam1_* with a matching whisk csv ----
sessions = struct('folder',{},'animal',{},'k',{});
ad = dir(char(dataRoot));
for a = 1:numel(ad)
    if ~ad(a).isdir || ~all(isstrprop(ad(a).name,'digit')), continue; end
    rr = dir(fullfile(char(dataRoot), ad(a).name, 'cam1_*'));
    [~,oo] = sort({rr.name}); rr = rr(oo);
    for kk = 1:numel(rr)
        if ~isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', ad(a).name, kk))))
            sessions(end+1) = struct('folder',fullfile(rr(kk).folder,rr(kk).name), ...
                                     'animal',ad(a).name,'k',kk); %#ok<AGROW>
        end
    end
end
assert(~isempty(sessions),'no sessions with whisk csv under %s',dataRoot);
sIdx = find(strcmp({sessions.folder}, char(startFolder)),1); if isempty(sIdx), sIdx=1; end

% ---- shared state (set by load_session) ----
[t,env,xL,xR,thr,breath_f,tBr] = deal([]);
epochs=zeros(0,2); cur=1; animal=""; k=0; wcsv=""; runFolderCur=""; ph=gobjects(0);

% ---- figure ----
f = figure('Color','w','Position',[60 70 1320 720],'KeyPressFcn',@onkey);
axF = axes('Parent',f,'Position',[0.06 0.69 0.90 0.27]); hold(axF,'on'); grid(axF,'on');
axC = axes('Parent',f,'Position',[0.06 0.43 0.90 0.22]); hold(axC,'on'); grid(axC,'on');
axP = axes('Parent',f,'Position',[0.06 0.17 0.90 0.22]); hold(axP,'on'); grid(axP,'on');
linkaxes([axF axC axP],'x');
B = {'add',@addep; 'adjust',@adjustep; 'delete',@delep; '<- prev',@()nav(-1); 'next ->',@()nav(1); ...
     'full',@fullview; 'SAVE',@saveep; '[ sess',@()jumpsess(-1); 'sess ]',@()jumpsess(1)};
for bi=1:size(B,1)
    uicontrol(f,'Style','pushbutton','String',B{bi,1},'Units','normalized', ...
        'Position',[0.07+0.098*(bi-1) 0.05 0.092 0.06],'Callback',@(~,~)B{bi,2}());
end
load_session(sIdx);

% ================= nested =================
    function load_session(i)
        sIdx = max(1,min(i,numel(sessions)));
        s = sessions(sIdx);
        runFolderCur = s.folder; animal = string(s.animal); k = s.k;
        wcsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, k));
        M = dlc_gate_interp(wcsv, 0.6);    % lik<0.6 -> linear interp
        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
        t  = (0:numel(La)-1)'/fpsW;
        xL = filtfilt(bb,aa, fillmissing(La-mean(La,'omitnan'),'linear'));
        xR = filtfilt(bb,aa, fillmissing(Ra-mean(Ra,'omitnan'),'linear'));
        env = (abs(hilbert(xL)) + abs(hilbert(xR)))/2;
        thr = THR_FRAC * prctile(env,95);
        sav = fullfile(char(runFolderCur), sprintf('%s_%s_n%d_whisk_epochs.mat',dstr,animal,k));
        if isfile(sav), Sv=load(sav); epochs=Sv.epochs;          % resume saved selection
        else, epochs=zeros(0,2); end                             % else start EMPTY (no auto-detect)
        cur=1;
        breath_f=[]; tBr=[];
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d',animal,k)), dataRoot);
            if isfile(Pn.breath)
                Bs=load(Pn.breath); br=Bs.breath(:); fb=double(Bs.fps); tBr=(0:numel(br)-1)'/fb;
                if BR_BP(2) < fb/2
                    [b2,a2] = butter(2, BR_BP/(fb/2), 'bandpass');
                    breath_f = filtfilt(b2, a2, fillmissing(br,'linear'));
                else
                    breath_f = br - mean(br,'omitnan');
                end
            end
        catch
        end
        cla(axF); plot(axF,t,env,'Color',[0.2 0.4 0.8]); yline(axF,thr,'r--');
        ylabel(axF,'whisk envelope'); xlim(axF,[t(1) t(end)]);
        cla(axC); plot(axC,t,xL,'-'); plot(axC,t,xR,'-');
        leg={sprintf('L %g-%g',BP(1),BP(2)), sprintf('R %g-%g',BP(1),BP(2))};
        if ~isempty(breath_f)
            plot(axC,tBr, zsc(breath_f)*std([xL;xR],'omitnan'),'k-');
            leg{end+1}=sprintf('breath %g-%g',BR_BP(1),BR_BP(2));
        end
        ylabel(axC,'whisk(BP) / breath(BP)'); legend(axC,leg,'Location','northeastoutside');
        % Hilbert instantaneous phase of the three filtered traces
        cla(axP);
        hp = plot(axP, t, angle(hilbert(xL)), '.'); hp(2) = plot(axP, t, angle(hilbert(xR)), '.');
        lp = leg(1:2);
        if ~isempty(breath_f)
            hp(3) = plot(axP, tBr, angle(hilbert(breath_f)), 'k.'); lp{3} = leg{end};
        end
        ylim(axP,[-pi pi]); yticks(axP,[-pi 0 pi]); yticklabels(axP,{'-\pi','0','\pi'});
        ylabel(axP,'Hilbert phase'); xlabel(axP,'s'); legend(axP, hp, lp, 'Location','northeastoutside');
        ph=gobjects(0); redraw(); fullview();
    end
    function redraw()
        delete(ph); ph=gobjects(size(epochs,1),1); yl=ylim(axF);
        for ii=1:size(epochs,1)
            lw=0.5; if ii==cur, lw=2.5; end
            ph(ii)=patch(axF, epochs(ii,[1 2 2 1]), yl([1 1 2 2]), [0.3 0.75 0.3], ...
                'FaceAlpha',0.25,'EdgeColor',[0.15 0.5 0.15],'LineWidth',lw, ...
                'ButtonDownFcn',@(~,~)pick(ii));
        end
        title(axF, sprintf('%s n%d  [%d/%d sessions]  |  %d epochs selected', ...
            animal,k,sIdx,numel(sessions),size(epochs,1)));
    end
    function zoomto(i)
        if isempty(epochs), return; end
        cur=max(1,min(i,size(epochs,1))); xlim(axF,[epochs(cur,1)-PAD epochs(cur,2)+PAD]); redraw();
    end
    function nav(d), zoomto(cur+d); end
    function pick(i), cur=i; redraw(); end
    function fullview(), xlim(axF,[t(1) t(end)]); redraw(); end
    function delep()
        if isempty(epochs), return; end
        epochs(cur,:)=[]; cur=max(1,min(cur,size(epochs,1))); redraw();
    end
    function addep()
        axes(axF); [x,~]=ginput(2); if numel(x)<2, return; end
        epochs(end+1,:)=sort([x(1) x(2)]); cur=size(epochs,1); redraw();
    end
    function adjustep()
        if isempty(epochs), return; end
        axes(axF); [x,~]=ginput(2); if numel(x)<2, return; end
        epochs(cur,:)=sort([x(1) x(2)]); redraw();
    end
    function saveep()
        out = fullfile(char(runFolderCur), sprintf('%s_%s_n%d_whisk_epochs.mat',dstr,animal,k));
        S = struct('epochs',epochs,'bp',BP,'br_bp',BR_BP,'fps',fpsW, ...
            'animal',char(animal),'run',k,'src_whisk',wcsv);
        save(out,'-struct','S');
        fprintf('saved %s  (%d epochs)\n', out, size(epochs,1));
    end
    function jumpsess(d), saveep(); load_session(sIdx+d); end   % auto-save then move
    function onkey(~,ev)
        switch ev.Key
            case 'a', addep();   case 'e', adjustep();   case 'd', delep();
            case 'leftarrow',  nav(-1);
            case 'rightarrow', nav(1);
            case 'f', fullview();   case 's', saveep();
            case 'leftbracket',  jumpsess(-1);
            case 'rightbracket', jumpsess(1);
        end
    end
end

% ================= local helpers =================
function z = zsc(x), z=(x-mean(x,'omitnan'))/std(x,'omitnan'); end
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
