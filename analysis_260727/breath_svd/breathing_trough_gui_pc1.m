function breathing_trough_gui_pc1()
%BREATHING_TROUGH_GUI_PC1  Inspiration-start (foot of rising flank) detector on
%                          the PC1 breath waveform produced by
%                          breathing_peak_gui_pc1.m.
%
%   Same as breathing_trough_gui.m but the source waveform comes from
%   breath_peak_pc1.mat (SVD PC1) instead of the DLC-based peak .mat. It loads
%   BP.breath, BP.t_breath and BP.insp_onset_idx verbatim; you only edit the
%   insp-start thresholds + the marks themselves.
%
%   Sessions = folders containing a breath_peak_pc1.mat file.
%
%   Output (per session):
%       breath_insp_start_pc1.mat
%           insp_start_idx, insp_starts_t, insp_starts_amp, insp_starts_train,
%           breath, t_breath, findstart_params, peak_file
%       breath_insp_start_pc1_figure.fig / .png
%
%   Requires: MATLAB R2021a+

%% ── Shared state ─────────────────────────────────────────────────────────
S.breath=[]; S.t=[]; S.fps=NaN; S.peakIdx=[]; S.troughs=[];
S.editMode='none'; S.masterFolder=''; S.sessionList={}; S.sessionIdx=0;

% Detector is PARAMETER-FREE (find_feet_deriv). Only a click-snap window remains.

%% ── Figure ──────────────────────────────────────────────────────────────
fig = uifigure('Name','Breathing Inspiration-Start Sorter (SVD PC1)', ...
               'Position',[40 30 1560 930], 'WindowState','maximized');
rootGL = uigridlayout(fig, [1 2]); rootGL.ColumnWidth={340,'1x'};
rootGL.Padding=[6 6 6 6]; rootGL.ColumnSpacing=8;

ctrlPan = uipanel(rootGL,'Title','Controls','FontSize',12,'FontWeight','bold'); ctrlPan.Layout.Column=1;
nCtrlRows=26; cGL=uigridlayout(ctrlPan,[nCtrlRows 2]);
cGL.RowHeight=repmat({'fit'},1,nCtrlRows); cGL.ColumnWidth={'fit','1x'};
cGL.Padding=[8 8 8 8]; cGL.RowSpacing=4; r=0;

r=r+1; h=uilabel(cGL,'Text','── Session Navigation ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnLoad=uibutton(cGL,'Text','Browse Master Folder','ButtonPushedFcn',@cb_loadMaster); btnLoad.Layout.Row=r; btnLoad.Layout.Column=[1 2];
r=r+1; lblMasterPath=uilabel(cGL,'Text','No folder selected','WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]); lblMasterPath.Layout.Row=r; lblMasterPath.Layout.Column=[1 2];
r=r+1; lblCurrentFolder=uilabel(cGL,'Text','—','FontWeight','bold','WordWrap','on','FontSize',10); lblCurrentFolder.Layout.Row=r; lblCurrentFolder.Layout.Column=[1 2];
r=r+1; lblNavStatus=uilabel(cGL,'Text','0 / 0  folders  |  0 complete','FontSize',10,'FontColor',[0.25 0.65 0.25]); lblNavStatus.Layout.Row=r; lblNavStatus.Layout.Column=[1 2];
r=r+1; btnPrev=uibutton(cGL,'Text','← Prev','ButtonPushedFcn',@cb_prevSession); btnPrev.Layout.Row=r; btnPrev.Layout.Column=1;
       btnNext=uibutton(cGL,'Text','Next →','ButtonPushedFcn',@cb_nextSession); btnNext.Layout.Row=r; btnNext.Layout.Column=2;

r=r+1; h=uilabel(cGL,'Text','── Source waveform (locked) ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; taInfo=uitextarea(cGL,'Value',{'(no session loaded)'},'Editable','off','FontSize',9,'FontName','Consolas'); taInfo.Layout.Row=r; taInfo.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','── Inspiration Start Detection ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; h=uilabel(cGL,'Text','Auto: derivative foot — no tuning.','FontSize',9,'FontColor',[0.4 0.4 0.4],'WordWrap','on'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; h=uilabel(cGL,'Text','Click-snap window (s):'); h.Layout.Row=r; h.Layout.Column=1;
       efSnap=uieditfield(cGL,'numeric','Value',0.30,'Limits',[0 10]); efSnap.Layout.Row=r; efSnap.Layout.Column=2;
r=r+1; btnDetect=uibutton(cGL,'Text','Detect Insp Starts','ButtonPushedFcn',@cb_detect,'BackgroundColor',[0.18 0.40 0.85],'FontColor','white'); btnDetect.Layout.Row=r; btnDetect.Layout.Column=[1 2];
r=r+1; btnGap=uibutton(cGL,'Text','Find Next Gap →','ButtonPushedFcn',@cb_find_gap,'BackgroundColor',[0.85 0.55 0.10],'FontColor','white'); btnGap.Layout.Row=r; btnGap.Layout.Column=[1 2];
r=r+1; btnExtra=uibutton(cGL,'Text','Find Next Extra →','ButtonPushedFcn',@cb_find_extra,'BackgroundColor',[0.55 0.30 0.70],'FontColor','white'); btnExtra.Layout.Row=r; btnExtra.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','── Manual Editing ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblMode=uilabel(cGL,'Text','Mode: Navigate','FontWeight','bold','FontSize',11,'FontColor',[0.90 0.80 0.10]); lblMode.Layout.Row=r; lblMode.Layout.Column=[1 2];
r=r+1; btnAdd=uibutton(cGL,'Text','+ Add Mode','ButtonPushedFcn',@cb_setAdd,'BackgroundColor',[0.20 0.42 0.85],'FontColor','white'); btnAdd.Layout.Row=r; btnAdd.Layout.Column=1;
       btnDel=uibutton(cGL,'Text','- Delete Mode','ButtonPushedFcn',@cb_setDelete,'BackgroundColor',[0.82 0.22 0.18],'FontColor','white'); btnDel.Layout.Row=r; btnDel.Layout.Column=2;
r=r+1; btnNav=uibutton(cGL,'Text','Navigate (N)','ButtonPushedFcn',@cb_setNav); btnNav.Layout.Row=r; btnNav.Layout.Column=1;
       btnClear=uibutton(cGL,'Text','Clear All Starts','ButtonPushedFcn',@cb_clearTroughs,'BackgroundColor',[0.40 0.10 0.10],'FontColor','white'); btnClear.Layout.Row=r; btnClear.Layout.Column=2;
r=r+1; lblCount=uilabel(cGL,'Text','Insp Starts: 0  |  Peaks: 0','FontWeight','bold','FontSize',12); lblCount.Layout.Row=r; lblCount.Layout.Column=[1 2];

r=r+1; h=uilabel(cGL,'Text','── Save ──','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; btnSave=uibutton(cGL,'Text','Save breath_insp_start_pc1.mat','ButtonPushedFcn',@cb_save,'BackgroundColor',[0.05 0.45 0.65],'FontColor','white'); btnSave.Layout.Row=r; btnSave.Layout.Column=[1 2];
% (Reset-Params button removed: the detector has no parameters to reset)

if ispref('BreathingTroughGUI_PC1','workingParams')
    try, rmpref('BreathingTroughGUI_PC1','workingParams'); catch, end
end

%% ── Right: plot ─────────────────────────────────────────────────────────
rightGL=uigridlayout(rootGL,[1 1]); rightGL.Layout.Column=2; rightGL.Padding=[0 0 0 0];
axMain=uiaxes(rightGL);
axMain.Title.String='PC1 breath waveform  |  ▲ peak ref  |  ▼ insp start (foot of rising flank)';
axMain.Title.Interpreter='none'; axMain.XLabel.String='Time (s)'; axMain.XLabel.Interpreter='none';
axMain.YLabel.String='PC1 (z-score)'; axMain.YLabel.Interpreter='none';
hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');
hTroughSc=gobjects(0); hPeakSc=gobjects(0);
fig.KeyPressFcn=@cb_key;

%% ═══════════════════════════════ CALLBACKS ═══════════════════════════════
    function cb_loadMaster(~,~)
        master=uigetdir(S.masterFolder,'Select Master Folder'); if isequal(master,0), return; end
        S.masterFolder=master; sessions=scanForSessions(master);
        if isempty(sessions), uialert(fig,'No breath_peak_pc1.mat found in any subfolder.','No Sessions'); return; end
        S.sessionList=sessions; S.sessionIdx=1; [~,mName]=fileparts(master); lblMasterPath.Text=mName;
        loadCurrentSession();
    end
    function cb_prevSession(~,~)
        if S.sessionIdx<=1, return; end; S.sessionIdx=S.sessionIdx-1; loadCurrentSession();
    end
    function cb_nextSession(~,~)
        if S.sessionIdx>=numel(S.sessionList), return; end; S.sessionIdx=S.sessionIdx+1; loadCurrentSession();
    end

    function loadCurrentSession()
        sess=S.sessionList{S.sessionIdx}; S.troughs=[];
        lblCurrentFolder.Text=sess.label; updateNavDisplay();
        try, PK=load(sess.peakMat); catch ME, uialert(fig,ME.message,'Load Error'); return; end
        assert(isfield(PK,'breath') && isfield(PK,'t_breath'),'Peak .mat missing breath/t_breath.');
        S.breath=double(PK.breath(:)); S.t=double(PK.t_breath(:));
        if isfield(PK,'findpeak_params') && isfield(PK.findpeak_params,'camera_fps')
            S.fps=PK.findpeak_params.camera_fps;
        elseif numel(S.t)>1
            S.fps=1/median(diff(S.t));
        else
            S.fps=NaN;
        end
        if isfield(PK,'insp_onset_idx'), S.peakIdx=sort(double(PK.insp_onset_idx(:))); else, S.peakIdx=[]; end

        info={sprintf('src  : %s', sess.label)};
        if isfield(PK,'findpeak_params')
            fp=PK.findpeak_params;
            if isfield(fp,'camera_fps'), info{end+1}=sprintf('fps  : %.3f',fp.camera_fps); end
            if isfield(fp,'source'),     info{end+1}=sprintf('signal: %s',fp.source);     end
            if isfield(fp,'inverted'),   info{end+1}=sprintf('inverted: %d',fp.inverted); end
        end
        info{end+1}=sprintf('peaks: %d',numel(S.peakIdx));
        taInfo.Value=info;

        if isfile(sess.matFile)
            try
                saved=load(sess.matFile,'insp_start_idx','findstart_params');
                if isfield(saved,'findstart_params')
                    fsp=saved.findstart_params;
                    if isfield(fsp,'snap_win_s'),     efSnap.Value=fsp.snap_win_s;         end
                end
                if isfield(saved,'insp_start_idx'), S.troughs=double(saved.insp_start_idx(:)); end
                lblCurrentFolder.Text=[sess.label '  ✓']; updateNavDisplay();
            catch
            end
        end
        refresh_mainPlot();
    end

    function updateNavDisplay()
        n=numel(S.sessionList); nDone=sum(cellfun(@(s) isfile(s.matFile), S.sessionList));
        lblNavStatus.Text=sprintf('%d / %d  folders  |  %d complete',S.sessionIdx,n,nDone);
    end

    function cb_detect(~,~)
        if isempty(S.breath) || isempty(S.peakIdx)
            uialert(fig,'Need breath waveform AND validated peaks.','No Data'); return;
        end
        S.troughs=find_feet_deriv(S.breath,S.fps,S.peakIdx);
        redraw_marks(); updateCount();
    end

    function cb_find_gap(~,~),   jump_to_issue('gap');   end
    function cb_find_extra(~,~), jump_to_issue('extra'); end

    function jump_to_issue(kind)
        % Zoom to the next inter-peak interval that breaks 1-start-per-cycle.
        %   'gap'   -- no insp start between two consecutive peaks (one missing)
        %   'extra' -- more than one (a false positive to delete)
        % Both are already counted by updateCount; this navigates to them.
        if numel(S.peakIdx)<2 || isempty(S.t), lblMode.Text='no peaks loaded'; return; end
        hits=[]; nIn=[];
        for k=1:numel(S.peakIdx)-1
            c=sum(S.troughs>S.peakIdx(k) & S.troughs<S.peakIdx(k+1));
            if (strcmp(kind,'gap') && c==0) || (strcmp(kind,'extra') && c>1)
                hits(end+1)=k; nIn(end+1)=c; %#ok<AGROW>
            end
        end
        if isempty(hits)
            lblMode.Text=sprintf('No %ss -- 1/cycle ✓',kind);
            lblMode.FontColor=[0.10 0.55 0.10]; return;
        end
        cur_center=mean(axMain.XLim);
        mid=arrayfun(@(g) (S.t(S.peakIdx(g))+S.t(S.peakIdx(g+1)))/2, hits);
        nxt=find(mid>cur_center,1,'first'); if isempty(nxt), nxt=1; end
        g=hits(nxt); t_lo=S.t(S.peakIdx(g)); t_hi=S.t(S.peakIdx(g+1));
        pad=max(0.5,(t_hi-t_lo)*1.5);
        axMain.XLim=[max(S.t(1),t_lo-pad), min(S.t(end),t_hi+pad)];
        if strcmp(kind,'gap')
            lblMode.Text=sprintf('Gap %d / %d -- between peaks %d & %d',nxt,numel(hits),g,g+1);
        else
            lblMode.Text=sprintf('Extra %d / %d -- %d starts between peaks %d & %d', ...
                nxt,numel(hits),nIn(nxt),g,g+1);
        end
        lblMode.FontColor=[0.85 0.45 0.10];
    end

    function cb_axClick(~,~)
        if isempty(S.breath), return; end
        cp=axMain.CurrentPoint; clickT=cp(1,1);
        if clickT<S.t(1) || clickT>S.t(end), return; end
        [~,nearIdx]=min(abs(S.t-clickT));
        switch S.editMode
            case 'add'
                % param-free snap: steepest rise in the click window -> walk back to foot
                snapFr=max(2,round(efSnap.Value*S.fps));
                lo=max(1,nearIdx-snapFr); hi=min(length(S.breath),nearIdx+snapFr);
                seg=(lo:hi)'; s=S.breath(seg); [~,mi]=max([diff(s);0]);
                cur=mi; while cur>1 && s(cur-1)<s(cur), cur=cur-1; end
                trIdx=seg(cur);
                if isempty(S.troughs) || min(abs(S.troughs-trIdx))>2, S.troughs=sort([S.troughs; trIdx]); end
            case 'delete'
                if isempty(S.troughs), return; end
                [~,ci]=min(abs(S.troughs-nearIdx)); S.troughs(ci)=[];
        end
        redraw_marks(); updateCount();
    end

    function cb_key(~,event)
        switch lower(event.Key)
            case 'n', cb_setNav();
            case 'a', cb_setAdd();
            case 'd', cb_setDelete();
        end
    end
    function cb_setAdd(~,~)
        S.editMode='add'; lblMode.Text='Mode: + ADD  (click → snap to foot)'; lblMode.FontColor=[0.30 0.62 1.00];
        disableDefaultInteractivity(axMain); axMain.ButtonDownFcn=@cb_axClick;
    end
    function cb_setDelete(~,~)
        S.editMode='delete'; lblMode.Text='Mode: - DELETE  (click → remove nearest)'; lblMode.FontColor=[1.00 0.35 0.25];
        disableDefaultInteractivity(axMain); axMain.ButtonDownFcn=@cb_axClick;
    end
    function cb_setNav(~,~)
        S.editMode='none'; lblMode.Text='Mode: Navigate'; lblMode.FontColor=[0.90 0.80 0.10];
        axMain.ButtonDownFcn=[]; enableDefaultInteractivity(axMain);
    end
    function cb_clearTroughs(~,~)
        S.troughs=[]; redraw_marks(); updateCount();
    end

    function cb_save(~,~)
        if isempty(S.sessionList), uialert(fig,'No session loaded.','Error'); return; end
        sess=S.sessionList{S.sessionIdx};
        insp_start_idx=S.troughs; insp_starts_t=[]; insp_starts_amp=[];
        insp_starts_train=zeros(size(S.breath));
        if ~isempty(S.troughs) && ~isempty(S.t)
            insp_starts_t=S.t(S.troughs); insp_starts_amp=S.breath(S.troughs);
            insp_starts_train(S.troughs)=1;
        end
        breath=S.breath; t_breath=S.t; %#ok<NASGU>
        findstart_params.snap_win_s=efSnap.Value; findstart_params.method='derivative_foot';
        findstart_params.camera_fps=S.fps; findstart_params.peak_file=sess.peakMat;
        findstart_params.saved_at=datestr(now); %#ok<TNOW1,DATST>
        peak_file=sess.peakMat; %#ok<NASGU>
        save(sess.matFile,'insp_start_idx','insp_starts_t','insp_starts_amp', ...
            'insp_starts_train','breath','t_breath','findstart_params','peak_file');
        saveScreenshots(sess.figBase);
        uialert(fig,sprintf('Saved %d insp-starts →\n%s\n+ .fig and .png',numel(S.troughs),sess.matFile),'Saved','Icon','success');
        updateNavDisplay();
    end

    function saveScreenshots(outBase)
        hFig=figure('Visible','off','Color','white','Position',[100 100 1400 460]);
        ax=axes(hFig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        copyobj(axMain.Children,ax); ax.XLim=axMain.XLim; ax.YLim=axMain.YLim;
        ax.XLabel.String=axMain.XLabel.String; ax.YLabel.String=axMain.YLabel.String;
        ax.Title.String=axMain.Title.String;
        ax.Title.Interpreter='none'; ax.XLabel.Interpreter='none'; ax.YLabel.Interpreter='none';
        savefig(hFig,[outBase '.fig']); exportgraphics(hFig,[outBase '.png'],'Resolution',150); close(hFig);
    end

%% ═══════════════════════════════ PLOT HELPERS ═══════════════════════════
    function refresh_mainPlot()
        cla(axMain); hold(axMain,'on');
        if isempty(S.t) || isempty(S.breath), return; end
        h=plot(axMain,S.t,S.breath,'Color',[0 0 0],'LineWidth',0.9); h.HitTest='off';
        h2=yline(axMain,0,'--','Color',[0.5 0.5 0.5],'LineWidth',0.6); h2.HitTest='off';
        hPeakSc=scatter(axMain,[],[],25,'^','filled','MarkerFaceColor',[0.95 0.50 0.05],'MarkerEdgeColor','none','Tag','peak_marker'); hPeakSc.HitTest='off';
        hTroughSc=scatter(axMain,[],[],75,'v','filled','MarkerFaceColor',[0 0.55 1],'MarkerEdgeColor','none','Tag','trough_marker'); hTroughSc.HitTest='off';
        xlim(axMain,[S.t(1) S.t(end)]); redraw_marks(); updateCount();
    end
    function redraw_marks()
        if isgraphics(hPeakSc)
            v=S.peakIdx; v(v<1 | v>length(S.breath))=[];
            if isempty(v), set(hPeakSc,'XData',[],'YData',[]); else, set(hPeakSc,'XData',S.t(v),'YData',S.breath(v)); end
        end
        if isgraphics(hTroughSc)
            v=S.troughs; v(v<1 | v>length(S.breath))=[];
            if isempty(v), set(hTroughSc,'XData',[],'YData',[]); else, set(hTroughSc,'XData',S.t(v),'YData',S.breath(v)); end
        end
        updateCount();
    end
    function updateCount()
        nGaps=0; nExtra=0;
        if numel(S.peakIdx)>=2 && ~isempty(S.troughs)
            for k=1:numel(S.peakIdx)-1
                n_in=sum(S.troughs>S.peakIdx(k) & S.troughs<S.peakIdx(k+1));
                if n_in==0, nGaps=nGaps+1; end
                if n_in>1,  nExtra=nExtra+1; end
            end
        end
        if nGaps+nExtra==0
            lblCount.Text=sprintf('Insp Starts: %d  |  Peaks: %d  |  ✓ 1/cycle',numel(S.troughs),numel(S.peakIdx));
            lblCount.FontColor=[0.10 0.55 0.10];
        else
            lblCount.Text=sprintf('Insp Starts: %d  |  Peaks: %d  |  gaps: %d  extra: %d',numel(S.troughs),numel(S.peakIdx),nGaps,nExtra);
            lblCount.FontColor=[0.85 0.20 0.10];
        end
    end
end  % breathing_trough_gui_pc1

%% ═══════════════════════════════ LOCAL HELPERS ══════════════════════════
function feet = find_feet_deriv(breath, fps, peak_idx)
% Derivative foot detector, anchored to the (already-validated) peaks. For each
% peak it searches ONLY that peak's inter-peak gap (prev peak -> this peak):
% take the steepest rising sample, then walk back to the baseline foot -- keep
% going through any local min that is still above baseline+20%*amp (so doublet
% NOTCHES don't trap it), bounded by the previous peak. No per-trace tuning.
% Assumes inspiration is a sharp UP deflection (peaks are maxima).
b = breath(:); peaks = sort(peak_idx(:)); npk = numel(peaks);
if npk==0, feet = []; return; end
base = prctile(b, 10);                                  % robust low baseline
if npk>1, medIPI = median(diff(peaks)); else, medIPI = round(fps); end
feet = zeros(npk,1);
for k = 1:npk
    p1 = peaks(k);
    if k==1, p0 = max(1, p1-round(1.5*medIPI)); else, p0 = peaks(k-1); end
    if k==1 && min(b(1:p1)) > base + 0.2*max(b(p1)-base, eps)
        feet(k) = NaN; continue;        % recording started mid-inspiration: onset off-screen
    end
    if p1 <= p0+2, feet(k) = p0; continue; end
    seg = (p0:p1)'; s = b(seg); [~, mi] = max([diff(s); 0]);   % steepest rising sample
    delta = 0.2*max(b(p1)-base, eps);                   % notch tolerance above baseline
    cur = mi;
    while cur>1 && (s(cur-1) < s(cur) || s(cur) > base+delta), cur = cur-1; end
    feet(k) = seg(cur);
end
feet = sort(unique(feet(~isnan(feet))));     % drop truncated first-inspiration (NaN)
end

function sessions = scanForSessions(masterFolder)
%SCANFORSESSIONS  Find every folder with a breath_peak_pc1.mat (PC1 peak output).
sessions = {};
hits = dir(fullfile(masterFolder,'**','breath_peak_pc1.mat'));
if isempty(hits), return; end
for i = 1:numel(hits)
    h = hits(i); [~, leaf] = fileparts(h.folder);
    sessions{end+1} = struct( ...                                   %#ok<AGROW>
        'folder',  h.folder, ...
        'label',   leaf, ...
        'peakMat', fullfile(h.folder, h.name), ...
        'matFile', fullfile(h.folder, 'breath_insp_start_pc1.mat'), ...
        'figBase', fullfile(h.folder, 'breath_insp_start_pc1_figure'));
end
[~, ord] = sort(cellfun(@(s) s.folder, sessions, 'UniformOutput',false));
sessions = sessions(ord);
end
