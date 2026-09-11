function breath_failure_gui_260825()
%BREATH_FAILURE_GUI_260825  Classify breath cycles by amplitude: GASP / normal / FAIL.
%
%   Sessions = folders holding both breath_peak_pc1.mat and
%   breath_insp_start_pc1.mat. Point it at a session root (or at the archive's
%   <Genotype>\<MMDD>) and it finds every recording underneath.
%
%   WHAT A "CYCLE" IS HERE
%       foot_i  ->  peak_i  ->  foot_{i+1}
%   one inspiration onset to the next, with the inspiratory peak between them. A
%   foot with no peak before the next foot is not a cycle and is not listed.
%
%   TWO SIGNALS TO THRESHOLD ON  (the "signal" dropdown)
%
%   1. fb / PC1  (breath_pc1.mat, robust-z)     threshold in Z
%      amplitude = trace(peak) - trace(foot)
%      The trace the landmarks were detected on. Its units are a projection onto
%      a basis, so they are arbitrary; robust-z puts the threshold in units of
%      that recording's own breath variability. NOT comparable in absolute terms
%      across recordings, though a z threshold roughly is.
%
%   2. rigid displacement (breath_fixedmetric.mat, `disp`)   threshold in PIXELS
%      amplitude = max - min WITHIN the cycle
%      Physical pixels on a common anatomical grid: the extractor applies each
%      run's own registration shift before cropping, so this IS comparable
%      across runs of a session. Measured on 260824: trace IQR varies with
%      CV 33% across runs for disp vs 66% for fb, i.e. about twice as
%      consistent. Still not perfect -- pixels are not millimetres, so across
%      SESSIONS it needs calibration (camera distance and zoom), and the ROI
%      follows anatomy per run but not posture changes.
%
%   WHY THE TWO AMPLITUDES ARE DEFINED DIFFERENTLY -- THIS MATTERS
%   `disp` reaches its extremes at a DIFFERENT PHASE than fb does. Sampling it
%   at fb's foot and peak measures an arbitrary phase difference, not an
%   amplitude: measured over 17 recordings that gives Spearman -0.29 against the
%   fb amplitude, and says gasps are 0.79x SMALLER than normal breaths. Taking
%   peak-to-peak within the cycle instead gives +0.90 and 2.17x LARGER, which
%   agrees with the frame-difference maps. So rigid displacement is peak-to-peak
%   by necessity, not by preference.
%
%   fb keeps peak-minus-foot, unchanged, so labels made before this dropdown
%   existed stay valid.
%
%   TWO THRESHOLDS, THREE CLASSES
%       amp >= high   ->  GASP    (+1)
%       amp <= low    ->  FAIL    (-1)
%       otherwise     ->  normal   (0)
%   Each signal keeps its OWN pair of thresholds, since z and px are not the
%   same scale. Both axes show the WHOLE recording on the same x range.
%
%   Clicking a cycle cycles it normal -> FAIL -> GASP -> normal. Overrides are
%   stored SEPARATELY from the thresholds, so retyping a threshold does not wipe
%   them.
%
%   OUTPUT (per recording)   breath_cycle_class_pc1.mat
%       cyc_foot_idx, cyc_peak_idx, cyc_next_foot_idx
%                       RAW frame indices, exactly as they appear in the pc1
%                       files -- NOT nDrop-corrected. A consumer that tosses
%                       nDrop frames must subtract nDrop itself, the same way it
%                       already does for insp_onset_idx. nDrop is in the params.
%       cyc_amp         amplitude in the CHOSEN signal (what was thresholded)
%       cyc_amp_fb      fb amplitude, z          } both are always written, so a
%       cyc_amp_disp    rigid displacement, px   } consumer can use either
%       cyc_ibi_s, cyc_t
%       class_auto      int8, the two thresholds alone (-1 fail, 0 normal, +1 gasp)
%       override_idx, override_class     the cycles you set by hand
%       class_final     int8  <-- THE ONE DOWNSTREAM USES
%                       = class_auto with override_idx replaced by override_class
%       is_fail, is_gasp                 logical copies of class_final
%       class_params    signal_source, amp_definition, thrLow, thrHigh, units,
%                       nDrop, fps, source files, saved_at
%
%   Keeping auto and the overrides separately is the peak editor's
%   insp_onsets_final convention: the hand-corrected set wins, and the pieces
%   are kept so the edit stays auditable and reversible instead of baked in.
%
%   Runqi Zhang / 2026-08-25

FAIL = int8(-1); NORM = int8(0); GASP = int8(1);

%% ---- shared state -----------------------------------------------------
S.sessions = {};
S.idx      = 0;
S.t        = [];
S.fps      = NaN;
S.cyc      = struct('foot',[],'peak',[],'nextFoot',[],'ibi',[],'t',[]);
S.trFB     = [];  S.ampFB = [];      % robust-z fb trace + peak-minus-foot amp
S.trDS     = [];  S.ampDS = [];      % rigid displacement (px) + peak-to-peak amp
S.hasDS    = false;
S.src      = 'fb';
S.thr      = struct('fb',[1.0 4.0], 'disp',[0.3 1.0]);   % [low high] per signal
S.ovr      = int8([]);
S.dirty    = false;
S.nDrop    = 30;

COL = struct('fail',[0.85 0.15 0.15], 'gasp',[0.10 0.35 0.85], 'norm',[0.45 0.45 0.45]);

%% ---- figure -----------------------------------------------------------
fig = uifigure('Name','Breath cycle classifier','Position',[40 40 1600 880]);
gl = uigridlayout(fig,[1 2]);
gl.ColumnWidth = {300,'1x'}; gl.Padding = [6 6 6 6]; gl.ColumnSpacing = 8;

cp  = uipanel(gl,'Title','Controls','FontWeight','bold'); cp.Layout.Column = 1;
cgl = uigridlayout(cp,[20 2]);
cgl.RowHeight = repmat({'fit'},1,20); cgl.ColumnWidth = {'fit','1x'};
cgl.Padding = [8 8 8 8]; cgl.RowSpacing = 4;
r = 0;

r=r+1; h=uilabel(cgl,'Text','-- Sessions --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; b=uibutton(cgl,'Text','Browse folder...','ButtonPushedFcn',@cb_browse); b.Layout.Row=r; b.Layout.Column=[1 2];
r=r+1; lblRoot=uilabel(cgl,'Text','no folder','WordWrap','on','FontSize',9,'FontColor',[0.5 0.5 0.5]); lblRoot.Layout.Row=r; lblRoot.Layout.Column=[1 2];
r=r+1; ddSess=uidropdown(cgl,'Items',{'(none)'},'ValueChangedFcn',@cb_pickSession); ddSess.Layout.Row=r; ddSess.Layout.Column=[1 2];
r=r+1; b=uibutton(cgl,'Text','<< Prev','ButtonPushedFcn',@(s,e)cb_step(-1)); b.Layout.Row=r; b.Layout.Column=1;
       b=uibutton(cgl,'Text','Next >>','ButtonPushedFcn',@(s,e)cb_step(1));  b.Layout.Row=r; b.Layout.Column=2;

r=r+1; h=uilabel(cgl,'Text','-- Signal --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; ddSrc=uidropdown(cgl,'Items',{'fb / PC1  (z)','rigid displacement  (px)'}, ...
           'ItemsData',{'fb','disp'},'Value','fb','ValueChangedFcn',@cb_src, ...
           'Tooltip',['fb: amplitude = peak - foot on the robust-z trace. ' ...
                      'rigid displacement: amplitude = max-min within the cycle, in pixels.']);
       ddSrc.Layout.Row=r; ddSrc.Layout.Column=[1 2];
r=r+1; lblSrc=uilabel(cgl,'Text','','WordWrap','on','FontSize',9,'FontColor',[0.35 0.35 0.35]); lblSrc.Layout.Row=r; lblSrc.Layout.Column=[1 2];

r=r+1; h=uilabel(cgl,'Text','-- Thresholds --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblHi=uilabel(cgl,'Text','GASP  amp >='); lblHi.Layout.Row=r; lblHi.Layout.Column=1;
       efHigh=uieditfield(cgl,'numeric','Value',S.thr.fb(2),'ValueChangedFcn',@cb_thr);
       efHigh.Layout.Row=r; efHigh.Layout.Column=2;
r=r+1; lblLo=uilabel(cgl,'Text','FAIL  amp <='); lblLo.Layout.Row=r; lblLo.Layout.Column=1;
       efLow=uieditfield(cgl,'numeric','Value',S.thr.fb(1),'ValueChangedFcn',@cb_thr);
       efLow.Layout.Row=r; efLow.Layout.Column=2;
r=r+1; lblThr=uilabel(cgl,'Text','','WordWrap','on','FontSize',9); lblThr.Layout.Row=r; lblThr.Layout.Column=[1 2];

r=r+1; h=uilabel(cgl,'Text','-- Manual overrides --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; lblOv=uilabel(cgl,'Text','none','WordWrap','on','FontSize',9); lblOv.Layout.Row=r; lblOv.Layout.Column=[1 2];
r=r+1; b=uibutton(cgl,'Text','Clear overrides','ButtonPushedFcn',@cb_clearOv); b.Layout.Row=r; b.Layout.Column=[1 2];

r=r+1; h=uilabel(cgl,'Text','-- Save --','FontWeight','bold'); h.Layout.Row=r; h.Layout.Column=[1 2];
r=r+1; b=uibutton(cgl,'Text','SAVE','FontWeight','bold','ButtonPushedFcn',@cb_save); b.Layout.Row=r; b.Layout.Column=[1 2];
r=r+1; lblStat=uilabel(cgl,'Text','','WordWrap','on','FontSize',10); lblStat.Layout.Row=r; lblStat.Layout.Column=[1 2];
r=r+1; hHelp=uilabel(cgl,'WordWrap','on','FontSize',9,'FontColor',[0.35 0.35 0.35], ...
    'Text','red = FAIL, blue = GASP, grey = normal. Click a cycle to step it normal -> FAIL -> GASP -> normal.');
       hHelp.Layout.Row=r; hHelp.Layout.Column=[1 2];

pp  = uipanel(gl); pp.Layout.Column = 2;
pgl = uigridlayout(pp,[2 1]); pgl.RowHeight = {'2x','1x'}; pgl.Padding = [4 4 4 4];
axT = uiaxes(pgl); axT.Layout.Row = 1;
axA = uiaxes(pgl); axA.Layout.Row = 2;
title(axT,'load a session'); xlabel(axA,'time (s)');
axT.ButtonDownFcn = @(s,e)cb_click(axT);
axA.ButtonDownFcn = @(s,e)cb_click(axA);

%% =======================================================================
    function cb_browse(~,~)
        d = uigetdir(pwd,'Pick a folder containing recordings');
        if ~ischar(d) || ~isfolder(d), return; end
        lblRoot.Text = d;
        f = dir(fullfile(d,'**','breath_peak_pc1.mat'));
        keep = {};
        for k = 1:numel(f)
            if isfile(fullfile(f(k).folder,'breath_insp_start_pc1.mat'))
                keep{end+1} = f(k).folder; %#ok<AGROW>
            end
        end
        if isempty(keep)
            uialert(fig,'No folder under there has both breath_peak_pc1.mat and breath_insp_start_pc1.mat.','Nothing found');
            return;
        end
        S.sessions = keep;
        ddSess.Items = sessionNames();  ddSess.ItemsData = 1:numel(keep);
        S.idx = 1; ddSess.Value = 1;
        loadSession();
    end

    function nm = sessionNames()
        nm = cell(size(S.sessions));
        for k = 1:numel(S.sessions)
            [~,a,b] = fileparts(regexprep(S.sessions{k},'[\\/]+$',''));
            done = isfile(fullfile(S.sessions{k},'breath_cycle_class_pc1.mat'));
            nm{k} = sprintf('%d/%d %s%s', k, numel(S.sessions), [a b], ternary(done,'  [saved]',''));
        end
    end

    function cb_pickSession(~,~)
        if askDiscard(), ddSess.Value = S.idx; return; end
        S.idx = ddSess.Value; loadSession();
    end
    function cb_step(d)
        if isempty(S.sessions), return; end
        if askDiscard(), return; end
        S.idx = min(max(S.idx + d,1), numel(S.sessions));
        ddSess.Value = S.idx; loadSession();
    end
    function tf = askDiscard()
        tf = false;
        if ~S.dirty, return; end
        c = uiconfirm(fig,'Unsaved changes on this recording. Leave anyway?', ...
                      'Unsaved','Options',{'Leave','Stay'},'DefaultOption',2);
        tf = strcmp(c,'Stay');
    end

    function loadSession()
        fp = S.sessions{S.idx};
        BP = load(fullfile(fp,'breath_peak_pc1.mat'));
        IP = load(fullfile(fp,'breath_insp_start_pc1.mat'));

        S.fps = 30;
        if isfield(BP,'findpeak_params') && isfield(BP.findpeak_params,'camera_fps')
            S.fps = BP.findpeak_params.camera_fps;
        end
        try
            addpath(fileparts(fileparts(mfilename('fullpath'))));
            S.fps = detect_session_fps(fp, S.fps);
        catch
        end

        bw = detrend(double(BP.breath(:)));
        S.trFB = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
        S.t = (0:numel(S.trFB)-1)'/S.fps;

        pk = sort(round(BP.insp_onset_idx(:)));
        ft = sort(round(IP.insp_start_idx(:)));
        pk = pk(pk>=1 & pk<=numel(S.trFB));
        ft = ft(ft>=1 & ft<=numel(S.trFB));
        foot=[]; peak=[]; nxt=[];
        for i = 1:numel(ft)-1
            q = pk(pk>ft(i) & pk<ft(i+1));
            if isempty(q), continue; end
            foot(end+1,1)=ft(i); peak(end+1,1)=q(1); nxt(end+1,1)=ft(i+1); %#ok<AGROW>
        end
        S.cyc.foot=foot; S.cyc.peak=peak; S.cyc.nextFoot=nxt;
        S.cyc.ibi = (nxt - foot)/S.fps;
        S.cyc.t   = S.t(foot);
        S.ampFB   = S.trFB(peak) - S.trFB(foot);        % peak MINUS foot, in z
        S.ovr     = repmat(int8(-2), numel(foot), 1);

        % ---- rigid displacement, if this session has a fixedmetric file -----
        [S.trDS, S.hasDS, whyNo] = load_disp(fp, numel(S.trFB));
        S.ampDS = nan(numel(foot),1);
        if S.hasDS
            % peak-to-peak WITHIN the cycle. disp peaks at a different phase
            % than fb, so sampling it at fb's landmarks measures phase, not
            % amplitude (Spearman -0.29 vs fb; peak-to-peak gives +0.90).
            for i = 1:numel(foot)
                s = foot(i):nxt(i);
                S.ampDS(i) = max(S.trDS(s)) - min(S.trDS(s));
            end
        end
        if S.hasDS
            ddSrc.Enable = 'on';  lblSrc.Text = 'rigid displacement available (px)';
        else
            ddSrc.Enable = 'off'; ddSrc.Value = 'fb'; S.src = 'fb';
            lblSrc.Text = ['rigid displacement unavailable: ' whyNo];
        end

        % a previous save wins over the defaults -- that is why it was saved
        prior = fullfile(fp,'breath_cycle_class_pc1.mat');
        if isfile(prior)
            P = load(prior);
            if isfield(P,'class_params')
                cp2 = P.class_params;
                src = 'fb';
                if isfield(cp2,'signal_source'), src = cp2.signal_source; end
                if strcmp(src,'disp') && ~S.hasDS, src = 'fb'; end
                S.src = src;  ddSrc.Value = src;
                if isfield(cp2,'thrLow') && isfield(cp2,'thrHigh')
                    S.thr.(src) = [cp2.thrLow cp2.thrHigh];
                elseif isfield(cp2,'thrLowZ') && isfield(cp2,'thrHighZ')
                    S.thr.fb = [cp2.thrLowZ cp2.thrHighZ];      % pre-dropdown files
                end
            end
            if isfield(P,'override_idx') && isfield(P,'override_class') && ...
               isfield(P,'class_final') && numel(P.class_final) == numel(foot)
                S.ovr = repmat(int8(-2), numel(foot), 1);
                S.ovr(P.override_idx) = int8(P.override_class);
            end
        end
        syncThrFields();
        S.dirty = false;
        redraw();
    end

    function [d, ok, why] = load_disp(fp, nExpect)
    %LOAD_DISP  This recording's rigid-displacement trace, in pixels.
    %  Walks up for the session-level breath_fixedmetric.mat and takes the
    %  column whose run_names entry matches this folder. A run absent from
    %  run_names has no measured registration shift behind it, so it is refused
    %  rather than silently matched to the wrong column.
        d = []; ok = false; why = '';
        [~,nm,ex] = fileparts(regexprep(fp,'[\\/]+$','')); recName = [nm ex];
        p = fp; mf = '';
        for up = 1:4
            p = fileparts(p); if isempty(p), break; end
            cand = fullfile(p,'breath_fixedmetric.mat');
            if isfile(cand), mf = cand; break; end
        end
        if isempty(mf), why = 'no breath_fixedmetric.mat above this folder'; return; end
        try
            F = load(mf,'disp','run_names');
        catch ME
            why = ME.message; return;
        end
        names = cellfun(@(c) string(c), F.run_names(:));
        j = find(names == string(recName), 1);
        if isempty(j), why = 'this run is not in breath_fixedmetric run_names'; return; end
        v = F.disp(:,j);
        v = v(1:min(numel(v), nExpect));
        if numel(v) < nExpect, v(end+1:nExpect,1) = NaN; end
        d = v; ok = true;
    end

    function a = ampNow()
        if strcmp(S.src,'disp'), a = S.ampDS; else, a = S.ampFB; end
    end
    function u = unitNow()
        if strcmp(S.src,'disp'), u = 'px'; else, u = 'z'; end
    end
    function tr = traceNow()
        if strcmp(S.src,'disp'), tr = S.trDS; else, tr = S.trFB; end
    end

    function c = classAuto()
        a = ampNow(); th = S.thr.(S.src);
        c = repmat(NORM, numel(a), 1);
        c(a <= th(1)) = FAIL;
        c(a >= th(2)) = GASP;
    end
    function c = classNow()
        c = classAuto();
        m = S.ovr ~= -2;
        c(m) = S.ovr(m);
    end

    function syncThrFields()
        th = S.thr.(S.src);
        efLow.Value = th(1); efHigh.Value = th(2);
        lblLo.Text = sprintf('FAIL  amp <= (%s)', unitNow());
        lblHi.Text = sprintf('GASP  amp >= (%s)', unitNow());
    end

    function cb_src(~,~)
        S.src = ddSrc.Value;
        syncThrFields();
        redraw();
    end

    function cb_thr(~,~)
        if efHigh.Value <= efLow.Value
            uialert(fig,'GASP threshold must be above the FAIL threshold.','Thresholds crossed');
            syncThrFields(); return;
        end
        S.thr.(S.src) = [efLow.Value efHigh.Value];
        S.dirty = true; redraw();
    end
    function cb_clearOv(~,~)
        S.ovr(:) = int8(-2); S.dirty = true; redraw();
    end

    function redraw()
        if isempty(S.trFB), return; end
        c  = classNow();  a = ampNow();  tr = traceNow();
        th = S.thr.(S.src);  u = unitNow();
        xr = [0 S.t(end)];

        cla(axT); hold(axT,'on');
        fin = tr(~isnan(tr));
        if isempty(fin), fin = 0; end
        yl = [min(fin)-0.05*range0(fin), max(fin)+0.05*range0(fin)];
        shade(axT, find(c==FAIL), yl, COL.fail);
        shade(axT, find(c==GASP), yl, COL.gasp);
        plot(axT, S.t, tr, 'k-','LineWidth',0.6,'HitTest','off');
        plot(axT, S.cyc.t, tr(S.cyc.foot), 'v','Color',[0.1 0.4 0.9], ...
             'MarkerFaceColor',[0.1 0.4 0.9],'MarkerSize',3,'HitTest','off');
        plot(axT, S.t(S.cyc.peak), tr(S.cyc.peak), '^','Color',[0.9 0.3 0.1], ...
             'MarkerFaceColor',[0.9 0.3 0.1],'MarkerSize',3,'HitTest','off');
        xlim(axT,xr); ylim(axT,yl);
        if strcmp(S.src,'disp'), ylabel(axT,'rigid displacement (px)');
        else,                    ylabel(axT,'breath fb / PC1 (z)'); end
        [~,aa,bb] = fileparts(regexprep(S.sessions{S.idx},'[\\/]+$',''));
        title(axT, sprintf('%s   |   %d cycles:  %d FAIL,  %d normal,  %d GASP%s', [aa bb], ...
              numel(c), nnz(c==FAIL), nnz(c==NORM), nnz(c==GASP), ...
              ternary(S.dirty,'   *unsaved*','')), 'Interpreter','none');

        cla(axA); hold(axA,'on');
        drawClass(axA, a, c==NORM, COL.norm, 4);
        drawClass(axA, a, c==FAIL, COL.fail, 5);
        drawClass(axA, a, c==GASP, COL.gasp, 5);
        ov = S.ovr ~= -2;
        if any(ov)
            plot(axA, S.cyc.t(ov), a(ov), 'o','MarkerEdgeColor','k', ...
                 'MarkerFaceColor','none','MarkerSize',9,'LineWidth',1,'HitTest','off');
        end
        yline(axA, th(1), '--','Color',COL.fail,'LineWidth',1,'HitTest','off');
        yline(axA, th(2), '--','Color',COL.gasp,'LineWidth',1,'HitTest','off');
        xlim(axA,xr); xlabel(axA,'time (s)');
        if strcmp(S.src,'disp'), ylabel(axA,'cycle amplitude, peak-to-peak (px)');
        else,                    ylabel(axA,'cycle amplitude, peak-foot (z)'); end

        med = median(a,'omitnan');
        lblThr.Text = sprintf(['median amp %.3g %s\nFAIL <= %.3g  (%.2f x median)\n' ...
            'GASP >= %.3g  (%.2f x median)'], med, u, th(1), th(1)/max(med,eps), ...
            th(2), th(2)/max(med,eps));
        lblOv.Text  = sprintf('%d cycle(s) overridden', nnz(ov));
        lblStat.Text = sprintf('%d/%d  |  %d FAIL (%.0f%%), %d GASP (%.0f%%)  [%s]', ...
            S.idx, numel(S.sessions), nnz(c==FAIL), 100*mean(c==FAIL), ...
            nnz(c==GASP), 100*mean(c==GASP), S.src);
    end

    function v = range0(x)
        v = max(x) - min(x); if v <= 0 || ~isfinite(v), v = 1; end
    end

    function shade(ax, ii, yl, col)
        if isempty(ii), return; end
        X = nan(5,numel(ii)); Y = nan(5,numel(ii));
        for k = 1:numel(ii)
            i = ii(k);
            aa = S.t(S.cyc.foot(i)); bb = S.t(S.cyc.nextFoot(i));
            X(:,k) = [aa bb bb aa aa]'; Y(:,k) = [yl(1) yl(1) yl(2) yl(2) yl(1)]';
        end
        patch(ax, X, Y, col, 'EdgeColor','none','FaceAlpha',0.40,'HitTest','off');
    end

    function drawClass(ax, a, m, col, sz)
        if ~any(m), return; end
        plot(ax, S.cyc.t(m), a(m), 'o','MarkerFaceColor',col, ...
             'MarkerEdgeColor','none','MarkerSize',sz,'HitTest','off');
    end

    function cb_click(ax)
        if isempty(S.cyc.foot), return; end
        p = ax.CurrentPoint(1,1);
        [~,i] = min(abs(S.cyc.t + S.cyc.ibi/2 - p));
        c = classNow();
        switch c(i)
            case NORM, nxtC = FAIL;
            case FAIL, nxtC = GASP;
            otherwise, nxtC = NORM;
        end
        auto = classAuto();
        % Record an override only where it disagrees with the thresholds, so
        % retyping a threshold cannot leave stale overrides behind.
        if nxtC == auto(i), S.ovr(i) = int8(-2); else, S.ovr(i) = nxtC; end
        S.dirty = true; redraw();
    end

    function cb_save(~,~)
        if isempty(S.cyc.foot), uialert(fig,'Nothing loaded.','Error'); return; end
        fp = S.sessions{S.idx};
        th = S.thr.(S.src);
        cyc_foot_idx      = S.cyc.foot;       %#ok<NASGU>
        cyc_peak_idx      = S.cyc.peak;       %#ok<NASGU>
        cyc_next_foot_idx = S.cyc.nextFoot;   %#ok<NASGU>
        cyc_amp      = ampNow();              %#ok<NASGU>
        cyc_amp_fb   = S.ampFB;               %#ok<NASGU>
        cyc_amp_disp = S.ampDS;               %#ok<NASGU>
        cyc_ibi_s    = S.cyc.ibi;             %#ok<NASGU>
        cyc_t        = S.cyc.t;               %#ok<NASGU>

        class_auto     = classAuto();          %#ok<NASGU>
        class_final    = classNow();
        override_idx   = find(S.ovr ~= -2);
        override_class = S.ovr(override_idx);  %#ok<NASGU>
        is_fail = class_final == FAIL;         %#ok<NASGU>
        is_gasp = class_final == GASP;         %#ok<NASGU>

        if strcmp(S.src,'disp')
            adef = 'max-min within the cycle';
        else
            adef = 'trace(peak) - trace(foot)';
        end
        class_params = struct( ...
            'signal_source', S.src, ...
            'amp_definition',adef, ...
            'units',         unitNow(), ...
            'thrLow',        th(1), ...
            'thrHigh',       th(2), ...
            'nDrop',         S.nDrop, ...
            'fps',           S.fps, ...
            'code',         'int8: -1 = FAIL, 0 = normal, +1 = GASP', ...
            'index_base',   'RAW pc1 frames, NOT nDrop-corrected', ...
            'peak_file',     fullfile(fp,'breath_peak_pc1.mat'), ...
            'start_file',    fullfile(fp,'breath_insp_start_pc1.mat'), ...
            'gui',          'breath_failure_gui_260825.m', ...
            'saved_at',      datestr(now)); %#ok<TNOW1,DATST,NASGU>

        save(fullfile(fp,'breath_cycle_class_pc1.mat'), ...
            'cyc_foot_idx','cyc_peak_idx','cyc_next_foot_idx', ...
            'cyc_amp','cyc_amp_fb','cyc_amp_disp','cyc_ibi_s','cyc_t', ...
            'class_auto','override_idx','override_class','class_final', ...
            'is_fail','is_gasp','class_params');
        S.dirty = false;
        ddSess.Items = sessionNames();
        redraw();
        uialert(fig, sprintf('%d cycles: %d FAIL, %d normal, %d GASP  [%s, %s]\n-> %s', ...
            numel(class_final), nnz(class_final==FAIL), nnz(class_final==NORM), ...
            nnz(class_final==GASP), S.src, unitNow(), ...
            fullfile(fp,'breath_cycle_class_pc1.mat')), 'Saved','Icon','success');
    end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
