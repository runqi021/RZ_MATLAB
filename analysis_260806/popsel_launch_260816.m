function popsel_launch_260816(C, tauN, ctrN, GROUP, TRIG, outDir, altFigDirs)
%POPSEL_LAUNCH_260816  The include/exclude window. See popsel_gui_260816.m.
%
%  Everything the population panels need is already in C (cached by
%  popsel_precompute_260816), so a click only re-averages arrays that are
%  already in memory. Nothing is recomputed and nothing is read from disk except
%  the per-cell PNG.

nC = numel(C);

% figures moved out of the main folder (gate-failing, merged away) still resolve.
% NOTE the loop variables here are deliberately not 'i': nested functions share
% this workspace, so a loop index reused inside one of them would clobber this.
for kk = 1:nC
    if isfile(C(kk).png), continue; end
    [~,b,e] = fileparts(C(kk).png);
    for a = 1:numel(altFigDirs)
        cand = fullfile(altFigDirs{a}, [b e]);
        if isfile(cand), C(kk).png = cand; break; end
    end
end

decF = fullfile(outDir, sprintf('popsel_decisions_%s.csv', GROUP));
dec  = repmat("undecided", nC, 1);
if isfile(decF)
    Tp = readtable(decF,'TextType','string');
    for kk = 1:height(Tp)
        j = find([C.cell] == Tp.cell(kk), 1);
        if ~isempty(j), dec(j) = Tp.decision(kk); end
    end
    fprintf('resumed %d earlier decision(s) from %s\n', nnz(dec~="undecided"), decF);
end

cur = 1;
incCol = [0.15 0.55 0.20];
excCol = [0.75 0.15 0.15];

%% ---------------- window ----------------
f = figure('Color','w','Units','normalized','Position',[0.03 0.05 0.94 0.86], ...
           'Name',sprintf('population select - %s', GROUP), 'NumberTitle','off', ...
           'KeyPressFcn',@onKey, 'CloseRequestFcn',@onClose);

axFig = axes(f,'Position',[0.015 0.10 0.50 0.86]);  axis(axFig,'off');
axDff = axes(f,'Position',[0.575 0.62 0.39 0.33]);
axHst = axes(f,'Position',[0.575 0.24 0.39 0.31]);
txt   = uicontrol(f,'Style','text','Units','normalized', ...
        'Position',[0.575 0.015 0.39 0.20],'HorizontalAlignment','left', ...
        'BackgroundColor','w','FontName','Consolas','FontSize',9);

uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.02 0.035 0.075 0.045], ...
    'String','<< Prev','FontSize',10,'Callback',@(s,e) step(-1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.10 0.035 0.075 0.045], ...
    'String','Next >>','FontSize',10,'Callback',@(s,e) step(1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.195 0.035 0.095 0.045], ...
    'String','INCLUDE (i)','FontSize',10,'FontWeight','bold', ...
    'BackgroundColor',[0.80 0.94 0.80],'Callback',@(s,e) mark("include"));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.295 0.035 0.095 0.045], ...
    'String','EXCLUDE (e)','FontSize',10,'FontWeight','bold', ...
    'BackgroundColor',[0.97 0.82 0.82],'Callback',@(s,e) mark("exclude"));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.395 0.035 0.075 0.045], ...
    'String','undo (u)','FontSize',10,'Callback',@(s,e) mark("undecided"));
hJump = uicontrol(f,'Style','edit','Units','normalized','Position',[0.485 0.035 0.05 0.045], ...
    'String','','FontSize',10,'Callback',@jumpTo, ...
    'TooltipString','type a cell id and press enter');
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.02 0.005 0.10 0.028], ...
    'String','save CSV','FontSize',9,'Callback',@(s,e) saveDecisions(true));
hStatus = uicontrol(f,'Style','text','Units','normalized','Position',[0.13 0.003 0.40 0.028], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);

show();

%% ---------------- callbacks ----------------
    function onKey(~, ev)
        switch ev.Key
            case 'i',          mark("include");
            case 'e',          mark("exclude");
            case 'u',          mark("undecided");
            case 'rightarrow', step(1);
            case 'leftarrow',  step(-1);
        end
    end

    function step(d)
        cur = max(1, min(nC, cur + d));
        show();
    end

    function jumpTo(~,~)
        v = str2double(get(hJump,'String'));
        j = find([C.cell] == v, 1);
        if isempty(j)
            set(hStatus,'String',sprintf('cell %g is not in this %s set', v, GROUP));
        else
            cur = j;  show();
        end
        set(hJump,'String','');
    end

    function mark(what)
        dec(cur) = what;
        saveDecisions(false);          % never lose a decision to a crash
        if cur < nC, cur = cur + 1; end
        show();
    end

    function saveDecisions(verbose)
        T = table([C.cell].', string({C.group}).', string({C.date}).', dec, ...
                  [C.nSpikes].', [C.IBI].', [C.logZ].', ...
                  'VariableNames',{'cell','group','date','decision','nSpikes','IBI','logZ'});
        writetable(T, decF);
        if verbose, set(hStatus,'String',sprintf('saved -> %s', decF)); end
    end

    function onClose(~,~)
        saveDecisions(false);
        fprintf('decisions saved -> %s\n', decF);
        delete(f);
    end

%% ---------------- drawing ----------------
    function show()
        % ---- left: the per-cell figure ----
        cla(axFig); axis(axFig,'off');
        if isfile(C(cur).png)
            I = imread(C(cur).png);
            image(axFig, I); axis(axFig,'image','off');
        else
            text(axFig,0.5,0.5,sprintf('no rendered figure for cell %d', C(cur).cell), ...
                 'Horizontal','center','FontSize',11);
            axis(axFig,[0 1 0 1]); axis(axFig,'off');
        end
        switch dec(cur)
            case "include", tc = incCol;  tag = 'INCLUDED';
            case "exclude", tc = excCol;  tag = 'EXCLUDED';
            otherwise,      tc = [0 0 0]; tag = 'undecided';
        end
        title(axFig, sprintf('[%d/%d]  cell %d  -  %s', cur, nC, C(cur).cell, tag), ...
              'Color',tc,'FontSize',12,'FontWeight','bold','Interpreter','none');
        drawPop();
    end

    function drawPop()
        inc = find(dec == "include");
        if strcmpi(TRIG,'onset')
            dffFld = 'dffOnsetN';  hstFld = 'histOnsetN';  trigLab = 'insp onset';
        else
            dffFld = 'dffPeakN';   hstFld = 'histPeakN';   trigLab = 'insp peak';
        end

        % ---- population dF/F ----
        cla(axDff); hold(axDff,'on');
        if ~isempty(inc)
            M = cell2mat(arrayfun(@(x) x.(dffFld), C(inc), 'uni',0).');
            mu = mean(M,1,'omitnan');
            n  = sum(isfinite(M),1);
            se = std(M,0,1,'omitnan') ./ max(sqrt(n),1);
            fill(axDff,[tauN fliplr(tauN)],[mu+se fliplr(mu-se)], [0.2 0.7 0.2], ...
                 'FaceAlpha',0.25,'EdgeColor','none');
            plot(axDff, tauN, mu, '-','Color',[0.15 0.55 0.20],'LineWidth',2);
        end
        xline(axDff,0,'-','Color',[0.35 0.75 1.00],'LineWidth',1.2);
        xlim(axDff,[-1 1]); grid(axDff,'on');
        xlabel(axDff, sprintf('time from %s (IBI)', trigLab));
        ylabel(axDff,'\DeltaF/F');
        title(axDff, sprintf('population mean \\pm SEM  (n = %d cells)', numel(inc)), ...
              'FontWeight','normal');

        % ---- population event histogram ----
        cla(axHst); hold(axHst,'on');
        if ~isempty(inc)
            H = cell2mat(arrayfun(@(x) x.(hstFld), C(inc), 'uni',0).');
            mh = mean(H,1,'omitnan');
            nh = sum(isfinite(H),1);
            sh = std(H,0,1,'omitnan') ./ max(sqrt(nh),1);
            bar(axHst, ctrN, mh, 1, 'FaceColor',[0.25 0.25 0.25],'EdgeColor','none');
            errorbar(axHst, ctrN, mh, sh, 'LineStyle','none','Color',[0.6 0.6 0.6], ...
                     'CapSize',0,'LineWidth',0.6);
        end
        xline(axHst,0,'-','Color',[0.35 0.75 1.00],'LineWidth',1.2);
        xlim(axHst,[-1 1]); grid(axHst,'on');
        xlabel(axHst, sprintf('time from %s (IBI)', trigLab));
        ylabel(axHst,'spk/cyc %');
        title(axHst,'population event histogram','FontWeight','normal');

        % ---- live statistics ----
        set(txt,'String', statText(inc));
    end

    function s = statText(inc)
        nInc = numel(inc);
        nExc = nnz(dec=="exclude");
        nUnd = nnz(dec=="undecided");

        L = strings(0,1);
        L(end+1) = sprintf('%-9s %3d included   %3d excluded   %3d undecided', ...
                           GROUP, nInc, nExc, nUnd);
        if nInc == 0
            L(end+1) = "";
            L(end+1) = "nothing included yet";
            s = char(strjoin(L, newline));  return;
        end

        ev  = sum([C(inc).nSpikes]);
        ibi = [C(inc).IBI];
        pOn = [C(inc).pOnset];  pPk = [C(inc).pPeak];
        nSig = nnz((pOn<0.05) | (pPk<0.05));

        % pooled occupancy-weighted Rayleigh. Each event keeps the weight of its
        % own recording, so pooling is a sum over the selected cells.
        S1 = 0; S2 = 0; R = 0;
        for q = inc(:).'
            w  = 1 ./ C(q).occ;                 % NaN where the bin was never visited
            e  = C(q).evSum;
            ok = isfinite(w) & isfinite(e) & e > 0;
            if ~any(ok), continue; end
            S1 = S1 + sum(e(ok).*w(ok));
            S2 = S2 + sum(e(ok).*w(ok).^2);
            R  = R  + sum(e(ok).*w(ok).*exp(1i*C(q).phaseCtrs(ok)));
        end
        if S1 > 0
            Rbar = abs(R)/S1;   nEff = S1^2/max(S2,eps);
            logZ = log(max(nEff*Rbar^2, eps));
            prefDeg = mod(rad2deg(angle(R)),360);
        else
            Rbar = NaN; nEff = NaN; logZ = NaN; prefDeg = NaN;
        end

        L(end+1) = sprintf('pooled events   %5d      IBI %.2f - %.2f s', ev, min(ibi), max(ibi));
        L(end+1) = sprintf('cells p<0.05    %3d of %3d  (%.0f%%)', nSig, nInc, 100*nSig/nInc);
        L(end+1) = "";
        L(end+1) = sprintf('pooled Rayleigh  logZ %6.2f   Rbar %.3f', logZ, Rbar);
        L(end+1) = sprintf('                 nEff %6.0f   pref %3.0f deg', nEff, prefDeg);
        L(end+1) = sprintf('                 alpha 0.05 -> 1.10 | 0.001 -> 1.93');
        L(end+1) = "";
        L(end+1) = sprintf('this cell: %d ev, IBI %.2f s, logZ %.2f', ...
                           C(cur).nSpikes, C(cur).IBI, C(cur).logZ);
        L(end+1) = sprintf('           pOnset %.4f  pPeak %.4f', C(cur).pOnset, C(cur).pPeak);
        s = char(strjoin(L, newline));
    end
end
