function breath_doublet_diag()
% breath_doublet_diag  Find DOUBLETS (two inspiratory pulses stacked within one
% expected cycle) across Vglut2+ChAT PC1 traces, and show how the peak/derivative
% detector handles them. A doublet = two prominent local maxima closer than
% ~0.5x the median cycle, where the dip between them does NOT return to baseline.

ROOTS = {'D:\Ventral_surface_summary\Vglut2', ...
         'D:\Ventral_surface_summary\ChAT'};
OUT   = 'D:\Ventral_surface_summary\_detect_derivative\DIAG_doublets.png';

hits = {};
total = 0;
for ri = 1:numel(ROOTS)
    h = dir(fullfile(ROOTS{ri},'**','breath_pc1.mat'));
    for j = 1:numel(h)
        folder = h(j).folder;
        P = load(fullfile(folder,'breath_pc1.mat'));
        x = zsc(double(P.breathTrace(:))); if mean(x.^3)<0, x=-x; end
        fps = 30; if isfield(P,'fps'), fps=double(P.fps); end
        n=numel(x); f0=dom(x,fps);
        % current (coarse) peaks vs fine peaks (small separation) to expose doublets
        sepC = max(1,round(0.6/f0*fps));
        pkF = find(islocalmax(x,'MinProminence',0.5*std(x),'MinSeparation',max(1,round(0.12*fps))));
        if numel(pkF)<3, continue; end
        base = prctile(x,10);
        medIPI = median(diff(pkF));
        for q = 1:numel(pkF)-1
            a=pkF(q); b=pkF(q+1); ipi=b-a;
            if ipi < 0.5*medIPI && ipi < sepC          % two pulses very close
                dip = min(x(a:b));
                amp = min(x(a),x(b)) - base;
                if amp>0 && (dip-base)/amp > 0.30       % dip doesn't return to baseline
                    [~,leaf]=fileparts(folder);
                    hits{end+1} = struct('x',x,'fps',fps,'f0',f0,'a',a,'b',b, ...
                        'leaf',leaf,'depth',(dip-base)/amp,'ipi_ms',1000*ipi/fps); %#ok<AGROW>
                    total = total+1;
                end
            end
        end
    end
end
fprintf('Found %d doublet candidates across all sessions.\n', total);
if isempty(hits), fprintf('No doublets detected by this criterion.\n'); return; end

% sort by shallowest dip (most "stacked") and show all, with current vs fixed onset
dep = cellfun(@(s) s.depth, hits); [~,ord]=sort(dep,'descend'); hits=hits(ord);
N = numel(hits);
f=figure('Visible','off','Color','w','Position',[40 40 1500 420]);
for k=1:N
    s=hits{k}; x=s.x; fps=s.fps; n=numel(x); f0=s.f0;
    pk = s.a; if x(s.b)>x(s.a), pk=s.b; end                 % detected (taller) peak
    base = prctile(x,10); delta = 0.20*(x(pk)-base);
    lo2 = max(1, s.a-round(1.2/f0*fps)); seg=(lo2:pk)'; ss=x(seg); dpf=[diff(ss);0]; [~,mi]=max(dpf);
    c=mi; while c>1 && ss(c-1)<ss(c), c=c-1; end; oCur=seg(c);             % stops at notch
    c=mi; while c>1 && (ss(c-1)<ss(c) || ss(c)>base+delta), c=c-1; end; oFix=seg(c); % -> baseline foot
    lo=max(1, s.a-round(0.7/f0*fps)); hi=min(n, s.b+round(0.7/f0*fps));
    ax=subplot(1,N,k); hold(ax,'on'); grid(ax,'on'); w=lo:hi; t=(0:n-1)'/fps;
    plot(ax,t(w),x(w),'k','LineWidth',1);
    plot(ax,t(s.a),x(s.a),'^','MarkerFaceColor',[0.9 0.25 0.25],'MarkerEdgeColor','none','MarkerSize',8);
    plot(ax,t(s.b),x(s.b),'^','MarkerFaceColor',[0.9 0.55 0.1],'MarkerEdgeColor','none','MarkerSize',8);
    plot(ax,t(oCur),x(oCur),'o','MarkerFaceColor',[0.1 0.45 0.95],'MarkerEdgeColor','k','MarkerSize',9);
    plot(ax,t(oFix),x(oFix),'s','MarkerFaceColor',[0.1 0.7 0.2],'MarkerEdgeColor','k','MarkerSize',9);
    yline(ax, base, ':', 'baseline');
    xlim(ax,[t(lo) t(hi)]); set(ax,'XTickLabel',[]);
    title(ax,sprintf('%s  gap=%.0fms', strrep(s.leaf,'_','\_'), s.ipi_ms),'FontSize',8);
    if k==1, legend(ax,{'PC1','pulse1','pulse2','onset: current (notch)','onset: fixed (foot)'},'FontSize',7,'Location','northwest'); end
end
sgtitle('Doublet onset: current stops at notch (blue) vs fixed walks to baseline foot (green)','FontWeight','bold');
exportgraphics(f,OUT,'Resolution',130); close(f);
fprintf('Saved -> %s\n', OUT);
end

function y=zsc(x), x=double(x(:)); s=std(x,'omitnan'); if s==0,s=1; end; y=(x-mean(x,'omitnan'))/s; end
function f0=dom(x,fps)
x=x(:)-mean(x(:)); n=numel(x); Pp=abs(fft(x.*hann(n))).^2; Pp=Pp(1:floor(n/2));
fa=(0:floor(n/2)-1)'*(fps/n); m=fa>=0.3&fa<=3.5; [~,j]=max(Pp.*m); f0=fa(j); if ~isfinite(f0)||f0<=0,f0=1; end
end
