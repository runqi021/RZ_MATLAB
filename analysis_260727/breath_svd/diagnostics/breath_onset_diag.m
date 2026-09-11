function breath_onset_diag()
% breath_onset_diag  Zoom into individual cycles of one slow-breather session to
% see mis-detected onsets and test the fix: local-min walk-back (current) vs a
% dx slope-threshold foot (stop where the rise slope first exceeds TAU per frame).

FOLDER = 'D:\Ventral_surface_summary\ChAT\0522\cell1\roi3_8x_x-1070y730z0_3000f_15lp_00001';
OUT    = 'D:\Ventral_surface_summary\_detect_derivative\DIAG_slow.png';
TAU    = 0.15;     % per-frame z-slope threshold for the dx foot
NSHOW  = 12;       % cycles to show

P = load(fullfile(FOLDER,'breath_pc1.mat'));
x = zsc(double(P.breathTrace(:))); if mean(x.^3)<0, x=-x; end
fps = 30; if isfield(P,'fps'), fps=double(P.fps); end
n=numel(x); t=(0:n-1)'/fps;
f0 = dom(x,fps);
pk = find(islocalmax(x,'MinProminence',0.5*std(x),'MinSeparation',max(1,round(0.6/f0*fps))));

oMin = zeros(size(pk)); oDx = zeros(size(pk));
for k=1:numel(pk)
    p1=pk(k); if k==1, p0=max(1,p1-round(1.5/f0*fps)); else, p0=pk(k-1); end
    if p1<=p0+2, oMin(k)=p0; oDx(k)=p0; continue; end
    seg=(p0:p1)'; s=x(seg); dpf=[diff(s);0]; [~,mi]=max(dpf);
    % current: walk back to local min
    c=mi; while c>1 && s(c-1)<s(c), c=c-1; end; oMin(k)=seg(c);
    % dx: walk back until rise slope drops below TAU (start of steep rise)
    c=mi; while c>1 && dpf(c-1)>TAU, c=c-1; end; oDx(k)=seg(c);
end

f=figure('Visible','off','Color','w','Position',[40 40 1500 820]);
rows=3; cols=ceil(NSHOW/rows);
for k=1:min(NSHOW,numel(pk))
    p1=pk(k); lo=max(1,p1-round(0.9/f0*fps)); hi=min(n,p1+round(0.3/f0*fps));
    ax=subplot(rows,cols,k); hold(ax,'on'); grid(ax,'on');
    w=lo:hi; plot(ax,t(w),x(w),'k','LineWidth',1);
    plot(ax,t(p1),x(p1),'^','MarkerFaceColor',[0.9 0.25 0.25],'MarkerEdgeColor','none','MarkerSize',9);
    plot(ax,t(oMin(k)),x(oMin(k)),'o','MarkerFaceColor',[0.1 0.45 0.95],'MarkerEdgeColor','k','MarkerSize',9);
    plot(ax,t(oDx(k)),x(oDx(k)),'s','MarkerFaceColor',[0.1 0.7 0.2],'MarkerEdgeColor','k','MarkerSize',9);
    xlim(ax,[t(lo) t(hi)]); set(ax,'XTickLabel',[]); title(ax,sprintf('cycle %d',k),'FontSize',8);
    if k==1, legend(ax,{'PC1','peak','onset: local-min','onset: dx>%g'},'FontSize',7,'Location','northwest'); end
end
sgtitle(sprintf('Slow breather onsets:  local-min (blue) vs dx-threshold TAU=%.2f (green)', TAU),'FontWeight','bold');
exportgraphics(f,OUT,'Resolution',130); close(f);
fprintf('median dt(local-min->dx) = %.0f ms\n', 1000*median((oMin-oDx)/fps));
fprintf('Saved -> %s\n', OUT);
end

function y=zsc(x), x=double(x(:)); s=std(x,'omitnan'); if s==0,s=1; end; y=(x-mean(x,'omitnan'))/s; end
function f0=dom(x,fps)
x=x(:)-mean(x(:)); n=numel(x); Pp=abs(fft(x.*hann(n))).^2; Pp=Pp(1:floor(n/2));
fa=(0:floor(n/2)-1)'*(fps/n); m=fa>=0.3&fa<=3.5; [~,j]=max(Pp.*m); f0=fa(j); if ~isfinite(f0)||f0<=0,f0=1; end
end
