function breath_detect_all()
% breath_detect_all  Run the DERIVATIVE onset detector on every PC1 trace
% (Vglut2 + ChAT), with a hard inspiration-duration cut [0.05, 1.0] s, and save
% one full-trace figure per session.
%
%   signal  : raw PC1, auto-inverted by skewness (inspiration = sharp UP spike)
%   peaks   : islocalmax (expiration onset)
%   onset   : per peak, anchor at steepest rising sample, walk back to the local
%             min (d/dt sign flip) = inspiration onset
%   cut     : keep a peak/onset pair only if (peak-onset) in [0.05, 1.0] s;
%             rejected pairs drawn in grey
%
% Output: one PNG per session in <root>\_detect_derivative\ + a console summary.

ROOTS  = {'D:\Ventral_surface_summary\Vglut2', ...
          'D:\Ventral_surface_summary\ChAT'};
OUTDIR = 'D:\Ventral_surface_summary\_detect_derivative';
if ~isfolder(OUTDIR), mkdir(OUTDIR); end
DUR = [0.05 1.0];        % hard inspiration-duration cut (s)

jobs = {};
for ri = 1:numel(ROOTS)
    h = dir(fullfile(ROOTS{ri}, '**', 'breath_pc1.mat'));
    [~, ds] = fileparts(ROOTS{ri});
    for j = 1:numel(h), jobs{end+1} = struct('folder',h(j).folder,'ds',ds); end %#ok<AGROW>
end
fprintf('%-7s %-34s %5s %6s %5s %5s %6s\n','ds','session','fps','f0','kept','rej','medDur');

for i = 1:numel(jobs)
    folder = jobs{i}.folder;
    P = load(fullfile(folder,'breath_pc1.mat'));
    x = zsc(double(P.breathTrace(:)));
    if mean(x.^3) < 0, x = -x; end                       % inspiration -> UP
    fps = 30; if isfield(P,'fps'), fps = double(P.fps); end
    n = numel(x); t = (0:n-1)'/fps;
    f0 = dominant_freq(x, fps, [0.3 3.5]);

    pk = find(islocalmax(x,'MinProminence',0.5*std(x),'MinSeparation',max(1,round(0.6/f0*fps))));
    base = prctile(x, 10);                               % robust low baseline
    on = zeros(size(pk));
    for k = 1:numel(pk)
        p1 = pk(k);
        if k==1, p0 = max(1, p1-round(1.5/f0*fps)); else, p0 = pk(k-1); end
        if k==1 && min(x(1:p1)) > base + 0.2*max(x(p1)-base,eps)
            on(k) = NaN; continue;            % recording started mid-inspiration
        end
        if p1 <= p0+2, on(k) = p0; continue; end
        seg = (p0:p1)'; s = x(seg); [~,mi] = max([diff(s);0]);    % steepest rising sample
        delta = 0.2*max(x(p1)-base, eps);                % "notch" tolerance above baseline
        % walk back to the baseline foot: keep going while descending OR still
        % above baseline+delta (so above-baseline notches/doublets don't trap it),
        % bounded by the previous peak (p0).
        cur = mi; while cur>1 && (s(cur-1) < s(cur) || s(cur) > base+delta), cur = cur-1; end
        on(k) = seg(cur);
    end
    dur = (pk - on)/fps;                                 % inspiration duration (NaN if truncated)
    keep = dur >= DUR(1) & dur <= DUR(2);                 % NaN -> false (dropped)

    [~, leaf] = fileparts(folder);
    f = figure('Visible','off','Color','w','Position',[40 40 1700 360]);
    ax = axes(f); hold(ax,'on'); grid(ax,'on');
    plot(ax, t, x, 'Color',[0 0 0], 'LineWidth',0.6);
    % rejected (grey)
    rj = ~keep;
    plot(ax, t(pk(rj)), x(pk(rj)), '^', 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',5);
    onRej = on(rj); onRej = onRej(~isnan(onRej));        % truncated first has NaN onset
    plot(ax, t(onRej), x(onRej), 'o', 'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor','none','MarkerSize',5);
    % kept
    plot(ax, t(pk(keep)), x(pk(keep)), '^', 'MarkerFaceColor',[0.9 0.25 0.25],'MarkerEdgeColor','none','MarkerSize',6);
    plot(ax, t(on(keep)), x(on(keep)), 'o', 'MarkerFaceColor',[0.1 0.45 0.95],'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',0.4);
    xlim(ax,[t(1) t(end)]); xlabel(ax,'time (s)'); ylabel(ax,'PC1 (z, inv)');
    title(ax, sprintf('%s / %s   fps=%g  f0=%.2f Hz   kept %d  rejected %d   median insp=%.0f ms', ...
        jobs{i}.ds, strrep(leaf,'_','\_'), fps, f0, sum(keep), sum(rj), 1000*median(dur(keep))), 'FontSize',9);
    exportgraphics(f, fullfile(OUTDIR, sprintf('%s_%s.png', jobs{i}.ds, leaf)), 'Resolution',130);
    close(f);

    fprintf('%-7s %-34s %5g %6.2f %5d %5d %5.0f\n', jobs{i}.ds, leaf(1:min(34,end)), ...
        fps, f0, sum(keep), sum(rj), 1000*median(dur(keep)));
end
fprintf('\nSaved %d full-trace figures -> %s\n', numel(jobs), OUTDIR);
end

%% helpers
function y = zsc(x)
x = double(x(:)); s = std(x,'omitnan'); if s==0, s=1; end
y = (x-mean(x,'omitnan'))/s;
end
function f0 = dominant_freq(x, fps, band)
x = x(:)-mean(x(:)); n = numel(x);
P = abs(fft(x.*hann(n))).^2; P = P(1:floor(n/2)); fa = (0:floor(n/2)-1)'*(fps/n);
m = fa>=band(1) & fa<=band(2); [~,j] = max(P.*m); f0 = fa(j);
if ~isfinite(f0) || f0<=0, f0 = 1; end
end
