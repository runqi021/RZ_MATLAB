function breath_onset_methods()
% breath_onset_methods  Compare PARAMETER-FREE inspiration-onset detectors on
% raw PC1, anchored to the (reliable) peaks so nothing needs per-trace tuning.
%
% Signal = raw PC1, auto-inverted by skewness (inspiration = sharp UP spike).
% NO amplitude filtering (raw is physically faithful; the sharp foot is the
% high-frequency content a filter would destroy).
%
% Methods (one FIGURE each, 3 sessions = rows, fast/noisy -> slow/clean):
%   M1 Trough        - baseline minimum of the cycle (no tuning)
%   M2 LineFitFoot   - fit the steep rise, extrapolate to baseline (your method)
%   M3 ChordFoot     - point of max gap below the baseline->peak chord (Kneedle;
%                      parameter-free, self-scaling per cycle)  <-- best bet
%   M4 SlopeWalkback - walk back from peak until slope < 20% of THIS cycle's max
%                      (fraction is scale-invariant, so one value fits all traces)
% Base/Stats/Image toolboxes only.

SESS = {
 'D:\Ventral_surface_summary\Vglut2\0224\cell3\pFN_roi1_z0_512x512_6x_2000f_00001'
 'D:\Ventral_surface_summary\ChAT\0522\IO\roi3_6x_x-700y600z-25_3000f_16lp_00001'
 'D:\Ventral_surface_summary\ChAT\0522\cell1\roi3_8x_x-1070y730z0_3000f_15lp_00001'};
OUTDIR = 'D:\Ventral_surface_summary\_onset_methods';
if ~isfolder(OUTDIR), mkdir(OUTDIR); end

methods   = {'M1_Trough','M2_LineFitFoot','M3_ChordFoot','M4_SlopeWalkback'};
methTitle = {'M1: Trough (cycle minimum)','M2: Line-fit foot (your method)', ...
             'M3: Chord-foot (Kneedle, parameter-free)','M4: Derivative sign-flip (local min at base of rise)'};

D = struct([]);
for s = 1:numel(SESS)
    P = load(fullfile(SESS{s},'breath_pc1.mat'));
    x = zsc(double(P.breathTrace(:)));
    if mean(x.^3) < 0, x = -x; end                 % auto-invert: sharp spikes UP
    fps = 30; if isfield(P,'fps'), fps = double(P.fps); end
    n = numel(x); t = (0:n-1)'/fps;
    f0 = dominant_freq(x, fps, [0.3 3.5]);
    pks = find_peaks_simple(x, fps, f0);           % peaks on RAW signal
    [~, leaf] = fileparts(SESS{s});
    win = t <= min(t(end), max(8, 10/f0));
    D(s).leaf=leaf; D(s).t=t; D(s).x=x; D(s).fps=fps; D(s).f0=f0; D(s).pks=pks; D(s).win=win;
end

for m = 1:numel(methods)
    f = figure('Color','w','Position',[60 60 1250 820]);
    for s = 1:numel(D)
        on = compute_onsets(D(s).x, D(s).fps, D(s).pks, D(s).f0, m);
        ax = subplot(numel(D),1,s); hold(ax,'on'); grid(ax,'on');
        w = D(s).win;
        plot(ax, D(s).t(w), D(s).x(w), 'Color',[0 0 0], 'LineWidth',0.9);
        pkw = D(s).pks(ismember(D(s).pks, find(w)));
        plot(ax, D(s).t(pkw), D(s).x(pkw), '^', 'MarkerFaceColor',[0.9 0.3 0.3], ...
            'MarkerEdgeColor','none', 'MarkerSize',8);
        onw = on(ismember(on, find(w)));
        plot(ax, D(s).t(onw), D(s).x(onw), 'o', 'MarkerFaceColor',[0.1 0.45 0.95], ...
            'MarkerEdgeColor','k', 'MarkerSize',8, 'LineWidth',0.5);
        xlim(ax,[D(s).t(find(w,1)) D(s).t(find(w,1,'last'))]);
        ylabel(ax,'PC1 (z, inverted)');
        title(ax, sprintf('%s   |   %.2f Hz   |   %d peaks / %d onsets', ...
            strrep(D(s).leaf,'_','\_'), D(s).f0, numel(D(s).pks), numel(on)), 'FontSize',9);
        if s==numel(D), xlabel(ax,'time (s)'); end
        if s==1, legend(ax,{'PC1','peak (expir. onset)','inspiration onset'},'Location','northeast','FontSize',8); end
    end
    sgtitle(sprintf('%s     [raw, auto-inverted]', methTitle{m}), 'FontWeight','bold');
    exportgraphics(f, fullfile(OUTDIR, [methods{m} '.png']), 'Resolution',130);
end
fprintf('Saved %d method figures -> %s\n', numel(methods), OUTDIR);
end

%% ---- onset methods (all anchored to peaks; per-cycle, self-scaling) ----
function on = compute_onsets(sig, fps, peaks, f0, method)
sig = sig(:); n = numel(sig); peaks = sort(peaks(:)); on = [];
lookback = round(1.5/f0*fps);
for k = 1:numel(peaks)
    p1 = peaks(k);
    if k==1, p0 = max(1, p1-lookback); else, p0 = peaks(k-1); end
    if p1 <= p0+2, continue; end
    s = sig(p0:p1); [~, it] = min(s); B = p0 + it - 1;     % baseline min before peak
    rise = (B:p1)';
    switch method
        case 1                                              % trough
            o = B;
        case 2                                              % line-fit foot
            o = linefit_foot(sig, fps, p1, p0);
        case 3                                              % chord-foot (Kneedle)
            if numel(rise) < 3, o = B;
            else
                y = sig(rise); chord = y(1) + (y(end)-y(1))*(0:numel(rise)-1)'/(numel(rise)-1);
                [~, jj] = max(chord - y);                   % max gap BELOW the chord
                o = rise(jj);
            end
        case 4                                              % derivative sign-flip:
            if numel(rise) < 4, o = B;                      % anchor at steepest
            else                                            % sample, walk back to the
                d = [diff(sig(rise)); 0]*fps; [~, mi] = max(d);   % LOCAL MIN (d: -/+)
                cur = mi;
                while cur > 1 && sig(rise(cur-1)) < sig(rise(cur)), cur = cur - 1; end
                o = rise(cur);
            end
    end
    on(end+1) = o; %#ok<AGROW>
end
on = unique(min(max(on(:),1),n));
end

function o = linefit_foot(breath, fps, p, lo)
breath = breath(:); idx_w = (lo:p)'; seg = breath(idx_w);
d = (seg(2:end)-seg(1:end-1))*fps; d = movmean(d, max(1,round(0.05*fps)));
max_d = max(d);
if ~isfinite(max_d) || max_d<=0, o = lo; return; end
steep = d > 0.30*max_d; s_idx = find(steep);
if isempty(s_idx), o = lo; return; end
g = find(diff(s_idx)>1);
if isempty(g), run = s_idx; else, run = s_idx(g(end)+1:end); end
if numel(run) < 2, o = idx_w(run(1)); return; end
ff = idx_w(run)+0.5; sf = (seg(run)+seg(run+1))/2;
pf = polyfit(ff, sf, 1); m = pf(1); c = pf(2);
if abs(m) < eps, o = lo; return; end
pre = run(1)-1; if pre>=1, base = median(seg(1:pre)); else, base = seg(1); end
xf = (base - c)/m; o = round(max(lo, min(p-1, xf)));
end

%% ---- shared helpers ----
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

function pk = find_peaks_simple(sig, fps, f0)
sep = max(1, round(0.6/f0*fps));
pk = find(islocalmax(sig, 'MinProminence', 0.5*std(sig), 'MinSeparation', sep));
end
