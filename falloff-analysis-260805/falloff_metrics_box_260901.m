function falloff_metrics_box_260901()
%FALLOFF_METRICS_BOX_260901  top5-bot5 falloff inside the rendered box.
%
% Two panels per dataset, WT and shiverer side by side:
%     1  top - bot                 raw, ADU
%     2  (top - bot) / mW^2        normalised by power
%
% TWO FITS, and they are different things:
%   * GLOBAL (red, power-normalised panel only): one exponential over the whole
%     fitted depth range. This is where l comes from. It is not drawn on the raw
%     panel because the power steps break that curve and a single exponential
%     through it would be meaningless.
%   * PER POWER ROUND (orange, both panels): the SAME l, HELD FIXED, with only
%     the amplitude free, fitted separately inside each round. Because l is not
%     refitted these are a TEST rather than a fit -- if a round's points do not
%     follow a line of the global slope, that round really does decay
%     differently. Each segment is annotated with its R^2 at that fixed l. On
%     the raw panel the segments also make the power jumps legible: parallel
%     runs offset in amplitude.
%
% TWO THINGS DIFFER FROM falloff_260804.m, both because this describes the block
% that was actually rendered rather than the whole acquisition:
%
%  1. THE METRICS COME FROM THE 200x200 um SUB-VOLUME, not the full 456 um field.
%     top/bot are the mean of the brightest/dimmest PCT_FRAC of pixels IN THE
%     BOX, so a large surface vessel outside the box no longer sets `top` -- the
%     render and the falloff curve describe one object.
%
%  2. ONE POINT PER DEPTH, NO OVERLAP DOTS. falloff_260804 plots every plane of
%     every power round, so where rounds overlap a depth carries two or three
%     dots at different powers. Here the metrics are read off the ASSEMBLED
%     volume (falloff_proj_<tag>_vol*.mat), where build_volume already resolved
%     each overlap at its MIDPOINT -- the shallower round keeps the top half, the
%     deeper round takes the bottom half. Every depth appears exactly once, from
%     the round that was exposed for it, i.e. the same planes the 3D render uses,
%     and pctPlane carries that round's power for the normalisation.
%
% Definitions are falloff_260804.m's, unchanged:
%     top    = mean of the brightest PCT_FRAC of pixels in the frame-averaged plane
%     bot    = mean of the dimmest  PCT_FRAC        (detector offset + background)
%     diff   = top - bot            -- the additive floor removed, because bot does
%                                      not scale with P^2 and would otherwise
%                                      dominate the deep end
%     diff_n = diff / mW^P_EXPONENT
% l is refit here on the boxed data, so it will not match the whole-field l in
% falloff_<tag>.mat exactly -- that is the intended comparison, not an error.
%
% Runqi Zhang / 2026-09-01

%% --------------------------- USER PARAMETERS -------------------------------
% {dataDir, roiTag, label, ML_RANGE, AP_RANGE, powerLut, FIT_RANGE}
% Same boxes as falloff_render3d_260831.
% FIT_RANGE overrides the FIT_DEPTH stored in falloff_<tag>.mat; [] keeps the
% stored one. The WT is extended to 200 from its stored 150, which also brings
% its deepest power round (238.5 mW, 175-265 um) inside the fit for the first
% time -- at 150 that round had no points in range and got no segment at all.
% Sensitivity, for reference: shiverer 60.1 um to 500 vs 58.2 um to 450.
DS = { ...
  'C:\fall-off\260728_vglut2_vessel',   'roi1', 'Vglut2 (WT)', [100 300], [50  250], 'post_260723', [10 200]
  'C:\fall-off\260804_shiver_dbh_vessel','roi2','shiverer',    [0   200], [100 300], 'post_260723', [] };

PIA_MODE   = 'off';       % re-zero depth on the MEASURED pia rather than the
                           % filename z0:
                           %  'auto'  pia = peak of the 7 um smoothed top-bot
                           %          curve, per dataset. Measured, not guessed:
                           %          +3 um (WT) and -3 um (shiverer) on this
                           %          data, i.e. filename z0 is a few um off in
                           %          OPPOSITE directions, so the two were not
                           %          referenced to the same landmark.
                           %  'off'   keep the filename z0
                           % l IS UNAFFECTED by this: a constant depth shift
                           % changes the fitted AMPLITUDE, not the slope. It
                           % moves the fit-range edges by a few um, which is the
                           % only way it can touch l at all.
PIA_SMOOTH = 7;            % um, smoothing before the peak search. The shiverer
                           % plateaus over ~15 um near the surface, so the raw
                           % argmax there is set by one bright plane.
XLIM       = [-10 500];    % SHARED depth axis for every panel, so the WT and the
                           % shiverer are read at one scale. Starts at -10 to
                           % keep the above-pia points the extended cache holds:
                           % the rise into the surface is what locates the pia,
                           % and clipping at 0 hides it.
VOL_FROM   = -10;          % which cache: 0 = pia-start, -10 = extended
PCT_FRAC   = 0.05;         % top/bot fraction -- falloff_260804's value
P_EXPONENT = 2;
PIX_TARGET = 1.0;
OUT        = 'C:\fall-off\falloff_metrics_box_260901';
FS         = 0.75;
%% ---------------------------------------------------------------------------

here = fileparts(mfilename('fullpath')); addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

n = size(DS,1);  M = cell(n,1);
for k = 1:n
    M{k} = box_metrics(DS{k,1}, DS{k,2}, DS{k,4}, DS{k,5}, DS{k,6}, ...
                       VOL_FROM, PCT_FRAC, P_EXPONENT, PIX_TARGET, DS{k,7}, ...
                       PIA_MODE, PIA_SMOOTH);
end

%% --- figure: 2 metrics down, 2 datasets across -------------------------------
% {field, ylabel, subtitle, y-scale}. The raw row is LINEAR: it is there to show
% the power steps as amplitude jumps, and a log axis flattens those into equal
% offsets. The power-normalised row stays log, where a single exponential is a
% straight line and the fit can be judged by eye.
PAN = { 'diff',   'top 5% - bot 5%  (ADU)',                   'raw',  'linear'
        'diff_n', sprintf('(top - bot) / mW^%g', P_EXPONENT), 'normalised by power', 'log' };

f = figure('Color','w','Units','pixels','Position',[30 30 round(1200*FS) round(900*FS)]);
set(f,'DefaultAxesFontSize',max(5,10*FS),'DefaultTextFontSize',max(5,10*FS));
tl = tiledlayout(f, size(PAN,1), n, 'TileSpacing','compact','Padding','compact');

for p = 1:size(PAN,1)
    for k = 1:n
        ax = nexttile(tl, (p-1)*n + k); hold(ax,'on');
        B = M{k};  y = B.(PAN{p,1});  L = B.L;
        ok = isfinite(y) & y > 0;

        % Colour by power round, so a step in the raw panel is visibly a change
        % of laser power rather than a change in the tissue.
        cols = lines(max(numel(B.rounds),3));
        for q = 1:numel(B.rounds)
            m = false(size(y)); m(B.rounds{q}) = true; m = m & ok;
            plot(ax, B.depth(m), y(m), '.', 'Color', cols(q,:), 'MarkerSize',6, ...
                 'DisplayName', sprintf('%.1f mW', B.roundMW(q)));
        end
        set(ax,'YScale',PAN{p,4});
        yl0 = ylim(ax);
        if strcmpi(PAN{p,4},'linear'), yl0(1) = 0; ylim(ax, yl0); end
        fitLo = max(B.fitD(1), min(B.depth));  fitHi = min(B.fitD(2), max(B.depth));

        % POWER-NORMALISED row: ONE fit, drawn only across the fitted range, so
        % the line marks its own extent and no shading is needed to say where it
        % applies. No per-round segments here -- with the power divided out the
        % rounds already lie on one curve, and overlaying four more lines on it
        % only obscures where they depart from it.
        if strcmp(PAN{p,1},'diff_n') && isfinite(L)
            dd = linspace(fitLo, fitHi, 200);
            plot(ax, dd, B.A*exp(-dd/L), '-', 'Color',[0.85 0.15 0.15], ...
                 'LineWidth',1.8, 'HandleVisibility','off');
        end

        % RAW row: per-round, l fixed at the global value, amplitude free. This
        % is where the segments earn their place -- they make the power steps
        % legible as parallel runs offset in amplitude.
        if isfinite(L) && strcmp(PAN{p,1},'diff')
            inFit = B.depth >= B.fitD(1) & B.depth <= B.fitD(2);
            Aq = nan(numel(B.rounds),1); d0q = Aq; d1q = Aq; r2q = Aq;
            for q = 1:numel(B.rounds)
                idx = B.rounds{q};
                % Only points inside the GLOBAL fit range: outside it the curve
                % is either the above-pia rise or the noise floor, and neither
                % is the exponential these segments are testing.
                m = idx(ok(idx) & inFit(idx));
                if numel(m) < 3, continue; end
                dq = B.depth(m);  yq = y(m);
                Aq(q)  = exp(mean(log(yq) + dq/L));
                rr     = log(yq) - (log(Aq(q)) - dq/L);
                r2q(q) = 1 - sum(rr.^2)/max(sum((log(yq)-mean(log(yq))).^2), eps);
                d0q(q) = min(dq);  d1q(q) = max(dq);
            end
            % Extend consecutive segments to meet at the round BOUNDARY, then
            % join them with a vertical riser, so the four fits read as one
            % sawtooth. The riser is the power step: same tissue, same l, only
            % the amplitude changes, so a vertical line is literally what
            % happens between the two models at that depth.
            good = find(isfinite(Aq));
            e0 = d0q; e1 = d1q;
            for j = 1:numel(good)-1
                a = good(j); b = good(j+1);
                db = (d1q(a) + d0q(b))/2;
                e1(a) = db;  e0(b) = db;
            end
            for q = good(:).'
                dd = linspace(e0(q), e1(q), 100);
                plot(ax, dd, Aq(q)*exp(-dd/L), '-', 'Color',[0.95 0.55 0.10], ...
                     'LineWidth',1.4, 'HandleVisibility','off');
                text(ax, d0q(q), Aq(q)*exp(-d0q(q)/L), sprintf('%.2f ', r2q(q)), ...
                     'Color',[0.70 0.40 0.05], 'FontSize',max(5,7.5*FS), ...
                     'HorizontalAlignment','right','VerticalAlignment','middle');
            end
            for j = 1:numel(good)-1
                a = good(j); b = good(j+1);
                db = e1(a);
                plot(ax, [db db], [Aq(a)*exp(-db/L), Aq(b)*exp(-db/L)], '-', ...
                     'Color',[0.95 0.55 0.10], 'LineWidth',1.4, 'HandleVisibility','off');
            end
        end
        ylim(ax, yl0);
        xline(ax, 0, ':', 'Color',[0.5 0.5 0.5], 'HandleVisibility','off');
        grid(ax,'on'); ax.GridAlpha = 0.12; box(ax,'off');
        xlim(ax, XLIM);
        if p == size(PAN,1), xlabel(ax,'depth (\mum)'); end
        if k == 1, ylabel(ax, PAN{p,2}); end

        if strcmp(PAN{p,1},'diff_n') && isfinite(L)
            text(ax, 0.04, 0.06, sprintf('global l = %.1f \\mum   R^2 = %.3f', L, B.r2), ...
                 'Units','normalized','Color',[0.85 0.15 0.15], ...
                 'FontSize',max(6,9*FS),'VerticalAlignment','bottom');
        else
            text(ax, 0.04, 0.06, sprintf('orange: l = %.1f \\mum fixed, amplitude free', L), ...
                 'Units','normalized','Color',[0.70 0.40 0.05], ...
                 'FontSize',max(6,8.5*FS),'VerticalAlignment','bottom');
        end
        if p == 1
            title(ax, sprintf('%s  --  %s\nML %g-%g, AP %g-%g \\mum   (%s)', ...
                  DS{k,3}, DS{k,2}, DS{k,4}(1), DS{k,4}(2), DS{k,5}(1), DS{k,5}(2), PAN{p,3}), ...
                  'FontWeight','normal');
            legend(ax,'Location','northeast','Box','off','FontSize',max(5,8*FS));
        else
            title(ax, PAN{p,3}, 'FontWeight','normal','FontSize',max(6,10*FS));
        end
    end
end

title(tl, {sprintf(['top %g%% - bot %g%% inside the RENDERED %g x %g \\mum box'], ...
        100*PCT_FRAC, 100*PCT_FRAC, diff(DS{1,4}), diff(DS{1,5})), ...
        ['one point per depth -- read off the assembled volume, so overlapping rounds ' ...
         'are already split at their midpoint   |   grey band = global fit range   |   ' ...
         'orange = global l held fixed, amplitude free, per round (number = R^2)']}, ...
        'FontWeight','normal','FontSize',max(6,10*FS));

exportgraphics(f,[OUT '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(f,[OUT '.pdf'],'ContentType','vector','BackgroundColor','white');
d = dir([OUT '.png']);
fprintf('\nsaved %s.png (%.0f KB) + .pdf\n', OUT, d.bytes/1024);

% Save the boxed fits so falloff_summary_scatter_260805 can pick them up instead
% of the whole-field ones. Written as a table rather than hardcoded into that
% script, so changing a box or a fit range here updates the scatter too and the
% two figures cannot quietly disagree about l.
BOX = struct('session',{},'tag',{},'ell',{},'r2',{},'fit',{},'ml',{},'ap',{},'volFrom',{});
for k = 1:n
    [~, sess] = fileparts(DS{k,1});
    BOX(k).session = sess;      BOX(k).tag = DS{k,2};
    BOX(k).ell = M{k}.L;        BOX(k).r2  = M{k}.r2;
    BOX(k).fit = M{k}.fitD;     BOX(k).ml  = DS{k,4};
    BOX(k).ap  = DS{k,5};       BOX(k).volFrom = VOL_FROM;
    BOX(k).pia = M{k}.pia;      BOX(k).piaMode = PIA_MODE;
end
save([OUT '.mat'], 'BOX', 'PCT_FRAC', 'P_EXPONENT');
fprintf('boxed fits -> %s.mat  (read by falloff_summary_scatter_260805)\n', OUT);

fprintf('\n%-14s %9s %8s %9s   per-round R2 at that fixed l\n', ...
        'dataset','l (box)','R2','l (full)');
for k = 1:n
    B = M{k};  L = B.L;  r2s = [];
    inFit = B.depth >= B.fitD(1) & B.depth <= B.fitD(2);
    for q = 1:numel(B.rounds)
        idx = B.rounds{q};  yq = B.diff_n(idx);
        m = idx(isfinite(yq) & yq > 0 & inFit(idx));
        if numel(m) < 3, r2s(end+1) = NaN; continue; end %#ok<AGROW>
        dq = B.depth(m);  v = B.diff_n(m);
        A  = exp(mean(log(v) + dq/L));
        rr = log(v) - (log(A) - dq/L);
        r2s(end+1) = 1 - sum(rr.^2)/max(sum((log(v)-mean(log(v))).^2), eps); %#ok<AGROW>
    end
    fprintf('%-14s %9.1f %8.3f %9.1f   %s\n', DS{k,3}, L, B.r2, B.ellFull, ...
            strjoin(compose('%.2f', r2s), '  '));
end
end

% =========================================================================
function B = box_metrics(dataDir, tag, mlR, apR, lut, volFrom, pctFrac, pExp, pixT, ...
                         fitOverride, piaMode, piaSmooth)
if volFrom >= 0
    vf = fullfile(dataDir, sprintf('falloff_proj_%s_vol.mat', tag));
else
    vf = fullfile(dataDir, sprintf('falloff_proj_%s_vol_from%d.mat', tag, volFrom));
end
assert(isfile(vf), 'no cache: %s', vf);
S  = load(vf);
Ff = load(fullfile(dataDir, sprintf('falloff_%s.mat', tag)));
B.fitD = Ff.FIT_DEPTH;  B.ellFull = Ff.Lp;
if nargin >= 10 && ~isempty(fitOverride)
    fprintf('%-12s FIT_DEPTH overridden %s -> %s\n', tag, mat2str(B.fitD), mat2str(fitOverride));
    B.fitD = fitOverride;
end

[nR, nC, ~] = size(S.V);
rI = sort(round(mlR/pixT)+1); rI = max(1,min(nR,rI));
cI = sort(round(apR/pixT)+1); cI = max(1,min(nC,cI));
Vol = S.V(rI(1):rI(2), cI(1):cI(2), :);
nz  = size(Vol,3);
nK  = max(1, round(pctFrac * numel(Vol(:,:,1))));

top = nan(nz,1); bot = nan(nz,1);
for z = 1:nz
    v = double(reshape(Vol(:,:,z), [], 1));
    top(z) = mean(maxk(v, nK));
    bot(z) = mean(mink(v, nK));
end
B.depth  = S.depth(:);
B.top    = top;  B.bot = bot;
B.diff   = top - bot;

% Re-zero on the MEASURED pia: the depth at which top-bot peaks. Above the
% surface there is no tissue to excite, so the signal rises into the pia and
% falls below it -- the maximum is the surface. Smoothed first, because a
% plateau (the shiverer holds maximum over ~15 um) makes the raw argmax depend
% on a single bright plane.
B.pia = 0;
if strcmpi(piaMode,'auto')
    ds = movmean(B.diff, max(1, round(piaSmooth/pixT)));
    [~, ip] = max(ds);
    B.pia = B.depth(ip);
    B.depth = B.depth - B.pia;
    fprintf('%-12s pia detected at filename z %+g um -> depth re-zeroed there\n', tag, B.pia);
end
mW       = laser_power_calibration(S.pctPlane, lut);
B.diff_n = B.diff ./ (mW(:).^pExp);

% Power rounds: contiguous runs of equal per-plane power. The assembled volume
% carries one plane per depth from its owning round, so a change in pctPlane IS
% a round boundary.
pp  = S.pctPlane(:);
bnd = [1; find(diff(pp) ~= 0) + 1; numel(pp)+1];
B.rounds   = arrayfun(@(q) (bnd(q):bnd(q+1)-1).', (1:numel(bnd)-1).', 'uni', 0);
B.roundPct = arrayfun(@(q) pp(bnd(q)), (1:numel(bnd)-1).');
% Label rounds in mW, not in the SI percentage. Percent is a Pockels setpoint,
% not a physical quantity, and mW is a non-linear function of it (~ pct^1.9 on
% this rig), so "48%" is not twice "21%" in any sense that matters. The
% normalisation already divides by mW^2, so the legend should be in the same
% units the maths uses.
B.roundMW = laser_power_calibration(B.roundPct, lut);

% Global exponential over FIT_DEPTH, on the power-normalised metric.
m = B.depth >= B.fitD(1) & B.depth <= B.fitD(2) & isfinite(B.diff_n) & B.diff_n > 0;
if nnz(m) >= 3
    pfit = polyfit(B.depth(m), log(B.diff_n(m)), 1);
    B.L  = -1/pfit(1);   B.A = exp(pfit(2));
    r    = log(B.diff_n(m)) - polyval(pfit, B.depth(m));
    B.r2 = 1 - sum(r.^2)/sum((log(B.diff_n(m)) - mean(log(B.diff_n(m)))).^2);
else
    B.L = NaN; B.A = NaN; B.r2 = NaN;
end

fprintf('%-12s box %dx%d px, %d planes %g..%g um, %d rounds (%s mW), l(box) %.1f (full field %.1f)\n', ...
    tag, diff(rI)+1, diff(cI)+1, nz, B.depth(1), B.depth(end), numel(B.rounds), ...
    strjoin(compose('%.1f', B.roundMW(:).'), ','), B.L, B.ellFull);
end
