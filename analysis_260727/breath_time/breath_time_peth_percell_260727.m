function breath_time_peth_percell_260727(triggerIn)
%% breath_time_peth_percell_260727  One QC panel per cell for the time-domain PETH.
% -----------------------------------------------------------------------
% Every number in the summary table is shown here in the form it was derived from,
% so nothing has to be taken on trust.
%
%   1  observed PETH with the SHUFFLE ENVELOPE (mean +/- 2 SD of the null) drawn
%      behind it. If the observed curve stays inside the grey band, the cell is not
%      modulated -- and that is visible directly, without reading a p-value. The
%      preferred latency is marked.
%   2  the null distribution of T_exc with the observed value marked. At these
%      event counts this null is discrete and right-skewed, which is exactly why
%      significance is called on the rank-based p and not on the z.
%   3  event-latency histogram: every event assigned to its immediately preceding
%      accepted onset. This is where precision and latency MAD come from, and it
%      is conditioned on the cell firing, unlike panel 1.
%   4  the null distribution of the precision statistic, with the preferred latency
%      re-estimated inside every shuffle (otherwise precision is circular).
%   5  the numbers.
%
% Input : breath_time_peth_data.mat
% Output: <dataset>\analysis_260727\breath_time\percell\cell_###.png
%
% Runqi Zhang / 2026-07-27
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
which_cells = 'tested';   % 'sig' | 'tested' | 'all'
showFirst   = 6;
% =====================================================================

TRIGGER = 'onset';   % 'onset' | 'peak' -- must match the run you want to plot
if nargin >= 1 && ~isempty(triggerIn), TRIGGER = char(triggerIn); end
outDir = fullfile(cfg.outRoot,'breath_time',lower(TRIGGER));
dataFile = fullfile(outDir,'breath_time_peth_data.mat');
assert(isfile(dataFile), ['breath_time_peth_data.mat not found:\n  %s\n' ...
    'Run breath_time_peth_260727.m first.'], dataFile);
S = load(dataFile,'R','ctrs','sig_exc','q_exc','sig_pre','q_pre','params');
R = S.R; ctrs = S.ctrs; params = S.params;

switch lower(which_cells)
    case 'sig', sel = find(S.sig_exc);
    case 'all', sel = find([R.n_events] > 0);
    otherwise,  sel = find([R.tested]);
end
assert(~isempty(sel), 'No cells selected (which_cells = %s).', which_cells);
[~,o] = sort([R(sel).mod_exc_z],'descend','MissingPlacement','last'); sel = sel(o);

sub = fullfile(outDir,'percell');
if ~isfolder(sub), mkdir(sub);
else, old = dir(fullfile(sub,'cell_*.png')); for k=1:numel(old), delete(fullfile(old(k).folder,old(k).name)); end
end

col = cfg.genotype_color;
hOff = figure('Color','w','Visible','off');
fprintf('rendering %d cells (%s, ranked by M_exc)...\n', numel(sel), which_cells);

for si = 1:numel(sel)
    c = sel(si); C = R(c);
    if si <= showFirst
        hf = figure('Color','w','Visible','on');
    else
        hf = hOff; clf(hf); set(0,'CurrentFigure',hf);   % clf() does NOT make hf current;
    end                                                  % without this subplot() draws into
                                                         % the last VISIBLE figure and every
                                                         % export after the 6th comes out blank
    set(hf,'Units','centimeters','Position',[1 1 36 10.5]);   % tall enough for the text block

    % --- 1. PETH vs the shuffle envelope ---
    a1 = subplot(1,5,1); hold(a1,'on'); box(a1,'on'); grid(a1,'on');
    if ~isempty(C.peth_shuffle_mean)
        mu = C.peth_shuffle_mean; sd = C.peth_shuffle_sd;
        fill(a1, [ctrs fliplr(ctrs)], [mu+2*sd fliplr(max(mu-2*sd,0))], [.85 .85 .85], ...
             'EdgeColor','none');
        plot(a1, ctrs, mu, '-','Color',[.5 .5 .5],'LineWidth',1);
    end
    stairs(a1, ctrs, C.peth, '-','Color',col,'LineWidth',1.4);
    yline(a1, C.mean_peth_rate, 'k:','LineWidth',1);
    xline(a1, 0, 'k-','LineWidth',1);
    if isfinite(C.preferred_latency_s), xline(a1, C.preferred_latency_s, '--','Color',[.85 .1 .1],'LineWidth',1.5); end
    xlim(a1,[ctrs(1) ctrs(end)]);
    xlabel(a1,'time from inspiration onset (s)'); ylabel(a1,'events / s of exposure');
    title(a1,{'PETH vs shuffle envelope (mean \pm2SD)','red = preferred latency'},'FontSize',7);

    % --- 2. null of T_exc ---
    a2 = subplot(1,5,2); hold(a2,'on'); box(a2,'on'); grid(a2,'on');
    nl = double(C.null_exc);
    if ~isempty(nl)
        histogram(a2, nl, 30, 'FaceColor',[.72 .72 .72],'EdgeColor','none');
        yl = ylim(a2); plot(a2,[C.T_exc C.T_exc], yl, '-','Color',col,'LineWidth',2);
        plot(a2,[prctile(nl,95) prctile(nl,95)], yl,'k--','LineWidth',1);
    end
    xlabel(a2,'T_{exc} under the null'); ylabel(a2,'count');
    title(a2,sprintf('excitation null (%d shuffles)\ncoloured = observed', params.nShuffle),'FontSize',7);

    % --- 3. event latencies ---
    a3 = subplot(1,5,3); hold(a3,'on'); box(a3,'on'); grid(a3,'on');
    dt = double(C.dt_obs);
    if ~isempty(dt)
        h = params.medCycle_s/2;                 % nearest-trigger latency is bounded by
        eb = -h:params.binWidth_s:h;              % half the inter-trigger interval
        histogram(a3, dt, eb, 'FaceColor',[.20 .45 .80],'EdgeColor','none');
        yl = ylim(a3);
        tp = C.pref_latency_nearest_s;            % preferred latency, wrapped into +/-half cycle
        if isfinite(tp)
            fill(a3, tp + params.relWin_s*[-1 1 1 -1], [0 0 yl(2) yl(2)], ...
                 [.85 .1 .1], 'FaceAlpha',0.12, 'EdgeColor','none');
            plot(a3, [tp tp], yl, '--', 'Color',[.85 .1 .1],'LineWidth',1.5);
        end
        xline(a3, 0, 'k-','LineWidth',1);
        plot(a3, [C.latency_median_s C.latency_median_s], yl, '-','Color',[.1 .1 .5],'LineWidth',1.2);
        xlim(a3,[min(eb) max(eb)]); ylim(a3, yl);
        ylabel(a3,'# events');
    end
    xlabel(a3,'signed latency to NEAREST trigger (s)');
    title(a3,sprintf('median %.0f ms, MAD %.0f ms\nshaded = \\pm%.0f ms precision window', ...
        1000*C.latency_median_s, 1000*C.latency_mad_s, 1000*params.relWin_s),'FontSize',7);

    % --- 4. null of precision ---
    a4 = subplot(1,5,4); hold(a4,'on'); box(a4,'on'); grid(a4,'on');
    np = double(C.null_pre);
    if ~isempty(np)
        histogram(a4, np, 30, 'FaceColor',[.72 .72 .72],'EdgeColor','none');
        yl = ylim(a4); plot(a4,[C.precision_fraction C.precision_fraction], yl,'-','Color',col,'LineWidth',2);
    end
    xlabel(a4,'precision fraction under the null'); ylabel(a4,'count');
    title(a4,{'precision null','(preferred latency re-estimated each shuffle)'},'FontSize',7);

    % --- 5. numbers ---
    a5 = subplot(1,5,5); axis(a5,'off');
    txt = {
        sprintf('\\bfcell %d\\rm   %d recording(s)', C.cell_id, C.n_obs)
        ''
        sprintf('events          %d', C.n_events)
        sprintf('duration        %.0f s', C.duration_s)
        sprintf('rate            %.2f /min', C.rate_per_min)
        sprintf('accepted breaths %d', C.n_accepted_breaths)
        ''
        sprintf('mean rate       %.3f ev/s', C.mean_peth_rate)
        sprintf('peak            %.3f ev/s', C.peak_rate)
        sprintf('peak/mean       %.2f', C.peak_over_mean)
        sprintf('pref latency    %.0f ms', 1000*C.preferred_latency_s)
        ''
        sprintf('\\bfM_exc          %.2f\\rm', C.mod_exc_z)
        sprintf('p_exc           %.4f', C.p_exc)
        sprintf('q_exc           %.4f', S.q_exc(c))
        sprintf('%s', ternary(S.sig_exc(c),'\bfEXCITED\rm','not sig (exc)'))
        ''
        sprintf('M_sup           %s', numstr(C.mod_sup_z))
        sprintf('p_sup           %s', numstr(C.p_sup))
        sprintf('sup tested      %d  (exp/win %.1f)', C.sup_tested, C.sup_expected_per_window)
        ''
        sprintf('recruitment     %.3f  (%d breaths)', C.recruitment, C.n_breaths_recruited)
        sprintf('precision       %.3f', C.precision_fraction)
        sprintf('precision_z     %.2f', C.precision_z)
        sprintf('q_precision     %.4f', S.q_pre(c))
        };
    text(a5, 0, 1, txt, 'Units','normalized','VerticalAlignment','top', ...
         'FontName','FixedWidth','FontSize',6.2,'Interpreter','tex');

    sgtitle(sprintf('cell %d  |  %s', C.cell_id, strjoin(cellstr(C.rec_names),'  ')), ...
            'Interpreter','none','FontSize',8);
    exportgraphics(hf, fullfile(sub, sprintf('cell_%03d_Mexc%05.1f.png', C.cell_id, ...
        max(C.mod_exc_z,0))), 'Resolution',150);
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d per-cell panels to %s\n', numel(sel), sub);
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function s = numstr(v)
if isnan(v), s = 'n/a'; else, s = sprintf('%.4f', v); end
end
