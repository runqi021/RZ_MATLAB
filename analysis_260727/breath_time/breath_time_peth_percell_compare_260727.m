function breath_time_peth_percell_compare_260727()
%% breath_time_peth_percell_compare_260727  One figure per cell, BOTH alignments.
% -----------------------------------------------------------------------
% Replaces flipping between breath_time\onset\percell\ and \peak\percell\.
% Each cell gets ONE figure, 2 rows x 5 panels:
%
%     row 1  ONSET-triggered   (zero = inspiration onset, foot_idx)
%     row 2  PEAK-triggered    (zero = inspiratory peak,  peak_idx)
%
% and the columns are the same quantity in both rows, so a column is a direct
% comparison:
%
%   1  PETH with the shuffle envelope (mean +/- 2 SD)
%   2  null distribution of T_exc, observed marked
%   3  signed event latency to the NEAREST trigger
%   4  null distribution of the precision statistic
%   5  the numbers
%
% The title states which alignment is TIGHTER for that cell, by latency MAD --
% which is the landmark that cell is actually locked to. MAD, not M_exc: M_exc
% mixes dispersion with firing rate and with how much evidence the cell provided,
% so it barely moves between alignments.
%
% Axes are created with an explicit 'Parent' rather than subplot(), because
% clf(hf) does NOT make an off-screen figure current and subplot() would then draw
% into the last visible window -- which silently produced 28 blank exports before.
%
% Input : breath_time\{onset,peak}\breath_time_peth_data.mat
% Output: breath_time\percell_compare\cell_###.png
%
% Runqi Zhang / 2026-07-28
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
which_cells = 'tested';   % 'sig' | 'tested' | 'all'
rankBy      = 'onset';    % which alignment's M_exc orders the output files
showFirst   = 4;
% =====================================================================

base = fullfile(cfg.outRoot,'breath_time');
A = load_side(base,'onset');  B = load_side(base,'peak');
nC = min(numel(A.R), numel(B.R));

switch lower(which_cells)
    case 'sig', sel = find(A.sig_exc(1:nC) | B.sig_exc(1:nC));
    case 'all', sel = find([A.R(1:nC).n_events] > 0);
    otherwise,  sel = find([A.R(1:nC).tested] & [B.R(1:nC).tested]);
end
assert(~isempty(sel), 'No cells selected (which_cells = %s).', which_cells);
key = ternary(strcmpi(rankBy,'peak'), [B.R(sel).mod_exc_z], [A.R(sel).mod_exc_z]);
[~,o] = sort(key,'descend','MissingPlacement','last');  sel = sel(o);

sub = fullfile(base,'percell_compare');
if ~isfolder(sub), mkdir(sub);
else, old = dir(fullfile(sub,'cell_*.png')); for k=1:numel(old), delete(fullfile(old(k).folder,old(k).name)); end
end

hOff = figure('Color','w','Visible','off');
fprintf('rendering %d cells, both alignments...\n', numel(sel));

for si = 1:numel(sel)
    c = sel(si);
    if si <= showFirst, hf = figure('Color','w','Visible','on'); else, hf = hOff; clf(hf); end
    set(hf,'Units','centimeters','Position',[1 1 36 19]);

    for row = 1:2
        if row == 1, S = A; tg = 'ONSET'; else, S = B; tg = 'PEAK'; end
        yb = 0.56 - 0.50*(row-1);                    % row base, top row higher
        draw_row(hf, yb, S, c, tg, cfg);
    end

    mA = A.R(c).latency_mad_s;  mB = B.R(c).latency_mad_s;
    if isfinite(mA) && isfinite(mB)
        if abs(mA-mB) < 1e-9
            verdict = sprintf('TIED on MAD (%.0f ms both)', 1000*mA);
        elseif mA < mB
            verdict = sprintf('tighter ONSET-aligned (MAD %.0f vs %.0f ms)', 1000*mA, 1000*mB);
        else
            verdict = sprintf('tighter PEAK-aligned (MAD %.0f vs %.0f ms)',  1000*mB, 1000*mA);
        end
    else
        verdict = 'MAD not available in one alignment';
    end
    sgtitle(sprintf('cell %d   |   %s   |   %s', A.R(c).cell_id, ...
            strjoin(cellstr(A.R(c).rec_names),'  '), verdict), ...
            'Interpreter','none','FontSize',9);

    rank = max(ternary(strcmpi(rankBy,'peak'), B.R(c).mod_exc_z, A.R(c).mod_exc_z), 0);
    exportgraphics(hf, fullfile(sub, sprintf('cell_%03d_Mexc%05.1f.png', ...
        A.R(c).cell_id, rank)), 'Resolution',140);
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d combined panels to %s\n', numel(sel), sub);
end

%% ========================= helpers =========================
function S = load_side(base, tg)
f = fullfile(base, tg, 'breath_time_peth_data.mat');
assert(isfile(f), ['Missing %s\nRun breath_time_peth_260727(''%s'') first.'], f, tg);
S = load(f,'R','ctrs','sig_exc','q_exc','sig_pre','q_pre','params');
end

function draw_row(hf, yb, S, c, tg, cfg)
C = S.R(c); ctrs = S.ctrs; params = S.params; col = cfg.genotype_color;
h = 0.33;  w = 0.145;  x0 = 0.055;  gap = 0.037;
ax = @(k) axes('Parent',hf,'Position',[x0+(k-1)*(w+gap), yb, w, h]); %#ok<LAXES>

% --- 1. PETH + shuffle envelope ---
a1 = ax(1); hold(a1,'on'); box(a1,'on'); grid(a1,'on');
if ~isempty(C.peth_shuffle_mean)
    mu = C.peth_shuffle_mean; sd = C.peth_shuffle_sd;
    fill(a1,[ctrs fliplr(ctrs)],[mu+2*sd fliplr(max(mu-2*sd,0))],[.85 .85 .85],'EdgeColor','none');
    plot(a1, ctrs, mu, '-','Color',[.5 .5 .5],'LineWidth',1);
end
if ~isempty(C.peth), stairs(a1, ctrs, C.peth, '-','Color',col,'LineWidth',1.4); end
yline(a1, C.mean_peth_rate,'k:'); xline(a1,0,'k-','LineWidth',1);
if isfinite(C.preferred_latency_s), xline(a1, C.preferred_latency_s,'--','Color',[.85 .1 .1],'LineWidth',1.5); end
xlim(a1,[ctrs(1) ctrs(end)]);
ylabel(a1,'ev / s exposure');
xlabel(a1, sprintf('time from %s (s)', lower(tg)));
title(a1, sprintf('\\bf%s\\rm  PETH vs null envelope', tg), 'FontSize',7.5);

% --- 2. T_exc null ---
a2 = ax(2); hold(a2,'on'); box(a2,'on'); grid(a2,'on');
nl = double(C.null_exc);
if ~isempty(nl)
    histogram(a2, nl, 30, 'FaceColor',[.72 .72 .72],'EdgeColor','none');
    yl = ylim(a2);
    plot(a2,[C.T_exc C.T_exc], yl, '-','Color',col,'LineWidth',2);
    plot(a2,[prctile(nl,95) prctile(nl,95)], yl,'k--','LineWidth',1);
end
xlabel(a2,'T_{exc} under null'); ylabel(a2,'count');
title(a2, sprintf('M_{exc} = %.2f,  p = %.4f', C.mod_exc_z, C.p_exc),'FontSize',7.5);

% --- 3. signed nearest-trigger latency ---
a3 = ax(3); hold(a3,'on'); box(a3,'on'); grid(a3,'on');
dt = double(C.dt_obs);
if ~isempty(dt)
    hc = params.medCycle_s/2;  eb = -hc:params.binWidth_s:hc;
    histogram(a3, dt, eb, 'FaceColor',[.20 .45 .80],'EdgeColor','none');
    yl = ylim(a3);
    tp = C.pref_latency_nearest_s;
    if isfinite(tp)
        fill(a3, tp + params.relWin_s*[-1 1 1 -1], [0 0 yl(2) yl(2)], [.85 .1 .1], ...
             'FaceAlpha',0.12,'EdgeColor','none');
        plot(a3,[tp tp], yl,'--','Color',[.85 .1 .1],'LineWidth',1.5);
    end
    plot(a3,[C.latency_median_s C.latency_median_s], yl,'-','Color',[.1 .1 .5],'LineWidth',1.2);
    xline(a3,0,'k-','LineWidth',1);
    xlim(a3,[min(eb) max(eb)]); ylim(a3,yl);
end
xlabel(a3,'signed latency to nearest (s)'); ylabel(a3,'# events');
title(a3, sprintf('median %.0f ms,  MAD %.0f ms', 1000*C.latency_median_s, 1000*C.latency_mad_s),'FontSize',7.5);

% --- 4. precision null ---
a4 = ax(4); hold(a4,'on'); box(a4,'on'); grid(a4,'on');
np = double(C.null_pre);
if ~isempty(np)
    histogram(a4, np, 30, 'FaceColor',[.72 .72 .72],'EdgeColor','none');
    yl = ylim(a4);
    plot(a4,[C.precision_fraction C.precision_fraction], yl,'-','Color',col,'LineWidth',2);
end
xlabel(a4,'precision under null'); ylabel(a4,'count');
title(a4, sprintf('precision %.3f,  z = %.2f', C.precision_fraction, C.precision_z),'FontSize',7.5);

% --- 5. numbers ---
a5 = ax(5); axis(a5,'off');
txt = {
    sprintf('\\bf%s-triggered\\rm', tg)
    sprintf('events        %d', C.n_events)
    sprintf('breaths       %d', C.n_accepted_breaths)
    sprintf('mean rate     %.3f ev/s', C.mean_peth_rate)
    sprintf('peak          %.3f ev/s', C.peak_rate)
    sprintf('peak/mean     %.2f', C.peak_over_mean)
    ''
    sprintf('pref latency  %.0f ms', 1000*C.preferred_latency_s)
    sprintf('  (nearest)   %.0f ms', 1000*C.pref_latency_nearest_s)
    ''
    sprintf('\\bfM_exc        %.2f\\rm', C.mod_exc_z)
    sprintf('p_exc         %.4f', C.p_exc)
    sprintf('q_exc         %.4f', S.q_exc(c))
    sprintf('%s', ternary(S.sig_exc(c),'\bfEXCITED\rm','not sig'))
    ''
    sprintf('recruitment   %.3f', C.recruitment)
    sprintf('precision     %.3f', C.precision_fraction)
    sprintf('latency MAD   %.0f ms', 1000*C.latency_mad_s)
    sprintf('M_sup         %s', numstr(C.mod_sup_z))
    };
text(a5, 0, 1, txt, 'Units','normalized','VerticalAlignment','top', ...
     'FontName','FixedWidth','FontSize',6.4,'Interpreter','tex');
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

function s = numstr(v)
if isnan(v), s = 'n/a'; else, s = sprintf('%.3f', v); end
end
