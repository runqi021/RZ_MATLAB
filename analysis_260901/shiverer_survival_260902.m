% shiverer_survival_260902.m
% -----------------------------------------------------------------------
%  Kaplan-Meier survival curve for the shiverer colony, from birth.
%
%  Eight animals, every one of them observed to death, so there is NO CENSORING:
%  the Kaplan-Meier estimator reduces to the plain empirical survival function
%  and the curve necessarily reaches 0%%. If any animal is ever still alive or
%  lost to follow-up, add it to the table with a NaT death date and the estimator
%  below has to change -- a censored animal must leave the risk set without
%  causing a step down, which is exactly what a naive "fraction still alive"
%  calculation gets wrong.
%
%  The curve starts at 100%% on day 0 = the animal's own birthday, so the x axis
%  is AGE, not calendar time. The animals were born on five different dates, so a
%  calendar-time plot would be a different figure answering a different question.
%
%  Steps are drawn with `stairs` and a post-step convention: survival holds flat
%  until a death and drops on the day of that death. Deaths on the same day
%  produce one step of the combined size, which is why the day-119..133 tail
%  looks like single drops.
%
%  Runqi Zhang / 2026-09-02
% -----------------------------------------------------------------------
clear; clc;

%% ===================== USER-EDITABLE =====================
outDir = 'C:\fall-off\shiverer_survival_260902';

% birth date, death date -- one row per animal
REC = { '5/28/2026', '8/27/2026'
        '4/27/2026', '8/26/2026'
        '3/16/2026', '7/27/2026'
        '3/16/2026', '7/13/2026'
        '3/16/2026', '7/20/2026'
        '3/16/2026', '6/23/2026'
        '5/28/2026', '9/2/2026'
        '6/11/2026', '9/3/2026' };
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end

born = datetime(REC(:,1), 'InputFormat','M/d/yyyy');
died = datetime(REC(:,2), 'InputFormat','M/d/yyyy');
life = days(died - born);
assert(all(life > 0), 'a death date is not after its birth date');

n = numel(life);
[ls, ord] = sort(life);

% survival after each death: with no censoring this is just (n - k)/n
tt = [0; ls(:)];
ss = [1; (n - (1:n).')/n] * 100;

% median survival = the first day at which survival has dropped to <= 50%
medSurv = ls(find(ss(2:end) <= 50, 1));

fprintf('shiverer survival, n = %d, no censoring\n\n', n);
fprintf('%-4s %-12s %-12s %8s\n','#','born','died','days');
for k = 1:n
    j = ord(k);
    fprintf('%-4d %-12s %-12s %8d\n', k, datestr(born(j),'mm/dd/yyyy'), ...
            datestr(died(j),'mm/dd/yyyy'), life(j));
end
fprintf('\nmedian %d d   mean %.1f d   range %d-%d d   IQR %.0f-%.0f d\n', ...
        medSurv, mean(life), min(life), max(life), prctile(life,25), prctile(life,75));

%% ===================== FIGURE =====================
fig = figure('Color','w','Units','centimeters','Position',[3 3 15 11], ...
             'Name','shiverer survival');
ax = axes(fig,'Position',[0.135 0.145 0.815 0.72]);
hold(ax,'on'); box(ax,'on');

stairs(ax, tt, ss, '-', 'Color',[0.10 0.10 0.10], 'LineWidth',2.0);
% a marker on each death, at the survival it steps DOWN to
plot(ax, ls, ss(2:end), 'o', 'MarkerFaceColor',[0.10 0.10 0.10], ...
     'MarkerEdgeColor','w', 'MarkerSize',5.5, 'LineStyle','none');
% median guide
plot(ax, [0 medSurv], [50 50], ':', 'Color',[0.70 0.15 0.15], 'LineWidth',1.0);
plot(ax, [medSurv medSurv], [0 50], ':', 'Color',[0.70 0.15 0.15], 'LineWidth',1.0);
text(ax, medSurv+2, 53, sprintf('median %d d', medSurv), ...
     'Color',[0.70 0.15 0.15], 'FontSize',9);

xlim(ax, [0 max(ls)*1.06]);  ylim(ax, [0 103]);
set(ax,'YTick',0:20:100,'TickDir','out','Layer','top');
xlabel(ax,'age (days)');  ylabel(ax,'survival (%)');
title(ax, {sprintf('shiverer survival   |   n = %d, all deaths observed (no censoring)', n), ...
           sprintf('median %d d, mean %.1f d, range %d-%d d', ...
                   medSurv, mean(life), min(life), max(life))}, ...
      'FontSize',10,'FontWeight','normal');

stem = fullfile(outDir,'shiverer_survival');
exportgraphics(fig,[stem '.png'],'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');

T = table((1:n).', born(ord), died(ord), ls, ss(2:end), ...
    'VariableNames',{'rank','born','died','lifespan_days','survival_pct_after'});
writetable(T, [stem '.csv']);
fprintf('\nwrote %s.png / .pdf / .csv\n', stem);
