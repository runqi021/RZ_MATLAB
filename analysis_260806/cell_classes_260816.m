function [T, INFO] = cell_classes_260816(sertMode)
%CELL_CLASSES_260816  Hand assignment of cells to respiratory classes.
%
%   T = cell_classes_260816()            % default Sert rule: 'included'
%   T = cell_classes_260816('all')       % see SERT RULE below
%
% Returns a table: cell, group, date, class, source.
%
% ASSIGNED BY RUNQI, 2026-08-16. Every id below was checked against the registry:
% the genotype filed here matches the genotype the registry reports, and no cell
% appears in two classes.
%
% SERT RULE -- this one is not a fixed list, it is a criterion, and the three
% readings differ by a lot, so it is explicit:
%   'included' (default)  cells INCLUDED in the popsel GUI AND p<0.01   -> 10 cells
%   'notexcluded'         all p<0.01 except those explicitly EXCLUDED   -> 17 cells
%   'all'                 every Sert cell with p<0.01, curation ignored -> 22 cells
% The gap is the 7 significant cells left UNDECIDED in the GUI (242 244 248 252
% 255 257 259) -- not rejected, just not reached. 'included' is the strict reading
% and matches the instruction to keep the included set for the rest of the
% analysis; the others are here so the choice can be changed in one place rather
% than re-derived by hand.
%
% p is the circular-shift permutation p on the PSTH (1200 shuffles), "either
% trigger" = min(pOnset, pPeak) < 0.01. NOTE the p floor at 1200 shuffles is
% 1/1201 = 0.00083, so p<0.01 is a real threshold here, unlike p<0.001.

if nargin < 1 || isempty(sertMode), sertMode = 'included'; end

sumRoot = 'D:\Ventral_surface_summary';

%% ---- explicit assignments ----
% Vglut2 I edits, 2026-08-16 (RZ):
%   -26  removed. Weakest of the group by a distance: logZ 1.09 against 2.1-5.1
%        for the rest, IQR 129 deg against 30-87 -- barely modulated and barely
%        phase-concentrated.
%   +49  added. Vglut2/0810, 71 events, logZ 3.34, permutation p at the 1/1201
%        floor on BOTH triggers -- well inside the group's range.
% NOTE this makes Vglut2 I a two-session group (0224 at IBI 0.47 s plus 0810 at
% 1.47 s), where it was 0224-only. The phase axis is cycle-normalised so that is
% legitimate, but the group is no longer one animal.
% 176 moved pre-I -> I, 2026-08-17 (RZ). Its median phase is 126 deg, the
% shallowest of the four pre-I cells and closest to the I group's 180-202 deg;
% logZ 3.89, 136 events, so it is a well-modulated cell either way. This makes
% Vglut2 pre-I n=3 and Vglut2 I n=7.
A = { 'pre-I',                  [178 54 48]
      % +53 added 2026-08-17 (RZ). Vglut2/0810 roi1_5x_z-51, 78 events, and its
      % permutation p sits at the 1/1201 floor on BOTH triggers -- it was simply
      % never assigned to anything, not deliberately left out.
      % -30, -34 removed 2026-08-17 (RZ). Both were at 180 deg like the rest;
      % they are now in NO class, so they drop out of every figure entirely
      % rather than moving somewhere else.
      'I',                      [33 32 25 49 176 53]
      'post-I',                 [20 18]
      'null',                   [224 200 223 225 222 214 226 227]
      % Vgat I edit, 2026-08-17 (RZ): -192 removed. Reason not recorded here --
      % fill it in. Note it was the only member of this group from a SURFACE
      % recording (Vgat/0730 cell1, roi1_2x_x1100y950_z13, 81 events); the three
      % that remain are all from the deep z250/z265 FOVs.
      'I',                      [277 279 278]
      'tonic pre-I suppressed', [185 190 195]
      'null',                   [186 194 193]
      % Vglut2 tonic/rhythmic null, added 2026-08-17 (RZ) from the per-cell
      % summary active folder. All five verified Vglut2 against the registry:
      % 24/26 from 0224, 35/36/37 from 0728.
      % 26: dropped from Vglut2 'I' on 08-16 (weakest of that group, logZ 1.09),
      % then in NO class at all, which meant it was silently absent from the map.
      % Filed as null 2026-08-17 so it is represented rather than missing.
      'null',                   [24 26 35 36 37]
      % ChAT tonic/rhythmic null, added 2026-08-17 (RZ). ChAT/0523, 219 events --
      % the most active cell in any null group here, so 'null' is a statement
      % about phase locking, not about firing rate.
      'null',                   [21]
      % Sert tonic/rhythmic null, added 2026-08-17 (RZ). All five Sert/0721 and
      % all p > 0.01, so none of them collides with the post-I criterion set
      % below -- these are the included-but-not-significant Sert, given a class
      % of their own rather than left out.
      % NOTE 250 is the odd one: it is UNDECIDED in the popsel GUI, not included
      % (the other four are 'include'). Classing it here overrides that, which is
      % fine as a deliberate call but is not something the curation supports.
      'null',                   [267 271 239 250 235] };

cells = []; class = strings(0,1); src = strings(0,1);
for k = 1:size(A,1)
    cells = [cells, A{k,2}];                                   %#ok<AGROW>
    class = [class; repmat(string(A{k,1}), numel(A{k,2}), 1)]; %#ok<AGROW>
    src   = [src;   repmat("explicit",     numel(A{k,2}), 1)]; %#ok<AGROW>
end

%% ---- Sert post-I, by criterion ----
L = load(fullfile(sumRoot,'pop_analysis_260813','pop_features.mat'),'T');
S = L.T(strcmp(L.T.group,'Sert'),:);
sig = S.cell( (S.pOnset < 0.01) | (S.pPeak < 0.01) );

decF = fullfile(sumRoot,'popsel_260816','popsel_decisions_Sert.csv');
inc = []; exc = [];
if isfile(decF)
    Td = readtable(decF,'TextType','string');
    inc = Td.cell(Td.decision == "include");
    exc = Td.cell(Td.decision == "exclude");
end

switch lower(sertMode)
    case 'included',    sertCells = intersect(sig, inc);
    case 'notexcluded', sertCells = setdiff(sig, exc);
    case 'all',         sertCells = sig;
    otherwise, error('sertMode must be included | notexcluded | all');
end
sertCells = sort(sertCells(:)).';

cells = [cells, sertCells];
class = [class; repmat("post-I", numel(sertCells), 1)];
src   = [src;   repmat("Sert p<0.01 (" + string(sertMode) + ")", numel(sertCells), 1)];

%% ---- attach genotype/date from the registry, and check ----
addpath(fileparts(mfilename('fullpath')));
D = load(fullfile(sumRoot,'event_latency_260811','event_latency_data.mat'),'CELL','OBS','REC');
obsOf = pooled_obs_260814(D.CELL, D.OBS);
grp = strings(numel(cells),1); dat = strings(numel(cells),1);
for k = 1:numel(cells)
    c = cells(k);
    assert(c <= numel(obsOf) && ~isempty(obsOf{c}), 'cell %d is not live in the registry', c);
    p = regexp(D.OBS(obsOf{c}(1)).label,'/','split');
    grp(k) = string(p{1});  dat(k) = string(p{2});
end
assert(numel(unique(cells)) == numel(cells), 'a cell was assigned to two classes');

T = table(cells(:), grp, dat, class, src, ...
          'VariableNames',{'cell','group','date','class','source'});
T = sortrows(T, {'group','class','cell'});

INFO = struct('sertMode',sertMode,'sertSig',sig(:).','sertIncluded',inc(:).', ...
              'sertExcluded',exc(:).','sertChosen',sertCells);
end
