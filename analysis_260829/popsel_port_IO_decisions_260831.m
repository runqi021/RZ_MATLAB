% popsel_port_IO_decisions_260831.m
% -----------------------------------------------------------------------
%  Carry the 16-Aug IO curation (107 include / 20 exclude / 13 undecided) onto
%  the CURRENT cell numbering, so it does not have to be done again.
%
%  THE PROBLEM. popsel_decisions_IO.csv is keyed on the cell NUMBER, and the
%  28-Aug registry rebuild renumbered nearly every cell (the number is a rank in
%  dir() scan order, not an identity -- see feedback_cell_ids_not_stable). Read
%  as-is, only 116 of its 140 rows still land on an IO cell at all, and those
%  that do are mostly the WRONG IO cell.
%
%  THE BRIDGE. popsel_cache_IO.mat, written in the same session as the CSV,
%  stores each cell's OBS label -- '<genotype>/<date>/<recording>/<roi>' -- which
%  is an identity, not a position. So:
%       old id --(old cache)--> label --(current registry)--> new id
%
%  IT IS VERIFIED, NOT TRUSTED. A silent mis-map would attach a decision to the
%  wrong neuron, which is worse than redoing the curation. Three checks, and the
%  script refuses to write if any fails:
%    1. labels must be unique on both sides (a duplicate label = ambiguous join)
%    2. every ported row's nSpikes must match the value the CSV recorded
%    3. every ported row's IBI must match to 1e-6
%  nSpikes and IBI were computed independently of the numbering, so agreement on
%  both is strong evidence the join found the same neuron.
%
%  Unmatched cells stay 'undecided' rather than being guessed at, and new cells
%  that did not exist in August are added as 'undecided'.
%
%  Runqi Zhang / 2026-08-31
% -----------------------------------------------------------------------

clear; clc;

%% ===================== USER-EDITABLE =====================
sumRoot  = 'D:\Ventral_surface_summary';
GROUP    = 'IO';
oldDir   = fullfile(sumRoot,'popsel_260816','_stale_pre260828');
newDir   = fullfile(sumRoot,'popsel_260816');
dryRun   = false;      % true = report the mapping, write nothing
%% =========================================================

oldCache = fullfile(oldDir, sprintf('popsel_cache_%s.mat', GROUP));
oldCsv   = fullfile(oldDir, sprintf('popsel_decisions_%s.csv', GROUP));
newCache = fullfile(newDir, sprintf('popsel_cache_%s.mat', GROUP));
newCsv   = fullfile(newDir, sprintf('popsel_decisions_%s.csv', GROUP));
for f = {oldCache, oldCsv, newCache}
    assert(isfile(f{1}), 'missing: %s', f{1});
end

O = load(oldCache,'C');  Cold = O.C;
N = load(newCache,'C');  Cnew = N.C;
T = readtable(oldCsv,'TextType','string');

fprintf('=========== popsel_port_%s_decisions_260831 ===========\n', GROUP);
fprintf('old cache : %d cells   old csv: %d rows\n', numel(Cold), height(T));
fprintf('new cache : %d cells\n', numel(Cnew));

labOld = string({Cold.label}.');
labNew = string({Cnew.label}.');
assert(numel(unique(labOld)) == numel(labOld), ...
    'old labels are not unique -- the join would be ambiguous');
assert(numel(unique(labNew)) == numel(labNew), ...
    'new labels are not unique -- the join would be ambiguous');

% old cell number -> label
oldIdOf = [Cold.cell].';
% csv row -> label
[tfR, locR] = ismember(T.cell, oldIdOf);
if any(~tfR)
    fprintf(2,'%d csv row(s) are not in the old cache and cannot be ported\n', nnz(~tfR));
end
rowLab = strings(height(T),1);
rowLab(tfR) = labOld(locR(tfR));

% label -> new cell number
[tfN, locN] = ismember(rowLab, labNew);
ok = tfR & tfN & rowLab ~= "";
fprintf('\nmapped %d of %d rows by label\n', nnz(ok), height(T));

%% ---- verification against quantities that do NOT depend on the numbering ----
bad = {};
for i = find(ok).'
    cn = Cnew(locN(i));
    if double(cn.nSpikes) ~= double(T.nSpikes(i))
        bad{end+1} = sprintf('%s: nSpikes %g -> %g', rowLab(i), T.nSpikes(i), cn.nSpikes); %#ok<AGROW>
    end
    if abs(double(cn.IBI) - double(T.IBI(i))) > 1e-6
        bad{end+1} = sprintf('%s: IBI %.6f -> %.6f', rowLab(i), T.IBI(i), cn.IBI); %#ok<AGROW>
    end
end
if ~isempty(bad)
    fprintf(2,'VERIFICATION FAILED on %d field(s):\n', numel(bad));
    fprintf(2,'   %s\n', bad{:});
    error('popsel_port:verify','refusing to write a mapping that does not check out');
end
fprintf('verified  : nSpikes and IBI agree for all %d ported rows\n', nnz(ok));

%% ---- how far did the numbers actually move ----
movedN = nnz(T.cell(ok) ~= [Cnew(locN(ok)).cell].');
fprintf('id drift  : %d of %d ported cells changed number\n', movedN, nnz(ok));

%% ---- build the new decisions table ----
dec = repmat("undecided", numel(Cnew), 1);
for i = find(ok).'
    dec(locN(i)) = T.decision(i);
end
nAdd = numel(Cnew) - nnz(ok);

Tn = table([Cnew.cell].', string({Cnew.group}.'), string({Cnew.date}.'), dec, ...
           [Cnew.nSpikes].', [Cnew.IBI].', [Cnew.logZ].', ...
    'VariableNames', {'cell','group','date','decision','nSpikes','IBI','logZ'});

fprintf('\nported    : %d include / %d exclude / %d undecided\n', ...
    nnz(dec=="include"), nnz(dec=="exclude"), nnz(dec=="undecided"));
fprintf('            (%d of those undecided are cells with no 16-Aug decision)\n', nAdd);

if dryRun
    fprintf('\nDRY RUN -- nothing written.\n');
    return;
end
if isfile(newCsv)
    bk = [newCsv '.bak_' datestr(now,'yymmdd_HHMMSS')];
    copyfile(newCsv, bk);
    fprintf('existing csv backed up -> %s\n', bk);
end
writetable(Tn, newCsv);
fprintf('wrote %s\n', newCsv);
fprintf('\nopen it with:  popsel_run_260831   (OPEN = ''%s'')\n', GROUP);
