% per_cell_summary_ibi_split_260814.m
% -----------------------------------------------------------------------
%  Split the rendered per-cell summaries into breathing-rate subfolders:
%
%      per-cell-summary_active_260812\IBI_le1\      IBI <= 1 s        (fast)
%                                    \IBI_1to2\     1 < IBI <= 2 s    (middle)
%                                    \IBI_2to3\     2 < IBI <= 3 s    (slow)
%
%  IBI is the cell's POOLED MEDIAN inter-onset interval -- the number printed on
%  each figure. Median, not mean: the cycle distribution has a tail out to 32 s
%  from missed onsets, and a mean would follow it.
%
%  COPIES the existing renders rather than re-rendering. Nothing that is DRAWN has
%  changed since that batch ran (the later edits added stats outputs and split the
%  histogram into compute/draw halves without altering a number), so a re-render
%  would produce byte-identical panels and cost ~30 min. Copying also guarantees
%  the three folders and the parent hold the same figure for the same cell.
%
%  The bins are MUTUALLY EXCLUSIVE: every cell lands in exactly one folder and no
%  figure is duplicated, so the three folders partition the population and can be
%  compared against each other directly. The script asserts that the three counts
%  sum to the number of cells.
%
%  Runqi Zhang / 2026-08-14
% -----------------------------------------------------------------------

clear; clc;
rootPath = 'D:\Ventral_surface_summary';
srcDir   = fullfile(rootPath,'per-cell-summary_active_260812');
L = load(fullfile(rootPath,'pop_analysis_260813','pop_features.mat'),'T');
T = L.T;
% lo < IBI <= hi, so the three bins are disjoint and cover everything
BINS = { 'IBI_le1',  -inf, 1
         'IBI_1to2',    1, 2
         'IBI_2to3',    2, 3 };

assert(isfolder(srcDir), 'source folder not found: %s', srcDir);
srcPng = dir(fullfile(srcDir,'*.png'));
srcPng = srcPng(~contains({srcPng.name},'_avgproj'));
fprintf('source: %d figures in %s\n\n', numel(srcPng), srcDir);

nAssigned = 0;
for bi = 1:size(BINS,1)
    lo = BINS{bi,2}; hi = BINS{bi,3};
    sel = find(T.IBI > lo & T.IBI <= hi);
    nAssigned = nAssigned + numel(sel);
    dst = fullfile(srcDir, BINS{bi,1});
    % Clear the destination first. These folders are pure copies, rebuilt from
    % srcDir every time, so nothing unique lives here. Without this, a figure left
    % over from an earlier run of the parent batch survives alongside the new one
    % -- and after the 260814 cell merges seven stems no longer exist upstream, so
    % the stale copies would be the ONLY version of cells that have been merged
    % away, silently inflating the folder and breaking the partition count.
    if isfolder(dst)
        old = dir(fullfile(dst,'*'));
        old = old(~[old.isdir]);
        for f = 1:numel(old), delete(fullfile(dst, old(f).name)); end
        if ~isempty(old), fprintf('%-12s cleared %d stale files\n', BINS{bi,1}, numel(old)); end
    else
        mkdir(dst);
    end
    nCopy = 0; missing = {};
    for k = 1:numel(sel)
        c = T.cell(sel(k));
        hit = find(~cellfun(@isempty, regexp({srcPng.name}, ...
                    sprintf('_cell%03d\\.png$', c), 'once')), 1);
        if isempty(hit), missing{end+1} = sprintf('cell%03d',c); continue; end %#ok<SAGROW>
        stem = erase(srcPng(hit).name, '.png');
        for ext = {'.png','.pdf','_avgproj.png'}
            f = fullfile(srcDir, [stem ext{1}]);
            if isfile(f), copyfile(f, fullfile(dst, [stem ext{1}])); end
        end
        nCopy = nCopy + 1;
    end
    Tc = T(sel, {'cell','key','group','session','nRec','nSpikes','IBI','logZ','pOnset','pPeak'});
    Tc = sortrows(Tc,'IBI');
    writetable(Tc, fullfile(dst, ['cells_' BINS{bi,1} '.csv']));
    if isinf(lo), rng_txt = sprintf('IBI <= %g', hi);
    else,         rng_txt = sprintf('%g < IBI <= %g', lo, hi); end
    fprintf('%-12s (%s) : %3d cells, %3d copied  ->  %s\n', ...
            BINS{bi,1}, rng_txt, numel(sel), nCopy, dst);
    if ~isempty(missing)
        fprintf(2,'   no render found for: %s\n', strjoin(missing,', '));
    end
    if ~isempty(sel)
        fprintf('     IBI %.2f - %.2f s | groups: ', min(Tc.IBI), max(Tc.IBI));
        g = categorical(Tc.group);
        for u = categories(g)', fprintf('%s %d  ', u{1}, nnz(g==u{1})); end
        fprintf('\n');
    end
end
assert(nAssigned == height(T), ...
    'bins are not a partition: %d assigned vs %d cells', nAssigned, height(T));
fprintf('\npartition checked: %d cells, each in exactly one folder\n', nAssigned);
