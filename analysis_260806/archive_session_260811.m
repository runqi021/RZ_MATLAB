function archive_session_260811(srcRoot, genotype, dateStr, mode)
%ARCHIVE_SESSION_260811  Copy one flat session into Ventral_surface_summary.
% -----------------------------------------------------------------------
%   archive_session_260811(src, 'Vglut2', '0810')          dry run (default)
%   archive_session_260811(src, 'Vglut2', '0810', 'copy')  do it
%
% Layout: <Genotype>\<MMDD>\<site>\<recording>\  with the archive's own rule --
%   SITE COMES FROM motorPosition, single-linkage clustered at mergeSiteUm.
% Not from the filename prefix: on this rig every recording of a session is often
% called roi1_*, which would collapse genuinely different fields into one site
% (archive_into_summary_260729's filename guess got Vgat/0730 wrong the same way).
%
% Only the small files the analyses read are copied -- no TIFF stacks.
% COPY, NOT MOVE. Nothing already in the archive is overwritten.
%
% Runqi Zhang / 2026-08-11

if nargin < 4, mode = 'dry'; end
doCopy = strcmpi(mode,'copy');
ARCHIVE    = 'D:\Ventral_surface_summary';
mergeSiteUm = 100;

PATTERNS = {'*_ch1_dFF.mat', '*_ch1_meta.mat', '*cpSAM_output.mat', ...
            '*_AVG_for_CP.tif', '*_AVG_ROIlabel.tif', '*_AVG_ROImask.tif', ...
            'breath_pc1.mat', 'breath_peak_pc1.mat', 'breath_insp_start_pc1.mat', ...
            'ca_spike_data.mat'};

d = dir(srcRoot); d = d([d.isdir]);
d = d(~ismember({d.name},{'.','..'}) & ~startsWith({d.name},'.'));
keep = false(numel(d),1); XY = nan(numel(d),2);
for i = 1:numel(d)
    fp = fullfile(srcRoot, d(i).name);
    m  = dir(fullfile(fp,'*_ch1_meta.mat'));
    if isempty(m) || isempty(dir(fullfile(fp,'*_ch1_dFF.mat'))), continue; end
    M = load(fullfile(m(1).folder, m(1).name),'motorPosition');
    XY(i,:) = M.motorPosition(1:2);  keep(i) = true;
end
d = d(keep); XY = XY(keep,:);
assert(~isempty(d), 'no usable recordings under %s', srcRoot);

% single-linkage clustering on stage XY -> site index
if size(XY,1) > 1
    site = cluster(linkage(pdist(XY),'single'), 'cutoff', mergeSiteUm, 'criterion','distance');
else
    site = 1;
end
[~,~,site] = unique(site, 'stable');       % renumber in first-seen order

fprintf('\n== %s %s: %d recordings -> %d site(s) at %d um ==\n', ...
        genotype, dateStr, numel(d), max(site), mergeSiteUm);
totF = 0; totB = 0;
for s = 1:max(site)
    ix = find(site == s);
    fprintf(' cell%d  (%d rec, centre %.0f, %.0f)\n', s, numel(ix), mean(XY(ix,1)), mean(XY(ix,2)));
    for k = ix.'
        rp  = fullfile(srcRoot, d(k).name);
        dst = fullfile(ARCHIVE, genotype, dateStr, sprintf('cell%d',s), d(k).name);
        nf = 0; nb = 0;
        for p = 1:numel(PATTERNS)
            src = dir(fullfile(rp, PATTERNS{p}));
            for j = 1:numel(src)
                tgt = fullfile(dst, src(j).name);
                if isfile(tgt), continue; end
                nf = nf + 1; nb = nb + src(j).bytes;
                if doCopy
                    if ~isfolder(dst), mkdir(dst); end
                    copyfile(fullfile(src(j).folder, src(j).name), tgt);
                end
            end
        end
        totF = totF + nf; totB = totB + nb;
        fprintf('    %-42s %2d files %6.1f MB\n', d(k).name(1:min(42,end)), nf, nb/1e6);
    end
end
fprintf('%s: %d files, %.2f GB\n', upper(mode), totF, totB/1e9);
if ~doCopy, fprintf('Dry run. Re-run with ''copy'' to write.\n'); end
end
