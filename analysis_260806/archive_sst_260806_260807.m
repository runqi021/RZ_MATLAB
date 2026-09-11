function archive_sst_260806_260807(mode)
%ARCHIVE_SST_260806_260807  Put the two new Sst sessions into Ventral_surface_summary.
% -----------------------------------------------------------------------
%   archive_sst_260806_260807('dry')    list what would be copied (default)
%   archive_sst_260806_260807('copy')   copy
%
% Layout follows the archive convention  <Genotype>\<MMDD>\<site>\<recording>\ .
% The site is taken from the recording's own roiN prefix (roi1_2x_... -> roi1),
% which is what the experimenter used to mark a field; it is NOT a claim that two
% recordings under the same roiN are the same cell.
%
% Only the small files the map and the breath analysis read are copied -- masks,
% average projection, meta, dFF, breath, spikes.  No TIFF stacks: the two sessions
% together are ~150 GB of movies and none of it is needed here.
%
% COPY, NOT MOVE.  C: is left untouched, and nothing already in the archive is
% overwritten.
%
% Runqi Zhang / 2026-08-06

if nargin < 1, mode = 'dry'; end
doCopy = strcmpi(mode,'copy');

SRC = { 'C:\260806_sst-soma-g8s\phys', '0806'
        'C:\260807_sst-soma-g8s\phys', '0807' };
ARCHIVE = 'D:\Ventral_surface_summary\Sst';

PATTERNS = {'*_ch1_dFF.mat', '*_ch1_meta.mat', '*cpSAM_output.mat', ...
            '*_AVG_for_CP.tif', '*_AVG_ROIlabel.tif', '*_AVG_ROImask.tif', ...
            'breath_pc1.mat', 'breath_peak_pc1.mat', 'breath_insp_start_pc1.mat', ...
            'ca_spike_data.mat'};

totB = 0; totF = 0; nRec = 0;
for s = 1:size(SRC,1)
    root = SRC{s,1}; dateStr = SRC{s,2};
    d = dir(root); d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}) & ~startsWith({d.name},'.'));
    fprintf('\n== %s -> %s\\%s ==\n', root, ARCHIVE, dateStr);
    for i = 1:numel(d)
        rp = fullfile(root, d(i).name);
        if isempty(dir(fullfile(rp,'*_ch1_dFF.mat'))), continue; end
        tok  = regexp(d(i).name, '^(roi\d+)', 'tokens', 'once');
        site = 'roi0'; if ~isempty(tok), site = tok{1}; end
        dst  = fullfile(ARCHIVE, dateStr, site, d(i).name);

        nf = 0; nb = 0;
        for p = 1:numel(PATTERNS)
            src = dir(fullfile(rp, PATTERNS{p}));
            for k = 1:numel(src)
                tgt = fullfile(dst, src(k).name);
                if isfile(tgt), continue; end
                nf = nf + 1; nb = nb + src(k).bytes;
                if doCopy
                    if ~isfolder(dst), mkdir(dst); end
                    copyfile(fullfile(src(k).folder, src(k).name), tgt);
                end
            end
        end
        nRec = nRec + 1; totF = totF + nf; totB = totB + nb;
        fprintf('  %-40s -> %s\\%-6s  %2d files  %5.1f MB\n', ...
                d(i).name(1:min(40,end)), dateStr, site, nf, nb/1e6);
    end
end
fprintf('\n%s: %d recordings, %d files, %.2f GB\n', upper(mode), nRec, totF, totB/1e9);
if ~doCopy, fprintf('Dry run. Re-run with ''copy'' to write.\n'); end
end
