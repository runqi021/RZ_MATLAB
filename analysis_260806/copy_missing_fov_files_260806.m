function copy_missing_fov_files_260806(mode)
%COPY_MISSING_FOV_FILES_260806  Fill in the archive's missing analysis files.
% -----------------------------------------------------------------------
%   copy_missing_fov_files_260806('dry')    list what would be copied  (default)
%   copy_missing_fov_files_260806('copy')   actually copy
%
% WHY. Ventral_surface_summary holds the FOV folders, but for many of them the
% analysis files were never copied in -- e.g. Vglut2\0728 has 16 folders and only
% 3 carry a dFF. The raw sessions on E: have dFF + breath for those same FOVs, so
% the breath x Ca analysis is missing ROIs purely because the small files are not
% in the archive. This copies ONLY those small files (a few MB per FOV), never the
% TIFF stacks.
%
% COPY, NOT MOVE -- the raw session on E: is left untouched.
% Nothing already in the archive is overwritten; existing files are skipped.
%
% A FOV whose folder does not exist in the archive is NOT created here: that needs
% a site (cellN) assignment, which is a judgement call, so it is only reported.
%
% Runqi Zhang / 2026-08-06

if nargin < 1, mode = 'dry'; end
doCopy = strcmpi(mode,'copy');

SRC = { 'Sert',   '0721', 'E:\260721_Sert_soma_G8s\phys\baseline'
        'Vglut2', '0728', 'E:\260728_vglut2_soma-g8s\phys'
        'Vgat',   '0730', 'E:\260730_vgat-g8m_shiverer\phys' };
ARCHIVE = 'D:\Ventral_surface_summary';

% the small files this analysis needs
PATTERNS = {'*_ch1_dFF.mat', '*_ch1_meta.mat', 'breath_pc1.mat', ...
            'breath_peak_pc1.mat', 'breath_insp_start_pc1.mat', 'ca_spike_data.mat'};

totBytes = 0; totFiles = 0; nFilled = 0; noFolder = {};

for s = 1:size(SRC,1)
    gen = SRC{s,1}; dt = SRC{s,2}; raw = SRC{s,3};
    if ~isfolder(raw), warning('raw session missing: %s', raw); continue; end

    % map recording name -> archive folder (any site level under <gen>)
    ad = dir(fullfile(ARCHIVE, gen, '**')); ad = ad([ad.isdir]);
    arcMap = containers.Map('KeyType','char','ValueType','char');
    for i = 1:numel(ad)
        if ~isKey(arcMap, ad(i).name), arcMap(ad(i).name) = fullfile(ad(i).folder, ad(i).name); end
    end

    d = dir(raw); d = d([d.isdir]);
    d = d(~ismember({d.name},{'.','..'}));
    d = d(~startsWith({d.name}, {'analysis_','roi_','coherence_','_','.','cell_'}));

    fprintf('\n== %s %s ==\n', gen, dt);
    for i = 1:numel(d)
        name = d(i).name;  rp = fullfile(raw, name);
        if isempty(dir(fullfile(rp,'*_ch1_dFF.mat'))), continue; end
        if ~isfile(fullfile(rp,'breath_peak_pc1.mat')), continue; end   % no trigger, unusable

        if ~isKey(arcMap, name)
            noFolder{end+1} = sprintf('%s/%s/%s', gen, dt, name); %#ok<AGROW>
            continue;
        end
        dst = arcMap(name);
        % already complete?
        if ~isempty(dir(fullfile(dst,'*_ch1_dFF.mat'))) && isfile(fullfile(dst,'breath_peak_pc1.mat'))
            continue;
        end

        nf = 0; nb = 0;
        for p = 1:numel(PATTERNS)
            src = dir(fullfile(rp, PATTERNS{p}));
            for k = 1:numel(src)
                tgt = fullfile(dst, src(k).name);
                if isfile(tgt), continue; end
                nf = nf + 1; nb = nb + src(k).bytes;
                if doCopy, copyfile(fullfile(src(k).folder, src(k).name), tgt); end
            end
        end
        nFilled = nFilled + 1; totFiles = totFiles + nf; totBytes = totBytes + nb;
        fprintf('  %-46s %2d files  %6.1f MB  -> %s\n', name(1:min(46,end)), nf, nb/1e6, ...
                strrep(dst, ARCHIVE, '<archive>'));
    end
end

fprintf('\n%s: %d FOV folders filled, %d files, %.2f GB\n', upper(mode), nFilled, totFiles, totBytes/1e9);
if ~isempty(noFolder)
    fprintf('NO archive folder (needs a site assignment, not created here): %d\n', numel(noFolder));
    fprintf('   %s\n', noFolder{:});
end
if ~doCopy, fprintf('\nDry run only. Re-run with ''copy'' to write.\n'); end
end
