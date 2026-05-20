% pair_behavior_to_phys.m
%  Pairs behavior AVIs + DLC CSVs with calcium FOV folders by timestamp.
%  Copies files into FOV folders and trims if behavior has extra frames.
%
%  Timestamp matching:
%    AVI  -> start time parsed from Basler filename (YYYYMMDD_HHMMSSmmm)
%    TIF  -> start time parsed from ScanImage epoch in TIFF header
%    Matched by nearest |TIF_epoch - AVI_start| within tolerance.
%
%  Usage: edit physDir and behavDir below, then run.

clear; clc;

%% ========================= USER PARAMETERS =========================
physDir     = "C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys";
behavDir    = "C:\Users\Admin\Desktop\260330_sst_soma_g8s_falling_edge";
matchTolSec = 120;    % max allowed time diff for a valid match (seconds)
nChFallback = 1;      % fallback nChannels if auto-detection fails
doTrim      = true;   % trim AVI + CSV to match TIF frame count
dryRun      = false;  % true = show pairing table only, no file operations
% ====================================================================

%% 1 -- Discover calcium FOV folders (raw TIF basename == parent folder)
allTifs = dir(fullfile(physDir, '**', '*.tif'));
fov = struct('folder',{},'tifPath',{},'tifName',{},'epoch',{},'nCh',{},'rawFramesPerCh',{});

for i = 1:numel(allTifs)
    [~, tifBase]    = fileparts(allTifs(i).name);
    [~, parentName] = fileparts(allTifs(i).folder);
    if strcmp(tifBase, parentName)
        n = numel(fov) + 1;
        fov(n).folder  = allTifs(i).folder;
        fov(n).tifPath = fullfile(allTifs(i).folder, allTifs(i).name);
        fov(n).tifName = allTifs(i).name;
        fov(n).epoch          = NaT;
        fov(n).nCh            = nChFallback;
        fov(n).rawFramesPerCh = [];
    end
end

nFov = numel(fov);
fprintf('Found %d FOV folders.\n', nFov);
assert(nFov > 0, 'No FOV folders found under %s', physDir);

%% 2 -- Parse TIF start times + metadata (epoch, nChannels, framesPerSlice)
for i = 1:nFov
    desc = '';

    % Try Tiff class first (fast, reads first IFD only)
    try
        t = Tiff(fov(i).tifPath, 'r');
        try, desc = t.getTag('ImageDescription'); catch, end
        if isempty(desc)
            try, desc = t.getTag('Software'); catch, end
        end
        t.close();
    catch
        % Tiff class fails on some BigTIFFs — fall back to imfinfo(1)
        try
            info1 = imfinfo(fov(i).tifPath);
            info1 = info1(1);
            if isfield(info1,'ImageDescription'), desc = info1.ImageDescription; end
            if isempty(desc) && isfield(info1,'Software'), desc = info1.Software; end
        catch, end
    end

    % Parse epoch
    tok = regexp(desc, 'epoch\s*=\s*\[([^\]]+)\]', 'tokens', 'once');
    if ~isempty(tok)
        nums = str2double(strsplit(strtrim(tok{1})));
        fov(i).epoch = datetime(nums(1),nums(2),nums(3),nums(4),nums(5),nums(6));
    else
        d = dir(fov(i).tifPath);
        fov(i).epoch = datetime(d.datenum, 'ConvertFrom', 'datenum');
        fprintf('  WARNING: no epoch in %s, using file mod time.\n', fov(i).tifName);
    end

    % nChannels from channelSave
    tok2 = regexp(desc, 'SI\.hChannels\.channelSave\s*=\s*([^\r\n]+)', 'tokens', 'once');
    if ~isempty(tok2)
        chNums = regexp(strtrim(tok2{1}), '\d+', 'match');
        if ~isempty(chNums), fov(i).nCh = numel(chNums); end
    end

    % If _meta.mat exists, use it for nCh override and cache framesPerSlice
    metaHits = dir(fullfile(fov(i).folder, '*_meta.mat'));
    if ~isempty(metaHits)
        m = load(fullfile(metaHits(1).folder, metaHits(1).name));
        if isfield(m,'channelSave'), fov(i).nCh = numel(m.channelSave); end
        if isfield(m,'framesPerSlice') && isfield(m,'numSlices')
            fov(i).rawFramesPerCh = m.framesPerSlice * max(m.numSlices, 1);
        end
    end
end

nValidEpoch = sum(~isnat([fov.epoch]));
fprintf('  %d/%d FOVs have valid timestamps.\n', nValidEpoch, nFov);

%% 3 -- Discover behavior AVIs + parse timestamps from filenames
aviHits = dir(fullfile(behavDir, '*.avi'));
avi = struct('path',{},'name',{},'time',{},'csvPath',{},'csvName',{});

for i = 1:numel(aviHits)
    nm = aviHits(i).name;
    tok = regexp(nm, '__(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})(\d{3})\.avi$', 'tokens', 'once');
    if isempty(tok)
        warning('Cannot parse timestamp from %s, skipping.', nm);
        continue;
    end
    v = cellfun(@str2double, tok);
    n = numel(avi) + 1;
    avi(n).path    = fullfile(aviHits(i).folder, nm);
    avi(n).name    = nm;
    avi(n).time    = datetime(v(1),v(2),v(3),v(4),v(5),v(6)+v(7)/1000);
    avi(n).csvPath = "";
    avi(n).csvName = "";
end

nAvi = numel(avi);
fprintf('Found %d behavior AVIs.\n', nAvi);
assert(nAvi > 0, 'No AVI files found in %s', behavDir);

%% 4 -- Find DLC CSVs and match to AVIs
csvHits = dir(fullfile(behavDir, '**', '*DLC*.csv'));
fprintf('Found %d DLC CSV files.\n', numel(csvHits));

for i = 1:nAvi
    aviBase = erase(avi(i).name, '.avi');
    for j = 1:numel(csvHits)
        if startsWith(csvHits(j).name, aviBase)
            avi(i).csvPath = fullfile(csvHits(j).folder, csvHits(j).name);
            avi(i).csvName = csvHits(j).name;
            break
        end
    end
    if strlength(avi(i).csvPath) == 0
        warning('No DLC CSV for %s', avi(i).name);
    end
end

%% 5 -- Timestamp matching: pair each AVI to nearest FOV
[~, ord] = sort([avi.time]);   avi = avi(ord);
[~, ord] = sort([fov.epoch]);  fov = fov(ord);

fovTimes = [fov.epoch];          % 1×nFov row
fovUsed  = false(1, nFov);       % 1×nFov row (must match fovTimes)

matchFov   = nan(nAvi, 1);
matchDelta = nan(nAvi, 1);

for i = 1:nAvi
    diffs = abs(seconds(avi(i).time - fovTimes));
    diffs(fovUsed | isnat(fovTimes)) = Inf;
    [bestDiff, bestJ] = min(diffs);

    if bestDiff <= matchTolSec
        matchFov(i)   = bestJ;
        matchDelta(i) = bestDiff;
        fovUsed(bestJ) = true;
    else
        warning('AVI %s: nearest FOV is %.1f sec away (tolerance=%d). No match.', ...
            avi(i).name, bestDiff, matchTolSec);
    end
end

nPaired = sum(~isnan(matchFov));
fprintf('\nMatched %d/%d AVIs to FOVs.  Unmatched FOVs: %d\n', ...
    nPaired, nAvi, sum(~fovUsed));

%% 6 -- Determine frame counts and build pairing table
pairs = struct('ai',{},'fi',{},'delta',{},'tifFr',{},'aviFr',{},...
               'trim',{},'skip',{},'hasCSV',{});

for i = 1:nAvi
    fi = matchFov(i);
    if isnan(fi), continue; end

    p.ai    = i;
    p.fi    = fi;
    p.delta = matchDelta(i);
    p.skip  = isfile(fullfile(fov(fi).folder, avi(i).name));
    p.hasCSV = strlength(avi(i).csvPath) > 0;

    % TIF frames per channel (raw, before any frame dropping)
    if ~isempty(fov(fi).rawFramesPerCh)
        p.tifFr = fov(fi).rawFramesPerCh;
    else
        fprintf('  Reading imfinfo for %s (no _meta.mat) ...\n', fov(fi).tifName);
        info = imfinfo(fov(fi).tifPath);
        p.tifFr = numel(info) / fov(fi).nCh;
    end

    % AVI frames
    vr = VideoReader(avi(i).path);
    p.aviFr = vr.NumFrames;

    p.trim = p.aviFr > p.tifFr;
    pairs(end+1) = p; %#ok<SAGROW>
end

nPairs = numel(pairs);

%% 7 -- Display pairing table
fprintf('\n');
fprintf('%-4s  %-12s  %-45s  %-12s  %7s  %8s  %8s  %-5s  %-5s\n', ...
    '#','AVI time','FOV folder','TIF epoch','dt(s)','AVI_fr','TIF_fr','Trim','Skip');
fprintf('%s\n', repmat('-', 1, 120));

nSkip = 0;  nTrim = 0;
for k = 1:nPairs
    p  = pairs(k);
    ai = p.ai;  fi = p.fi;
    [~, fovName] = fileparts(fov(fi).folder);

    trimStr = ''; skipStr = '';
    if p.skip,  skipStr = 'yes'; nSkip = nSkip+1; end
    if p.trim,  trimStr = 'yes'; nTrim = nTrim+1; end

    fprintf('%-4d  %s  %-45s  %s  %7.1f  %8d  %8d  %-5s  %-5s\n', ...
        k, datestr(avi(ai).time,'HH:MM:SS'), fovName, ...
        datestr(fov(fi).epoch,'HH:MM:SS'), p.delta, ...
        p.aviFr, p.tifFr, trimStr, skipStr);
end

fprintf('\nSummary: %d pairs, %d to trim, %d already exist (skip)\n', nPairs, nTrim, nSkip);

if nPairs - nSkip == 0
    fprintf('Nothing to copy. Done.\n');
    return
end

if dryRun
    fprintf('Dry run — no files copied.\n');
    return
end

%% 8 -- User confirmation
reply = input('Proceed with copy? (y/n): ', 's');
if ~strcmpi(reply, 'y')
    fprintf('Aborted.\n');
    return
end

%% 9 -- Copy + trim
nCopied = 0;  nTrimmed = 0;
for k = 1:nPairs
    p  = pairs(k);
    ai = p.ai;  fi = p.fi;
    if p.skip
        fprintf('[%d/%d] SKIP %s (already in FOV)\n', k, nPairs, avi(ai).name);
        continue
    end

    fovFolder = fov(fi).folder;
    aviDst    = fullfile(fovFolder, avi(ai).name);
    nKeep     = p.tifFr;

    % ---- AVI ----
    if p.trim && doTrim
        fprintf('[%d/%d] Trim-copy AVI (%d -> %d frames) ...', k, nPairs, p.aviFr, nKeep);
        vr = VideoReader(avi(ai).path);
        if vr.BitsPerPixel <= 8
            prof = 'Grayscale AVI';
        else
            prof = 'Uncompressed AVI';
        end
        vw = VideoWriter(aviDst, prof);
        vw.FrameRate = vr.FrameRate;
        open(vw);
        for f = 1:nKeep
            writeVideo(vw, readFrame(vr));
        end
        close(vw);
        nTrimmed = nTrimmed + 1;
        fprintf(' done\n');
    else
        fprintf('[%d/%d] Copy AVI ...', k, nPairs);
        copyfile(avi(ai).path, aviDst);
        fprintf(' done\n');
    end

    % ---- DLC CSV ----
    if p.hasCSV
        csvDst = fullfile(fovFolder, avi(ai).csvName);
        if p.trim && doTrim
            fprintf('         Trim-copy CSV (%d -> %d rows) ...', p.aviFr, nKeep);
            lines = readlines(avi(ai).csvPath);
            % 3 header lines + nKeep data rows
            nLines = min(3 + nKeep, numel(lines));
            writelines(lines(1:nLines), csvDst);
            fprintf(' done\n');
        else
            fprintf('         Copy CSV ...');
            copyfile(avi(ai).csvPath, csvDst);
            fprintf(' done\n');
        end
    end

    nCopied = nCopied + 1;
end

%% 10 -- Summary
fprintf('\n========== Complete ==========\n');
fprintf('Pairs:    %d\n', nPairs);
fprintf('Copied:   %d\n', nCopied);
fprintf('Trimmed:  %d\n', nTrimmed);
fprintf('Skipped:  %d (already existed)\n', nSkip);
fprintf('Unmatched FOVs (no behavior): %d\n', sum(~fovUsed));
fprintf('Done.\n');
