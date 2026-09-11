% pair_behavior_to_phys.m
%  Pairs behavior AVIs (+ DLC CSVs) with calcium FOV folders by timestamp,
%  then copies each video into its matching FOV folder.
%
%  Handles BOTH Basler naming/formats automatically:
%    legacy Pylon-Viewer :  ...__YYYYMMDD_HHMMSSmmm.avi     (flat in behavDir)
%    new dual-cam GUI    :  <name>_YYYYMMDD_HHMMSS_runNNN.avi
%                           nested as  <name>/<name>_<ts>_run<NNN>/*.avi
%                           with a sibling timestamps.csv
%
%  Timestamp matching:
%    AVI  -> start time parsed from the filename (both formats above)
%    TIF  -> start time parsed from ScanImage epoch in TIFF header
%    Matched by nearest |TIF_epoch - AVI_start| within tolerance.
%
%  Trimming (both formats trim to the TIF frame count):
%    legacy (VideoReader-readable) : trimmed via VideoReader->VideoWriter (re-encode)
%    new (FFV1, lossless)          : MATLAB VideoReader CANNOT decode FFV1, so trim
%                                    via ffmpeg "-frames:v N -c copy" -- lossless and
%                                    frame-exact (FFV1 is all-intra). Frame count for
%                                    the new format is read from the sibling
%                                    timestamps.csv (no video decode needed).
%    The new cam is hardware-triggered by the 2P frames so counts usually already
%    match; trim is the safety net. See feedback_ffv1_video_codec.
%
%  Usage: edit physDir and behavDir below, then run.

clear; clc;

%% ========================= USER PARAMETERS =========================
physDir     = "C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\phys";

% "" = find it automatically: look beside the session folder for a sibling
% whose name starts with the session name (..._falling, ..._falling_edge,
% ..._behavior), then fall back to the session folder itself. Set it
% explicitly to override. The layout INSIDE it does not matter -- the AVI
% search is recursive, so cam1/<run>/*.avi, a flat pile of *.avi, or anything
% else all work, and non-video files (timesheet .xlsx, notes) are ignored
% because only *.avi is globbed.
behavDir    = "C:\Users\Admin\Desktop\260909_ChAT-g8m-shiverer_falling";

matchTolSec = 60;    % max allowed time diff for a valid match (seconds)
nChFallback = 1;      % fallback nChannels if auto-detection fails
doTrim      = true;   % trim AVI + CSV to match TIF frame count
dryRun      = false;  % true = show pairing table only, no file operations
ffmpegExe   = "";     % "" = auto-detect. Needed only to trim new-format FFV1 videos.
% ====================================================================

assert(isfolder(physDir), 'physDir does not exist: %s', physDir);
if strlength(behavDir) == 0
    behavDir = find_behav_dir(physDir);
end
assert(isfolder(behavDir), 'behavDir does not exist: %s', behavDir);
fprintf('  phys : %s\n  behav: %s\n', physDir, behavDir);

% Auto-detect ffmpeg (used for lossless FFV1 trimming of new-format videos).
% ORDER MATTERS: the conda ffmpeg on this machine (both on PATH and under
% envs\*\Library\bin) dies with 0xC0000139, and `ffmpeg -version` returns 0
% anyway, so a PATH-first search picks a binary that then fails every trim
% with only "ffmpeg trim failed -> copying untrimmed" to show for it. The
% imageio_ffmpeg binary works, so it is searched first.
if strlength(ffmpegExe) == 0
    cand = strings(0,1);
    for env = ["dlc310", "cellpose-gpu"]
        binDir = fullfile(getenv("USERPROFILE"), ".conda", "envs", env, ...
                          "lib", "site-packages", "imageio_ffmpeg", "binaries");
        hits = dir(fullfile(binDir, 'ffmpeg*.exe'));      % version is in the name
        for h = 1:numel(hits)
            cand(end+1,1) = string(fullfile(hits(h).folder, hits(h).name)); %#ok<SAGROW>
        end
    end
    cand = [cand
            fullfile(getenv("USERPROFILE"), ".conda", "envs", "dlc310",       "Library", "bin", "ffmpeg.exe")
            fullfile(getenv("USERPROFILE"), ".conda", "envs", "cellpose-gpu", "Library", "bin", "ffmpeg.exe")
            "ffmpeg"];   % on PATH, last resort
    ffmpegExe = "";
    for c = cand(:)'
        if c == "ffmpeg"
            [st,~] = system('ffmpeg -version');
            if st == 0, ffmpegExe = "ffmpeg"; break; end
        elseif isfile(c)
            ffmpegExe = c; break;
        end
    end
end
if strlength(ffmpegExe) == 0
    fprintf('  NOTE: ffmpeg not found -> new-format (FFV1) videos will be copied UNTRIMMED.\n');
else
    fprintf('  ffmpeg: %s\n', ffmpegExe);
end

%% 1 -- Discover calcium FOV folders (raw TIF basename == parent folder)
allTifs = dir(fullfile(physDir, '**', '*.tif'));
fov = struct('folder',{},'tifPath',{},'tifName',{},'epoch',{},'nCh',{},'rawFramesPerCh',{});

for i = 1:numel(allTifs)
    [~, tifBase]    = fileparts(allTifs(i).name);
    % NOT fileparts(): on a FOLDER path it treats the last dot as an extension,
    % so "fov1_1.7x_y1650x1320_..." comes back as "fov1_1", the name test below
    % fails, and the FOV is dropped with NO warning -- the script just reports a
    % smaller count. Any decimal zoom or power triggers it (1.7x, 2.1x, 2.5x,
    % 15.5lp, 17.5lp); on 260728 it silently lost 11 of 19 FOVs.
    parentName      = lastdir(allTifs(i).folder);
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

%% 3 -- Discover behavior AVIs (recursive) + parse timestamps from filenames
%  Recursive so it catches the new dual-cam GUI layout <name>/<name>_<ts>_runNNN/*.avi
%  as well as legacy AVIs sitting flat in behavDir.
aviHits = dir(fullfile(behavDir, '**', '*.avi'));
avi = struct('path',{},'name',{},'time',{},'csvPath',{},'csvName',{},'tsPath',{}, ...
             'isNew',{},'cam',{});

% Anything already sitting inside physDir is a video this script copied on an
% earlier run. Skipping those is what makes behavDir safe to point anywhere --
% including the session root, or physDir itself when there is no separate
% behaviour folder at all (260721_Sert is laid out that way).
physFull = string(fullfile(char(physDir)));

nSkipCopied = 0;
for i = 1:numel(aviHits)
    if aviHits(i).isdir, continue; end
    nm = aviHits(i).name;

    if startsWith(string(fullfile(aviHits(i).folder)), physFull)
        nSkipCopied = nSkipCopied + 1;
        continue;
    end

    % Both acquisition routes are equally valid and MAY BE MIXED IN ONE SESSION
    % (260721_Sert has one Pylon file among GUI ones), so the format is decided
    % per file, never per session. The leading capture is the camera label --
    % everything before the timestamp -- so a second camera is recognised
    % whatever it was named (basler_dual_acq.py takes --cam2-name freely, so it
    % is not necessarily "cam2").
    %   MATLAB GUI  : <cam>_YYYYMMDD_HHMMSS_runNNN.avi        (run suffix, no ms)
    %   Pylon Viewer: <model>__<serial>__YYYYMMDD_HHMMSSmmm.avi (3-digit ms)
    tokNew = regexp(nm, '^(.*)_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})_run\d+\.avi$', 'tokens', 'once');
    tokOld = regexp(nm, '^(.*)_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})(\d{3})\.avi$', 'tokens', 'once');

    if ~isempty(tokNew)
        camLbl = tokNew{1};
        v = cellfun(@str2double, tokNew(2:end));
        tt = datetime(v(1),v(2),v(3),v(4),v(5),v(6));
        isNewFmt = true;
    elseif ~isempty(tokOld)
        camLbl = tokOld{1};
        v = cellfun(@str2double, tokOld(2:end));
        tt = datetime(v(1),v(2),v(3),v(4),v(5),v(6)+v(7)/1000);
        isNewFmt = false;
    else
        warning('Cannot parse timestamp from %s, skipping.', nm);
        continue;
    end
    camLbl = regexprep(camLbl, '_+$', '');      % Pylon leaves a trailing '_'

    % new format: a per-run timestamps.csv sits beside the AVI (row count = frame count)
    tsPath = "";
    tsA = fullfile(aviHits(i).folder, 'timestamps.csv');
    tsB = fullfile(aviHits(i).folder, [erase(nm,'.avi') '_timestamps.csv']);
    if     isfile(tsA), tsPath = string(tsA);
    elseif isfile(tsB), tsPath = string(tsB);
    end

    n = numel(avi) + 1;
    avi(n).path    = fullfile(aviHits(i).folder, nm);
    avi(n).name    = nm;
    avi(n).time    = tt;
    avi(n).csvPath = "";
    avi(n).csvName = "";
    avi(n).tsPath  = tsPath;
    avi(n).isNew   = isNewFmt;
    avi(n).cam     = camLbl;
end
if nSkipCopied > 0
    fprintf('  Ignored %d AVI(s) already inside physDir (copied by an earlier run).\n', nSkipCopied);
end

nAvi = numel(avi);
fprintf('Found %d behavior AVIs (%d GUI-format, %d Pylon-format).\n', ...
    nAvi, sum([avi.isNew]), sum(~[avi.isNew]));
if nAvi == 0 && nSkipCopied > 0
    % Every video found was already inside physDir, i.e. a previous run paired
    % this session. That is a finished state, not an error.
    fprintf('All %d video(s) are already in their FOV folders. Nothing to do.\n', nSkipCopied);
    return
end
assert(nAvi > 0, 'No AVI files found in %s', behavDir);

camList = unique(string({avi.cam}), 'stable');
fprintf('  %d camera(s): %s\n', numel(camList), strjoin(cellstr(camList), ', '));

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
fovUsed  = false(1, nFov);       % any camera matched this FOV

matchFov   = nan(nAvi, 1);
matchDelta = nan(nAvi, 1);

% MATCH ONE CAMERA AT A TIME. Each camera records the SAME FOV, so a single
% shared "already used" mask would let cam1 claim a FOV and leave the second
% camera's video of it unmatched. The exclusivity that matters is one video per
% camera per FOV, not one video per FOV.
for c = camList(:)'
    idxC  = find(strcmp({avi.cam}, c));      % already in time order
    usedC = false(1, nFov);
    for i = idxC
        diffs = abs(seconds(avi(i).time - fovTimes));
        diffs(usedC | isnat(fovTimes)) = Inf;
        [bestDiff, bestJ] = min(diffs);

        if bestDiff <= matchTolSec
            matchFov(i)   = bestJ;
            matchDelta(i) = bestDiff;
            usedC(bestJ)  = true;
        else
            warning('AVI %s: nearest free FOV is %.1f sec away (tolerance=%d). No match.', ...
                avi(i).name, bestDiff, matchTolSec);
        end
    end
    fovUsed = fovUsed | usedC;
end

nPaired = sum(~isnan(matchFov));
fprintf('\nMatched %d/%d AVIs to FOVs.  FOVs with no behaviour: %d\n', ...
    nPaired, nAvi, sum(~fovUsed));

%% 6 -- Determine frame counts and build pairing table
pairs = struct('ai',{},'fi',{},'delta',{},'tifFr',{},'aviFr',{},...
               'trim',{},'skip',{},'hasCSV',{},'isNew',{},'readErr',{});

for i = 1:nAvi
    fi = matchFov(i);
    if isnan(fi), continue; end

    p.ai    = i;
    p.fi    = fi;
    p.delta = matchDelta(i);
    p.skip  = isfile(fullfile(fov(fi).folder, avi(i).name));
    p.hasCSV = strlength(avi(i).csvPath) > 0;

    % TIF frames per channel (raw, before any frame dropping).
    %
    % THE RAW FILE IS AUTHORITATIVE, not _meta.mat.  _meta.mat records the
    % CONFIGURED framesPerSlice; when an acquisition is stopped by hand -- which
    % is routine -- the configured value exceeds what was actually written
    % (6000 vs 3007 on 260728_vglut2), and trusting it makes aviFr > tifFr false,
    % so the trim is skipped in exactly the case that needs it.  Counting TIFF
    % directories costs ~1.0 s on a 3 GB / 6000-page file (imfinfo takes 4.3 s),
    % so there is no reason to prefer the header.
    nPages = count_tif_pages(fov(fi).tifPath);
    if isnan(nPages)
        assert(~isempty(fov(fi).rawFramesPerCh), ...
            'Cannot read %s and no _meta.mat to fall back on.', fov(fi).tifName);
        p.tifFr = fov(fi).rawFramesPerCh;
        fprintf(2, '  WARNING: cannot read %s -- falling back to _meta.mat (%d fr)\n', ...
            fov(fi).tifName, p.tifFr);
    else
        p.tifFr = nPages / fov(fi).nCh;
        if ~isempty(fov(fi).rawFramesPerCh) && fov(fi).rawFramesPerCh ~= p.tifFr
            fprintf('  %s: stopped early -- %d frames on disk vs %d configured; using disk.\n', ...
                fov(fi).tifName, p.tifFr, fov(fi).rawFramesPerCh);
        end
    end

    % AVI frame count.
    %   new format (FFV1): from sibling timestamps.csv (no video decode)
    %   legacy format    : from VideoReader
    p.isNew   = avi(i).isNew;
    p.readErr = false;
    p.aviFr   = NaN;
    if avi(i).isNew
        p.aviFr = count_ts_frames(avi(i).tsPath);   % NaN if timestamps.csv absent
    else
        try
            vr = VideoReader(avi(i).path);
            p.aviFr = vr.NumFrames;
        catch ME
            p.readErr = true;
            fprintf(2, '  VideoReader failed on %s (%s)\n', avi(i).name, ME.message);
        end
    end

    % Trim only if we know the AVI has extra frames (and, for new format, ffmpeg exists).
    p.trim = isfinite(p.aviFr) && p.aviFr > p.tifFr;
    if p.isNew && strlength(ffmpegExe) == 0, p.trim = false; end
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
    fovName = lastdir(fov(fi).folder);   % see the note at the discovery loop

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
    if p.trim && doTrim && p.isNew
        % FFV1 (lossless, all-intra): keep first nKeep frames with -c copy (no re-encode).
        fprintf('[%d/%d] ffmpeg trim-copy FFV1 (%d -> %d frames) ...', k, nPairs, p.aviFr, nKeep);
        if isfile(aviDst), delete(aviDst); end
        cmd = sprintf('"%s" -y -v error -i "%s" -frames:v %d -c copy "%s"', ...
            ffmpegExe, avi(ai).path, nKeep, aviDst);
        st = system(cmd);
        if st == 0 && isfile(aviDst)
            nTrimmed = nTrimmed + 1; fprintf(' done\n');
        else
            fprintf(2, ' ffmpeg trim failed -> copying untrimmed\n');
            copyfile(avi(ai).path, aviDst);
        end
    elseif p.trim && doTrim
        % legacy VideoReader-readable path (re-encode a shorter copy)
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
        if isfinite(p.aviFr) && p.aviFr ~= p.tifFr
            fprintf(2, ' [AVI %d fr vs TIF %d fr -- copied UNTRIMMED]', p.aviFr, p.tifFr);
        end
        fprintf(' done\n');
    end

    % ---- timestamps.csv (new format) ----
    if avi(ai).isNew && strlength(avi(ai).tsPath) > 0
        tsName = [erase(avi(ai).name, '.avi') '_timestamps.csv'];
        tsDst  = fullfile(fovFolder, tsName);
        if p.trim && doTrim
            % keep header + first nKeep rows so the CSV matches the trimmed video
            tl = readlines(avi(ai).tsPath);
            tl = tl(strlength(strtrim(tl)) > 0);
            write_lines(tl(1:min(1 + nKeep, numel(tl))), tsDst);
            fprintf('         Trim-copy timestamps.csv -> %s\n', tsName);
        else
            copyfile(avi(ai).tsPath, tsDst);
            fprintf('         Copy timestamps.csv -> %s\n', tsName);
        end
    end

    % ---- DLC CSV ----
    if p.hasCSV
        csvDst = fullfile(fovFolder, avi(ai).csvName);
        if p.trim && doTrim
            fprintf('         Trim-copy CSV (%d -> %d rows) ...', p.aviFr, nKeep);
            lines = readlines(avi(ai).csvPath);
            % 3 header lines + nKeep data rows
            nLines = min(3 + nKeep, numel(lines));
            write_lines(lines(1:nLines), csvDst);
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

%% =============================== LOCAL FUNCTIONS ===============================
function bd = find_behav_dir(physDir)
% Locate the behaviour folder from physDir alone.
%
% The convention is a sibling of the session folder sharing its name plus a
% suffix (260728_vglut2_soma-g8s -> 260728_vglut2_soma-g8s_falling), but it is
% not universal: 260721_Sert has no separate behaviour folder at all and the
% videos live beside the TIFFs. So: prefer a sibling that actually contains
% AVIs, and otherwise hand back the session folder and let the recursive search
% find them wherever they are.
    sessionRoot = fileparts(char(physDir));           % ...\<session>\phys -> ...\<session>
    [parentDir, sessionName] = fileparts(sessionRoot);

    hits = dir(fullfile(parentDir, [sessionName '*']));
    hits = hits([hits.isdir]);
    best = "";
    for k = 1:numel(hits)
        if any(strcmp(hits(k).name, {'.', '..'})), continue; end
        cand = fullfile(hits(k).folder, hits(k).name);
        if strcmpi(cand, sessionRoot), continue; end   % the session folder itself
        if ~isempty(dir(fullfile(cand, '**', '*.avi')))
            best = string(cand);
            break
        end
    end

    if strlength(best) > 0
        bd = best;
        fprintf('  behavDir auto-detected: %s\n', bd);
    else
        bd = string(sessionRoot);
        fprintf(['  No sibling behaviour folder with AVIs found; searching the session\n' ...
                 '  folder itself (%s). Videos already inside physDir are ignored.\n'], bd);
    end
end

function n = count_tif_pages(tifPath)
% Number of IFDs actually present in a TIFF, by walking the directory chain.
% ~4x faster than imfinfo (1.0 s vs 4.3 s on a 3 GB / 6000-page BigTIFF)
% because it never builds the per-page metadata structs.  Returns NaN if the
% file cannot be opened, so the caller can fall back.
    n = NaN;
    try
        t = Tiff(tifPath, 'r');
        c = onCleanup(@() t.close());
        n = 1;
        while ~t.lastDirectory()
            t.nextDirectory();
            n = n + 1;
        end
    catch
        n = NaN;
    end
end

function nm = lastdir(p)
% Last component of a folder path, dots and all. fileparts() cannot be used on
% a folder: it splits on the final dot as if the name had a file extension.
    parts = split(string(p), filesep);
    parts(strlength(parts) == 0) = [];      % trailing separator
    nm = char(parts(end));
end

function write_lines(L, dst)
% writelines() only exists from R2022a; this machine runs R2021b, where the
% trim path died with "Unrecognized function or variable 'writelines'".
    fid = fopen(dst, 'w');
    assert(fid > 0, 'Cannot open %s for writing', dst);
    c = onCleanup(@() fclose(fid));
    for i = 1:numel(L)
        fprintf(fid, '%s\n', char(L(i)));
    end
end

function n = count_ts_frames(tsPath)
% Number of data rows in a per-run timestamps.csv (= AVI frame count).
% Returns NaN if the file is missing/unreadable. Assumes 1 header row.
    n = NaN;
    if strlength(tsPath) == 0 || ~isfile(tsPath), return; end
    try
        L = readlines(tsPath);
        L = L(strlength(strtrim(L)) > 0);   % ignore blank lines
        n = max(numel(L) - 1, 0);
    catch
        n = NaN;
    end
end
