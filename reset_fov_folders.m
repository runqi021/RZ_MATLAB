%% reset_fov_folders.m
%  Reset FOV subfolders to raw state: keep only the raw TIFF, delete
%  everything else (processed TIFFs, MATs, PNGs, PDFs, subfolders, etc.).
%
%  Raw TIFF = {FOV_NAME}/{FOV_NAME}.tif  (filename matches parent folder)
%
%  Two input modes:
%    masterFolder: scan all experiments (auto-dives phys/processed/)
%    folderPath:   single folder
%
%  Set dryRun = true to preview without deleting.
%
%  RZ 2026-04

clear; clc;

%% ============================== USER ==============================
masterFolder = "";
folderPath = "C:\Users\Admin\Desktop\260224_vglut2_soma_g8s\phys\processed\breathing";
%folderPath   = "D:\251124_live_vglut2_soma_g8s+cy5\phys\IO";

keepMC = false;    % true = keep raw TIFF + preproc + MC + meta; false = keep raw TIFF only
dryRun = false;    % true = preview only; false = actually delete

%% ==================================================================
%  DISCOVER FOV FOLDERS
%  ==================================================================
fovDirs = discover_fov_folders(masterFolder, folderPath);
nFOV = numel(fovDirs);
fprintf('\n=== Found %d FOV folder(s) ===\n\n', nFOV);

if nFOV == 0
    fprintf('Nothing to do.\n');
    return
end

%% ==================================================================
%  SCAN & REPORT
%  ==================================================================
totalFiles   = 0;
totalFolders = 0;
totalBytes   = 0;

for ff = 1:nFOV
    fovDir  = fovDirs(ff);
    parts   = split(fovDir, filesep);
    fovName = parts(end);
    rawTif  = fullfile(fovDir, fovName + ".tif");

    fprintf('[%d/%d] %s\n', ff, nFOV, fovDir);

    if ~isfile(rawTif)
        fprintf('  WARNING: raw TIFF not found (%s.tif) — skipping\n', fovName);
        continue
    end

    % Build keep set
    keepSet = string(rawTif);
    if keepMC
        keepSet = [keepSet; find_mc_files(fovDir, fovName)]; %#ok<AGROW>
    end

    % List everything in this folder
    allItems = dir(fullfile(fovDir, '*'));
    allItems = allItems(~ismember({allItems.name}, {'.', '..'}));

    nDel = 0;
    nKeep = 0;
    for ii = 1:numel(allItems)
        itemPath = string(fullfile(allItems(ii).folder, allItems(ii).name));
        if ~allItems(ii).isdir && any(itemPath == keepSet)
            nKeep = nKeep + 1;
            fprintf('  KEEP:  %s\n', allItems(ii).name);
            continue
        end

        if allItems(ii).isdir
            subFiles = dir(fullfile(char(itemPath), '**', '*'));
            subFiles = subFiles(~[subFiles.isdir]);
            subBytes = sum([subFiles.bytes]);
            fprintf('  DELETE dir:  %s/  (%d files, %.1f MB)\n', ...
                allItems(ii).name, numel(subFiles), subBytes / 1e6);
            totalFolders = totalFolders + 1;
            totalBytes   = totalBytes + subBytes;
        else
            fprintf('  DELETE file: %s  (%.1f MB)\n', ...
                allItems(ii).name, allItems(ii).bytes / 1e6);
            totalBytes = totalBytes + allItems(ii).bytes;
        end
        nDel = nDel + 1;
        totalFiles = totalFiles + 1;
    end

    if nDel == 0
        fprintf('  (already clean)\n');
    end
end

fprintf('\n--- Summary: %d items to delete across %d FOVs (%.1f GB) ---\n', ...
    totalFiles, nFOV, totalBytes / 1e9);

if dryRun
    fprintf('\n*** DRY RUN — nothing deleted. Set dryRun = false to execute. ***\n');
    return
end

%% ==================================================================
%  CONFIRM & DELETE
%  ==================================================================
resp = input('\nType YES to confirm deletion: ', 's');
if ~strcmp(strtrim(resp), 'YES')
    fprintf('Aborted.\n');
    return
end

nDeleted = 0;
nErrors  = 0;

for ff = 1:nFOV
    fovDir  = fovDirs(ff);
    parts   = split(fovDir, filesep);
    fovName = parts(end);
    rawTif  = fullfile(fovDir, fovName + ".tif");

    if ~isfile(rawTif), continue; end

    keepSet = string(rawTif);
    if keepMC
        keepSet = [keepSet; find_mc_files(fovDir, fovName)]; %#ok<AGROW>
    end

    allItems = dir(fullfile(fovDir, '*'));
    allItems = allItems(~ismember({allItems.name}, {'.', '..'}));

    for ii = 1:numel(allItems)
        itemPath = string(fullfile(allItems(ii).folder, allItems(ii).name));
        if ~allItems(ii).isdir && any(itemPath == keepSet), continue; end

        try
            if allItems(ii).isdir
                rmdir(char(itemPath), 's');
            else
                delete(char(itemPath));
            end
            nDeleted = nDeleted + 1;
        catch ME
            fprintf(2, '  ERROR deleting %s: %s\n', itemPath, ME.message);
            nErrors = nErrors + 1;
        end
    end
end

fprintf('\nDone. Deleted: %d  |  Errors: %d\n', nDeleted, nErrors);


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function fovDirs = discover_fov_folders(masterFolder, folderPath)
%DISCOVER_FOV_FOLDERS  Find all FOV subfolders containing a raw TIFF.
%  A FOV folder is any folder where {folderName}/{folderName}.tif exists.

fovDirs = string.empty;

if masterFolder ~= ""
    scanDirs = resolve_all_scan_dirs(string(masterFolder));
elseif folderPath ~= ""
    scanDirs = resolve_scan_dirs(string(folderPath));
else
    error('Set either masterFolder or folderPath.');
end

for i = 1:numel(scanDirs)
    fovDirs = [fovDirs, find_fov_dirs(scanDirs(i))]; %#ok<AGROW>
end
fovDirs = unique(fovDirs, 'stable');
end


function scanDirs = resolve_all_scan_dirs(masterDir)
%RESOLVE_ALL_SCAN_DIRS  From a master directory, find all phys/ dirs and resolve.
assert(isfolder(masterDir), 'masterFolder not found: %s', masterDir);
allDirs = dir(fullfile(masterDir, '**', 'phys'));
scanDirs = string.empty;
for i = 1:numel(allDirs)
    if allDirs(i).isdir
        physDir = string(fullfile(allDirs(i).folder, allDirs(i).name));
        scanDirs = [scanDirs, resolve_scan_dirs(physDir)]; %#ok<AGROW>
    end
end
if isempty(scanDirs)
    scanDirs = resolve_scan_dirs(masterDir);
end
end


function scanDirs = resolve_scan_dirs(rootDir)
%RESOLVE_SCAN_DIRS  Smart dive: phys/ -> processed/ -> category subfolders.
scanDirs = string.empty;

physDir = fullfile(rootDir, "phys");
if isfolder(physDir)
    procDir = fullfile(physDir, "processed");
    if isfolder(procDir)
        scanDirs(end+1) = procDir;
    else
        scanDirs(end+1) = physDir;
    end
    return
end

procDir = fullfile(rootDir, "processed");
if isfolder(procDir)
    scanDirs(end+1) = procDir;
    return
end

scanDirs(end+1) = rootDir;
end


function keepFiles = find_mc_files(fovDir, fovName)
%FIND_MC_FILES  Return list of MC-related files to keep.
%  Keeps: _preproc.tif, _preproc_MC.tif, _preproc_MC_MC.tif,
%         _MC_output.mat, _MC_QC.png, _MC_shifts.png,
%         _MCinfo.mat, _MCshift.png, _meta.mat
keepSuffixes = [
    "_preproc.tif"
    "_preproc_MC.tif"
    "_preproc_MC_MC.tif"
    "_preproc_MC_output.mat"
    "_preproc_MC_QC.png"
    "_preproc_MC_shifts.png"
    "_preproc_MC_MC_output.mat"
    "_preproc_MC_MC_QC.png"
    "_preproc_MC_MC_shifts.png"
    "_MCinfo.mat"
    "_MCshift.png"
    "_meta.mat"
];
stem = fovName + "_ch1";
keepFiles = string.empty;
for s = 1:numel(keepSuffixes)
    f = fullfile(fovDir, stem + keepSuffixes(s));
    if isfile(f)
        keepFiles(end+1,1) = string(f); %#ok<AGROW>
    end
end
end


function fovDirs = find_fov_dirs(scanDir)
%FIND_FOV_DIRS  Recursively find folders where {name}/{name}.tif exists.
fovDirs = string.empty;

dd = dir(fullfile(scanDir, '**', '*.tif'));
for j = 1:numel(dd)
    if dd(j).isdir, continue; end

    parts  = split(string(dd(j).folder), filesep);
    parName = parts(end);
    tifBase = erase(string(dd(j).name), ".tif");

    if tifBase == parName
        fovDirs(end+1) = string(dd(j).folder); %#ok<AGROW>
    end
end
fovDirs = unique(fovDirs, 'stable');
end
