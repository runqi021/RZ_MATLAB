%% tiff_timeline_260722.m
% Standalone: discover all raw recording TIFFs (same logic as Batch_dffQC_260325.m),
% read each file's LAST-MODIFIED time, and print a chronological timeline.
%
% Read-only: this script NEVER moves, deletes, or writes any data. It only
% lists files and prints their modification times.
%
% Uses dir() for mod time (fast) instead of imfinfo() so it does not open the
% multi-GB TIFFs.
%
% Run:  run tiff_timeline_260722.m

clear; clc;

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);

%% ============================== USER ==============================
masterFolder = "";                          % multi-experiment master dir (leave "" for single)
folderPath   = "D:\260721_Sert_soma_G8s\phys";   % single folder (ignored if masterFolder set)

SortBy   = "time";   % "time" = chronological | "name" = discovery order
SaveCSV  = true;     % also write a CSV timeline next to nothing-destructive (into folderPath root)

%% ============================== DISCOVER TIFS ==============================
tifList = discover_tifs_readonly(masterFolder, folderPath);
nF = numel(tifList);
assert(nF > 0, "No recording TIFFs found. Check folderPath: %s", folderPath);
fprintf("\n=== Found %d TIF files ===\n\n", nF);

%% ============================== READ MOD TIMES ==============================
modNum  = zeros(nF,1);           % datenum, last modified
sizeGB  = zeros(nF,1);
names   = strings(nF,1);
for i = 1:nF
    d = dir(char(tifList(i)));
    if isempty(d)
        modNum(i) = NaN; sizeGB(i) = NaN;
    else
        modNum(i) = d.datenum;
        sizeGB(i) = d.bytes / 1e9;
    end
    [~, nm, ext] = fileparts(char(tifList(i)));
    names(i) = string([nm ext]);
end
modDT = datetime(modNum, 'ConvertFrom','datenum');

%% ============================== SORT ==============================
if SortBy == "time"
    [~, ord] = sort(modNum);
else
    ord = (1:nF)';
end

%% ============================== PRINT TIMELINE ==============================
fprintf("%-4s  %-19s  %8s  %8s  %s\n", "#", "last modified", "Δ(min)", "size GB", "file");
fprintf("%s\n", repmat('-',1,90));

prev = NaN;
for k = 1:nF
    i = ord(k);
    if isnan(modNum(i))
        fprintf("%-4d  %-19s  %8s  %8s  %s\n", k, "(missing)", "-", "-", names(i));
        continue
    end
    if isnan(prev)
        dMin = 0;
    else
        dMin = (modNum(i) - prev) * 24 * 60;   % datenum days -> minutes
    end
    fprintf("%-4d  %-19s  %8.1f  %8.2f  %s\n", ...
        k, string(modDT(i), 'yyyy-MM-dd HH:mm:ss'), dMin, sizeGB(i), names(i));
    prev = modNum(i);
end

%% ============================== SUMMARY ==============================
valid = ~isnan(modNum);
if any(valid)
    tspanMin = (max(modNum(valid)) - min(modNum(valid))) * 24 * 60;
    fprintf("%s\n", repmat('-',1,90));
    fprintf("First: %s\n", string(min(modDT(valid)), 'yyyy-MM-dd HH:mm:ss'));
    fprintf("Last : %s\n", string(max(modDT(valid)), 'yyyy-MM-dd HH:mm:ss'));
    fprintf("Span : %.1f min  (%.2f h)   |   files: %d   |   total: %.1f GB\n", ...
        tspanMin, tspanMin/60, nF, sum(sizeGB(valid)));
end

%% ============================== OPTIONAL CSV ==============================
if SaveCSV
    T = table((1:nF)', names(ord), modDT(ord), sizeGB(ord), tifList(ord).', ...
        'VariableNames', {'order','file','last_modified','size_GB','full_path'});
    csvPath = fullfile(char(folderPath), sprintf('tiff_timeline_%s.csv', datestr(now,'yymmdd_HHMMSS')));
    try
        writetable(T, csvPath);
        fprintf("\nSaved timeline CSV:\n%s\n", csvPath);
    catch ME
        fprintf("\n[warn] could not write CSV (%s): %s\n", csvPath, ME.message);
    end
end

%% ======================================================================
%% ==================== DISCOVERY (READ-ONLY) ===========================
%% Same resolution logic as Batch_dffQC_260325.m, but find_tifs never moves
%% loose files — it only lists recording TIFFs.
function tifList = discover_tifs_readonly(masterFolder, folderPath)
tifList = string.empty;

if masterFolder ~= ""
    masterFolder = string(masterFolder);
    assert(isfolder(masterFolder), "masterFolder not found: %s", masterFolder);
    allDirs = dir(fullfile(masterFolder, '**', 'phys'));
    physDirs = string.empty;
    for i = 1:numel(allDirs)
        if allDirs(i).isdir
            physDirs(end+1) = string(fullfile(allDirs(i).folder, allDirs(i).name)); %#ok<AGROW>
        end
    end
    scanDirs = string.empty;
    for i = 1:numel(physDirs)
        scanDirs = [scanDirs, resolve_scan_dirs(physDirs(i))]; %#ok<AGROW>
    end
    for i = 1:numel(scanDirs)
        tifList = [tifList, find_tifs_readonly(scanDirs(i))]; %#ok<AGROW>
    end
elseif folderPath ~= ""
    folderPath = string(folderPath);
    assert(isfolder(folderPath), "folderPath not found: %s", folderPath);
    scanDirs = resolve_scan_dirs(folderPath);
    for i = 1:numel(scanDirs)
        tifList = [tifList, find_tifs_readonly(scanDirs(i))]; %#ok<AGROW>
    end
else
    error("Set either masterFolder or folderPath.");
end

tifList = unique(tifList, 'stable');
fprintf("[discover] Found %d TIF files\n", numel(tifList));
end

function scanDirs = resolve_scan_dirs(rootDir)
scanDirs = string.empty;
physDir = fullfile(rootDir, "phys");
if isfolder(physDir)
    procDir = fullfile(physDir, "processed");
    if isfolder(procDir), scanDirs(end+1) = procDir; else, scanDirs(end+1) = physDir; end
    return
end
procDir = fullfile(rootDir, "processed");
if isfolder(procDir), scanDirs(end+1) = procDir; return; end
scanDirs(end+1) = rootDir;
end

function tifList = find_tifs_readonly(scanDir)
% Recursively find recording TIFFs (recording_name/recording_name.tif).
% READ-ONLY: loose .tif files are listed in place, never moved.
tifList = string.empty;
dd = dir(fullfile(scanDir, '**', '*.tif'));
for j = 1:numel(dd)
    if dd(j).isdir, continue; end
    parts = split(string(dd(j).folder), filesep);
    parName = parts(end);
    tifBase = erase(string(dd(j).name), ".tif");
    if string(tifBase) == string(parName)
        tifList(end+1) = string(fullfile(dd(j).folder, dd(j).name)); %#ok<AGROW>
    end
    % else: .tif whose name doesn't match its parent folder -> skip (not a recording)
end
end
