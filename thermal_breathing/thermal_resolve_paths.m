function P = thermal_resolve_paths(dlcCsv, dataRoot)
% THERMAL_RESOLVE_PATHS  Resolve all thermal paths from JUST the DLC csv.
%   P = thermal_resolve_paths(dlcCsv, dataRoot)
%
% Mapping (matches thermal_nostril_breath_single.py):
%   animal = csv name prefix (before first '_')
%   _n#    = run index into chronologically-sorted DATA_ROOT/<animal>/cam1_* folders
%            (n1 = 1st cam1 run, n2 = 2nd, ...)
% Returns struct P with: animal, k, folder, ats, nostrilC, nostrilROI.

    [~, base] = fileparts(char(dlcCsv));
    animal = extractBefore(base, '_');
    tok = regexp(base, '_n(\d+)', 'tokens');
    assert(~isempty(tok), 'cannot parse run index (_n#) from "%s"', base);
    k = str2double(tok{end}{1});

    runs = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    runs = runs([runs.isdir]);
    assert(~isempty(runs), 'no cam1_* folders under %s', fullfile(char(dataRoot), animal));
    [~, ord] = sort({runs.name});           % chronological (timestamp in name)
    runs = runs(ord);
    assert(k <= numel(runs), 'csv says _n%d but only %d cam1 runs for %s', k, numel(runs), animal);
    folder = fullfile(runs(k).folder, runs(k).name);

    a = dir(fullfile(folder, 'Rec-*.ats'));
    assert(~isempty(a), 'no Rec-*.ats in %s', folder);
    stem = fullfile(folder, erase(a(1).name, '.ats'));

    P.animal     = animal;
    P.k          = k;
    P.folder     = folder;
    P.ats        = [stem '.ats'];
    P.nostrilC   = [stem '_nostrilC.mat'];
    P.nostrilROI = [stem '_nostrilROI.mat'];
    P.breath     = [stem '_breath.mat'];     % final averaged breathing signals only
end
