% breath_whisk_overlay_RZ.m
% SIMPLE raw overlay of breathing + LEFT/RIGHT whisking for one session.
% Paired by animal + run. Breathing = _breath.mat (inhale-up). Whisking = raw
% base->tip angle (deg), LEFT mirrored so protraction is +up on both sides.
% Both z-scored onto one axis (units differ); no fancy sync -- each on its own
% time-from-0 (cameras ~co-triggered at 400 Hz).
close all; clc; clear;

%
animal   = "5916300";    % needs a _breath.mat: have 5916296 n4, 5840027 n3
runIdx   = 1;
dataRoot = "D:\260615_thermalNbasler";
WIN      = [];      % display window (s); [] = full clip

whiskDir = fullfile(char(dataRoot),'whisk','260615_whisk-RZ-2026-06-17','videos');
noseDir  = fullfile(char(dataRoot),'nose_thermal','260615_thermal_nose-RZ-2026-06-17','videos');
addpath(fullfile(fileparts(mfilename('fullpath')),'thermal_breathing'));   % thermal_resolve_paths

% ---- breathing ----
noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', char(animal), runIdx));
P = thermal_resolve_paths(noseCsv, dataRoot);
assert(isfile(P.breath), 'no _breath.mat for %s n%d (run drawROI_N_lpsub first)', char(animal), runIdx);
B = load(P.breath);  breath = B.breath(:);  tB = (0:numel(breath)-1)'/double(B.fps);

% ---- whisking (raw angle; y up, LEFT x mirrored -> protraction + on both) ----
wcsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', char(animal), runIdx));
M = dlc_gate_interp(wcsv, 0.6);   % lik<0.6 -> linear interp
La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));   % LEFT
Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));   % RIGHT
tW = (0:numel(La)-1)'/400;

%%
La(size(breath)+1:size(La)) = [];

%%
figure;
plot(breath, La);

%%
% ---- simple z-scored overlay ----
z = @(x) (x - mean(x,'omitnan')) / std(x,'omitnan');
figure('Color','w','Position',[100 120 1250 450]); hold on;
plot(tB, z(breath), 'k-', 'LineWidth',1.3);
plot(tW, z(La), '-');
plot(tW, z(Ra), '-');
grid on; xlabel('s'); ylabel('z-score');
legend({'breath (inhale up)','L whisk','R whisk'}, 'Orientation','horizontal','Location','northoutside');
title(sprintf('%s n%d — breathing vs whisking (raw)', char(animal), runIdx), 'Interpreter','none');
if ~isempty(WIN), xlim(WIN); end

% ---- helper ----
function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(dirPath, [prefix '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bestN = arrayfun(@(x) localbest(x.name), d);
    [~,ix] = max(bestN);
    csv = fullfile(d(ix).folder, d(ix).name);
end
function n = localbest(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
