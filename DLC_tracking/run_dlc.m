function run_dlc(step, configPath)
clear all; close all; clc;
% run_dlc  —  Run DeepLabCut pipeline steps from MATLAB.
%
% WORKFLOW:
%   run_dlc('full')       — create project + extract frames + open labeling GUI
%   run_dlc('check')      — check labels after labeling
%   run_dlc('train')      — create training set + train network
%   run_dlc('evaluate')   — evaluate trained network
%   run_dlc('analyze')    — run inference on all videos
%   run_dlc('labeled')    — create labeled overlay videos for QC
%   run_dlc('refine')     — extract outlier frames → then run_dlc('label') → retrain
%
% USAGE:
%   run_dlc('full')                           % first time — sets up everything
%   run_dlc('train', 'C:\...\config.yaml')    % after labeling
%   run_dlc('analyze', 'C:\...\config.yaml')  % after training

%% ------------------- USER SETTINGS --------------------------------------
videoDir     = "C:\260407_KA_electro_dorsal_whisking";
projName     = "whisker_dorsal";
experimenter = "RZ";
bodyparts    = "vL1,vL2,vL3,vL4,vR1,vR2,vR3";
%bodyparts    = "vL1,vL2,vL3,vL4,vL5,vL6,vR1,vR2,vR3,vR4,vR5,vR6";

net          = "resnet_101";          % "resnet_50", "resnet_101", "hrnet_w32"
maxiters     = 200000;                % training iterations
saveiters    = 10000;                 % save snapshot every N iters
numframes    = 20;                    % frames to extract per video (kmeans)
algo         = "uniform";             % "kmeans" or "uniform"

%% ------------------- PATHS ----------------------------------------------
condaEnv  = "C:\Users\Admin\.conda\envs\dlc310\python.exe";
scriptDir = fileparts(mfilename('fullpath'));
pyScript  = fullfile(scriptDir, "dlc_pipeline.py");

if nargin < 1, step = 'full'; end

%% ------------------- RUN ------------------------------------------------
switch lower(step)

    case 'full'
        % Create project + extract frames + label GUI
        cmd = sprintf('"%s" "%s" full --videoDir "%s" --projName %s --experimenter %s --bodyparts %s --numframes %d --net %s --algo %s', ...
            condaEnv, pyScript, videoDir, projName, experimenter, bodyparts, numframes, net, algo);
        fprintf("Creating project, extracting frames, opening labeling GUI...\n");
        run_cmd(cmd);
        fprintf("\n>>> After labeling, run:  run_dlc(''train'', ''<config.yaml path>'')\n");

    case 'label'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" label --config "%s"', condaEnv, pyScript, cfg);
        run_cmd(cmd);

    case 'check'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" check --config "%s"', condaEnv, pyScript, cfg);
        run_cmd(cmd);

    case 'train'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" train --config "%s" --net %s --maxiters %d --saveiters %d', ...
            condaEnv, pyScript, cfg, net, maxiters, saveiters);
        fprintf("Training network (this will take a while)...\n");
        run_cmd(cmd);

    case 'evaluate'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" evaluate --config "%s"', condaEnv, pyScript, cfg);
        run_cmd(cmd);

    case 'analyze'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" analyze --config "%s" --videoDir "%s"', ...
            condaEnv, pyScript, cfg, videoDir);
        run_cmd(cmd);

    case 'labeled'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" labeled --config "%s" --videoDir "%s"', ...
            condaEnv, pyScript, cfg, videoDir);
        run_cmd(cmd);

    case 'refine'
        cfg = get_config(configPath, nargin);
        cmd = sprintf('"%s" "%s" refine --config "%s" --videoDir "%s"', ...
            condaEnv, pyScript, cfg, videoDir);
        run_cmd(cmd);

    otherwise
        error("Unknown step '%s'. Use: full, label, check, train, evaluate, analyze, labeled, refine", step);
end

end

%% =========================================================================
function run_cmd(cmd)
    fprintf(">> %s\n\n", cmd);
    status = system(cmd);
    if status ~= 0
        warning("Command exited with status %d", status);
    end
end

function cfg = get_config(configPath, nargs)
    if nargs < 2 || isempty(configPath)
        error("This step requires a config path. Usage: run_dlc('%s', 'C:\\...\\config.yaml')", 'step');
    end
    cfg = configPath;
end
