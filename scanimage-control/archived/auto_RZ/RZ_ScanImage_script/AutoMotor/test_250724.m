%%
addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
%%
evalin('base', 'whos')

%%
hCtrl = evalin('base', 'hSICtl');  % Get the ScanImage controller

%% MOVING MOTOR
% Set step size for X axis to 100 µm
hCtrl.motorStepSize(1) = 100;
%%
% Move +X three times, with 0.5s pause between
for i = 1:10
    hCtrl.motorStepPosition(1, 1);  % Step +1 on X axis
    pause(1);                     % Wait for move to complete
end

%%
hCtrl.motorStepPosition(1, 3);  % Step +1 on X axis
% (step, axis) axis - 1:x, 2:y, 3:z


%% Actually moiving Z to set Z start and stop
%% hacking motor control v5 GUI
hCtrl = evalin('base', 'hSICtl');  % Get controller

%% Step 1: Move to stack START position
hCtrl.hModel.hStackManager.stackZStartPos = 100;   % µm
hCtrl.hModel.hStackManager.stackZEndPos   = 200;   % µm

%%
hCtrl.hModel.hStackManager.stackReturnHome = true; % return to start Z after stack
hCtrl.hModel.hStackManager.stackStartCentered = false; % use absolute startZ
%% Step 3: Set Z step and number of slices
hCtrl.hModel.hStackManager.stackZStepSize = 5;
hCtrl.hModel.hStackManager.numSlices      = 21;

%%
hCtrl.hModel.hStackManager.stackReturnHome = true; % return to start Z after stack
hCtrl.hModel.hStackManager.stackStartCentered = false; % use absolute startZ

%% Power control, hacking power control v4 gui
hBeams = hCtrl.hModel.hBeams;

hBeams.stackStartPower = 20;
hBeams.stackEndPower   = 40;

hBeams.stackUserOverrideLz = true;
hBeams.stackUseStartPower  = true;

hCtrl.pzGUI()
hBeams.pzAdjust = true;


%% manipulate Chanel control GUI v4
hChannels = hCtrl.hModel.hChannels;

% Enable display and save for Channel 1 only
hChannels.channelDisplay = 1;
hChannels.channelSave = 1;


%% Main Control
hCtrl.hModel.hStackManager.framesPerSlice = 500;


%% set DIR and file name
hScan2D = hCtrl.hModel.hScan2D;

% Set logging folder and filename
hScan2D.logFilePath = 'C:\Users\dklab\Desktop\RZ_ScanImage_script\test_250724';
hScan2D.logFileStem = 'Trial1_MouseX';

%%
hCtrl.hModel.startGrab();

%%
pos = hCtrl.hModel.hMotors.motorPosition;
fprintf('Current position: X = %.1f µm, Y = %.1f µm, Z = %.1f µm\n', pos(1), pos(2), pos(3));


