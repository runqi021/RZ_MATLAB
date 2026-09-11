%% Auto Acq Script - Snake Pattern Tiling (Simple Step-Per-Tile Logic)
% Moves stage by 10 steps (X per image, Y per row) and acquires images.

%% CONFIGURATION
numCols = 4;     % Number of tiles along X (columns)
numRows = 8;     % Number of tiles along Y (rows)
pauseTime = 1;   % Pause after each motor step (s)
filePrefix = 'dbx_g8s_1.2x_930nm_6-20lp_z210_-8_400f';  % Base filename

%% INIT
addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
hCtrl = evalin('base', 'hSICtl');
hMotors = hCtrl.hModel.hMotors;
hScan2D = hCtrl.hModel.hScan2D;

xStepSign = -1;  % -1 for left, +1 for right
yStepSign = 1;

for row = 0:(numRows-1)
    if mod(row,2) == 0
        xDir = xStepSign;
    else
        xDir = -xStepSign;
    end

    for col = 0:(numCols-1)
        % Move 10 X steps
        for i = 1:10
            hCtrl.motorStepPosition(1, xDir);
            pause(1);
        end

        % Set filename
        filename = sprintf('%s_x%02d_y%02d', filePrefix, col, row);
        hScan2D.logFileStem = filename;

        % Start acquisition
        hCtrl.hModel.startGrab();
        while ~strcmp(hCtrl.hModel.acqState, 'idle')
            pause(1);
        end
    end

    % Move 10 Y steps after each row
    for j = 1:10
        hCtrl.motorStepPosition(2, yStepSign);
        pause(1);
    end
end
