%% Auto Acq Script - Snake Pattern Tiling
% Automatically moves stage and acquires tiled images with ScanImage

%% CONFIGURATION - EDIT THESE
numCols = 4;     % Number of tiles along X (columns)
numRows = 8;     % Number of tiles along Y (rows)
xStep = -644.8;   % Step size in X (µm)
yStep = 644.8;   % Step size in Y (µm)
filePrefix = 'cleared_vgat_g8s_-450um_1.2x_930nm_8-32LP_z1500-15_200f';  % Base name

%filePrefix = 'dbx_g8s_1.2x_930nm_6-20lp_z210_-8_400f';  % Base name

%% INIT
addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
hCtrl = evalin('base', 'hSICtl');
hMotors = hCtrl.hModel.hMotors;
hScan2D = hCtrl.hModel.hScan2D;

% Get initial position (X,Y,Z)
origin = hMotors.motorPosition;
x0 = origin(1);
y0 = origin(2);
z0 = origin(3);

%% TILING LOOP
for row = 0:(numRows-1)
    if mod(row,2) == 0
        colRange = 0:(numCols-1);   % left to right
    else
        colRange = (numCols-1):-1:0; % right to left (snake pattern)
    end

    for col = colRange
        % Compute target coordinates
        xTarget = x0 + col * xStep;
        yTarget = y0 + row * yStep;

        % Get current position
        currPos = hMotors.motorPosition;
        xCurr = currPos(1);
        yCurr = currPos(2);

        % Compute delta and step
        dx = (xTarget - xCurr) / 10;
        dy = (yTarget - yCurr) / 10;

        % Move in 10 smaller steps
        for step = 1:10
            intermediateX = xCurr + dx * step;
            intermediateY = yCurr + dy * step;
            hMotors.moveStartRelative([intermediateX, intermediateY, z0]);
            hMotors.moveWaitForFinish([1 2 3]);
            pause(1);
        end

        % Set filename with coordinate info
        pos = hMotors.motorPosition;
        filename = sprintf('%s_col%02d_row%02d_x%.0f_y%.0f', filePrefix, col, row, pos(1), pos(2));
        hScan2D.logFileStem = filename;

        % Start acquisition
        hCtrl.hModel.startGrab();
        while ~strcmp(hCtrl.hModel.acqState, 'idle')  % Wait while acquisition is active
            pause(1);
        end
    end
end
