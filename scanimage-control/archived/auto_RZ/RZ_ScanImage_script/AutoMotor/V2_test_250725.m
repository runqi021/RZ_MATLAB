%% Auto Acq Script - Snake Pattern Tiling with 10-Step Moves
numCols = 3;
numRows = 3;
xStep = -644.8;
yStep = 644.8;
pauseTime = 0.5;  % short pause for motor stabilization
filePrefix = 'dbx_g8s_1.2x_930nm_6-20lp_z210_-8_400f';

addpath(genpath('C:\Users\dklab\Desktop\RZ_ScanImage_script'));
hCtrl = evalin('base', 'hSICtl');
hMotors = hCtrl.hModel.hMotors;
hScan2D = hCtrl.hModel.hScan2D;

% Initial position
origin = hMotors.motorPosition;
x = origin(1);
y = origin(2);
z = origin(3);

%% Tiling Loop
for row = 0:(numRows-1)
    if mod(row,2) == 0
        colRange = 0:(numCols-1);
        xDir = sign(xStep);
    else
        colRange = (numCols-1):-1:0;
        xDir = -sign(xStep);
    end

    for col = colRange
        % Move in X direction by total step (divided into 10 steps)
        stepSize = xDir * abs(xStep)/10;
        for i = 1:10
            hMotors.moveStartRelative([stepSize 0 0]);
            hMotors.moveWaitForFinish([1 2 3]);
            pause(pauseTime);
        end
        x = x + xDir * abs(xStep);  % Update logical X position

        % Set filename with tile indices + stage position
        filename = sprintf('%s_col%02d_row%02d_x%.0f_y%.0f', filePrefix, col, row, x, y);
        hScan2D.logFileStem = filename;

        % Acquire image
        hCtrl.hModel.startGrab();
        while ~strcmp(hCtrl.hModel.acqState, 'idle')
            pause(0.1);
        end
    end

    % After each row, move Y by 10 steps
    stepSizeY = sign(yStep) * abs(yStep)/10;
    for j = 1:10
        hMotors.moveStartRelative([0 stepSizeY 0]);
        hMotors.moveWaitForFinish([1 2 3]);
        pause(pauseTime);
    end
    y = y + sign(yStep) * abs(yStep);  % Update logical Y position
end
