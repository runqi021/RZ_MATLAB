%%
%% Auto Acq Script - Snake Pattern Tiling
% Automatically moves stage and acquires tiled images with ScanImage
% ONE ROUND OF ACQ. SHUTTER OFF AFTER
%% CONFIGURATION - EDIT THESE
numCols = 9;     % Number of tiles along X (columns)
numRows = 9;     % Number of tiles along Y (rows)
xStep = 400;   % Step size in X (µm)
yStep = 400;   % Step size in Y (µm)
%filePrefix = 'autoOFF_test_01';
%filePrefix = 'cleared_vgat_g8s_-450um_1.2x_930nm_8-32LP_z1500-15_200f';  % Base name
filePrefix = '260909_ChAT_g8m_shiverer';  % Base name

%%
addpath('C:\Users\dklab\Desktop\LaserControl-master\code');
rehash toolboxcache;

laserControl.settings.readSettings;  
C = laserControl.chameleon('COM4');  %# your laser port
C.turnOn();
C.openShutter();
C.setWavelength(930);
waitForLaser(C);
C.setGDD(12500)
waitForLaser(C);
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

%while strcmp(hCtrl.hModel.acqState,'idle'), C.closeShutter(); C.turnOff(); end

%%
function waitForLaser(C)
    while true
        [ready, msg] = C.isReady();
        tuning       = C.isTuning();
        lambda = C.readWavelength();
        if ismethod(C,'readGDD')
            gdd = C.readGDD();
        else
            gdd = NaN;
        end
        
        if ready && ~tuning
            fprintf('Laser is ready! Wavelength = %.1f nm', lambda);
            if ~isnan(gdd)
                fprintf(', GDD = %.0f fs^2', gdd);
            end
            fprintf('\n');
            return;
        end
        
        if tuning
            fprintf('  …waiting: tuning Wavelength/GDD (now %.1f nm', lambda);
            if ~isnan(gdd)
                fprintf(', GDD %.0f fs^2', gdd);
            end
            fprintf(')...\n');
        else
            fprintf('  …waiting: %s (%.1f nm', msg, lambda);
            if ~isnan(gdd)
                fprintf(', GDD %.0f fs^2', gdd);
            end
            fprintf(')\n');
        end
        pause(1);
    end
end
