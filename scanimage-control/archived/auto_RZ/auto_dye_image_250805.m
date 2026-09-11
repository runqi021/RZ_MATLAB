%%
%M=[];

%% -------------------------------
%% 0) Prepare LaserControl & ScanImage
%% -------------------------------
zoom = 4;
addpath('C:\Users\dklab\Desktop\LaserControl-master\code');
rehash toolboxcache;

laserControl.settings.readSettings;  
C = laserControl.chameleon('COM4');  %# your laser port

C.turnOn();
C.openShutter();

hCtrl   = evalin('base','hSICtl');
hScan2D = hCtrl.hModel.hScan2D;

%% -------------------------------
%% 1) Define your wavelength??power lookup
wls  = M(:,1);
pcts = M(:,2);
maxPower_mW = 20;  
beamIdx = 1;  

%% -------------------------------
%% 3) ALL Wavelengths @ GDD = 8000
%% -------------------------------
fixedGDD = 8000;  % fs^2

for idx = 1:numel(wls)
    wl  = wls(idx);
    pct = pcts(idx);
    
    C.setWavelength(wl);
    waitForLaser(C);
    %C.setGDDMode(true);
    %waitForLaser(C);
    C.setGDD(fixedGDD);
    waitForLaser(C);
    
    hCtrl.hModel.hBeams.powers(beamIdx) = pct;
    
    fname = sprintf('n4_%dnm_%dlp_%dmW_%dx_GDD-%d', wl, round(pct), maxPower_mW, zoom, fixedGDD);
    hScan2D.logFileStem = fname;
    
    hCtrl.hModel.startGrab();
    while ~strcmp(hCtrl.hModel.acqState,'idle')
        pause(0.1);
    end
end

%% -------------------------------
%% 4) 930 nm @ GDD = 0:2000:14000
%% -------------------------------
wl = 930;
idx930 = find(wls==wl,1);
pct930 = pcts(idx930);

for GDD = 0:2000:14000
    C.setWavelength(wl);
    waitForLaser(C);
    %C.setGDDMode(true);
    C.setGDD(GDD);
    waitForLaser(C);
    
    hCtrl.hModel.hBeams.powers(beamIdx) = pct930;
    
    fname = sprintf('n4_%dnm_%dlp_%dmW_%dx_GDD-%d', wl, round(pct930), maxPower_mW, zoom, GDD);
    hScan2D.logFileStem = fname;
    
    hCtrl.hModel.startGrab();
    while ~strcmp(hCtrl.hModel.acqState,'idle')
        pause(0.1);
    end
end

%% -------------------------------
%% 5) Clean up LaserControl
%% -------------------------------
%C.closeShutter();
%C.turnOff();

% ------------- PUT THIS AT THE END! -------------
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
