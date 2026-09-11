addpath('C:\Users\dklab\Desktop\LaserControl-master\code')
rehash toolboxcache
which laserControl.loghandler -all


%%
% 2) (Re)generate the settings file if you haven’t yet—
%    this lets you point each class at the right COM port.
laserControl.settings.readSettings;

%%
% 3) Create the Chameleon object (replace 'COM4' with your port)
C = laserControl.chameleon('COM4');

% 4) Make sure the laser is on
C.turnOn;

% 5) Open its shutter
C.openShutter;

%%
C.setWavelength(930);   % Set wavelength to 920 nm

%%
%C.setGDDMode(true);     % Switch to manual GDD mode
%%
C.setGDD(8000);         % Set GDD to 2000 fs^2

%%
while true
    [ready, msg] = C.isReady();
    tuning = C.isTuning();
    if ready && ~tuning
        disp('Laser is ready!');
        break;
    else
        if tuning
            fprintf('Waiting for laser: currently tuning wavelength...\n');
        else
            fprintf('Waiting for laser: %s\n', msg);
        end
        pause(2);
    end
end

% --- Now print tuned parameters ---
lambda = C.readWavelength();
fprintf('Current wavelength: %.1f nm\n', lambda);

if ismethod(C, 'readGDD')
    gdd = C.readGDD();
    fprintf('Current GDD: %.1f fs^2\n', gdd);
else
    disp('GDD query not implemented for this laser class.');
end

