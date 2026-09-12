function galil_encoder_probe_260911()
%GALIL_ENCODER_PROBE_260911  Ask the Galil directly whether x/y have encoders.
%
% READ-ONLY. Sends only query commands (ID, TP, TD, RP, and MG of internal
% variables). NOTHING here commands motion, changes a setting, or writes to the
% controller. No acquisition, no laser, no scanner.
%
% WHY: hSI.hMotors.motorPosition is Galil TD (a step-pulse count), and the
% driver only ever queries axes A,B,C (DMC4040.m:68-69) even though a DMC-4040
% has four. A stage_backlash run showed TP moving on C (z, 0.1 um encoder) and
% sitting at exactly 0 on A and B. That zero has at least four explanations and
% they are NOT equivalent:
%     1. no encoder on those axes at all
%     2. encoder on the motor but not cabled to this controller
%     3. cabled but the controller axis is not configured to count it
%     4. wired to axis D, which the driver never reads
% This distinguishes them.
%
% WHAT TO LOOK AT
%   _MTx  motor type.  +-2 / +-2.5 = stepper, +-1 = servo.  A stepper axis can
%         still have an encoder; motor type alone does not decide it.
%   _CEx  configure-encoder.  Encodes main+aux encoder type (quadrature, pulse-
%         dir, reversed). Non-default values imply someone configured one.
%   TP    main encoder position.  TD on a stepper axis = step count.
%   Move the stage BY HAND (or with the GUI) between the two passes: an encoder
%   that is wired and counting will change even with the motor idle. That is
%   the single most decisive test, and it needs no commanded motion at all.
%
% CAUTION: this shares the serial link with ScanImage. Each query follows the
% driver's own send/read/confirm pattern so the buffer stays in sync, but do not
% run it during an acquisition.
%
% Runqi Zhang / 2026-09-11.  NOT YET RUN.

hSI = evalin('base','hSI');
try
    hLSC = hSI.hMotors.hMotor(1).hLSC;
catch
    error('galil_probe:noLSC','could not reach hSI.hMotors.hMotor(1).hLSC');
end

fprintf('\n===== DRIVER CONFIGURATION =====\n');
fprintf('numDeviceDimensions : %d\n', hLSC.numDeviceDimensions);
fprintf('activeAxes          : %s   <- axes the driver ever queries\n', ...
        strjoin(hLSC.activeAxes,','));
fprintf('positionDeviceUnits : [%g %g %g] m  = [%.5f %.5f %.5f] um/count\n', ...
        hLSC.positionDeviceUnits, hLSC.positionDeviceUnits*1e6);

fprintf('\n===== CONTROLLER IDENTITY =====\n');
try
    fprintf('%s\n', strtrim(hLSC.infoHardware));
catch ME
    fprintf('could not read ID: %s\n', ME.message);
end

fprintf('\n===== PER-AXIS INTERNALS (A B C D) =====\n');
fprintf('%-6s %12s %12s %12s %12s\n','var','A','B','C','D');
for v = {'MT','CE','TP','TD','RP','TE'}
    row = cell(1,4);
    for k = 1:4
        ax = char('A'+k-1);
        row{k} = q(hLSC, sprintf('MG _%s%c', v{1}, ax));
    end
    fprintf('%-6s %12s %12s %12s %12s\n', v{1}, row{:});
end

fprintf('\n===== RAW TP / TD ON ALL FOUR AXES =====\n');
fprintf('TP ABCD : %s\n', q(hLSC,'TP ABCD'));
fprintf('TD ABCD : %s\n', q(hLSC,'TD ABCD'));

fprintf('\n===== DECISIVE TEST =====\n');
fprintf(['Now push the stage in x and y BY HAND (or jog with the Motor Controls\n' ...
         'arrows), then run this again.\n' ...
         '  TP on A/B CHANGES  -> an encoder IS wired and counting. It is only\n' ...
         '                        that ScanImage never reads it.\n' ...
         '  TP on A/B STAYS 0  -> nothing is counting on those inputs.\n' ...
         'Hand motion is the cleanest version of this: the motor is idle, so\n' ...
         'anything that moves can only be a real position sensor.\n']);
fprintf('\n_MT legend: +-1 servo, +-2 stepper(active low), +-2.5 stepper(active high)\n');
fprintf('_CE legend: 0 is the default; non-zero means an encoder type was set.\n');
end

function s = q(hLSC, cmd)
% One query, using the driver's own send/read/confirm sequence so the serial
% buffer is left exactly as it was found.
s = '?';
try
    hLSC.hRS232.sendCommand([cmd ';']);
    r = hLSC.hRS232.readStringRaw;
    try, hLSC.readCommandConfirmation(); catch, end
    r = strrep(r, [cmd ';'], '');      % strip any command echo
    r = strtrim(strrep(r, ':', ''));
    if ~isempty(r), s = r; end
catch
    s = 'ERR';
end
end
