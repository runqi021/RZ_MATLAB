classdef chameleon < laserControl.laser
%%  chameleon - control class for Coherent Chameleon lasers
%
% Example
% C = chameleon('COM1');
%
% Rob Campbell - Basel 2017
% No longer loghandler, RUNQI as 250805
    properties
       faultMessage = ...
            {'Laser head interlock', 'External interlock', ...
              'PS cover interlock', 'LBO temperature', ...
              'LBO not locked at set temperature', 'Vanadate temperature', ...
              'Etalon temperature', 'Diode 1 temperature', ...
              'Diode 2 temperature', 'Baseplate temperature', ...
              'Heatsink 1 temperature', 'Heatsink 2 temperature', ...
              '', '', '', ...
              'Diode 1 over-current', 'Diode 2 over-current', ...
              'Over-current', 'Diode 1 under-voltage', ...
              'Diode 2 under-voltage', 'Diode 1 over-voltage', ...
              'Diode 2 over-voltage', '', '', 'Diode 1 EEPROM', ...
              'Diode 2 EEPROM', 'Laser head EEPROM', 'PS EEPROM',...
              'PS-head mismatch', 'LBO battery', 'Shutter state mismatch', ...
              'CPU PROM checksum', 'Head PROM checksum', ...
              'Diode 1 PROM checksum', 'Diode 2 PROM checksum', ...
              'CPU PROM range','Head PROM range', 'Diode 1PROM range', ...
              'Diode 2 PROM range', 'Head-diode mismatch', '', '', ...
              'Lost modelock', '', '', '', 'Ti-Sapph temperature', '', ...
              'PZT X', 'Cavity humidity', 'Tuning stepper motor homing', ...
              'Lasing', 'Laser failed to begin modelocking', 'Headboard comms', ...
              'System lasing', 'PS-head EEPROM mismatch', ...
              'Modelock slit stepper motor homing', 'Chameleon-verdi EEPROM', ...
              'Chameleon precompensator homing', 'Chameleon curve EEPROM'};
    end

    methods

        % Constructor (second arg is ignored, for backward compatibility)
        function obj = chameleon(serialComms, varargin)
            if nargin<1
                error('chameleon requires at least one input argument: you must supply the laser COM port as a string')
            end

            obj.maxWavelength=1100;
            obj.minWavelength=700;
            obj.friendlyName = 'Chameleon';

            fprintf('\nSetting up Chameleon laser communication on serial port %s\n', serialComms);
            laserControl.clearSerial(serialComms)
            obj.controllerID=serialComms;
            success = obj.connect;

            if ~success
                fprintf('Component chameleon failed to connect to laser over the serial port.\n')
                return
            end

            obj.targetWavelength=obj.currentWavelength;
            fprintf('Connected to Chameleon laser on serial port %s\n\n', serialComms)
        end

        function delete(obj)
            fprintf('Disconnecting from Chameleon laser\n')
            if ~isempty(obj.hC) && isa(obj.hC,'serial') && isvalid(obj.hC)
                fprintf('Closing serial communications with Chameleon laser\n')
                flushinput(obj.hC)
                fclose(obj.hC);
                delete(obj.hC);
            end
        end

        function success = connect(obj)
            obj.hC=serial(obj.controllerID,'BaudRate',19200, ...
                        'TimeOut',5, ...
                        'Terminator', 'CR/LF');

            try
                fopen(obj.hC);
            catch ME
                fprintf(' * ERROR: Failed to connect to Chameleon:\n%s\n\n', ME.message)
                success=false;
                return
            end

            flushinput(obj.hC)
            success = false;
            if ~isempty(obj.hC)
                s1 = obj.sendAndReceiveSerial('ECHO=0');
                s2 = obj.sendAndReceiveSerial('PROMPT=0');
                s3 = obj.setWatchDogTimer(0);

                if s1==1 && s2==1 && s3==1
                    success=true;
                else
                    [~,s] = obj.isShutterOpen;
                    if s==true
                        success=true;
                    else
                        fprintf('Failed to communicate with Chameleon laser\n')
                        success=false;
                    end
                end
            end
            obj.isLaserConnected=success;
        end

        function success = isControllerConnected(obj)
            if strcmp(obj.hC.Status,'closed')
                success=false;
            else
                [~,success] = obj.isShutterOpen;
            end
            obj.isLaserConnected=success;
        end

        function success = turnOn(obj)
            success=false;

            if ~obj.readKeySwitch
                fprintf('Key switch is set to "STANDBY". Can not turn on Chameleon laser.\n')
                obj.isLaserOn=success;
                return
            end

            faultInd = obj.readFaultState;
            if length(faultInd)==1 && faultInd==0
                success=obj.sendAndReceiveSerial('L=1');
            end

            obj.isLaserOn=success;
        end

        function success = turnOff(obj)
            success=obj.sendAndReceiveSerial('L=0');
            if success
                obj.isLaserOn=false;
            end
        end
        
        function success = setGDDMode(obj, useManual)
            % useManual = true for manual GDD, false for lookup/auto GDD
            if useManual
                mode = 1;
            else
                mode = 0;
            end
            cmd = sprintf('GDD_MODE=%d', mode);   % <-- Check your manual! 
            [success, ~] = obj.sendAndReceiveSerial(cmd, false);
        end

        function success = setGDD(obj, value)
            % Set manual GDD value (in fs^2 or as required by your hardware)
            cmd = sprintf('GDD=%d', round(value));   % <-- Check your manual for correct syntax!
            [success, ~] = obj.sendAndReceiveSerial(cmd, false);
        end

        function gdd = readGDD(obj)
            [success, gdd] = obj.sendAndReceiveSerial('?GDD');
            if ~success
                gdd = [];
                return
            end
            gdd = str2double(gdd);
        end
                
        function [powerOnState,reply] = isPoweredOn(obj)
            [success,reply]=obj.sendAndReceiveSerial('?L');
            if ~success
                powerOnState=0;
                return
            end
            powerOnState = (str2double(reply)==1);
            obj.isLaserOn=powerOnState;
        end

        function [laserReady,msg] = isReady(obj)
            laserReady = false;
            msg='';
            [shutterState,success] = obj.isShutterOpen;

            if ~success
                msg='No connection to laser';
                obj.isLaserReady=false;
                return
            end
            if ~obj.isPoweredOn
                msg='Laser not powered on';
                obj.isLaserReady=false;
                return
            end
            if shutterState==0
                msg='Laser shutter is closed';
                obj.isLaserReady=false;
                return
            end
            if ~obj.isModeLocked
                msg='Laser not modelocked';
                obj.isLaserReady=false;
                return
            end

            laserReady=true;
            obj.isLaserReady=laserReady;
        end

        function modelockState = isModeLocked(obj)
            [success,reply]=obj.sendAndReceiveSerial('?MDLK');
            if ~success
                modelockState=false;
                obj.isLaserModeLocked=modelock;
                return
            end

            modelockState = str2double(reply);
            modelockState = (modelockState==1);
            obj.isLaserModeLocked=modelockState;
        end

        function success = openShutter(obj)
            success=obj.sendAndReceiveSerial('SHUTTER=1');
            pause(0.75)
            if success
                obj.isLaserShutterOpen=true;
            end
        end

        function success = closeShutter(obj)
            success=obj.sendAndReceiveSerial('SHUTTER=0');
            pause(0.75)
            if success
                obj.isLaserShutterOpen=false;
            end
        end

        function [shutterState,success] = isShutterOpen(obj)
            [success,reply]=obj.sendAndReceiveSerial('?S');
            if ~success
                shutterState=[];
                return
            end
            shutterState = str2double(reply);
            obj.isLaserShutterOpen=shutterState;
        end

        function wavelength = readWavelength(obj)
            [success,wavelength]=obj.sendAndReceiveSerial('?VW');
            if ~success
                wavelength=[];
                return
            end
            wavelength = str2double(wavelength);
            if ~isnan(wavelength)
                obj.currentWavelength=wavelength;
            else
                fprintf('Failed to read wavelength from Chameleon. Likely laser is tuning.\n')
            end
        end

        function success = setWavelength(obj,wavelengthInNM)
            success=false;
            if length(wavelengthInNM)>1
                fprintf('wavelength should be a scalar')
                return
            end
            if ~obj.isTargetWavelengthInRange(wavelengthInNM)
                return
            end
            cmd = sprintf('WAVELENGTH=%d', round(wavelengthInNM));
            [success,wavelength]=obj.sendAndReceiveSerial(cmd,false);
            if ~success
                return
            end
            obj.currentWavelength=wavelength;
            obj.targetWavelength=wavelengthInNM;
        end

        function tuning = isTuning(obj)
            [success,reply]=obj.sendAndReceiveSerial('?TS');
            if ~success
                tuning=nan;
                return
            end

            reply = str2double(reply);
            tuning = reply > 0;
        end

        function laserPower = readPower(obj)
            [success,laserPower]=obj.sendAndReceiveSerial('?UF');
            if ~success
                laserPower=[];
                return
            end
            laserPower = str2double(laserPower);
        end

        function laserID = readLaserID(obj)
            [success,laserID]=obj.sendAndReceiveSerial('?SN');
            if ~success
                laserID=[];
                return
            end
            laserID = ['Chameleon, Serial Number: ', laserID];
        end

        function laserStats = returnLaserStats(obj)
            lambda = obj.readWavelength;
            outputPower = obj.readPower;
            humidity = obj.readHumidity;

            laserStats=sprintf('wavelength=%dnm,outputPower=%dmW,humidity=%0.1f', ...
                lambda,outputPower,humidity);
        end

        function success=setWatchDogTimer(obj,value)
            if value <= 0
                [success,~] = obj.sendAndReceiveSerial('HB=0');
                return
            else
                [success,~] = obj.sendAndReceiveSerial('HB=1');
                if ~success
                    return
                end
                if value>100
                    value=100;
                elseif value<1
                    value=1;
                end
                value = num2str(round(value));
                [success,~] = obj.sendAndReceiveSerial(['HBR=',value]);
            end
        end

        function laserHumidity = readHumidity(obj)
            [success,laserHumidity]=obj.sendAndReceiveSerial('?RH');
            if ~success
                laserHumidity=[];
                return
            end
            laserHumidity = str2double(laserHumidity);
        end

        function warmedUpValue = readWarmedUp(obj)
            [success,warmedUpValue]=obj.sendAndReceiveSerial('?ST');
            if ~success
                warmedUpValue=[];
                return
            end
            if strfind(warmedUpValue,'OK')
                warmedUpValue=true;
            else
                warmedUpValue=false;
            end
        end

        function keyState = readKeySwitch(obj)
            [success,reply] = obj.sendAndReceiveSerial('?K');
            if ~success
                keyState=[];
                return
            end
            keyState = str2double(reply);
        end

        function [faultNumbers,faultStateString] = readFaultState(obj)
            [success,faultNumbers]=obj.sendAndReceiveSerial('?F');
            if ~success
                faultNumbers=[];
                faultStateString = '';
                return
            end
            faultNumbers = cellfun(@str2double, strsplit(faultNumbers, '&'));
            if length(faultNumbers)==1 && faultNumbers(1) == 0
                faultStateString = 'no faults';
                return
            end
            faultMSG = cell(1,length(faultNumbers));
            for ii=1:length(faultNumbers)
                if faultNumbers(ii) > length(obj.faultMessage) || isempty(obj.faultMessage{faultNumbers(ii)})
                    faultMSG{ii} = sprintf('Unknown fault state: %d. ', faultNumbers(ii));
                else
                    faultMSG{ii} = sprintf('Fault code %d: %s fault. ', ...
                        faultNumbers(ii), obj.faultMessage{faultNumbers(ii)});
                end
            end
            if length(faultNumbers) == 1
                fprintf('Laser reports 1 error:\n')
            elseif length(faultNumbers) > 1
                fprintf('Laser reports %d errors:\n', length(faultNumbers))
            end
            cellfun(@(x) fprintf('%s\n',x), faultMSG)
            faultStateString = [faultMSG{:}];
        end

        function [success,reply]=sendAndReceiveSerial(obj,commandString,waitForReply)
            if nargin<3
                waitForReply=true;
            end

            if isempty(commandString) || ~ischar(commandString)
                reply='';
                success=false;
                % logMessage removed
                return
            end

            fprintf(obj.hC,commandString);

            if ~waitForReply
                reply=[];
                success=true;
                if obj.hC.BytesAvailable>0
                    fprintf('Not waiting for reply by there are %d BytesAvailable\n',obj.hC.BytesAvailable)
                end
                return
            end

            reply=fgets(obj.hC);
            doFlush=1;
            if obj.hC.BytesAvailable>0
                if doFlush
                    fprintf('Read in from the Chameleon buffer using command "%s" but there are still %d BytesAvailable. Flushing.\n', ...
                        commandString, obj.hC.BytesAvailable)
                    flushinput(obj.hC)
                else
                    fprintf('Read in from the Chameleon buffer using command "%s" but there are still %d BytesAvailable. NOT FLUSHING.\n', ...
                        commandString, obj.hC.BytesAvailable)
                end
            end

            if ~isempty(reply)
                reply(end)=[];
            else
                % logMessage removed
                success=false;
                return
            end

            reply = strrep(reply,commandString,'');

            success=true;
        end

    end % methods

end % classdef
