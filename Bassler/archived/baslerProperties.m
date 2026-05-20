classdef baslerProperties < handle
    properties
        mdl;
        triggers;
        maxFrameSize;
    end
    properties (Dependent)
        AcquisitionFrameRate;
        ExposureTime;
        
        FramesPerTrigger;
        ROIPosition;
        currtriggertype;
        
        BytesPerFrame;
    end
    
    methods
        function obj = baslerProperties(mdl)
            obj.mdl = mdl;
            
            %Load parameters from camera
            obj.AcquisitionFrameRate = mdl.src.AcquisitionFrameRate;
            
            obj.triggers = triggerinfo(mdl.vid);
            obj.FramesPerTrigger = mdl.vid.FramesPerTrigger;
            obj.ROIPosition = mdl.vid.ROIPosition;
            
            obj.maxFrameSize = [obj.mdl.src.AutoFunctionROIWidth,obj.mdl.src.AutoFunctionROIHeight];
            
            %Default parameter set for fast acquisition
            triggerconfig(mdl.vid,obj.triggers(1)); %Change this to 3 for hardware trigger
            mdl.src.AcquisitionFrameRateEnable = 'True';
            mdl.src.SensorReadoutMode = 'Fast';
            mdl.src.ExposureTime = 600;
            mdl.src.ExposureAuto = 'Off';
            mdl.src.ExposureMode = 'Timed';
            mdl.src.DeviceLinkThroughputLimit = 419430400;
        end
        
        %AcquisitionFrameRate
        function val = get.AcquisitionFrameRate(obj)
            val = obj.mdl.src.AcquisitionFrameRate;
        end
        function set.AcquisitionFrameRate(obj,val)
            try
                obj.mdl.src.AcquisitionFrameRate = val;
            catch
                obj.propSetError();
                obj.mdl.src.AcquistionFrameRate = 60;
            end
        end
        
        %ExposureTime
        function val = get.ExposureTime(obj)
            val = obj.mdl.src.ExposureTime;
        end
        function set.ExposureTime(obj,val)
            try
                obj.mdl.src.ExposureTime = val;
            catch
                obj.propSetError();
                obj.mdl.src.ExposureTime = 5000;
            end
        end
        
        %FramesPerTrigger
        function val = get.FramesPerTrigger(obj)
            val = obj.mdl.vid.FramesPerTrigger;
        end
        function set.FramesPerTrigger(obj,val)
            try
                obj.mdl.vid.FramesPerTrigger = val;
            catch
                obj.propSetError();
                obj.mdl.vid.FramesPerTrigger = 1;
            end
        end
        
        %ROIPosition
        function val = get.ROIPosition(obj)
            val = obj.mdl.vid.ROIPosition;
        end
        function set.ROIPosition(obj,val)
            try
                obj.mdl.vid.ROIPosition = double(val);
            catch
                obj.propSetError();
            end
        end
        
        %currtriggertype
        function val = get.currtriggertype(obj)
            val = triggerconfig(obj.mdl.vid);
        end
        function set.currtriggertype(obj,val)
            try
                triggerconfig(obj.mdl.vid,obj.triggers(val));
            catch
                obj.propSetError();
            end
        end
        
        %BytesPerFrame
        function val = get.BytesPerFrame(obj)
            val = obj.ROIPosition(3)*obj.ROIPosition(4);
        end
        function set.BytesPerFrame(obj,val)
            warning('Cannot set this! Change ROIPosition instead.');
        end
        
        function propSetError(obj)
            warning('Error setting proprety; setting to default value');
        end
    end
end
        