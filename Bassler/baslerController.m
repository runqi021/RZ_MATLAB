classdef baslerController < handle
    properties
        vid;
        src;
        params;
    end
    properties (Dependent)
        logger;
    end
    
    methods
        function obj = baslerController(vid_input)
            %Handles to camera objects
            obj.vid = vid_input;
            obj.src = getselectedsource(obj.vid);
          
            %Parameter handler
            obj.params = baslerProperties(obj);
        end
        
        function val = get.logger(obj)
            try
                val = obj.vid.DiskLogger;
            catch
                val = [];
            end
        end
        function set.logger(obj,val)
            if isa(val,'VideoWriter')
                obj.vid.LoggingMode = 'disk';
                obj.vid.DiskLogger = val;
            else
                warn('Invalid videowriter, defaulting to memory.');
                obj.vid.LoggingMode = 'memory';
            end
        end
        
        %Handle disconnecting camera when controller is destroyed
        function delete(obj)
            delete(obj.vid);
        end
    end
end