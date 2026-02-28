classdef baslerGUI < handle
    properties
        mdl;
        gui;
        isPreviewing;
    end
    
    methods
        function obj = baslerGUI(ctl)
            obj.mdl = ctl;
            obj.isPreviewing = 0;
            
            obj.gui.f = figure('Toolbar','none','Menubar','none','NumberTItle','Off','Position',[100,50,1000,750]);
            
            obj.gui.h = axes('Parent',obj.gui.f,'Position',[0.1,0.29,0.5,0.66]);
            obj.gui.img = image(zeros(obj.mdl.src.AutoFunctionROIWidth,obj.mdl.src.AutoFunctionROIHeight));
            obj.gui.h.XTick = [];
            obj.gui.h.YTick = [];
            
            obj.gui.preview = uicontrol('Style','pushbutton','String','Preview','Position',[150,100,100,25]);
            obj.gui.grab = uicontrol('Style','pushbutton','String','Grab','Position',[450,100,100,25]);
            obj.gui.abort = uicontrol('Style','pushbutton','String','Abort','Position',[450,60,100,25]);
            obj.gui.cdir = uicontrol('Style','pushbutton','String','Dir..','Position',[465,150,75,25]);
            obj.gui.fname = uicontrol('Style','edit','String','file','Position',[150,150,300,25]);
            obj.gui.setROI = uicontrol('Style','pushbutton','String','Set imaging ROI','Position',[225,185,100,25]);
            obj.gui.resetROI = uicontrol('Style','pushbutton','String','Reset imaging ROI','Position',[375,185,100,25]);
            
            obj.gui.params = uipanel(obj.gui.f,'Title','Parameters','Position',[0.65,0.05,0.3,0.9]);
            obj.gui.FRateTxt = uicontrol(obj.gui.params,'Style','text','String','Frame Rate','Units','normalized','Position',[0.1,0.8,0.8,0.025]);
            obj.gui.FRate = uicontrol(obj.gui.params,'Style','edit','String',num2str(obj.mdl.params.AcquisitionFrameRate),'Units','normalized','Position',[0.1,0.775,0.8,0.025]);
            obj.gui.ExpTxt = uicontrol(obj.gui.params,'Style','text','String','Exposure Time','Units','normalized','Position',[0.1,0.7,0.8,0.025]);
            obj.gui.Exp = uicontrol(obj.gui.params,'Style','edit','String',num2str(obj.mdl.params.ExposureTime),'Units','normalized','Position',[0.1,0.675,0.8,0.025]);
            obj.gui.FramesTxt = uicontrol(obj.gui.params,'Style','text','String','Frames/trigger','Units','normalized','Position',[0.1,0.6,0.8,0.025]);
            obj.gui.Frames = uicontrol(obj.gui.params,'Style','edit','String',num2str(obj.mdl.params.FramesPerTrigger),'Units','normalized','Position',[0.1,0.575,0.8,0.025]);
            obj.gui.TriggerTxt = uicontrol(obj.gui.params,'Style','text','String','Trigger type','Units','normalized','Position',[0.1,0.5,0.8,0.025]);
            obj.gui.Trigger = uicontrol(obj.gui.params,'Style','popupmenu','String',{obj.mdl.params.triggers.TriggerType},'Units','normalized','Position',[0.1,0.475,0.8,0.025]);
            obj.gui.loadparms = uicontrol(obj.gui.params,'Style','pushbutton','String','Load params','Units','normalized','Position',[0.1,0.05,0.35,0.05]);
            obj.gui.saveparms = uicontrol(obj.gui.params,'Style','pushbutton','String','Save params','Units','normalized','Position',[0.55,0.05,0.35,0.05]);
            
            obj.gui.bandTtl = uicontrol(obj.gui.f,'Style','text','String','Estimated bandwidth (MB/s)','Position',[275,100,150,25]);
            bw = obj.mdl.params.ROIPosition(3)*obj.mdl.params.ROIPosition(4)*obj.mdl.params.AcquisitionFrameRate/1e6;
            obj.gui.bandwidth = uicontrol(obj.gui.f,'Style','text','String',[num2str(bw),' MB/s'],'Position',[275,60,150,25]);
            
            obj.gui.status = uicontrol(obj.gui.f,'Style','text','String','Not acquiring.','Position',[400,25,200,25]);
            
            obj.gui.f.CloseRequestFcn = {@obj.fClose};
            obj.gui.preview.Callback = @obj.preview;
            obj.gui.grab.Callback = @obj.grab;
            obj.gui.abort.Callback = @obj.abort;
            obj.gui.cdir.Callback = @obj.changedir;
            obj.gui.fname.Callback = @obj.changename;
            obj.gui.FRate.Callback = @obj.changeFRate;
            obj.gui.Exp.Callback = @obj.changeExp;
            obj.gui.Frames.Callback = @obj.changeFrames;
            obj.gui.Trigger.Callback = @obj.changeTrigger;
            obj.gui.setROI.Callback = @obj.setROI;
            obj.gui.resetROI.Callback = @obj.resetROI;
            obj.gui.loadparms.Callback = @obj.loadparams;
            obj.gui.saveparms.Callback = @obj.saveparams;
            
            obj.mdl.vid.TimerPeriod = 2;
            obj.mdl.vid.TimerFcn = @obj.statusUpd;
            obj.mdl.vid.StartFcn = @obj.startNot;
            obj.mdl.vid.StopFcn = @obj.stopNot;
            obj.mdl.vid.UserData = obj.gui.status;
            
            logger = VideoWriter('file.avi','Grayscale Avi');
            assignin('base','logger',logger);
            obj.mdl.logger = evalin('base','logger');
        end
        
        %Callbacks
        function preview(obj,src,evt)
            if obj.isPreviewing
                stoppreview(obj.mdl.vid);
                obj.isPreviewing = 0;
            else
                preview(obj.mdl.vid,obj.gui.img);
                obj.isPreviewing = 1;
            end
        end
        
        function grab(obj,src,evt)
            stoppreview(obj.mdl.vid);          
            obj.changename(obj.gui.fname);
            ss = strsplit(obj.mdl.logger.Filename,'.');
            obj.gui.fname.String = ss{1};
            
            start(obj.mdl.vid);
        end
        
        function abort(obj,src,evt)
            stop(obj.mdl.vid);
        end
        
        %Having the videoinput in an object somehow reduces speed by a
        %lot. This has to do with the way the videowriter is
        %created...we have to assign it in the base workspace...
        function changedir(obj,src,evt)
            path = uigetdir(obj.mdl.logger.Path);
            name = [path,'\',obj.gui.fname.String];
            if exist([name '.avi'],'file') == 2
                name = [name '_new'];
            end
            logger = VideoWriter([name,'.avi'],'Grayscale Avi');
            assignin('base','logger',logger);
            obj.mdl.logger = evalin('base','logger');
        end
        
        function changename(obj,src,evt)
            path = obj.mdl.logger.Path;
            name = [path,'\',src.String];
            if exist([name '.avi'],'file') == 2
                name = [name '_new'];
            end
            logger = VideoWriter([name,'.avi'],'Grayscale Avi');
            assignin('base','logger',logger);
            obj.mdl.logger = evalin('base','logger');
        end
        
        function changeFRate(obj,src,evt)
            try
                obj.mdl.params.AcquisitionFrameRate = str2double(src.String);
                obj.updateBW();
            catch
                src.String = obj.mdl.params.AcquisitionFrameRate;
            end
        end
        
        function changeExp(obj,src,evt)
            try
                obj.mdl.params.ExposureTime = str2double(src.String);
            catch
                src.String = obj.mdl.params.ExposureTime;
            end
        end
        
        function changeFrames(obj,src,evt)
            try
                obj.mdl.params.FramesPerTrigger = str2double(src.String);
            catch
                src.String = obj.mdl.params.FramesPerTrigger;
            end
        end
        
        function changeTrigger(obj,src,evt)
            try
                if strcmp(src.Value.TriggerType,'hardware')
                    obj.mdl.src.TriggerMode = 'On';
                else
                    obj.mdl.src.TriggerMode = 'Off';
                end
                obj.mdl.params.currtriggertype = src.Value;
            catch
            end
        end
        
        function setROI(obj,src,evt)
            try
                axes(obj.gui.h);
                mask = roipoly;
                info = regionprops(mask,'Boundingbox');
                bb = info.BoundingBox;
                bb(2:-1:1) = bb(1:2);
                bb(4:-1:3) = bb(3:4);
                obj.mdl.params.ROIPosition = info.BoundingBox;
                obj.updateBW();
            catch
            end
        end
        
        function resetROI(obj,src,evt)
            obj.mdl.params.ROIPosition = [0,0,obj.mdl.params.maxFrameSize(1),obj.mdl.params.maxFrameSize(2)];
            obj.updateBW();
        end
        
        %Note that I'm not saving the trigger type. So change it yourself.
        function loadparams(obj,src,evt)
            [n,p] = uigetfile('*.mat');
            S = load(fullfile(p,n));
            
            obj.mdl.params.AcquisitionFrameRate = S.props.frate;
            obj.mdl.params.ExposureTime = S.props.exp;
            obj.mdl.params.FramesPerTrigger = S.props.frames;
            obj.mdl.params.ROIPosition = S.props.ROIPosition;
            
            obj.gui.FRate.String = num2str(S.props.frate);
            obj.gui.Exp.String = num2str(S.props.exp);
            obj.gui.Frames.String = num2str(S.props.frames);
            obj.updateBW();
        end
         
        function saveparams(obj,src,evt)
            [n,p] = uiputfile('*.mat','Save as...');
            mfile = matfile(fullfile(p,n),'Writable',true);
            props.frate = obj.mdl.params.AcquisitionFrameRate;
            props.exp = obj.mdl.params.ExposureTime;
            props.frames = obj.mdl.params.FramesPerTrigger;
            props.ROIPosition = obj.mdl.params.ROIPosition;
            mfile.props = props;
        end
        
        function startNot(obj,vid,~)
            z = vid.UserData;
            z.String = 'Acquiring : ';
        end
        
        function statusUpd(obj,vid,~)
            z = vid.UserData;
            z.String = ['Acquiring : ', num2str(vid.FramesAcquired - vid.DiskLoggerFrameCount), ' frames in buffer.'];
        end
        
        function stopNot(obj,vid,~)
            z = vid.UserData;
            z.String = 'Not acquiring.';
        end
        
        %----------------------------------
        function updateBW(obj)
            bw = obj.mdl.params.BytesPerFrame*obj.mdl.params.AcquisitionFrameRate/1e6;
            obj.gui.bandwidth.String = [num2str(bw),' MB/s'];
        end
        
        %Handle the closing of the figure;
        function fClose(obj,src,evt)
            try
                closepreview(obj.mdl.vid);
                delete(obj.gui.f);
            catch
                delete(obj.gui.f);
            end
        end
    end
end