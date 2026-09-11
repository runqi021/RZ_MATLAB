%% 0) Clean up any leftover port objects
try
    if exist('sp','var') && sp.IsOpen,  sp.Close();  end
    if exist('sp','var'),               sp.Dispose(); end
    clear sp
catch
end

%% 1) Open the port fresh
NET.addAssembly('System');
portName = 'COM4';  % ? your laser port
sp = System.IO.Ports.SerialPort();
sp.PortName    = portName;
sp.BaudRate    = 9600;
sp.Parity      = System.IO.Ports.Parity.None;
sp.DataBits    = 8;
sp.StopBits    = System.IO.Ports.StopBits.One;

% **Use ASCII encoding** so MATLAB decodes the bytes correctly
sp.Encoding    = System.Text.Encoding.ASCII;

% **Tell it that commands end in CR (char(13)) only**
sp.NewLine     = char(13);

% Timeouts so ReadLine errors instead of hanging
sp.ReadTimeout = 2000;
sp.WriteTimeout= 2000;

if ~sp.IsOpen, sp.Open(); end

%% 2) Flush any junk
if sp.BytesToRead > 0
    sp.DiscardInBuffer();
end

%% 3) Open the shutter
fprintf('>> SHUTTER OPEN\n');
sp.WriteLine('SHUTTER OPEN');  
pause(0.2);

%% 4) Query its state
fprintf('>> SHUTTER?\n');
sp.WriteLine('SHUTTER?');
pause(0.2);
try
    resp = strtrim(char(sp.ReadLine()));
    fprintf('Laser shutter is: %s\n', resp);
catch ME
    fprintf('No reply (timeout or bad command): %s\n', ME.message);
end

%% 5) Clean up
if sp.IsOpen, sp.Close(); end
sp.Dispose();
clear sp



%%
%% ———— Open & configure port ————
NET.addAssembly('System');
portName = 'COM4';
sp = System.IO.Ports.SerialPort(portName,9600, ...
         System.IO.Ports.Parity.None,8,System.IO.Ports.StopBits.One);
sp.Encoding    = System.Text.Encoding.ASCII;
sp.NewLine     = char(13);
sp.ReadTimeout = 1000;
sp.WriteTimeout= 1000;
if ~sp.IsOpen, sp.Open(); end
if sp.BytesToRead>0, sp.DiscardInBuffer(); end

%% ———— Ask “who are you?” ————
sp.WriteLine('*IDN?');
pause(0.2);
idn = char(sp.ReadExisting());   % read everything the laser sent
disp(['IDN response:  ' idn]);

%% ———— Clean up ————
sp.Close(); sp.Dispose(); clear sp
