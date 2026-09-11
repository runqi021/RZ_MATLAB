hCtrl = evalin('base','hSICtl');          % grab the ScanImage controller
hCtrl.hModel.hBeams.powers(1) = 7;      % set beam #1 to 20%

%%
M=[];
%%
wls  = M(:,1);   % wavelengths
pcts = M(:,2);   % percent values