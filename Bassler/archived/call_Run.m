addpath(genpath("C:\Users\zhang\Desktop\DK\MATLAB\Bassler"));
savepath;

setenv("GENICAM_GENTL64_PATH", "C:\Program Files\DK\TeachingScope\Runtime\x64");
imaqreset;

info = imaqhwinfo('gentl');
info.DeviceInfo

%%
vid = videoinput('gentl',1);
src = getselectedsource(vid);

% ---- set binning on camera ----
if isprop(src,'BinningHorizontal')
    src.BinningHorizontal = 2;
    src.BinningVertical   = 2;
else
    disp("No binning properties exposed via GenTL.");
end

src.ExposureAuto = 'Off';
src.ExposureMode = 'Timed';    % if this property exists
src.ExposureTime = 5000;       % microseconds

preview(vid);   % optional sanity check

%%
stoppreview(vid);

%%
run_basler_n_runs
