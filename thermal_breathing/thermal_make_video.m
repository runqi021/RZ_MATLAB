function thermal_make_video(matPath, playFps)
% thermal_make_video  Write a plain .mp4 of a thermal *_thermal.mat so you can
% watch the recording in any player (VLC etc.) WITHOUT the FLIR GUI.
%
% FIXED linear temperature window across the whole clip, NO gamma, NO per-frame
% auto-scaling -- so the brightness you see IS true relative temperature and the
% nostril airflow oscillation is visible as-is. (Same rule as the 2P videos.)
%
%   thermal_make_video()                      % glob D:\260611_thermal_breathing
%   thermal_make_video('C:\...\xxx_thermal.mat')
%   thermal_make_video(matPath, 50)           % set playback fps (default = real time)
%
% Output: <stem>_view.mp4 next to the .mat.

if nargin < 1 || isempty(matPath)
    d = dir(fullfile('D:\260611_thermal_breathing', '**', '*_thermal.mat'));
    assert(~isempty(d), 'No *_thermal.mat found.');
    matPath = fullfile(d(1).folder, d(1).name);
end
if nargin < 2 || isempty(playFps), playFps = []; end   % [] = use real fps

UPSCALE = 6;          % nearest-neighbour zoom (192x96 is tiny); honest, no blur
CLIP_PCT = [0.5 99.7];% robust global window so a few hot/cold px don't wash it out
CMAP = hot(256);      % thermal LUT (fixed); set to gray(256) for grayscale

S = load(matPath);
stack = single(S.stack);            % [T x H x W] deg C
fps = double(S.fps);
[T, H, W] = size(stack);
if isempty(playFps), playFps = fps; end

% --- fixed global window (no per-frame scaling) ---
lo = prctile(stack(:), CLIP_PCT(1));
hi = prctile(stack(:), CLIP_PCT(2));
fprintf('thermal_make_video: %d frames %dx%d, fps=%.2f, window [%.2f %.2f] C\n', ...
    T, H, W, fps, lo, hi);

[outDir, stem] = fileparts(matPath);
outMp4 = fullfile(outDir, [stem '_view.mp4']);
vw = VideoWriter(outMp4, 'MPEG-4');
vw.FrameRate = playFps;
vw.Quality = 95;
open(vw);

lut = uint8(round(CMAP * 255));     % 256x3
for k = 1:T
    fr = squeeze(stack(k, :, :));               % H x W deg C
    g = (fr - lo) / (hi - lo);                  % LINEAR map, fixed window
    g = min(max(g, 0), 1);
    idx = uint16(round(g * 255)) + 1;           % 1..256
    rgb = reshape(lut(idx, :), H, W, 3);        % H x W x 3 uint8
    rgb = imresize(rgb, UPSCALE, 'nearest');
    writeVideo(vw, rgb);
end
close(vw);
fprintf('  saved %s  (%.1f s at %.0f fps)\n', outMp4, T/playFps, playFps);
end
