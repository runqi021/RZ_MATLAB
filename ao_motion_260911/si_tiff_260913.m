function varargout = si_tiff_260913(mode, varargin)
%SI_TIFF_260913  Reading raw ScanImage-logged TIFFs. One implementation.
%
%   [nCh, chIdx, fps] = si_tiff_260913('header', path, channel)
%   si_tiff_260913('stream', path, channel, chunkFrames, feedFcn)
%   F = si_tiff_260913('readall', path, channel)
%
% Shared by dm_axial_analyze_260913 and dm_axial_loop_eval_260913 so the
% channel-interleave arithmetic exists in exactly one place. That arithmetic is
% the dangerous part: pages interleave channels PER FRAME, so
%       page = (frame-1)*nCh + chIdx
% and a stride that gets chIdx wrong still returns correctly-sized images and
% never errors -- it silently hands back the OTHER label.
%
% Runqi Zhang / 2026-09-13.

switch lower(mode)
    case 'header'
        [varargout{1}, varargout{2}, varargout{3}] = si_header(varargin{:});
    case 'stream'
        si_stream(varargin{:});
    case 'readall'
        varargout{1} = si_read_all(varargin{:});
    otherwise
        error('si_tiff:mode','unknown mode "%s"', mode);
end
end

function [nCh, chIdx, fps] = si_header(path, channel)
assert(exist(path,'file') == 2, 'missing file: %s', path);
tf = Tiff(path,'r');
c = onCleanup(@() close(tf)); %#ok<NASGU>
sw = '';
try, sw = tf.getTag('Software'); catch, end
if isempty(sw)
    try, sw = tf.getTag('ImageDescription'); catch, end
end
chSave = parse_vec(sw, 'SI\.hChannels\.channelSave');
if isempty(chSave), chSave = channel; end
nCh   = numel(chSave);
chIdx = find(chSave == channel, 1);
assert(~isempty(chIdx), ...
    ['channel %d was not saved in %s (saved: %s).\n' ...
     'Pages interleave channels per frame, so the wrong index would silently ' ...
     'return the other label instead of erroring.'], ...
    channel, path, mat2str(chSave));
fps = parse_num(sw, 'SI\.hRoiManager\.scanFrameRate');
if isempty(fps) || ~isfinite(fps), fps = NaN; end
end

function si_stream(path, channel, chunkFrames, feed)
% Walk the file ONCE with nextDirectory. Never setDirectory(t,k) in a loop --
% that is quadratic in page count and these files run to thousands of pages.
[nCh, chIdx] = si_header(path, channel);
tf = Tiff(path,'r');
c  = onCleanup(@() close(tf)); %#ok<NASGU>

buf = [];  tbuf = [];  nIn = 0;  page = 0;
while true
    page = page + 1;
    if (mod(page-1, nCh) + 1) == chIdx
        im = tf.read();
        ts = page_time(tf);
        if isempty(buf)
            buf  = zeros([size(im) chunkFrames], 'single');
            tbuf = nan(1, chunkFrames);
        end
        nIn = nIn + 1;
        buf(:,:,nIn) = single(im); %#ok<AGROW>
        tbuf(nIn) = ts;
        if nIn == chunkFrames
            feed(buf, tbuf);
            nIn = 0;
        end
    end
    if tf.lastDirectory(), break; end
    tf.nextDirectory();
end
if nIn > 0, feed(buf(:,:,1:nIn), tbuf(1:nIn)); end
end

function F = si_read_all(path, channel)
out = {};
    function grab(B, ~)
        out{end+1} = B; %#ok<AGROW>
    end
si_stream(path, channel, 64, @grab);
assert(~isempty(out), 'no frames read from %s', path);
F = cat(3, out{:});
end

function ts = page_time(tf)
% frameTimestamps_sec is in the PER-PAGE ImageDescription, not the header.
ts = NaN;
try
    v = parse_num(tf.getTag('ImageDescription'), 'frameTimestamps_sec');
    if ~isempty(v), ts = v; end
catch
end
end

function v = parse_num(s, key)
v = [];
if isempty(s), return; end
m = regexp(s, [key '\s*=\s*([-\d\.eE+]+)'], 'tokens', 'once');
if ~isempty(m), v = str2double(m{1}); end
end

function v = parse_vec(s, key)
v = [];
if isempty(s), return; end
m = regexp(s, [key '\s*=\s*\[([^\]]*)\]'], 'tokens', 'once');
if ~isempty(m)
    v = str2double(regexp(m{1}, '[-\d\.eE+]+', 'match'));
    v = v(isfinite(v));
    return
end
m = regexp(s, [key '\s*=\s*([-\d\.eE+]+)'], 'tokens', 'once');
if ~isempty(m), v = str2double(m{1}); end
end
