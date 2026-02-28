loadtiff("E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\raw\fov2_-200_L_30lp_1100nm_5x_00003.tif");

%%
fn = "C:\Users\Admin\Desktop\260210_chat_soma_g8s\roi1_-1050_400_z-25-115_IO_00001.tif";
info = imfinfo(fn);

info(2)   % look at first frame


%%
parse_scanimage_power_from_info(info)

%%
function out = parse_scanimage_power_from_info(info)
% out = parse_scanimage_power_from_info(info)
% Parse ScanImage metadata from imfinfo() output and compute power-vs-Z schedule.
%
% Usage:
%   info = imfinfo(fn);
%   out = parse_scanimage_power_from_info(info);
%   out

% --- get the metadata text ---
% Depending on ScanImage version, this may live in info(1).Software or ImageDescription.
meta = "";
if isfield(info(1), 'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1), 'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
else
    error('No Software/ImageDescription metadata found.');
end

% Your paste uses ↵ as line breaks. Normalize to \n:
meta = replace(meta, char(8593), newline);  % just in case (unlikely)
meta = replace(meta, "↵", newline);

% helper to grab numeric or string values from "KEY = VALUE" lines
getNum = @(key) local_get_num(meta, key);
getStr = @(key) local_get_str(meta, key);

% --- extract what you asked for ---
out.lengthConstant = getNum("SI.hBeams.lengthConstants");      % "length constant"
out.startPower     = getNum("SI.hBeams.stackStartPower");      % stackStartPower
out.endPower       = getNum("SI.hBeams.stackEndPower");        % stackEndPower
out.zStart         = getNum("SI.hStackManager.stackZStartPos");% first Z
out.zEnd           = getNum("SI.hStackManager.stackZEndPos");  % last Z
out.zStep          = getNum("SI.hStackManager.stackZStepSize");% step (can be negative)

% Optional: framesPerSlice, numSlices (nice sanity checks)
out.framesPerSlice = getNum("SI.hStackManager.framesPerSlice");
out.numSlices      = getNum("SI.hStackManager.numSlices");

% Try to parse explicit zs list if present (often best truth)
out.zs = local_get_vec(meta, "SI.hStackManager.zs");

% --- build Z vector ---
if ~isempty(out.zs)
    z = out.zs(:);
else
    % include endpoints robustly even with negative step
    z = out.zStart : out.zStep : out.zEnd;
    z = z(:);
end
out.z = z;

% --- compute exponential power schedule ---
% ScanImage stores "lengthConstants" (LC) for power vs Z. A common model is:
%   P(z) = P0 * exp((z - z0)/LC)
% But ScanImage also allows specifying start/end power across a Z range.
% We will compute a schedule that:
%   (1) matches startPower at zStart
%   (2) matches endPower at zEnd
% If LC is provided, we report the implied LC too, and also provide an LC-based curve.
%
% First: simple log-linear interpolation from (zStart,Pstart) to (zEnd,Pend):
Pstart = out.startPower;
Pend   = out.endPower;

% Guard against zeros/NaNs
if any(~isfinite([Pstart,Pend])) || Pstart<=0 || Pend<=0
    out.power_loginterp = nan(size(z));
else
    % log(P) linear in z
    out.power_loginterp = exp( interp1([out.zStart,out.zEnd], log([Pstart,Pend]), z, 'linear','extrap') );
end

% Second: LC-based curve anchored at startPower:
LC = out.lengthConstant;
if isfinite(LC) && LC ~= 0 && isfinite(Pstart) && Pstart>0
    out.power_LC = Pstart * exp((z - out.zStart)/LC);
else
    out.power_LC = nan(size(z));
end

% Implied LC from endpoints (what LC would be required to hit Pend exactly):
if isfinite(Pstart) && isfinite(Pend) && Pstart>0 && Pend>0 && out.zEnd ~= out.zStart
    out.lengthConstant_implied = (out.zEnd - out.zStart) / log(Pend / Pstart);
else
    out.lengthConstant_implied = NaN;
end

% Provide a recommended "schedule" (usually: loginterp that exactly hits endpoints)
out.power_schedule = out.power_loginterp;

end

% ---------------- local helpers ----------------
function v = local_get_num(meta, key)
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), v = NaN; return; end
% take the first occurrence
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
% strip quotes/brackets if any, then parse first numeric token
tok = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function s = local_get_str(meta, key)
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), s = ""; return; end
s2 = extractAfter(meta, i(1) + strlength(pat) - 1);
line = strtrim(extractBefore(s2, newline));
s = line;
end

function vec = local_get_vec(meta, key)
% parse "key = [ ... ]" as numeric vector if present
vec = [];
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
% find bracket content
m = regexp(line, '\[(.*)\]', 'tokens', 'once');
if isempty(m), return; end
nums = regexp(m{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end

