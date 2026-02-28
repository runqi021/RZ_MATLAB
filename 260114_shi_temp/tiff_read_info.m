tiffPath = "E:\live_251104_definitely_not_sst_tdtomato_fitc\raw\live_not_sst_fitc_1.2x_930nm_7-26lp_z220_-10_100f_col00_row00_x-1500_y-499_00001.tif";
%%
info   = imfinfo(tiffPath);

%% ---------------- READ ScanImage META (POWER ONLY) ----------------
% Robust metadata extraction (Software -> ImageDescription -> low-level tag)
meta = "";
if isfield(info(1),'Software') && ~isempty(info(1).Software)
    meta = string(info(1).Software);
elseif isfield(info(1),'ImageDescription') && ~isempty(info(1).ImageDescription)
    meta = string(info(1).ImageDescription);
else
    t0 = Tiff(tiffPath,'r');
    c0 = onCleanup(@() t0.close());
    try
        meta = string(t0.getTag('ImageDescription'));
    catch
        try
            meta = string(t0.getTag('Software'));
        catch
            error("No ScanImage metadata found in Software/ImageDescription tags.");
        end
    end
end
meta = replace(meta, "↵", newline);

getNum = @(key) local_get_num(meta, key);
getVec = @(key) local_get_vec(meta, key);

LC_SI = getNum("SI.hBeams.lengthConstants");
P0_SI = getNum("SI.hBeams.stackStartPower");   % percent
P1_SI = getNum("SI.hBeams.stackEndPower");     % percent

% Power schedule from ScanImage (aligned by SLICE INDEX; independent of your z_um labeling)
% Model: P(z_motor) = P0 * exp((z_motor - z0)/LC)
p_pct  = P0_SI .* exp((z_motor - z_motor(1)) ./ LC_SI);

% Optional: exact endpoint match (if you want strict P0->P1 regardless of LC rounding)
% p_pct = exp(interp1([z_motor(1), z_motor(end)], log([P0_SI,P1_SI]), z_motor, 'linear','extrap'));

p_frac = p_pct / 100;

% sanity print
depth_meta = z_motor(1) - z_motor;  % depth-like axis from SI (not used for your plot)
if any(diff(depth_meta) <= 0), depth_meta = -depth_meta; end

fprintf('[meta] SI power : %.2f%% -> %.2f%%, LC=%.4g um\n', ...
    p_pct(1), p_pct(end), LC_SI);

%%
function v = local_get_num(meta, key)
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), v = NaN; return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
tok = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = local_get_vec(meta, key)
vec = [];
pat = key + " = ";
i = strfind(meta, pat);
if isempty(i), return; end
s = extractAfter(meta, i(1) + strlength(pat) - 1);
line = extractBefore(s, newline);
m = regexp(line, '\[(.*)\]', 'tokens', 'once');
if isempty(m), return; end
nums = regexp(m{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
if isempty(nums), return; end
vec = str2double(nums(:));
end
