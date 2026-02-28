svd_reconstruct_firstK_perZ()

%%
function svd_reconstruct_firstK_perZ()
% svd_reconstruct_firstK_perZ
% ------------------------------------------------------------
% For Z-stack TIFF with 50 frames per Z:
%   For each Z group:
%     - load Y [pixels x T]
%     - (optional) demean per pixel
%     - SVD via econ on time-cov (fast): eig(Y'Y) -> V, S; then U = Y*V/S
%     - reconstruct with first K modes
%     - add mean back (if demeaned)
%     - write reconstructed frames to output TIFF (same class/tags as input)
%
% NOTE: This is per-Z reconstruction; output has same #pages as input.
% ------------------------------------------------------------

%% ---------------- USER PARAMS ----------------
tiffPath   = "C:\Users\Admin\Downloads\shi_250115\wt_fitc\live_not_sst_fitc_1.2x_930nm_7-26lp_z220_-10_100f_col01_row01_x-900_y101_00001_minusBgPerZ_drop1_cropLR20.tif";
framesPerZ = 99;
dropFirstN = 0;
demeanPerPixel = true;

Krecon = 8;                   % <-- reconstruct using first 10 modes

outSuffix = sprintf("_svdK%02d", Krecon);

%% ---------------- METADATA / TYPE ----------------
assert(isfile(tiffPath), "File not found: %s", tiffPath);
info   = imfinfo(tiffPath);
nPages = numel(info);
H = info(1).Height;
W = info(1).Width;

if mod(nPages, framesPerZ) ~= 0
    error("nPages (%d) not divisible by framesPerZ (%d).", nPages, framesPerZ);
end
nZ = nPages / framesPerZ;

T = framesPerZ - dropFirstN;
if T <= 2, error("framesPerZ-dropFirstN too small: %d", T); end
Krecon = min(Krecon, T);

[inClass, bits, sf] = detect_tiff_class_from_tags(tiffPath);
fprintf("[perZ-recon] class=%s Bits=%d SampleFormat=%s\n", inClass, bits, sf);

% output path
[folder, base, ext] = fileparts(char(tiffPath));
outPath = fullfile(folder, base + outSuffix + ext);
if isfile(outPath), delete(outPath); end

%% ---------------- IO OPEN ----------------
tIn  = Tiff(tiffPath, 'r');
cIn  = onCleanup(@() tIn.close());

tOut = Tiff(outPath, 'w');
cOut = onCleanup(@() tOut.close());

tag = make_tiff_tag_from_info(info(1), bits, sf);

%% ---------------- BUFFERS ----------------
Y  = zeros(H*W, T, 'single');

%% ---------------- LOOP Z: SVD + RECON + WRITE ----------------
pageOut = 0;
for iz = 1:nZ
    basePage = (iz-1)*framesPerZ;

    % load frames into Y
    kk = 0;
    for f = 1:framesPerZ
        pageIdx = basePage + f;
        setDirectory(tIn, pageIdx);
        fr = single(tIn.read());
        if f <= dropFirstN
            continue;
        end
        kk = kk + 1;
        Y(:,kk) = fr(:);
    end

    % mean + demean
    if demeanPerPixel
        mu = mean(Y, 2, 'omitnan');          % [P x 1]
        Yd = Y - mu;
    else
        mu = zeros(H*W,1,'single');
        Yd = Y;
    end

    % ---- FAST SVD via time-covariance (T x T) ----
    % C = Yd'Yd = V * diag(s^2) * V'
    C = double(Yd' * Yd);
    C = (C + C')/2;                          % symmetrize numeric noise
    [V, D] = eig(C);
    lam = real(diag(D));
    lam(lam < 0) = 0;

    % sort desc
    [lam, ord] = sort(lam, 'descend');
    V = V(:, ord);

    s = sqrt(lam);                           % singular values
    % Avoid divide-by-zero
    s_eps = s;
    s_eps(s_eps < 1e-12) = 1e-12;

    % U = Yd * V * diag(1/s)
    % We'll only need first Krecon
    Vk = V(:, 1:Krecon);
    sk = s_eps(1:Krecon);
    Uk = Yd * single(Vk) ./ single(sk');     % [P x K] (broadcast over columns)

    % Reconstruction: Yhat = Uk * diag(sk) * Vk'
    % Uk .* sk'  gives [P x K], then * Vk' -> [P x T]
    Yhat = (Uk .* single(sk')) * single(Vk'); % [P x T]
    Yhat = Yhat + mu;                         % add mean back

    % ---- WRITE reconstructed frames for this Z group ----
    for f = 1:framesPerZ
        if f <= dropFirstN
            % if you dropped frames, pass-through originals (optional)
            % Here: pass-through originals to keep page count identical.
            setDirectory(tIn, basePage + f);
            fr0 = tIn.read();
            frOut = fr0;
        else
            tt = f - dropFirstN;
            frOut = reshape(Yhat(:,tt), [H W]);
            frOut = cast_like(frOut, inClass);
        end

        pageOut = pageOut + 1;
        tOut.setTag(tag);
        tOut.write(frOut);

        if pageOut < nPages
            tOut.writeDirectory();
        end
    end

    if mod(iz, max(1,floor(nZ/10)))==0 || iz==nZ
        fprintf("[perZ-recon] z %d/%d done (K=%d)\n", iz, nZ, Krecon);
    end
end

fprintf("Wrote: %s\n", outPath);

end

%% ================= helpers =================

function [cls,bits,fmtStr] = detect_tiff_class_from_tags(fp)
t = Tiff(fp,'r');
bits = t.getTag('BitsPerSample');
sf   = t.getTag('SampleFormat');
t.close();
switch sf
    case Tiff.SampleFormat.UInt
        fmtStr='UInt';
        if bits==8, cls='uint8';
        elseif bits==16, cls='uint16';
        elseif bits==32, cls='uint32';
        else, cls='uint16';
        end
    case Tiff.SampleFormat.Int
        fmtStr='Int';
        if bits==16, cls='int16';
        elseif bits==32, cls='int32';
        else, cls='int16';
        end
    case Tiff.SampleFormat.IEEEFP
        fmtStr='IEEEFP';
        cls='single';
    otherwise
        fmtStr='Unknown';
        cls='int16';
end
end

function tag = make_tiff_tag_from_info(info1, bits, fmtStr)
tag.ImageLength = info1.Height;
tag.ImageWidth  = info1.Width;
tag.Photometric = Tiff.Photometric.MinIsBlack;
tag.SamplesPerPixel = 1;
tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tag.Compression = Tiff.Compression.None;
tag.RowsPerStrip = min(info1.Height, 256);
tag.BitsPerSample = bits;
switch fmtStr
    case 'UInt',   tag.SampleFormat = Tiff.SampleFormat.UInt;
    case 'Int',    tag.SampleFormat = Tiff.SampleFormat.Int;
    case 'IEEEFP', tag.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise,     tag.SampleFormat = Tiff.SampleFormat.Int;
end
end

function out = cast_like(x, cls)
x = double(x);
switch cls
    case 'int16'
        x = min(max(x, -32768), 32767);
        out = int16(round(x));
    case 'int32'
        x = min(max(x, double(intmin('int32'))), double(intmax('int32')));
        out = int32(round(x));
    case 'uint16'
        x = min(max(x, 0), 65535);
        out = uint16(round(x));
    case 'uint8'
        x = min(max(x, 0), 255);
        out = uint8(round(x));
    case 'single'
        out = single(x);
    otherwise
        x = min(max(x, -32768), 32767);
        out = int16(round(x));
end
end
