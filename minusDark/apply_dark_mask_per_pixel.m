function out = apply_dark_mask_per_pixel(inputArg, darkOut, varargin)
% apply_dark_mask_per_pixel
%   Use darkOut (from estimate_dark_from_movie) to:
%     1) subtract dark current from a movie
%     2) compute per-pixel mean (after dark subtraction)
%     3) mask pixels whose mean is NOT greater than Factor * dark_std
%        (i.e., mean <= Factor * dark_std -> set to 0 for all frames)
%
% Inputs:
%   inputArg : either
%              (1) tiffPath : string / char
%              (2) V        : preloaded 3D array [H x W x T]
%   darkOut  : struct from estimate_dark_from_movie (must have
%              .dark_mean and .dark_std)
%
% Optional name/value:
%   'Factor'        : threshold multiplier on dark_std (default: 3)
%   'ClipNegative'  : true -> clip after dark subtraction to >= 0 (default: true)
%   'WriteTiff'     : true -> if input is a file, write masked TIFF (default: true)
%   'OutSuffix'     : suffix for output TIFF name (default: '_darkMasked')
%
% Output struct out:
%   .V_masked   : masked movie (same class as input if array given)
%   .mask       : [H x W] logical; true = kept (mean > Factor * dark_std)
%   .meanImg    : [H x W] per-pixel mean after dark subtraction
%   .thr        : scalar threshold (= Factor * dark_std)
%   .dark_mean  : copied from darkOut
%   .dark_std   : copied from darkOut
%   .tiffPathIn : input path (if any)
%   .tiffPathOut: output TIFF path (if written)

    % -------------------- Parse options --------------------
    p = inputParser;
    addParameter(p, 'Factor',       3,    @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ClipNegative', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'WriteTiff',    true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'OutSuffix',    '_darkMasked', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    P = p.Results;

    if ~isfield(darkOut, 'dark_mean') || ~isfield(darkOut, 'dark_std')
        error('darkOut must contain fields .dark_mean and .dark_std');
    end

    dark_mean = darkOut.dark_mean;
    dark_std  = darkOut.dark_std;
    thr       = P.Factor * dark_std;

    fprintf('[darkMask] dark_mean = %.2f, dark_std = %.2f, Factor = %.1f, thr = %.2f\n', ...
        dark_mean, dark_std, P.Factor, thr);

    % ------------------ File or array? ---------------------
    tiffPathLocal = '';
    if ischar(inputArg) || isstring(inputArg)
        % Case 1: filename
        tiffPathLocal = char(inputArg);
        assert(isfile(tiffPathLocal), 'File not found: %s', tiffPathLocal);

        info = imfinfo(tiffPathLocal);
        fprintf('[darkMask] Loading TIFF from file: %s (%d frames)\n', ...
                tiffPathLocal, numel(info));

        V = tiffreadVolume(tiffPathLocal);  % preserve original class
        inClass = class(V);

        % cast to single for math
        V = single(V);

    elseif isnumeric(inputArg)
        % Case 2: already loaded array
        V = inputArg;
        fprintf('[darkMask] Using pre-loaded movie array: [%d x %d x %d]\n', ...
                size(V,1), size(V,2), size(V,3));

        inClass = class(V);
        % do math in single
        if ~isa(V, 'single')
            V = single(V);
        end

    else
        error('Input must be a file path or a numeric movie array.');
    end

    [H, W, T] = size(V);

    % ------------------ Dark subtraction -------------------
    fprintf('[darkMask] Subtracting dark_mean = %.2f\n', dark_mean);
    V_sub = V - dark_mean;     % single

    if P.ClipNegative
        V_sub(V_sub < 0) = 0;
    end

    % ------------------ Per-pixel mean ---------------------
    fprintf('[darkMask] Computing per-pixel mean over T = %d frames\n', T);
    meanImg = mean(V_sub, 3);   % [H x W] single

    % ------------------ Build mask -------------------------
    mask = meanImg > thr;       % logical
    fracKeep = nnz(mask) / numel(mask) * 100;
    fprintf('[darkMask] Keeping %.2f %% of pixels (mean > %.2f)\n', fracKeep, thr);

    % ------------------ Apply mask -------------------------
    V_masked = V_sub .* single(mask);   % broadcast mask (single)

    % ------------------ Cast back to original class --------
    switch inClass
        case {'uint8','uint16','uint32'}
            V_masked(V_masked < 0) = 0;
            maxType = double(intmax(inClass));
            V_masked(V_masked > maxType) = maxType;
            V_masked = cast(V_masked, inClass);

        case {'int8','int16','int32'}
            minType = double(intmin(inClass));
            maxType = double(intmax(inClass));
            V_masked(V_masked < minType) = minType;
            V_masked(V_masked > maxType) = maxType;
            V_masked = cast(V_masked, inClass);

        otherwise
            % leave as single (or whatever) if it was float
            % If you want strictly 'single', uncomment:
            % V_masked = single(V_masked);
    end

    % ------------------ Optional TIFF write ----------------
    tiffPathOut = '';
    if ~isempty(tiffPathLocal) && P.WriteTiff
        [folder, name, ext] = fileparts(tiffPathLocal);
        outName    = sprintf('%s%s%s', name, P.OutSuffix, ext);
        tiffPathOut = fullfile(folder, outName);

        fprintf('[darkMask] Writing masked TIFF: %s\n', tiffPathOut);

        % write frame by frame to avoid huge memory duplication
        for k = 1:T
            frOut = V_masked(:,:,k);

            if k == 1
                imwrite(frOut, tiffPathOut, 'Compression', 'none');
            else
                imwrite(frOut, tiffPathOut, 'WriteMode', 'append', 'Compression', 'none');
            end
        end
    end

    % ------------------ Pack outputs -----------------------
    out = struct();
    out.V_masked    = V_masked;
    out.mask        = mask;
    out.meanImg     = meanImg;
    out.thr         = thr;
    out.dark_mean   = dark_mean;
    out.dark_std    = dark_std;
    out.tiffPathIn  = tiffPathLocal;
    out.tiffPathOut = tiffPathOut;
end
