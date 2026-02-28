function darkOut = estimate_dark_from_movie(inputArg, varargin)
% estimate_dark_from_movie
%   Accepts either:
%      (1) tiffPath   : string / char
%      (2) V          : preloaded 3D array (HxWxT)
%
% Usage:
%   out = estimate_dark_from_movie("movie.tif");
%   out = estimate_dark_from_movie(V);   % preloaded array

    % -------------------- Parse options --------------------
    p = inputParser;
    addParameter(p, 'MaxFrames', 200);
    addParameter(p, 'Subsample', 4);
    parse(p, varargin{:});
    P = p.Results;

    % ------------------ File or array? ---------------------
    tiffPathLocal = '';
    if ischar(inputArg) || isstring(inputArg)
        % Case 1: filename
        tiffPathLocal = char(inputArg);
        assert(isfile(tiffPathLocal), 'File not found: %s', tiffPathLocal);

        info = imfinfo(tiffPathLocal);
        fprintf('[dark] Loading TIFF from file: %s (%d frames)\n', ...
                tiffPathLocal, numel(info));

        V = single(tiffreadVolume(tiffPathLocal));   % load whole stack

    elseif isnumeric(inputArg)
        % Case 2: already loaded array
        V = single(inputArg);
        fprintf('[dark] Using pre-loaded movie array: [%d x %d x %d]\n', ...
                size(V,1), size(V,2), size(V,3));

    else
        error('Input must be a file path or a numeric movie array.');
    end

    % ----------------- Frame skip + sampling -----------------
    skipN = 300;                          % skip initial frames
    if size(V,3) > skipN
        V = V(:,:,skipN+1:end);
    end

    Ttot = size(V,3);
    nFrames    = min(Ttot, P.MaxFrames);
    idx_frames = round(linspace(1, Ttot, nFrames));

    ss = max(1, P.Subsample);             % subsample >= 1

    fprintf('[dark] Sampling %d frames (skip first %d), subsample=%d\n', ...
            nFrames, skipN, ss);

    Vss  = V(1:ss:end, 1:ss:end, idx_frames);
    vals = Vss(:);
    vals = vals(isfinite(vals));

    fprintf('[dark] Collected %d pixel samples\n', numel(vals));

    % ------------------- Histogram -------------------
    lo = max(0, prctile(vals, 0));
    hi = prctile(vals, 100);
    if hi <= lo
        hi = lo + 1;
    end

    nBins = 512;
    edges = linspace(lo, hi, nBins+1);
    [counts, ~] = histcounts(vals, edges);
    bin_centers = (edges(1:end-1) + edges(2:end))/2;

    counts_smooth = smoothdata(counts, 'gaussian', 9);

    % ----------------- First peak (dark) -----------------
    [pks, locs] = findpeaks(counts_smooth, bin_centers, ...
                            'MinPeakProminence', max(counts_smooth)*0.01);

    if isempty(pks)
        warning('[dark] No peaks found; falling back to low intensity edge.');
        dark_peak_loc = lo;
    else
        % left-most peak
        [~, idx_min]  = min(locs);
        dark_peak_loc = locs(idx_min);
    end

    % Use first peak location as mean
    dark_mean = dark_peak_loc;

    % ---------------- One-sided std from left of mean ----------------
    mask_left      = vals <= dark_mean;
    dark_vals_left = vals(mask_left);

    if isempty(dark_vals_left)
        warning('[dark] No values left of peak; using vals <= peak as dark.');
        mask_left      = vals <= dark_peak_loc;
        dark_vals_left = vals(mask_left);
    end

    if isempty(dark_vals_left)
        warning('[dark] Still no dark values; falling back to overall std.');
        dark_std_trunc = std(vals);
    else
        dev_left       = dark_vals_left - dark_mean;  % mostly negative
        dark_std_trunc = std(dev_left);
    end

    % Half-normal correction: std_trunc = sigma * sqrt(1 - 2/pi)
    scale_half = sqrt(1 - 2/pi);          % ~0.6028
    dark_std   = dark_std_trunc / scale_half;

    fprintf('[dark] dark_mean = %.2f, dark_std (extrapolated) = %.2f (ADU)\n', ...
            dark_mean, dark_std);

    % ---------------- Gaussian curve for visualization ----------------
    gauss_pdf = (1/(dark_std * sqrt(2*pi))) * ...
                exp( -0.5 * ((bin_centers - dark_mean)./dark_std).^2 );

    % Scale Gaussian so peak matches smoothed histogram near dark peak
    mask_peak = (bin_centers >= dark_mean - 3*dark_std) & ...
                (bin_centers <= dark_mean + 3*dark_std);

    if any(mask_peak)
        peak_hist  = max(counts(mask_peak));
        peak_gauss = max(gauss_pdf(mask_peak));
    else
        peak_hist  = max(counts_smooth);
        peak_gauss = max(gauss_pdf);
    end

    if peak_gauss > 0
        scale = peak_hist / peak_gauss;
    else
        scale = 1;
    end

    gauss_curve = gauss_pdf * scale;

    % ----------------- Pack outputs -----------------
    darkOut = struct();
    darkOut.dark_mean        = dark_mean;
    darkOut.dark_std         = dark_std;
    darkOut.dark_std_trunc   = dark_std_trunc;      % one-sided / truncated std
    darkOut.bin_centers      = bin_centers;
    darkOut.counts           = counts;
    darkOut.counts_smooth    = counts_smooth;
    darkOut.dark_peak_loc    = dark_peak_loc;
    darkOut.gauss_curve      = gauss_curve;
    darkOut.n_samples        = numel(vals);
    darkOut.n_dark_left      = numel(dark_vals_left);
    if ~isempty(tiffPathLocal)
        darkOut.tiffPath     = tiffPathLocal;
    else
        darkOut.tiffPath     = '(preloaded array)';
    end

    % ----------------- QC plot -----------------
    figure('Color','w');
    plot(bin_centers, counts, 'Color', [0 0 0], 'LineWidth', 1.5); hold on;
    %plot(bin_centers, counts_smooth, 'k', 'LineWidth', 1.5);
    plot(bin_centers, gauss_curve, 'b-', 'LineWidth', 1);
    yl = ylim;
    plot([dark_mean dark_mean], yl, 'r--', 'LineWidth', 1.5);

    xlabel('Intensity (ADU)');
    ylabel('Count');
    title(sprintf('Dark histogram: mean = %.1f, std = %.1f (ADU)', ...
                  dark_mean, dark_std));
    legend({'Raw hist','Gaussian dark noise','Dark current'});
    grid on;
end

