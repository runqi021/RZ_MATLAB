function [] = Fhist(inputArg, varargin)
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
    addParameter(p, 'MaxFrames', 1000);
    addParameter(p, 'Subsample', 4);
    addParameter(p, "Ylim", []);
    addParameter(p, "sat", 65535);
    parse(p, varargin{:});
    P = p.Results;

    % ------------------ File or array? ---------------------
    tiffPathLocal = '';
    if ischar(inputArg) || isstring(inputArg)
        % Case 1: filename
        tiffPathLocal = char(inputArg);
        assert(isfile(tiffPathLocal), 'File not found: %s', tiffPathLocal);

        info = imfinfo(tiffPathLocal);
        fprintf('[Fhist] Loading TIFF from file: %s (%d frames)\n', ...
                tiffPathLocal, numel(info));

        V = single(tiffreadVolume(tiffPathLocal));   % load whole stack

    elseif isnumeric(inputArg)
        % Case 2: already loaded array
        V = single(inputArg);
        fprintf('[Fhist] Using pre-loaded movie array: [%d x %d x %d]\n', ...
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

    fprintf('[Fhist] Sampling %d frames (skip first %d), subsample=%d\n', ...
            nFrames, skipN, ss);

    Vss  = V(1:ss:end, 1:ss:end, idx_frames);
    vals = Vss(:);
    vals = vals(isfinite(vals));

    fprintf('[Fhist] Collected %d pixel samples\n', numel(vals));

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


       % ----------------- QC plot -----------------
    figure('Color','w');
    plot(bin_centers, counts, 'Color', [0 0 0], 'LineWidth', 1.5); hold on;
    %plot(bin_centers, counts_smooth, 'k', 'LineWidth', 1.5);
    %plot(bin_centers, gauss_curve, 'b-', 'LineWidth', 1);
   
    if ~isempty(P.Ylim) ylim(P.Ylim); end
    if ~isempty(P.sat) xline(P.sat); end
    
    xlabel('Intensity (ADU)');
    ylabel('Count');
    
    %legend({'Raw hist','Gaussian dark noise','Dark current'});
    grid on;

    figure('Color','w');
    loglog(bin_centers, counts, 'Color', [0 0 0], 'LineWidth', 1.5); hold on;
    %plot(bin_centers, counts_smooth, 'k', 'LineWidth', 1.5);
    %plot(bin_centers, gauss_curve, 'b-', 'LineWidth', 1);
   
    if ~isempty(P.Ylim) ylim(P.Ylim); end
    if ~isempty(P.sat) xline(P.sat); end

    xlabel('Intensity (ADU)');
    ylabel('Count');
    %legend({'Raw hist','Gaussian dark noise','Dark current'});
    grid on;

end