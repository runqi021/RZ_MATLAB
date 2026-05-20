function hmapDFF(dFF, fps, prcLim, cmap)
% hmapDFF
%   Heatmap of ΔF/F with correct y-axis orientation.
%
% Inputs:
%   dFF    : [T x N] ΔF/F matrix (ROIs as columns)
%   fps    : sampling rate (Hz). Default = 30
%   prcLim : [low high] percentiles for color scaling (e.g. [0 99])
%   cmap   : colormap (e.g. parula, turbo, gray). Default = parula

    if nargin < 2 || isempty(fps)
        fps = 30;
    end
    if nargin < 3 || isempty(prcLim)
        prcLim = [0 100];
    end
    if nargin < 4 || isempty(cmap)
        cmap = bluewhitered;
    end

    [T, ~] = size(dFF);
    t = (0:T-1)' / fps;

    % percentile-based color limits
    cLow  = prctile(dFF(:), prcLim(1));
    cHigh = prctile(dFF(:), prcLim(2));
    cAbs = max(abs(cLow), abs(cHigh));
    
    figure('Name','\DeltaF/F heatmap','Color','w');
    imagesc(t, 1:size(dFF,2), dFF');   % rows = ROIs
    set(gca, 'YDir', 'normal');

    xlabel('Time (s)');
    ylabel('ROI');
    title('\DeltaF/F heatmap');
    colormap(cmap);
    caxis([-cAbs cAbs]);

    cb = colorbar;
    ylabel(cb, '\DeltaF/F');
end
