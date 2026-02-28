function stackDFF(dFF, dFF_denoise, fps, colors, dFF_scale, spacingFrac, S, a)
% stackDFF  Quick stacked ΔF/F plot with optional denoised overlay and events,
% plus optional random subsampling of ROIs.
%
% Inputs:
%   ...
%   S   : [T x N] binary event matrix
%   a   : fraction of ROIs to plot (0–1). Default = 1.
%         Example: a=0.5 → randomly pick 50% of ROIs.

% -------------------- defaults --------------------
if nargin < 2 || isempty(dFF_denoise)
    dFF_denoise = [];
end
if nargin < 3 || isempty(fps)
    fps = 30;
end
if nargin < 4 || isempty(colors)
    colors = {'k','r'};  
end
if ischar(colors)
    ccell = cell(1, numel(colors));
    for i = 1:numel(colors)
        ccell{i} = colors(i);
    end
    colors = ccell;
end
if nargin < 5 || isempty(dFF_scale)
    dFF_scale = 1;
end
if nargin < 6 || isempty(spacingFrac)
    spacingFrac = 0.5;
end
spacingFrac = max(0, min(1, spacingFrac));

if nargin < 7
    S = [];
end
hasEvents = ~isempty(S);

if nargin < 8 || isempty(a)
    a = 1;            % default: plot all ROIs
end
a = max(0, min(1, a));   % clamp to [0,1]

[T, N] = size(dFF);

% -------------------- subsample ROIs if a < 1 --------------------
if a < 1
    N_keep = max(1, round(a * N));    % number of ROIs to keep
    idx = randperm(N, N_keep);        % random subset

    % apply to all matrices
    dFF = dFF(:, idx);

    if ~isempty(dFF_denoise)
        dFF_denoise = dFF_denoise(:, idx);
    end

    if hasEvents
        S = S(:, idx);
    end

    fprintf('[stackDFF] Randomly selected %d of %d ROIs (a = %.2f)\n', ...
            N_keep, N, a);

    N = N_keep;   % update ROI count
end

% -------------------- sanity check for S --------------------
if hasEvents && ~isequal(size(S), [T, N])
    error('stackDFF: S must be same size as dFF after ROI subsampling.');
end

% -------------------- time vector --------------------
t = (0:T-1)' / fps;

% -------------------- vertical spacing --------------------
amp = max(dFF(:)) - min(dFF(:));
if amp == 0, amp = 1; end
gap = spacingFrac * amp;
offsets = gap * (0:N-1);

% stacked traces
Y_main = dFF + offsets;
if ~isempty(dFF_denoise)
    Y_deno = dFF_denoise + offsets;
end

% -------------------- plotting --------------------
figure('Name','stackDFF','Color','w');
hold on;

plot(t, Y_main, 'Color', colors{1}, 'LineWidth', 0.5);

if ~isempty(dFF_denoise)
    c2 = colors{min(2, numel(colors))};
    plot(t, Y_deno, 'Color', c2, 'LineWidth', 0.5);
end

% -------------------- event markers --------------------
if hasEvents
    if numel(colors) >= 3
        evColor = colors{3};
    else
        evColor = 'g';
    end

    for i = 1:N
        idx_evt = find(S(:, i) ~= 0);
        if isempty(idx_evt), continue; end

        plot(t(idx_evt), dFF(idx_evt, i) + offsets(i), ...
             'o', 'MarkerFaceColor', evColor, ...
             'MarkerEdgeColor', evColor, 'MarkerSize', 3);
    end
end

% -------------------- scale bar --------------------
xRange = t(end) - t(1);
x0 = xRange;

yTop = offsets(end) + 0.5*amp;
y0   = yTop - dFF_scale/2;

plot([x0 x0], [y0 y0 + dFF_scale], 'b', 'LineWidth', 2);
text(x0, y0 + dFF_scale + 0.05*dFF_scale, ...
     sprintf('%.2f \\DeltaF/F', dFF_scale), ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 9);

ylim([-gap, offsets(end) + amp + gap]);
xlim([t(1) t(end)]);

yt  = offsets;
ytk = 1:N;
keep = (mod(ytk,5)==1) | (ytk==N);
set(gca, 'YTick', yt(keep), 'YTickLabel', string(ytk(keep)));

xlabel('Time (s)');
ylabel('ROI#');
title('Stacked \DeltaF/F traces');
box on;
hold off;

end