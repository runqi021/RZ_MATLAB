function plt_vector(X,f,plt,Xerr,c,w)
% plot_vector  Plot spectrum/coherence with optional shaded CI
% Usage:
%   plot_vector(X,f,plt,Xerr,c,w)
%
% Inputs:
%   X    – data (length N)
%   f    – frequency vector (length N)
%   plt  – 'l' to convert to 10·log10 (dB), anything else → linear
%   Xerr – either
%           • scalar → horizontal CI line at that level, or
%           • 2×N matrix → [lower; upper] bounds for fill
%   c    – color spec (e.g. 'k','r','g'), default 'b'
%   w    – line width, default 1

N = numel(X);
if nargin<2 || isempty(f), f = (1:N); end
if numel(f)~=N, error('Length of X and f must match'); end
if nargin<3 || isempty(plt), plt = 'l'; end
if nargin<4, Xerr = []; end
if nargin<5 || isempty(c), c = 'b'; end
if nargin<6 || isempty(w), w = 1; end

% apply dB if requested
if strcmp(plt,'l')
    X = 10*log10(X);
    if ~isempty(Xerr) && ~isscalar(Xerr)
        Xerr = 10*log10(Xerr);
    end
end

hold on

if isempty(Xerr)
    % no CI at all
    plot(f, X, c, 'LineWidth', w);

elseif isscalar(Xerr)
    % single‐level CI (dashed line)
    plot(f, X, c, 'LineWidth', w);
    yl = get(gca,'YLim');
    plot(get(gca,'XLim'), [Xerr Xerr], [c '--'], 'LineWidth', w);

else
    % CI shading with performance tweaks
if size(Xerr,1)==2 && size(Xerr,2)==N
    pos_idx = f > 0;
    if any(pos_idx)
        f_plot = f(pos_idx);
        X_plot = X(pos_idx);
        lower = Xerr(1,pos_idx);
        upper = Xerr(2,pos_idx);

        xpoly = [f_plot(:); flipud(f_plot(:))];
        ypoly = [lower(:); flipud(upper(:))];
        fill(xpoly, ypoly, c, 'EdgeColor','none', 'FaceAlpha',0.3, 'HandleVisibility','off');
        line('XData', f_plot, 'YData', X_plot, 'Color', c, 'LineWidth', w);
    end

    else
        error('Xerr must be scalar or 2×N');
    end
end

xlabel('f');
if strcmp(plt,'l')
    ylabel('10·log_{10}(X)');
else
    ylabel('X');
end

hold off
end
