%% Visual MATLAB Color Reference Map
%  Displays named colors as a grid with [R G B] values.
%  Click any swatch to copy its RGB to clipboard & command window.

colors = {
    % --- Blues ---
    'Navy',        [0.00 0.00 0.50];
    'DarkBlue',    [0.00 0.00 0.55];
    'MediumBlue',  [0.00 0.00 0.80];
    'Blue',        [0.00 0.00 1.00];
    'RoyalBlue',   [0.25 0.41 0.88];
    'DodgerBlue',  [0.12 0.56 1.00];
    'DeepSkyBlue', [0.00 0.75 1.00];
    'SteelBlue',   [0.27 0.51 0.71];
    'CornflowerBlue',[0.39 0.58 0.93];
    'LightBlue',   [0.68 0.85 0.90];

    % --- Greens ---
    'DarkGreen',   [0.00 0.39 0.00];
    'Green',       [0.00 0.50 0.00];
    'ForestGreen', [0.13 0.55 0.13];
    'SeaGreen',    [0.18 0.55 0.34];
    'LimeGreen',   [0.20 0.80 0.20];
    'Lime',        [0.00 1.00 0.00];
    'SpringGreen', [0.00 1.00 0.50];
    'MediumAquamarine',[0.40 0.80 0.67];
    'YellowGreen', [0.60 0.80 0.20];
    'PaleGreen',   [0.60 0.98 0.60];

    % --- Reds / Oranges ---
    'DarkRed',     [0.55 0.00 0.00];
    'Maroon',      [0.50 0.00 0.00];
    'FireBrick',   [0.70 0.13 0.13];
    'Crimson',     [0.86 0.08 0.24];
    'Red',         [1.00 0.00 0.00];
    'OrangeRed',   [1.00 0.27 0.00];
    'DarkOrange',  [1.00 0.55 0.00];
    'Orange',      [1.00 0.65 0.00];
    'Coral',       [1.00 0.50 0.31];
    'Tomato',      [1.00 0.39 0.28];

    % --- Yellows / Browns ---
    'Gold',        [1.00 0.84 0.00];
    'Yellow',      [1.00 1.00 0.00];
    'Khaki',       [0.94 0.90 0.55];
    'Peru',        [0.80 0.52 0.25];
    'Chocolate',   [0.82 0.41 0.12];
    'SaddleBrown', [0.55 0.27 0.07];
    'Sienna',      [0.63 0.32 0.18];
    'Tan',         [0.82 0.71 0.55];
    'Wheat',       [0.96 0.87 0.70];
    'BurlyWood',   [0.87 0.72 0.53];

    % --- Purples / Pinks ---
    'Indigo',      [0.29 0.00 0.51];
    'DarkViolet',  [0.58 0.00 0.83];
    'Purple',      [0.50 0.00 0.50];
    'MediumOrchid',[0.73 0.33 0.83];
    'Magenta',     [1.00 0.00 1.00];
    'Violet',      [0.93 0.51 0.93];
    'DeepPink',    [1.00 0.08 0.58];
    'HotPink',     [1.00 0.41 0.71];
    'Pink',        [1.00 0.75 0.80];
    'Plum',        [0.87 0.63 0.87];

    % --- Cyans / Teals ---
    'Teal',        [0.00 0.50 0.50];
    'DarkCyan',    [0.00 0.55 0.55];
    'Cyan',        [0.00 1.00 1.00];
    'Turquoise',   [0.25 0.88 0.82];
    'Aquamarine',  [0.50 1.00 0.83];

    % --- Grays ---
    'Black',       [0.00 0.00 0.00];
    'DimGray',     [0.41 0.41 0.41];
    'Gray',        [0.50 0.50 0.50];
    'DarkGray',    [0.66 0.66 0.66];
    'Silver',      [0.75 0.75 0.75];
    'LightGray',   [0.83 0.83 0.83];
    'White',       [1.00 1.00 1.00];
};

nColors = size(colors, 1);
nCols = 10;
nRows = ceil(nColors / nCols);

fig = figure('Name','MATLAB Color Map — click a swatch to copy RGB', ...
    'NumberTitle','off', 'MenuBar','none', 'Color','w', ...
    'Units','normalized', 'Position',[0.05 0.1 0.9 0.75]);

ax = axes(fig, 'Position',[0.01 0.01 0.98 0.98]);
hold(ax, 'on');
axis(ax, 'ij');  % row 1 at top
xlim(ax, [0 nCols]);
ylim(ax, [0 nRows]);
axis(ax, 'off');

patchHandles = gobjects(nColors, 1);
for k = 1:nColors
    col = mod(k-1, nCols);
    row = floor((k-1) / nCols);
    rgb = colors{k, 2};
    name = colors{k, 1};

    % color swatch
    patchHandles(k) = patch(ax, ...
        col + [0.05 0.95 0.95 0.05], ...
        row + [0.05 0.05 0.55 0.55], ...
        rgb, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);

    % pick text color for contrast
    lum = 0.299*rgb(1) + 0.587*rgb(2) + 0.114*rgb(3);
    if lum > 0.5, tc = [0 0 0]; else, tc = [1 1 1]; end

    % name label on swatch
    text(ax, col+0.5, row+0.20, name, ...
        'HorizontalAlignment','center', 'FontSize',7, ...
        'FontWeight','bold', 'Color', tc);

    % RGB label on swatch
    rgbStr = sprintf('[%.2f %.2f %.2f]', rgb);
    text(ax, col+0.5, row+0.42, rgbStr, ...
        'HorizontalAlignment','center', 'FontSize',6, 'Color', tc);

    % store RGB in patch UserData for click callback
    patchHandles(k).UserData = rgb;
end

% Click callback: print & copy RGB
set(patchHandles, 'ButtonDownFcn', @(src,~) copyRGB(src));

fprintf('Click any color swatch to copy its [R G B] value.\n');

function copyRGB(src)
    rgb = src.UserData;
    str = sprintf('[%.2f %.2f %.2f]', rgb);
    clipboard('copy', str);
    fprintf('Copied: %s\n', str);
end
