svd_eigspectra_tiled_perZ

%%
function svd_eigspectra_tiled_perZ()
% svd_eigspectra_tiled_perZ
% ------------------------------------------------------------
% You have a Z-stack TIFF with 50 frames per Z.
% For EACH Z slice-group:
%   - build Y = [pixels x time] (time = 50 frames)
%   - optional demean per pixel
%   - compute eigenvalues lambda_k of (Y'Y)  (size 50x50)
%   - plot ONE curve: eigenvalue (y) vs mode index (x)
% Then TILE all slices into one big figure (many small panels).
%
% NO background/dark correction. NO mode maps. Just the spectra.
% ------------------------------------------------------------

%% ---------------- USER PARAMS ----------------
tiffPath   = "C:\Users\Admin\Downloads\shi_250115\wt_fitc\live_not_sst_fitc_1.2x_930nm_7-26lp_z220_-10_100f_col01_row01_x-900_y101_00001_minusBgPerZ_drop1_cropLR20.tif";

framesPerZ = 99;
dropFirstN = 0;              % within each Z group
demeanPerPixel = true;       % recommended

% Z axis labeling (optional)
zStart_um = 0;
zEnd_um   = 220;
zStep_um  = 10;

% Plot
K = 20;                      % # modes to plot (<= 50-dropFirstN)
useSemilogy = true;          % semilogy makes "falloff" readable
normalizeMode = "none";      % "none" | "sum" | "first"
                             % "sum"   -> lambda/sum(lambda) (variance fraction)
                             % "first" -> lambda/lambda(1)   (shape only)

% Save
saveFig = true;

%% ---------------- METADATA ----------------
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
K = min(K, T);

z_expected = (zStart_um:zStep_um:zEnd_um);
if numel(z_expected) == nZ
    z_um = z_expected(:);
else
    z_um = linspace(zStart_um, zEnd_um, nZ).';
end

fprintf("[perZ-SVD] Pages=%d, HxW=%dx%d, framesPerZ=%d, nZ=%d, T=%d\n", ...
    nPages, H, W, framesPerZ, nZ, T);

%% ---------------- COMPUTE EIGENVALUES PER Z ----------------
lam_all = nan(nZ, T);               % eigenvalues of Y'Y, sorted desc

t = Tiff(tiffPath,'r');
cleanupObj = onCleanup(@() t.close());

Y = zeros(H*W, T, 'single');        % reuse buffer

for iz = 1:nZ
    basePage = (iz-1)*framesPerZ;

    kk = 0;
    for f = 1:framesPerZ
        setDirectory(t, basePage + f);
        fr = single(t.read());
        if f <= dropFirstN
            continue;
        end
        kk = kk + 1;
        Y(:,kk) = fr(:);
    end

    if demeanPerPixel
        Yd = Y - mean(Y, 2);        % single
    else
        Yd = Y;
    end

    % Fast spectrum: C = Y'Y (T x T)
    C = double(Yd' * Yd);           % keep multiply in single, cast result to double
    lam = eig(C);
    lam = sort(real(lam), 'descend');
    lam(lam < 0) = 0;               % numeric cleanup
    lam_all(iz,:) = lam(:)';

    if mod(iz, max(1,floor(nZ/10)))==0 || iz==nZ
        fprintf("[perZ-SVD] iz=%d/%d  lam1=%.3g\n", iz, nZ, lam(1));
    end
end

%% ---------------- PREPARE PLOT VALUES ----------------
M = lam_all(:,1:K);

switch lower(string(normalizeMode))
    case "sum"
        M = M ./ max(sum(M,2), eps);
        yLabelStr = "lambda_k / sum(lambda)";
    case "first"
        M = M ./ max(M(:,1), eps);
        yLabelStr = "lambda_k / lambda_1";
    otherwise
        yLabelStr = "lambda_k (eigenvalue of Y^T Y)";
end

% global y-lims for consistent panels
if useSemilogy
    ymin = max(min(M(M>0), [], 'all'), eps);
    ymax = max(M, [], 'all');
else
    ymin = min(M, [], 'all');
    ymax = max(M, [], 'all');
end

%% ---------------- TILED PLOTS (one per Z) ----------------
% Choose a near-square grid
nCols = ceil(sqrt(nZ));
nRows = ceil(nZ / nCols);

fig = figure('Color','w', 'Position',[50 50 1800 950]);
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

for iz = 1:nZ
    ax = nexttile(tl);
    y = M(iz,:);
    x = 1:K;

    if useSemilogy
        semilogy(ax, x, max(y, eps), 'k-', 'LineWidth', 0.7);
        ylim(ax, [ymin ymax]);
    else
        plot(ax, x, y, 'k-', 'LineWidth', 0.7);
        ylim(ax, [ymin ymax]);
    end
    xlim(ax, [1 K]);

    % minimal clutter
    ax.TickDir = 'out';
    ax.Box = 'on';
    ax.FontSize = 6;

    % Title per panel
    title(ax, sprintf('z=%.0f', z_um(iz)), 'FontSize', 7);

    % only label outer axes to reduce mess
    if iz <= (nCols)         % top row
        % no special
    end
    if iz ~= 1 && iz ~= nCols && iz ~= (nRows-1)*nCols+1
        ax.XTickLabel = {};
        ax.YTickLabel = {};
    end
end

sgtitle(tl, sprintf('Per-Z eigenvalue falloff (T=%d frames/Z, demean=%d, norm=%s)', ...
    T, demeanPerPixel, normalizeMode), 'FontSize', 12);

% Add one invisible axes label block (cleaner than labeling each)
axL = axes(fig, 'Visible','off');
axL.XLabel.Visible = 'on';
axL.YLabel.Visible = 'on';
xlabel(axL, 'mode index k');
ylabel(axL, yLabelStr);

%% ---------------- SAVE ----------------
[folder, base, ~] = fileparts(char(tiffPath));
outMat = fullfile(folder, base + "_eigPerZ.mat");
save(outMat, 'lam_all', 'z_um', 'framesPerZ', 'dropFirstN', 'demeanPerPixel', 'normalizeMode', '-v7.3');
fprintf("Saved MAT: %s\n", outMat);

if saveFig
    outPng = fullfile(folder, base + "_eigPerZ_tiled.png");
    exportgraphics(fig, outPng, 'Resolution', 250);
    fprintf("Saved FIG: %s\n", outPng);
end

end
