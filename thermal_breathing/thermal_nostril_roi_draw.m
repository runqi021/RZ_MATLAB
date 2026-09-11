% thermal_nostril_roi_draw.m
% DRAWING ONLY. Draw an ellipse/circle on each nostril's average projection,
% extract the RAW deg-C ROI trace from the tracking-aligned stack, and save the
% ROI params + raw traces to <stem>_nostrilROI.mat.
%
% All downstream analysis (bandpass, invert, Hilbert phase, L+R average, PSD)
% lives in thermal_nostril_breath_analyze.m -- keep THIS script drawing-only.
%
% RUN IN MATLAB YOURSELF (interactive: you draw the ROI; cannot run headless).
% No flips: imagesc shows the native .ats orientation (row 0 top).

% Give ONLY the DLC csv; the _nostrilC.mat is resolved via thermal_resolve_paths.
dlcCsv    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n4DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv";
dataRoot  = "D:\260615_thermalNbasler";
ROI_SHAPE = 'ellipse';   % 'ellipse' | 'circle'
STAT      = 'mean';     % per-frame ROI readout: 'mean' | 'median' | 'max'

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here),'mod','bluewhitered'));
P = thermal_resolve_paths(dlcCsv, dataRoot);
S = load(P.nostrilC);
sides = {'L','R'};
R = struct();

for i = 1:2
    sd  = sides{i};
    avg = double(S.(sd+"_avg"));
    stk = double(S.(sd+"_stack"));          % [T x win x win], deg C (decimated)
    sfps = double(S.(sd+"_stack_fps"));
    sidelbl = char(S.(sd+"_side"));
    win = size(avg,1);
    cc  = (win+1)/2;                          % nostril (dot) sits at window center

    figure('Color','w','Position',[120 120 620 560],'Name',sprintf('%s nostril ROI',sidelbl));
    ax = axes;
    imagesc(ax, avg); axis(ax,'image'); colormap(ax, bluewhitered(256)); colorbar(ax);
    hold(ax,'on'); plot(ax, cc, cc, 'c+','MarkerSize',10,'LineWidth',1.2);
    title(ax, sprintf('%s avg proj (\\circC) — DRAW ROI, double-click to confirm', sidelbl));

    if strcmpi(ROI_SHAPE,'circle'), h = drawcircle(ax); else, h = drawellipse(ax); end
    wait(h);                                  % blocks until you double-click the ROI
    mask = createMask(h);                     % [win x win] logical

    % apply ROI to every frame of the aligned stack -> raw deg-C trace
    T = size(stk,1);
    flat = reshape(stk, T, []);
    roipix = flat(:, mask(:));
    switch STAT
        case 'mean',   tr = mean(roipix, 2, 'omitnan');
        case 'median', tr = median(roipix, 2, 'omitnan');
        case 'max',    tr = max(roipix, [], 2, 'omitnan');
        otherwise,     error('STAT must be mean | median | max');
    end

    R.(sd).side  = sidelbl;
    R.(sd).mask  = mask;
    R.(sd).trace = tr;        % RAW ROI readout (deg C) -- analysis happens elsewhere
    R.(sd).stat  = STAT;
    R.(sd).fps   = sfps;
    R.(sd).npix  = nnz(mask);
    R.(sd).shape = ROI_SHAPE;
    if strcmpi(ROI_SHAPE,'circle')
        R.(sd).center = h.Center; R.(sd).radius = h.Radius;
    else
        R.(sd).center = h.Center; R.(sd).semiaxes = h.SemiAxes; R.(sd).angle = h.RotationAngle;
    end
    fprintf('%s: ROI %d px, raw mean %.2f C, swing %.2f C\n', ...
        sidelbl, nnz(mask), mean(tr,'omitnan'), max(tr)-min(tr));
end
R.src = char(dlcCsv);  R.nostrilC = P.nostrilC;

% quick sanity check: raw L/R traces (drawing confirmation only)
figure('Color','w','Position',[80 80 1000 360]);
for i = 1:2
    sd = sides{i}; tr = R.(sd).trace; t = (0:numel(tr)-1)/R.(sd).fps;
    subplot(2,1,i); plot(t, tr, '-'); grid on; ylabel('\circC');
    title(sprintf('%s RAW ROI (%d px)', R.(sd).side, R.(sd).npix));
end
xlabel('s'); sgtitle('drawn-ROI raw traces (run thermal\_nostril\_breath\_analyze for analysis)');

outPath = P.nostrilROI;
save(outPath, '-struct', 'R');
fprintf('saved %s\n', outPath);
