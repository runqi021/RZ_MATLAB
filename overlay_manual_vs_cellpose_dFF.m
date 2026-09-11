%% overlay_manual_vs_cellpose_dFF.m
% Overlay a hand-drawn (ImageJ) ROI against a Cellpose ROI for the same
% recording. One figure, two stacked panels:
%   top    = raw F     (manual vs cellpose)
%   bottom = dF/F      (manual vs cellpose)
%
% Manual dF/F is computed with helper.dFF_RZ using the SAME parameters the
% cellpose pipeline used (read from the *_dFF.mat 'params'), so the two
% traces are directly comparable.

clear; clc;

%% ---------------- user params ----------------
% Hand-drawn ROI fluorescence (ImageJ "Values": col1 = slice, col2 = Mean F)
manualCsv  = "D:\ChAT_analysis\0523\roi1_4x_x-900y700z-15_6000f_12lp_00001\ROI1_largeROI_Values.csv";

% Cellpose pipeline dF/F output for this recording
cpDffMat   = "D:\ChAT_analysis\0523\roi1_4x_x-900y700z-15_6000f_12lp_00001\roi1_4x_x-900y700z-15_6000f_12lp_00001_ch1_dFF.mat";

cpROI      = 1;     % which cellpose ROI column to compare against
FPS        = 30;    % imaging rate (Hz) -- used if params has no fps
% ----------------------------------------------

addpath(fileparts(mfilename('fullpath')));   % so helper. is on path

%% ---- load cellpose F + dFF ----
cp = load(cpDffMat);                 % expects F_roi [T x N], dFF [T x N], params
assert(cpROI >= 1 && cpROI <= size(cp.dFF,2), ...
    'cpROI=%d out of range (file has %d ROIs).', cpROI, size(cp.dFF,2));

F_cp   = cp.F_roi(:, cpROI);
dFF_cp = cp.dFF(:,  cpROI);

% pull dF/F params from the pipeline so the manual trace matches exactly
p = cp.params;
DropFirstSec   = pickfield(p, {'DropFirstSec','dropFirstSec','TossSec'}, 0);
BaselineWinSec = pickfield(p, {'BaselineWinSec','baselineWinSec','baseWinSec'}, 20);
fps            = pickfield(p, {'FPS','fps','frameRate'}, FPS);

%% ---- load manual F, compute dFF with matching params ----
M = readmatrix(manualCsv);
F_man_raw = M(:,2);

% align lengths (defensive: both should be the full recording length)
n = min(numel(F_man_raw), size(cp.F_roi,1));
F_man_raw = F_man_raw(1:n);
F_cp      = F_cp(1:n);
dFF_cp    = dFF_cp(1:n);

out = helper.dFF_RZ(F_man_raw, ...
    'FPS', fps, ...
    'DropFirstSec', DropFirstSec, ...
    'BaselineWinSec', BaselineWinSec);

t      = out.t_dff;        % time (s) after dropping frames
dFF_man = out.dFF;
F_man   = out.F_dff;

% cellpose traces: drop the same lead frames so the x-axes line up
nDrop = n - numel(t);
F_cp   = F_cp(nDrop+1:end);
dFF_cp = dFF_cp(nDrop+1:end);

%% ---- plot: F (top), dFF (bottom) ----
figure('Color','w','Position',[100 100 1150 680]);

ax1 = subplot(2,1,1); hold on;
plot(t, F_man, 'k',  'LineWidth', 0.5, 'DisplayName','manual ROI');
plot(t, F_cp,  'r',  'LineWidth', 0.5, 'DisplayName',sprintf('cellpose ROI %d',cpROI));
ylabel('F (a.u.)'); title('Raw fluorescence');
legend('Location','best'); box off; xlim([t(1) t(end)]);

ax2 = subplot(2,1,2); hold on;
plot(t, dFF_man, 'k', 'LineWidth', 0.5, 'DisplayName','manual ROI');
plot(t, dFF_cp,  'r', 'LineWidth', 0.5, 'DisplayName',sprintf('cellpose ROI %d',cpROI));
ylabel('\DeltaF/F'); xlabel('Time (s)'); title('\DeltaF/F');
legend('Location','best'); box off; xlim([t(1) t(end)]);

linkaxes([ax1 ax2],'x');
sgtitle('Manual ROI vs Cellpose ROI');

%% ---- local helper: first matching field or default ----
function v = pickfield(s, names, dflt)
    v = dflt;
    if ~isstruct(s); return; end
    for k = 1:numel(names)
        if isfield(s, names{k}) && ~isempty(s.(names{k}))
            v = s.(names{k}); return;
        end
    end
end

