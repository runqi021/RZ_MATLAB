function out = dFF_RZ(F, varargin)
% F_dff           = outDFF.F_dff;
% t_dff           = outDFF.t_dff;
% dFF             = outDFF.dFF;
% dFF_oasis_deconv = outDFF.dFF_oasis_deconv;
% spikes_oasis     = outDFF.spikes_oasis;

% dFF_RZ  Compute dF/F and (optionally) OASIS AR(1) deconvolution
%
% fps = 30;
% 
% outDFF = dFF_RZ(F, t, ...
%     'FPS', fps, ...
%     'DropFirstSec', 1, ...
%     'BaselineWinSec', 20, ...
%     'UseOASIS', true, ...
%     'TauDecay', 0.45, ...
%     'PythonExe', 'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe');
% 
% Inputs:
%   F  : T x N fluorescence (one column per ROI)
%   t  : T x 1 time vector (seconds). If empty, constructed from FPS.
%
% Params (name/value):
%   'FPS'            : imaging rate (Hz). Default 30.
%   'DropFirstSec'   : how many seconds to drop from start. Default 1.
%   'BaselineWinSec' : baseline window (seconds) for movmedian. Default 20.
%   'UseOASIS'       : true/false, run OASIS AR(1) deconv. Default false.
%   'TauDecay'       : calcium decay (s) for AR(1). Default 0.45.
%   'PythonExe'      : Python exe with oasis-deconv installed.
%
% Output struct fields:
%   out.F_dff, out.t_dff, out.dFF
%   out.params (all settings)
%   If UseOASIS:
%     out.F_oasis_deconv
%     out.dFF_oasis_deconv
%     out.spikes_oasis
%     out.baseline_oasis
%     out.g_AR1

    %% Input
    if nargin < 1
        error('Usage: out = dFF_RZ(F, ...);');
    end
    
    p = inputParser;
    addParameter(p, 'FPS', 30);
    addParameter(p, 'DropFirstSec', 0);        % seconds
    addParameter(p, 'BaselineWinSec', 20);     % seconds
    addParameter(p, 'UseOASIS', false);
    % GCaMP8s half decay time: 310ms --> Tau = ~450ms
    addParameter(p, 'TauDecay', 0.45);         % seconds
    addParameter(p, 'PythonExe', ...
        'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe');
    parse(p, varargin{:});
    P = p.Results;

    fps = P.FPS;
   
    [T, N] = size(F);
    t = (0:T-1)'/fps;
    %% dFF Math
    nDrop   = round(P.DropFirstSec * fps);          % frames to drop
    win     = round(P.BaselineWinSec * fps);        % baseline window (frames)
    if mod(win,2) == 0
        win = win + 1;  % enforce odd
    end

    if T <= nDrop
        error('Data too short: T=%d <= nDrop=%d.', T, nDrop);
    end

    F_dff = F;
    t_dff = t;
    F_dff(1:nDrop,:) = [];
    t_dff(1:nDrop)   = [];

    [T2, N2] = size(F_dff);
    fprintf('  [dFF_RZ] After dropping %d frames: T=%d, N=%d\n', nDrop, T2, N2);

    % sliding median baseline
    F0 = movmedian(F_dff, win, 1);
    F0(F0 <= 0) = eps;
    dFF = (F_dff - F0) ./ F0;

    %% OASIS AR(1) deconvulution
    F_oasis_deconv   = [];
    dFF_oasis_deconv = [];
    spikes_oasis     = [];
    baseline_oasis   = [];
    g_AR1            = [];

    if P.UseOASIS
        tau_decay = P.TauDecay;   % seconds
        % AR(1) coefficient; keep your rounding if you like
        g = round(exp(-(1/fps)/tau_decay), 2);
        g_AR1 = g;

        fprintf('  [dFF_RZ] Running OASIS AR(1), g = %.2f on %d ROIs...\n', g, N2);

        [F_oasis_deconv, dFF_oasis_deconv, spikes_oasis, baseline_oasis] = ...
            oasis_deconv_and_dff_AR1(F_dff, g, P.PythonExe);

        % enforce ~zero-mean per ROI, then detrend
        dFF_oasis_deconv = dFF_oasis_deconv - mean(dFF_oasis_deconv, 1);
        dFF_oasis_deconv = detrend(dFF_oasis_deconv);
    end

    %% Output
    out = struct();
    out.F_dff = F_dff;
    out.t_dff = t_dff;
    out.dFF   = dFF;

    out.params = P;
    out.params.nDropFrames    = nDrop;
    out.params.winFrames      = win;
    out.params.t_in           = t;

    if P.UseOASIS
        out.F_oasis_deconv   = F_oasis_deconv;
        out.dFF_oasis_deconv = dFF_oasis_deconv;
        out.spikes_oasis     = spikes_oasis;
        out.baseline_oasis   = baseline_oasis;
        out.g_AR1            = g_AR1;
    end
end
