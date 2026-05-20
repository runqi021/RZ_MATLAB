function [F_oasis_deconv, dFF_oasis_deconv, spikes_oasis, baseline_oasis] = ...
         oasis_deconv_and_dff_AR1(F, g, varargin)
% AR(1), instantaneous rise, not enforcing monotonicity. bound decay time.

% assume exponential decay: c(t)=c(0)e^(−t/τ)
% g = e^(−Δt/τ)
% τ=−Δt/log(g)
% at 30fps, g=0.95 --> τ=650ms​

% GCaMP8s half decay time: 310ms --> Tau = ~450ms

% oasis_deconv_and_dff_AR1
% Run OASIS AR(1) via Python and compute dF/F per ROI
%
% Inputs
%   F  : T x N fluorescence (one column per ROI)
%   g  : scalar AR(1) coefficient (0<g<1)
%
% Optional name/value:
%   'PythonExe' : path to Python with oasis-deconv installed
%                 (default: 'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe')
%
% Outputs
%   F_oasis_deconv   : T x N deconvolved calcium traces c(t)
%   dFF_oasis_deconv : T x N dF/F of c(t) using a single baseline b per ROI
%   spikes_oasis     : T x N inferred spike train
%   baseline_oasis   : 1 x N baseline b for each ROI

    [T, N] = size(F);

    % ---- parse optional PythonExe ----
    p = inputParser;
    addParameter(p, 'PythonExe', ...
        'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe');
    parse(p, varargin{:});
    pythonExe = p.Results.PythonExe;

    %% --- Set Python env for OASIS --- %%
    pe = pyenv;
    if pe.Status == "NotLoaded"
        try
            pyenv('Version', pythonExe, 'ExecutionMode', 'OutOfProcess');
        catch
            % fallback: try default mode
            pyenv('Version', pythonExe);
        end
    else
        if ~contains(string(pe.Version), pythonExe)
            warning("MATLAB Python already set to: %s\nTo change, restart MATLAB and rerun.", ...
                    pe.Version);
        end
    end

    try
        py.importlib.import_module('oasis');
        py.importlib.import_module('oasis.functions');
    catch ME
        error(['Python OASIS not available in this env.\n' ...
               'Make sure you ran:\n' ...
               '  conda activate <env>\n' ...
               '  conda install -c conda-forge oasis-deconv\n\n' ...
               'Original error:\n%s'], ME.message);
    end

    try
        np     = py.importlib.import_module('numpy');
        oasisF = py.importlib.import_module('oasis.functions');
    catch ME
        error(['Python OASIS not available.\n' ...
               'Make sure pyenv points to an environment where you ran:\n' ...
               '  conda install -c conda-forge oasis-deconv\n\n' ...
               'Original error:\n%s'], ME.message);
    end

    F_oasis_deconv   = zeros(T, N);
    dFF_oasis_deconv = zeros(T, N);
    spikes_oasis     = zeros(T, N);
    baseline_oasis   = zeros(1, N);

    % AR(1) tuple for Python
    g_tuple = py.tuple(num2cell(g(:).'));  % 1-element tuple

    for i = 1:N
        y = double(F(:, i));
        y_np = np.array(y);

        try
            % AR(1) with fixed g
            out = oasisF.deconvolve(y_np, g_tuple, pyargs('penalty', int32(0)));
            c_py = out{1};
            s_py = out{2};
            b_py = out{3};

            c = double(py.array.array('d', c_py.tolist()));
            s = double(py.array.array('d', s_py.tolist()));
            b = double(b_py);

        catch ME
            warning('OASIS AR(1) failed on ROI %d: %s\nUsing fallback.', i, ME.message);
            c = y;
            s = zeros(size(y));

            finiteY = y(isfinite(y));
            if isempty(finiteY)
                b = 1e-6;
            else
                b = prctile(finiteY, 10);
            end
        end

        if ~isfinite(b) || b <= 0
            b = max(1e-6, prctile(c, 10));
        end

        F_oasis_deconv(:, i)   = c(:);
        dFF_oasis_deconv(:, i) = c / b;
        spikes_oasis(:, i)     = s(:);
        baseline_oasis(1, i)   = b;
    end
end
