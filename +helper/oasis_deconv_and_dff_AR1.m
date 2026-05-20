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
%                 (auto-detected from conda envs 'oasis' or 'cellpose-gpu')
%
% Outputs
%   F_oasis_deconv   : T x N deconvolved calcium traces c(t)
%   dFF_oasis_deconv : T x N dF/F of c(t) using a single baseline b per ROI
%   spikes_oasis     : T x N inferred spike train
%   baseline_oasis   : 1 x N baseline b for each ROI

    [T, N] = size(F);

    % ---- parse optional PythonExe ----
    p = inputParser;
    addParameter(p, 'PythonExe', '');
    parse(p, varargin{:});
    pythonExe = p.Results.PythonExe;

    % Auto-detect if not provided
    if isempty(pythonExe)
        up = getenv('USERPROFILE');
        pyCands = { ...
            fullfile(up, '.conda',     'envs','cellpose-gpu','python.exe'), ...
            fullfile(up, 'anaconda3',  'envs','cellpose-gpu','python.exe'), ...
            fullfile(up, 'miniconda3', 'envs','cellpose-gpu','python.exe'), ...
            fullfile(up, '.conda',     'envs','oasis','python.exe'), ...
            fullfile(up, 'anaconda3',  'envs','oasis','python.exe'), ...
            fullfile(up, 'miniconda3', 'envs','oasis','python.exe') };
        for ii = 1:numel(pyCands)
            if isfile(pyCands{ii}), pythonExe = pyCands{ii}; break; end
        end
    end
    assert(~isempty(pythonExe) && isfile(pythonExe), ...
        'Python executable not found. Provide ''PythonExe'' or install oasis conda env.');

    %% --- Set Python env for OASIS --- %%
    pe = pyenv;
    needsSwitch = false;
    if pe.Status == "NotLoaded"
        % First time — set to our env
        try
            pyenv('Version', pythonExe, 'ExecutionMode', 'OutOfProcess');
        catch
            pyenv('Version', pythonExe);
        end
    elseif ~strcmpi(string(pe.Executable), pythonExe)
        needsSwitch = true;
        % Try to terminate OutOfProcess Python and reload
        try
            terminate(pyenv);
            pause(1);
            pyenv('Version', pythonExe, 'ExecutionMode', 'OutOfProcess');
            needsSwitch = false;
        catch
            % terminate failed (InProcess or other reason)
        end
        if needsSwitch
            error(['MATLAB Python is already loaded with:\n  %s\n\n' ...
                   'OASIS needs:\n  %s\n\n' ...
                   'Please restart MATLAB and run OASIS before any other ' ...
                   'Python calls (e.g. Cellpose).'], ...
                   string(pe.Executable), pythonExe);
        end
    end

    % Import oasis
    try
        py.importlib.import_module('oasis');
        oasisF = py.importlib.import_module('oasis.functions');
        np     = py.importlib.import_module('numpy');
    catch ME
        curPy = '';
        try curPy = string(pyenv().Executable); catch; end
        error(['Python OASIS not available in this env.\n' ...
               'Make sure you ran:\n' ...
               '  conda activate oasis\n' ...
               '  conda install -c conda-forge oasis-deconv\n\n' ...
               'Current pyenv: %s\n' ...
               'Original error:\n%s'], curPy, ME.message);
    end

    F_oasis_deconv   = zeros(T, N);
    dFF_oasis_deconv = zeros(T, N);
    spikes_oasis     = zeros(T, N);
    baseline_oasis   = zeros(1, N);

    % AR(1) tuple for Python
    g_tuple = py.tuple(num2cell(g(:).'));  % 1-element tuple

    for i = 1:N
        y = double(F(:, i)');
        y_np = np.array(y);

        % Estimate noise std in MATLAB (MAD of first differences)
        % This avoids calling estimate_parameters() inside OASIS,
        % which segfaults due to scipy.linalg.toeplitz + numpy 2.x incompatibility.
        sn = median(abs(diff(y))) / 0.6745;

        try
            % AR(1) with fixed g and pre-computed sn
            out = oasisF.deconvolve(y_np, g_tuple, pyargs('sn', sn, 'penalty', int32(0)));
            c_py = out{1};
            s_py = out{2};
            b_py = out{3};

            c = double(py.array.array('d', c_py.tolist()));
            s = double(py.array.array('d', s_py.tolist()));
            b = double(b_py);
            fprintf('ROI %d: b=%.2f, max(c)=%.2f, nnz(s)=%d\n', i, b, max(c), nnz(s));

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
