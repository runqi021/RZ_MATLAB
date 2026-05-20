%% Euler Polar Plot — Oscillation Demo
%  Animated phasor on polar plot with real-time amplitude & phase traces
%  Close the figure to stop.

clear; clc;

%% ---- Parameters (edit these) ----
freq_Hz   = 1;        % oscillation frequency (Hz)
duration  = 5;        % total duration (s)
dt        = 0.0005;     % time step (s) — controls animation speed
radius    = 1;        % phasor amplitude

%% ---- Precompute ----
t_all   = 0:dt:duration;
nFrames = numel(t_all);
omega   = 2*pi*freq_Hz;
phase_all = omega * t_all;                    % unwrapped phase (rad)
amp_all   = radius * cos(omega * t_all);     % real part = amplitude projection

%% ---- Set up figure ----
fig = figure('Name','Euler Phasor Demo','Color','w', ...
    'Position',[100 100 1100 500]);

% --- Left: polar phasor ---
ax1 = polaraxes('Position',[0.02 0.12 0.28 0.78]);
hold(ax1,'on');
ax1.RLim = [0 radius*1.15];
ax1.ThetaZeroLocation = 'top';
ax1.ThetaDir = 'clockwise';
title(ax1, '$A \cdot e^{i\omega t} \quad (r = A,\; \theta = \omega t)$', 'Interpreter','latex');

hLine  = polarplot(ax1, [0 0], [0 radius], 'b-', 'LineWidth', 2);
hDot   = polarplot(ax1, 0, radius, 'ro', 'MarkerSize', 10, ...
    'MarkerFaceColor','r');
hTrail = polarplot(ax1, NaN, NaN, 'b-', 'LineWidth', 0.5, ...
    'Color', [0.6 0.6 1 0.4]);

% --- Middle: amplitude (real projection) ---
ax2 = axes('Position',[0.38 0.15 0.27 0.72]);
hold(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'Time (s)');
ylabel(ax2,'$x(t)$', 'Interpreter','latex');
title(ax2, '$x(t) = \mathrm{Re}\{A \cdot e^{i\omega t}\} = A\cos(\omega t)$', 'Interpreter','latex');
xlim(ax2,[0 duration]); ylim(ax2,[-radius*1.2 radius*1.2]);
hAmp = plot(ax2, NaN, NaN, 'b-', 'LineWidth', 1.5);
hAmpDot = plot(ax2, NaN, NaN, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');

% --- Right: phase ---
ax3 = axes('Position',[0.72 0.15 0.27 0.72]);
hold(ax3,'on'); grid(ax3,'on');
xlabel(ax3,'Time (s)');
ylabel(ax3,'$\theta(t)$ (rad)', 'Interpreter','latex');
title(ax3, '$\theta(t) = \mathrm{Arg}(A \cdot e^{i\omega t}) = \omega t \;\mathrm{mod}\; 2\pi$', 'Interpreter','latex');
xlim(ax3,[0 duration]); ylim(ax3,[0 2*pi]);
yticks(ax3, [0 pi/2 pi 3*pi/2 2*pi]);
yticklabels(ax3, {'0','\pi/2','\pi','3\pi/2','2\pi'});
hPhase = plot(ax3, NaN, NaN, 'm-', 'LineWidth', 1.5);
hPhaseDot = plot(ax3, NaN, NaN, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');

%% ---- Precompute one full cycle for the polar trail (unit circle) ----
nCycle      = round(1/freq_Hz / dt);  % frames per cycle
cycle_theta = linspace(0, 2*pi, nCycle);
cycle_r     = radius * ones(1, nCycle);

%% ---- Animate (loops forever — close the figure to stop) ----
% Use a sliding window so the amplitude/phase traces scroll continuously
winSec  = 1/freq_Hz * 3;             % show 3 full cycles in the traces
winN    = round(winSec / dt);

% Draw the full unit-circle trail once
hTrail.ThetaData = cycle_theta;
hTrail.RData     = cycle_r;

% Preallocate rolling buffers
buf_t   = NaN(1, winN);
buf_amp = NaN(1, winN);
buf_ph  = NaN(1, winN);

k = 0;  % global frame counter
while isvalid(fig)
    k = k + 1;
    t_now   = (k-1) * dt;
    theta_unwrap = omega * t_now;              % unwrapped phase for trace
    theta_polar  = mod(theta_unwrap, 2*pi);   % wrapped for polar phasor
    amp_k        = radius * cos(theta_unwrap);

    % ---- polar phasor ----
    hLine.ThetaData = [0 theta_polar];
    hLine.RData     = [0 radius];
    hDot.ThetaData  = theta_polar;
    hDot.RData      = radius;

    % ---- shift rolling buffers ----
    buf_t   = [buf_t(2:end),   t_now];
    buf_amp = [buf_amp(2:end), amp_k];
    buf_ph  = [buf_ph(2:end),  theta_polar];

    % ---- amplitude trace (scrolling window) ----
    hAmp.XData    = buf_t;
    hAmp.YData    = buf_amp;
    hAmpDot.XData = t_now;
    hAmpDot.YData = amp_k;
    ax2.XLim = [t_now - winSec, t_now + dt];

    % ---- phase trace (scrolling window) ----
    hPhase.XData    = buf_t;
    hPhase.YData    = buf_ph;
    hPhaseDot.XData = t_now;
    hPhaseDot.YData = theta_polar;
    ax3.XLim = [t_now - winSec, t_now + dt];

    drawnow limitrate;
end
