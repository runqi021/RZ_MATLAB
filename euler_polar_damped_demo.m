%% Euler Polar Plot — Damped Harmonic Oscillator (correct physics)
%  ẍ + 2γẋ + ω₀²x = 0,  x(0) = A₀,  v(0) = 0
%  Three regimes: underdamped / critically damped / overdamped
%  Close the figure to stop. Restarts automatically.

clear; clc;

%% ---- Parameters (edit these) ----
omega0  = 2*pi*1;   % natural frequency ω₀ (rad/s) — try 2*pi*1 = 1 Hz
gamma   = 6.28;      % damping coefficient γ (1/s)
                     %   γ < ω₀  → underdamped  (oscillates + decays)
                     %   γ = ω₀  → critically damped (fastest no-oscillation decay)
                     %   γ > ω₀  → overdamped  (sluggish, slow return)
A0      = 1;         % initial displacement
duration = 10;       % seconds per run before restart
dt       = 0.0005;     % time step (s)

%% ---- Determine regime ----
if gamma < omega0 - 1e-6
    regime = 'underdamped';
    omega_d = sqrt(omega0^2 - gamma^2);   % damped frequency
elseif gamma > omega0 + 1e-6
    regime = 'overdamped';
    beta = sqrt(gamma^2 - omega0^2);
else
    regime = 'critical';
end

%% ---- Analytical solution: x(t) and dx/dt for x(0)=A₀, v(0)=0 ----
switch regime
    case 'underdamped'
        % x(t) = A₀·e^{-γt}·[cos(ω_d t) + (γ/ω_d)·sin(ω_d t)]
        x_fun = @(t) A0*exp(-gamma*t) .* (cos(omega_d*t) + (gamma/omega_d)*sin(omega_d*t));
        env_fun = @(t) A0*exp(-gamma*t) .* sqrt(1 + (gamma/omega_d)^2);
        regime_str = sprintf('Underdamped: \\gamma=%.2f < \\omega_0=%.2f', gamma, omega0);
    case 'critical'
        % x(t) = A₀·(1 + γt)·e^{-γt}
        x_fun = @(t) A0*(1 + gamma*t) .* exp(-gamma*t);
        env_fun = @(t) abs(x_fun(t));
        regime_str = sprintf('Critically damped: \\gamma=%.2f = \\omega_0=%.2f', gamma, omega0);
    case 'overdamped'
        % x(t) = A₀·e^{-γt}·[cosh(βt) + (γ/β)·sinh(βt)]
        x_fun = @(t) A0*exp(-gamma*t) .* (cosh(beta*t) + (gamma/beta)*sinh(beta*t));
        env_fun = @(t) abs(x_fun(t));
        regime_str = sprintf('Overdamped: \\gamma=%.2f > \\omega_0=%.2f', gamma, omega0);
end

fprintf('Regime: %s\n', regime);

%% ---- Set up figure ----
fig = figure('Name','Damped Harmonic Oscillator','Color','w', ...
    'Position',[100 100 1200 550]);

% --- Top-left: polar / phase-space view ---
ax1 = polaraxes('Position',[0.02 0.10 0.28 0.78]);
hold(ax1,'on');
ax1.RLim = [0 A0*1.3];
ax1.ThetaZeroLocation = 'top';
ax1.ThetaDir = 'clockwise';
if strcmp(regime,'underdamped')
    title(ax1, {'$A_0 e^{-\gamma t} \cdot e^{i\omega_d t}$', ...
        ['$(' regime_str ')$']}, 'Interpreter','latex');
else
    title(ax1, {'Phase space', ['$(' regime_str ')$']}, 'Interpreter','latex');
end

hLine   = polarplot(ax1, [0 0], [0 A0], 'b-', 'LineWidth', 2);
hDot    = polarplot(ax1, 0, A0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
hSpiral = polarplot(ax1, NaN, NaN, '-', 'LineWidth', 1, 'Color',[0.4 0.4 1 0.5]);

% --- Top-right: signal x(t) ---
ax2 = axes('Position',[0.38 0.58 0.58 0.34]);
hold(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'Time (s)');
ylabel(ax2,'$x(t)$', 'Interpreter','latex');
title(ax2, '$\ddot{x} + 2\gamma\dot{x} + \omega_0^2 x = 0$', 'Interpreter','latex');
ylim(ax2,[-A0*1.2 A0*1.2]);
hEnvUp  = plot(ax2, NaN, NaN, 'r--', 'LineWidth', 1);
hEnvDn  = plot(ax2, NaN, NaN, 'r--', 'LineWidth', 1);
hSig    = plot(ax2, NaN, NaN, 'b-',  'LineWidth', 1.5);
hSigDot = plot(ax2, NaN, NaN, 'ro',  'MarkerSize', 8, 'MarkerFaceColor','r');

% --- Bottom-left: |x(t)| or envelope ---
ax3 = axes('Position',[0.38 0.10 0.27 0.34]);
hold(ax3,'on'); grid(ax3,'on');
xlabel(ax3,'Time (s)');
ylabel(ax3,'Envelope', 'Interpreter','latex');
if strcmp(regime,'underdamped')
    title(ax3, '$A_0 e^{-\gamma t}\sqrt{1+(\gamma/\omega_d)^2}$', 'Interpreter','latex');
else
    title(ax3, '$|x(t)|$', 'Interpreter','latex');
end
ylim(ax3,[0 A0*1.3]);
hAbs    = plot(ax3, NaN, NaN, 'r-',  'LineWidth', 1.5);
hAbsDot = plot(ax3, NaN, NaN, 'ro',  'MarkerSize', 8, 'MarkerFaceColor','r');

% --- Bottom-right: phase (only meaningful for underdamped) ---
ax4 = axes('Position',[0.72 0.10 0.25 0.34]);
hold(ax4,'on'); grid(ax4,'on');
xlabel(ax4,'Time (s)');
if strcmp(regime,'underdamped')
    ylabel(ax4,'$\theta(t)$ (rad)', 'Interpreter','latex');
    title(ax4, '$\omega_d t \;\mathrm{mod}\; 2\pi$', 'Interpreter','latex');
    ylim(ax4,[0 2*pi]);
    yticks(ax4, [0 pi/2 pi 3*pi/2 2*pi]);
    yticklabels(ax4, {'0','\pi/2','\pi','3\pi/2','2\pi'});
else
    ylabel(ax4,'$x(t)$', 'Interpreter','latex');
    title(ax4, 'No oscillation — phase undefined', 'Interpreter','latex');
    ylim(ax4,[-A0*1.2 A0*1.2]);
end
hPh    = plot(ax4, NaN, NaN, 'm-',  'LineWidth', 1.5);
hPhDot = plot(ax4, NaN, NaN, 'ro',  'MarkerSize', 8, 'MarkerFaceColor','r');

%% ---- Animate ----
while isvalid(fig)

    nFrames  = round(duration / dt);
    buf_t    = NaN(1, nFrames);
    buf_sig  = NaN(1, nFrames);
    buf_env  = NaN(1, nFrames);
    buf_ph   = NaN(1, nFrames);
    sp_theta = NaN(1, nFrames);
    sp_r     = NaN(1, nFrames);

    hSpiral.ThetaData = NaN;
    hSpiral.RData     = NaN;

    for k = 1:nFrames
        if ~isvalid(fig), return; end

        t_now = (k-1) * dt;
        sig_k = x_fun(t_now);
        env_k = env_fun(t_now);

        % ---- polar phasor ----
        if strcmp(regime,'underdamped')
            theta_k = mod(omega_d * t_now, 2*pi);
            r_k = env_k;
        else
            % no rotation — phasor just shrinks along θ=0
            theta_k = 0;
            r_k = sig_k;
        end

        hLine.ThetaData = [0 theta_k];
        hLine.RData     = [0 r_k];
        hDot.ThetaData  = theta_k;
        hDot.RData      = r_k;

        sp_theta(k) = theta_k;
        sp_r(k)     = r_k;
        hSpiral.ThetaData = sp_theta(1:k);
        hSpiral.RData     = sp_r(1:k);

        % ---- buffers ----
        buf_t(k)   = t_now;
        buf_sig(k) = sig_k;
        buf_env(k) = env_k;

        % ---- signal + envelope ----
        hSig.XData    = buf_t(1:k);
        hSig.YData    = buf_sig(1:k);
        hSigDot.XData = t_now;
        hSigDot.YData = sig_k;
        if strcmp(regime,'underdamped')
            hEnvUp.XData = buf_t(1:k);  hEnvUp.YData =  buf_env(1:k);
            hEnvDn.XData = buf_t(1:k);  hEnvDn.YData = -buf_env(1:k);
        end
        ax2.XLim = [0 duration];

        % ---- envelope / |x| ----
        hAbs.XData    = buf_t(1:k);
        hAbs.YData    = buf_env(1:k);
        hAbsDot.XData = t_now;
        hAbsDot.YData = env_k;
        ax3.XLim = [0 duration];

        % ---- phase ----
        if strcmp(regime,'underdamped')
            buf_ph(k) = theta_k;
        else
            buf_ph(k) = sig_k;   % just mirror x(t) since phase is undefined
        end
        hPh.XData    = buf_t(1:k);
        hPh.YData    = buf_ph(1:k);
        hPhDot.XData = t_now;
        hPhDot.YData = buf_ph(k);
        ax4.XLim = [0 duration];

        drawnow limitrate;
    end

    pause(1);
end
