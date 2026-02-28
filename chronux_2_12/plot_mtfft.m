function [] = plot_mtfft(F, fs, fpass, TW);

TW_spec = TW;
params_spec.Fs     = fs;
params_spec.tapers = [TW_spec, 2*TW_spec - 1];
params_spec.pad    = 0;
params_spec.fpass  = fpass;
params_spec.err    = [2 0.05];

% Multitaper spectrum of derivative
%[Sk_vel, fk_vel, Sconfk_vel] = mtspectrumc(breath_vel, params_spec);

% Multitaper spectrum of *raw* breathing (no derivative) for inset
[Sk_raw, fk_raw, Sconfk_raw] = mtspectrumc(F, params_spec);
%%
figure;
plot(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
hold on;
plot(fk_raw, Sconfk_raw(1,:), 'k--','LineWidth', 0.4);
hold on;
plot(fk_raw, Sconfk_raw(2,:), 'k--','LineWidth', 0.4);
hold on;
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([fpass]); 
grid('on');

figure;
semilogx(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
hold on;
semilogx(fk_raw, Sconfk_raw(1,:), 'k--','LineWidth', 0.4);
hold on;
semilogx(fk_raw, Sconfk_raw(2,:), 'k--','LineWidth', 0.4);
hold on;
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([fpass]); 
grid('on');

figure;
loglog(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
hold on;
loglog(fk_raw, Sconfk_raw(1,:), 'k--','LineWidth', 0.4);
hold on;
loglog(fk_raw, Sconfk_raw(2,:), 'k--','LineWidth', 0.4);
hold on;
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([fpass]); 
grid('on');

end