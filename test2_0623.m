% ================= ZERO-PHASE FILTER BREATHING =================
% Assumes:
%   brw  = breathing temperature signal, already subtracted LP < 1 Hz
%   fpsW = sampling rate, e.g. 400 Hz

fs = fpsW;

x = brw(:);

% filtfilt cannot handle NaNs
x = fillmissing(x, 'linear', 'EndValues', 'nearest');

% temperature signal should be inverted
brw_inv = -x;

% ---------------- Band-pass 2-20 Hz ----------------
ord = 4;
BP_breath = [1 20];   % Hz

[b_br, a_br] = butter(ord, BP_breath/(fs/2), 'bandpass');
brw_bp2_20 = filtfilt(b_br, a_br, brw_inv);

% ---------------- Quick plot ----------------
t_brw = (0:numel(brw_inv)-1)' / fs;

figure('Color','w'); hold on;
plot(t_brw, brw_inv, 'Color', [0.7 0.7 0.7]);
plot(t_brw, brw_bp2_20, 'r', 'LineWidth', 1.2);

xlabel('time (s)');
ylabel('inverted breathing temperature');
legend({'inverted raw breathing', 'BP 2-20 Hz breathing'});
title('Breathing temperature signal');
box off;