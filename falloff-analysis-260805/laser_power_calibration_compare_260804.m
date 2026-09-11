function laser_power_calibration_compare_260804()
% LASER_POWER_CALIBRATION_COMPARE_260804  Analyse the %-to-mW curves alone.
%
% No imaging data.  This looks only at the two measured calibrations held in
% laser_power_calibration.m -- the pre-2026-07-23 one and the post-2026-07-23
% one -- and asks what they say about the rig and about anything already
% analysed with the wrong one.
%
% Two-photon signal goes as P^2, so every error here is squared downstream.
% That is the reason to look at the curves on their own before using them.
%
% PANELS
%   A  mW vs %, both tables, measured points marked
%   B  the same on log-log -- a straight line here means a power law, and the
%      low end is where the two tables disagree most
%   C  new/old ratio over the shared range, in mW and (squared) in signal
%   D  incremental efficiency dmW/d% -- superlinear low end, rollover high end
%   E  relative two-photon yield (mW^2), normalised at 20 %
%   F  inverse curve: what setpoint delivers a target mW
%
% USAGE
%   laser_power_calibration_compare_260804

%% --------------------------- USER PARAMETERS -------------------------------
MARK_PCT   = [11 21 35 48];      % setpoints to annotate (the 260728 roi1 rounds)
MARK_LABEL = '260728 roi1 rounds';
OUT_DIR    = 'C:\Users\Admin\Desktop\RZ_MATLAB';
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);              % laser_power_calibration lives here
addpath(fileparts(thisDir));   % repo root, unconditionally

% Pull both tables straight out of the single source of truth.
[~, calOld] = laser_power_calibration(10, datetime(2026,1,1));
[~, calNew] = laser_power_calibration(10, datetime(2026,8,1));

fprintf('\n%s  (%s)\n', calOld.name, 'used for data acquired BEFORE 2026-07-23');
fprintf('  %d points, %g - %g %%, %g - %g mW\n', numel(calOld.pct), ...
    min(calOld.pct), max(calOld.pct), min(calOld.mW), max(calOld.mW));
fprintf('%s  (%s)\n', calNew.name, 'used for data acquired ON/AFTER 2026-07-23');
fprintf('  %d points, %g - %g %%, %g - %g mW\n', numel(calNew.pct), ...
    min(calNew.pct), max(calNew.pct), min(calNew.mW), max(calNew.mW));

%% --- dense grids -------------------------------------------------------------
gOld = linspace(min(calOld.pct), max(calOld.pct), 2000).';
gNew = linspace(min(calNew.pct), max(calNew.pct), 2000).';
yOld = interp1(calOld.pct, calOld.mW, gOld, 'pchip');
yNew = interp1(calNew.pct, calNew.mW, gNew, 'pchip');

shLo = max(min(calOld.pct), min(calNew.pct));
shHi = min(max(calOld.pct), max(calNew.pct));
gSh  = linspace(shLo, shHi, 2000).';
rSh  = interp1(calNew.pct, calNew.mW, gSh, 'pchip') ./ ...
       interp1(calOld.pct, calOld.mW, gSh, 'pchip');

cO = [0.30 0.30 0.75];
cN = [0.85 0.25 0.15];

%% --- the numbers that matter -------------------------------------------------
mwO = interp1(calOld.pct, calOld.mW, MARK_PCT, 'pchip');
mwN = interp1(calNew.pct, calNew.mW, MARK_PCT, 'pchip');
fprintf('\n%s -- what the two tables say\n', MARK_LABEL);
fprintf('  %6s %10s %10s %8s %12s\n', '%', 'mW old', 'mW new', 'new/old', 'signal x');
fprintf('  %s\n', repmat('-', 1, 50));
for i = 1:numel(MARK_PCT)
    fprintf('  %6g %10.1f %10.1f %8.2f %12.2f\n', MARK_PCT(i), mwO(i), mwN(i), ...
        mwN(i)/mwO(i), (mwN(i)/mwO(i))^2);
end

fprintf('\nRATIOS BETWEEN CONSECUTIVE ROUNDS -- this is what a P^2 correction divides by\n');
fprintf('  %14s %10s %10s %10s\n', 'pair', 'pct^2', 'mW^2 old', 'mW^2 new');
fprintf('  %s\n', repmat('-', 1, 48));
for i = 1:numel(MARK_PCT)-1
    fprintf('  %6g -> %-5g %10.2f %10.2f %10.2f\n', MARK_PCT(i), MARK_PCT(i+1), ...
        (MARK_PCT(i+1)/MARK_PCT(i))^2, (mwO(i+1)/mwO(i))^2, (mwN(i+1)/mwN(i))^2);
end

% Rollover and floor -- the two places each curve stops being usable.
[pkO, iO] = max(calOld.mW);
fprintf('\nOLD table rolls over at %g %% (%.0f mW); above that more setpoint LOSES power\n', ...
    calOld.pct(iO), pkO);
fprintf('  and reads a flat %.1f mW at both %g %% and %g %% -- that is a floor\n', ...
    calOld.mW(1), calOld.pct(1), calOld.pct(2));
fprintf('  (meter noise / Pockels leakage), not a measurement.  Unusable below ~3 %%.\n');
fprintf('NEW table stops at %g %% (%.0f mW) -- nothing above that can be quoted in mW.\n', ...
    max(calNew.pct), max(calNew.mW));

% Local power-law exponent, mW ~ pct^m.
mO = gradient(log(yOld)) ./ gradient(log(gOld));
mN = gradient(log(yNew)) ./ gradient(log(gNew));
fprintf('\nLocal exponent m in mW ~ %%^m (pure linear optics would give m = 1):\n');
for q = [2 5 10 20 35 50]
    sO = interp1(gOld, mO, q, 'linear', NaN);
    sN = interp1(gNew, mN, q, 'linear', NaN);
    fprintf('   at %5g %%   old m = %5.2f    new m = %5.2f\n', q, sO, sN);
end

%% --- figure -------------------------------------------------------------------
f = figure('Color', 'w', 'Position', [50 50 1500 850], 'Name', 'power calibration');
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- A: mW vs % ---------------------------------------------------------------
ax = nexttile; hold(ax, 'on');
plot(ax, gOld, yOld, '-', 'Color', cO, 'LineWidth', 1.6, 'DisplayName', 'pre 2026-07-23');
plot(ax, gNew, yNew, '-', 'Color', cN, 'LineWidth', 1.6, 'DisplayName', 'post 2026-07-23');
plot(ax, calOld.pct, calOld.mW, 'o', 'Color', cO, 'MarkerFaceColor', 'w', ...
    'MarkerSize', 4, 'HandleVisibility', 'off');
plot(ax, calNew.pct, calNew.mW, 'o', 'Color', cN, 'MarkerFaceColor', 'w', ...
    'MarkerSize', 4, 'HandleVisibility', 'off');
plot(ax, MARK_PCT, mwN, 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 6, ...
    'DisplayName', MARK_LABEL);
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'ScanImage setpoint (%)'); ylabel(ax, 'power at sample (mW)');
title(ax, 'A  measured calibration', 'FontWeight', 'normal');
legend(ax, 'Location', 'northwest', 'Box', 'off');

% --- B: log-log ---------------------------------------------------------------
ax = nexttile; hold(ax, 'on');
plot(ax, gOld, yOld, '-', 'Color', cO, 'LineWidth', 1.6);
plot(ax, gNew, yNew, '-', 'Color', cN, 'LineWidth', 1.6);
plot(ax, calOld.pct, calOld.mW, 'o', 'Color', cO, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
plot(ax, calNew.pct, calNew.mW, 'o', 'Color', cN, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
set(ax, 'XScale', 'log', 'YScale', 'log');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'setpoint (%)'); ylabel(ax, 'mW');
title(ax, 'B  log-log -- straight = power law', 'FontWeight', 'normal');
text(ax, 0.12, 2.6, {'old table floors at 2.3 mW', 'below ~3 % -- not real'}, ...
    'Color', cO, 'FontSize', 8);

% --- C: new/old ---------------------------------------------------------------
ax = nexttile; hold(ax, 'on');
plot(ax, gSh, rSh, 'k-', 'LineWidth', 1.6, 'DisplayName', 'mW ratio');
plot(ax, gSh, rSh.^2, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.4, ...
    'DisplayName', 'signal ratio (\propto mW^2)');
yline(ax, 1, 'k:', 'HandleVisibility', 'off');
for i = 1:numel(MARK_PCT)
    if MARK_PCT(i) >= shLo && MARK_PCT(i) <= shHi
        xline(ax, MARK_PCT(i), 'k--', sprintf('%g%%', MARK_PCT(i)), ...
            'FontSize', 8, 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
    end
end
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'setpoint (%)'); ylabel(ax, 'new / old');
title(ax, 'C  how wrong the old table is now', 'FontWeight', 'normal');
legend(ax, 'Location', 'southeast', 'Box', 'off');

% --- D: incremental efficiency ------------------------------------------------
ax = nexttile; hold(ax, 'on');
plot(ax, gOld, gradient(yOld) ./ gradient(gOld), '-', 'Color', cO, 'LineWidth', 1.6);
plot(ax, gNew, gradient(yNew) ./ gradient(gNew), '-', 'Color', cN, 'LineWidth', 1.6);
yline(ax, 0, 'k:');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'setpoint (%)'); ylabel(ax, 'dmW / d%');
title(ax, 'D  incremental efficiency (<0 = rollover)', 'FontWeight', 'normal');

% --- E: two-photon yield ------------------------------------------------------
ax = nexttile; hold(ax, 'on');
nO = interp1(calOld.pct, calOld.mW, 20, 'pchip');
nN = interp1(calNew.pct, calNew.mW, 20, 'pchip');
plot(ax, gOld, (yOld / nO).^2, '-', 'Color', cO, 'LineWidth', 1.6, 'DisplayName', 'pre');
plot(ax, gNew, (yNew / nN).^2, '-', 'Color', cN, 'LineWidth', 1.6, 'DisplayName', 'post');
set(ax, 'YScale', 'log');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'setpoint (%)'); ylabel(ax, 'relative signal (mW^2), = 1 at 20 %');
title(ax, 'E  two-photon yield vs setpoint', 'FontWeight', 'normal');
legend(ax, 'Location', 'northwest', 'Box', 'off');

% --- F: inverse ---------------------------------------------------------------
ax = nexttile; hold(ax, 'on');
% invert only the monotonic rising part of each curve (the old one rolls over)
mo = [true; diff(yOld) > 0];  mn = [true; diff(yNew) > 0];
plot(ax, yOld(mo), gOld(mo), '-', 'Color', cO, 'LineWidth', 1.6, 'DisplayName', 'pre');
plot(ax, yNew(mn), gNew(mn), '-', 'Color', cN, 'LineWidth', 1.6, 'DisplayName', 'post');
grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'desired power at sample (mW)'); ylabel(ax, 'setpoint to dial in (%)');
title(ax, 'F  inverse -- rising branch only', 'FontWeight', 'normal');
legend(ax, 'Location', 'southeast', 'Box', 'off');

title(tl, 'ScanImage setpoint -> power at sample: the two measured calibrations', ...
    'FontWeight', 'bold');

if SAVE_FIG
    base = fullfile(OUT_DIR, 'laser_power_calibration_compare_260804');
    exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor', 'white');
    exportgraphics(f, [base '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
    fprintf('\nsaved %s.png/.pdf\n', base);
end
end
