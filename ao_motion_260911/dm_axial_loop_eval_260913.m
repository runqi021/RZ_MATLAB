function R = dm_axial_loop_eval_260913(sessDir)
%DM_AXIAL_LOOP_EVAL_260913  Did the closed loop actually stabilise the plane?
%
%   R = dm_axial_loop_eval_260913('D:\Data\RZ\ao_axial_260913\loop_260913_2140')
%
% THE LOOP CANNOT BE ITS OWN REFEREE. When correction works, the dz the loop
% measures goes to zero BY CONSTRUCTION -- that is what it was driving to. So
% the live log proves nothing on its own.
%
% What proves it is the LOGGED FRAMES. Every frame of the OPEN-CLOSED-OPEN run
% was written to disk. This re-measures dz on all of them, offline, against the
% same reference stack, with the estimator blind to which window was which.
% If the CLOSED window is flatter, the loop worked.
%
% Reported per window: RMS, 1-99% excursion, and the spectrum -- because a loop
% with too much delay SUPPRESSES the slow component and AMPLIFIES the fast one,
% and a single RMS can hide that by averaging the two together.
%
% Runqi Zhang / 2026-09-13.  NOT YET RUN ON REAL DATA.

assert(exist(sessDir,'dir') == 7, 'no such session: %s', sessDir);
M = load(fullfile(sessDir,'loop_meta.mat'));
fprintf('\n=== %s ===\n', sessDir);
fprintf('channel %d | %.2f Hz | loop gain %.2f | extrapolation %.2f%s\n', ...
        M.channel, M.fps, M.loopGain, M.extrapGain, ...
        ternary(M.dryRun, '   [DRY RUN -- nothing was sent]', ''));
fprintf('segments OPEN %g s | CLOSED %g s | OPEN %g s\n', M.segSeconds);

R = struct('sessDir',sessDir,'meta',M);

%% ---------- re-measure dz on the logged frames ----------
E = dm_axial_estimator_260913('prep', M.refVol, M.zRef, ...
        struct('nStrips',8,'pixFrac',0.08));
% 8 strips here, not the 1 the loop ran with: offline there is no time budget,
% and 8x the sampling resolves the motion the loop could never see.

acc = struct('dz',[],'t',[]);
zPrev = NaN;
    function feed(F, tF)
        o = dm_axial_estimator_260913('run', E, F, zPrev);
        zPrev = o.zEnd;
        acc.dz = [acc.dz, o.dz];
        acc.t  = [acc.t,  tF];
        fprintf('    %d frames\r', numel(acc.t));
    end
si_tiff_260913('stream', M.file, M.channel, 200, @feed);
fprintf('\n');

S  = E.S;
nF = numel(acc.t);
if all(isfinite(acc.t)) && range(acc.t) > 0
    tF = acc.t - acc.t(1);
    fps = (nF-1)/range(acc.t);
else
    fps = M.fps;
    tF = (0:nF-1)/fps;
end
t  = reshape(tF(:).' + ((1:S)'-0.5)/S/fps, [], 1);
dz = acc.dz(:);
good = isfinite(dz);
fsFine = (numel(t)-1)/range(t);

%% ---------- split into the three windows ----------
e = M.tEdges;
win = {'OPEN 1', t < e(1); ...
       'CLOSED', t >= e(1) & t < e(2); ...
       'OPEN 2', t >= e(2)};

% Drop a settling margin at each boundary: the loop needs a few frames to pull
% in after switching on, and counting that as "closed" understates it.
marg = 1.0;
for i = 1:3
    win{i,2} = win{i,2} & good;
end
win{2,2} = win{2,2} & t > e(1)+marg;
win{3,2} = win{3,2} & t > e(2)+marg;

fprintf('\n%-9s %7s %10s %10s %10s\n','window','n','RMS um','1-99%% um','peak Hz');
stats = struct();
for i = 1:3
    v = dz(win{i,2});
    v = v - mean(v);
    if numel(v) < 20, fprintf('%-9s  too few frames\n', win{i,1}); continue; end
    [fp, spec, fAx] = spec_peak(v, fsFine);
    fprintf('%-9s %7d %10.3f %10.3f %10.2f\n', win{i,1}, numel(v), ...
            sqrt(mean(v.^2)), qtile(v,0.99)-qtile(v,0.01), fp);
    stats(i).name = win{i,1};
    stats(i).rms  = sqrt(mean(v.^2));
    stats(i).exc  = qtile(v,0.99)-qtile(v,0.01);
    stats(i).fPeak= fp;
    stats(i).spec = spec;
    stats(i).specF= fAx;
    stats(i).v    = v;
end

%% ---------- the verdict ----------
if numel(stats) >= 3 && ~isempty(stats(2).rms)
    openRms = mean([stats(1).rms, stats(3).rms]);
    ratio   = stats(2).rms / openRms;
    fprintf('\n  open (mean of both) %.3f um  ->  closed %.3f um\n', openRms, stats(2).rms);
    fprintf('  RESIDUAL %.0f%% of the uncorrected motion\n', 100*ratio);
    if M.dryRun
        fprintf('  >> DRY RUN: nothing was sent, so this SHOULD be ~100%%.\n');
        fprintf('     If it is not, the windows differ for some other reason\n');
        fprintf('     (drift, bleaching, the animal) and that confounds the real run.\n');
    elseif ratio < 0.5
        fprintf('  >> IT WORKS. The plane is measurably more stable with the loop closed.\n');
    elseif ratio < 0.9
        fprintf('  >> PARTIAL. Real but modest. Most likely delay-limited --\n');
        fprintf('     check the band-by-band table below before adding loop gain.\n');
    elseif ratio <= 1.1
        fprintf('  >> NO EFFECT. Check: was DRY_RUN off? did the DM actually move\n');
        fprintf('     (cmd RMS in the live log)? is the sign right?\n');
    else
        fprintf('  >> IT MADE IT WORSE. Two candidates, and they are distinguishable:\n');
        fprintf('     (a) SIGN INVERTED -- residual would be ~2x, broadband;\n');
        fprintf('     (b) DELAY -- slow band improves, fast band gets worse. See below.\n');
    end

    % Band split. This is what separates "too slow" from "wrong sign", and a
    % single RMS cannot: a delayed loop helps below the break-even frequency
    % and hurts above it, so the two effects partly cancel in the total.
    fBreak = 1/(6*(1/M.fps));
    fprintf('\n  by band (break-even for this loop rate is %.1f Hz):\n', fBreak);
    fprintf('%14s %12s %12s %10s\n','band','open um','closed um','ratio');
    bands = [0.3 fBreak; fBreak min(4*fBreak, fsFine/2)];
    nm = {'below break','above break'};
    for b = 1:2
        o1 = band_rms(stats(1).v, fsFine, bands(b,:));
        o3 = band_rms(stats(3).v, fsFine, bands(b,:));
        cc = band_rms(stats(2).v, fsFine, bands(b,:));
        oo = mean([o1 o3]);
        fprintf('%14s %12.3f %12.3f %9.0f%%\n', nm{b}, oo, cc, 100*cc/max(oo,eps));
    end
    R.ratio = ratio;
end

R.stats = stats;  R.dz = dz;  R.t = t;  R.good = good;  R.fsFine = fsFine;
save(fullfile(sessDir,'loop_eval.mat'), '-struct', 'R', '-v7.3');
fprintf('\nsaved %s\n', fullfile(sessDir,'loop_eval.mat'));
end

%% ============================ HELPERS ============================
function [fp, P, f] = spec_peak(v, fs)
n = 2^nextpow2(min(numel(v), 2^15));
x = v(1:min(numel(v),n));
x = (x - mean(x)) .* (0.5 - 0.5*cos(2*pi*(0:numel(x)-1).'/(numel(x)-1)));
P = abs(fft(x, n)).^2;
f = (0:n-1).' * fs / n;
h = f <= fs/2;
P = P(h);  f = f(h);
sel = f > 0.3;
[~, i] = max(P(sel));
fsub = f(sel);
fp = fsub(i);
end

function r = band_rms(v, fs, band)
n = numel(v);
X = fft(v(:) - mean(v));
f = (0:n-1).' * fs / n;
f(f > fs/2) = f(f > fs/2) - fs;
X(abs(f) < band(1) | abs(f) > band(2)) = 0;
r = sqrt(mean(real(ifft(X)).^2));
end

function q = qtile(x, p)
x = sort(double(x(isfinite(x))));
if isempty(x), q = NaN; return; end
q = x(max(1, min(numel(x), round(p*(numel(x)-1)) + 1)));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
