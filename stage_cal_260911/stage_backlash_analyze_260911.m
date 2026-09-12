function stage_backlash_analyze_260911(matFile)
%STAGE_BACKLASH_ANALYZE_260911  Re-derive the backlash result from raw counts.
%
% Reads a stage_backlash_260911 .mat and redoes the analysis WITHOUT relying on
% the um conversion the acquiring script used, which was wrong: it scaled the
% encoder register TP by the STEPPER quantum (positionDeviceUnits, 0.7815 um on
% x/y and 0.3125 um on z). The encoder has its own, finer resolution, so every
% um number in the original printout is off by a constant factor.
%
% Nothing here touches hardware. Load and fit only.
%
% Method:
%   1. Fit tpCnt against the COMMANDED target (um), per axis, over all data.
%      slope = encoder counts per micron -> quantum = 1/slope. This is the
%      calibration the acquiring script should have used.
%   2. Backlash per (axis, step) = mean(tpCnt DOWN) - mean(tpCnt UP), in counts.
%      Up and down visit the SAME target set, so the frame offset and any scale
%      error cancel exactly in the difference. This is why the backlash column
%      survived the unit bug while the err UP/err DOWN columns did not.
%   3. The same difference on tdCnt is a CONTROL. TD is the commanded step
%      count, so its backlash must be ~0. If it is not, something other than
%      the stage is direction-dependent and the TP result cannot be trusted.
%
% Runqi Zhang / 2026-09-11

if nargin < 1 || isempty(matFile)
    d = dir(fullfile(fileparts(mfilename('fullpath')),'stage_backlash_*.mat'));
    d = d(~contains({d.name},'analyze'));
    assert(~isempty(d),'no stage_backlash_*.mat found next to this script');
    [~,i] = max([d.datenum]);
    matFile = fullfile(d(i).folder, d(i).name);
end
fprintf('reading %s\n', matFile);
S = load(matFile);
assert(isfield(S,'R') && ~isempty(S.R), 'no R struct in that file');
R = S.R;

assert(isfield(R,'tpCnt'), ...
    ['This .mat predates the raw-count fields (tdCnt/tpCnt/rpCnt). Only the ' ...
     'mis-scaled um values were saved, so the encoder quantum cannot be ' ...
     'recovered -- the sweep would need re-running.']);

stepQ = [0.7815 -0.7815 0.3125];   % stepper quantum actually used, um/count
nm    = {'x','y','z'};

%% ---- 1. encoder scale, per axis ----
fprintf('\n================ ENCODER SCALE ================\n');
fprintf('%-4s %7s %12s %12s %12s %9s\n', ...
        'axis','n','cnt/um','quantum um','stepper um','ratio');
encQ = nan(1,3);
for ax = 1:3
    m = [R.ax] == ax;
    if nnz(m) < 3, continue; end
    t = [R(m).target].';  c = [R(m).tpCnt].';
    ok = isfinite(t) & isfinite(c);
    if nnz(ok) < 3 || std(t(ok)) == 0, continue; end
    p = polyfit(t(ok), c(ok), 1);      % counts = p(1)*um + p(2)
    encQ(ax) = 1/p(1);
    fprintf('%-4s %7d %12.4f %12.5f %12.4f %9.2f\n', ...
            nm{ax}, nnz(ok), p(1), encQ(ax), stepQ(ax), stepQ(ax)/encQ(ax));
end
fprintf(['\nratio = how many times too large the original um numbers were.\n' ...
         'A ratio near 1 would mean the stepper quantum was right after all.\n']);

%% ---- 2. backlash per condition, in counts and corrected um ----
fprintf('\n================ BACKLASH ================\n');
fprintf('%-4s %8s %7s %11s %11s %11s %10s\n', ...
        'axis','step','n','TP counts','TP um','TD counts','sd cnt');
axes_ = unique([R.ax]);
out = struct('ax',{},'step',{},'blCnt',{},'blUm',{},'tdCnt',{},'sd',{});
for ax = axes_(:).'
    steps = unique([R([R.ax]==ax).step]);
    for s = steps(:).'
        m  = [R.ax]==ax & [R.step]==s;
        up = m & [R.dir] > 0;   dn = m & [R.dir] < 0;
        if ~any(up) || ~any(dn), continue; end

        blCnt = mean([R(dn).tpCnt]) - mean([R(up).tpCnt]);
        tdBl  = mean([R(dn).tdCnt]) - mean([R(up).tdCnt]);
        blUm  = blCnt * encQ(ax);
        % Repeatability = scatter about each direction's own LINE, not about its
        % mean. The raw std just measures how far apart the targets are (it came
        % out as exactly targetSpread x counts/um), which says nothing about the
        % stage.
        sd = pooled_resid(R, up, dn);

        fprintf('%-4s %8.2f %7d %11.2f %11.4f %11.2f %10.2f\n', ...
                nm{ax}, s, nnz(m), blCnt, blUm, tdBl, sd);
        out(end+1) = struct('ax',ax,'step',s,'blCnt',blCnt,'blUm',blUm, ...
                            'tdCnt',tdBl,'sd',sd); %#ok<AGROW>
    end
end

%% ---- 2b. is the STEP SIZE itself direction-dependent? ----
% Distinct from backlash. Backlash is a one-off offset AT the reversal; a
% direction-dependent step size means every +step differs in length from every
% -step, so the error ACCUMULATES along a row instead of being a constant
% offset. In a snake raster these have different signatures: backlash shifts a
% whole row, an asymmetric step size shears it progressively.
%
% Measured from consecutive increments WITHIN each direction run, with the
% FIRST increment after each reversal dropped -- that one carries the backlash
% and would otherwise contaminate the step size.
fprintf('\n========== STEP ERROR by STEP SIZE and DIRECTION ==========\n');
fprintf('%-5s %9s %11s %10s %9s %6s\n', ...
        'dir','cmd um','measured um','err um','err %','n');
for ax = axes_(:).'
    for dirSign = [+1 -1]
        % '+' is 43, '-' is 45 -- ADD 2, do not subtract (41 is ')').
        lbl = sprintf('%c%s', char(43 + (dirSign<0)*2), nm{ax});
        if ~isfinite(encQ(ax))
            fprintf(['%-5s %9s %11s %10s %9s %6s   NO ENCODER -- not ' ...
                     'measurable from registers\n'], lbl,'-','-','-','-','-');
            continue
        end
        steps = unique([R([R.ax]==ax).step]);
        for s = steps(:).'
            d = incr_for(R, ax, s, dirSign);
            if isempty(d), continue; end
            meas = abs(mean(d) * encQ(ax));
            fprintf('%-5s %9.2f %11.4f %10.4f %8.2f%% %6d\n', ...
                    lbl, s, meas, meas-s, 100*(meas-s)/s, numel(d));
        end
    end
end
fprintf(['\n"measured" is the length of ONE step in that direction, from\n' ...
         'consecutive encoder readings, with the first step after each\n' ...
         'reversal dropped (that one carries the backlash, not the step size).\n' ...
         'Compare +ax against -ax at the same cmd size for a direction effect.\n']);

%% ---- 3. verdict ----
fprintf('\n================ SUMMARY ================\n');
for ax = axes_(:).'
    o = out([out.ax] == ax);
    if isempty(o), continue; end

    % An axis whose TP never moves has no readout at all. Say so instead of
    % reporting a confident zero, which is the exact failure this whole script
    % was built to avoid.
    if ~isfinite(encQ(ax)) || all([o.blCnt] == 0)
        fprintf(['%s: TP NEVER CHANGED -- no encoder readout on this axis.\n' ...
                 '   The zeros above are absence of measurement, NOT absence of\n' ...
                 '   backlash. x/y backlash cannot be measured this way; it needs\n' ...
                 '   the image or a dial indicator.\n'], nm{ax});
        continue
    end

    % Small steps under-report: at a step of a few encoder counts the reversal
    % is not cleanly resolved. Quote the PLATEAU (upper half of step sizes) and
    % show the trend separately rather than averaging everything together.
    st   = [o.step];
    good = o(st >= median(st));
    b = [good.blCnt];  bu = [good.blUm];  td = abs([o.tdCnt]);
    fprintf('%s: encoder %.5f um/count\n', nm{ax}, abs(encQ(ax)));
    fprintf(['   backlash (steps >= %g um): %+.2f +- %.2f counts = ' ...
             '%+.3f +- %.3f um   [n = %d]\n'], median(st), ...
            mean(b), std(b), mean(bu), std(bu), numel(good));
    fprintf('   trend over all steps:');
    for q = 1:numel(o), fprintf(' %g:%+.2f', o(q).step, o(q).blCnt); end
    fprintf('  (counts)\n');
    fprintf('   TD control: max |backlash| %.3f counts', max(td));
    if max(td) > 0.5
        fprintf('  *** NOT ~0 -- the command itself is direction-dependent,\n');
        fprintf('       so the TP figure above is NOT purely mechanical ***\n');
    else
        fprintf('  (~0, as it must be: the command is symmetric)\n');
    end
    if std(b) < 0.25*abs(mean(b))
        fprintf(['   step-size INDEPENDENT across the plateau -> consistent ' ...
                 'with true mechanical backlash\n']);
    else
        fprintf(['   varies with step size even on the plateau -> NOT a simple ' ...
                 'fixed backlash; read the per-step rows\n']);
    end
end
fprintf(['\nREMINDER: if hSI.hMotors.backlashCompensation is non-zero, ' ...
         'Motors.moveStartRelative\n(Motors.m:725-746) already overshoots on ' ...
         'reversal, and these numbers are the\nRESIDUAL after that ' ...
         'compensation, not the raw stage.\n']);
end

function d = incr_for(R, ax, s, dirSign)
% Consecutive TP increments within each direction run, in counts per commanded
% step. The first increment of every run is DROPPED: immediately after a
% reversal the stage is still taking up the backlash, so that increment measures
% backlash, not step size.
d = [];
m0 = [R.ax]==ax & [R.step]==s & [R.dir]==dirSign;
if ~any(m0), return; end
for rep = unique([R(m0).rep])
    m = m0 & [R.rep]==rep;
    ks = [R(m).k];  cs = [R(m).tpCnt];
    if numel(ks) < 3, continue; end
    % visit order: ascending k going up, descending k going down
    if dirSign > 0, [~,o] = sort(ks,'ascend'); else, [~,o] = sort(ks,'descend'); end
    di = diff(cs(o));
    if numel(di) >= 2, di = di(2:end); end
    d = [d, di]; %#ok<AGROW>
end
end

function sd = pooled_resid(R, up, dn)
% Repeatability: RMS residual about a straight line through (target, tpCnt),
% pooled over the two directions. Detrending is what separates "the stage does
% not return to the same place" from "the targets are far apart".
r = [];
for m = {up, dn}
    t = [R(m{1}).target].';  c = [R(m{1}).tpCnt].';
    ok = isfinite(t) & isfinite(c);
    if nnz(ok) < 3 || std(t(ok)) == 0, continue; end
    p = polyfit(t(ok), c(ok), 1);
    r = [r; c(ok) - polyval(p, t(ok))]; %#ok<AGROW>
end
if isempty(r), sd = NaN; else, sd = sqrt(mean(r.^2)); end
end
