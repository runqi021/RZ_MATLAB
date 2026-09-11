function [mW, cal] = laser_power_calibration(pct, acqDate)
% LASER_POWER_CALIBRATION  ScanImage Pockels setpoint (%) -> mW at the sample.
%
%   mW        = laser_power_calibration(pct, acqDate)
%   [mW, cal] = laser_power_calibration(pct, acqDate)
%
% pct     : setpoint(s) in percent
% acqDate : datetime (or datenum, or 'yyyy-mm-dd') the data was ACQUIRED.
%           Required -- there is no safe default, see below.
%           May instead be a table NAME ('pre_260723' / 'post_260723') to
%           force one explicitly, which is honest about overriding the date
%           rather than passing a fake one.
%
% cal     : struct with the table actually used (.pct .mW .name .validFrom
%           .validTo) plus .clamped (true where pct fell outside the table).
%
% WHY THIS IS A SEPARATE FUNCTION
%   The rig was re-measured on 2026-07-23 and the curve MOVED, so "% -> mW" is
%   not a property of ScanImage, it is a property of (rig, date).  Picking the
%   wrong table is silent and quadratic: two-photon signal goes as P^2, so a
%   20%% power error is a 44%% intensity error, and the two tables differ by
%   nearly 4x at 1%%.  Every script that turns % into mW must call this and
%   pass the acquisition date, never paste a table inline.
%
%   ScanImage % is a Pockels setpoint and the curve is strongly superlinear at
%   the low end, so linear interpolation is wrong there; pchip is used.
%
% THE TWO TABLES
%   pre_260723  0.10-90 %, peaks at 392 mW near 75-80 %% then ROLLS OVER
%               (385 at 85 %, 370 at 90 %) -- above ~70 %% more setpoint buys
%               nothing, and the 2.3 mW reading at both 0.10 %% and 1 %% is a
%               floor (meter noise / Pockels leakage), not a measurement.
%   post_260723 0-85 %%, re-measured 2026-08-05.  Starts at a TRUE zero, so it
%               is usable at the low end where the old table is not.  Peaks at
%               371 mW near 81 %% and rolls over by 85 %%.
%
%   Requesting a pct outside a table's range CLAMPS to the endpoint and warns.
%   For post_260723 that means nothing above 50 %% can be quoted in mW at all.
%
% EXAMPLE
%   mW = laser_power_calibration([11 21 35 48], datetime(2026,7,28))
%
% See also LASER_POWER_VS_Z_260530 (which hard-codes the pre-260723 table).

narginchk(2, 2);

CAL(1).name      = 'pre_260723';
CAL(1).validFrom = datetime(1900,1,1);
CAL(1).validTo   = datetime(2026,7,23);        % exclusive
CAL(1).pct = [0.10 1 3 5 8 10 12 15 18 20 22 25 27 30 32 35 37 40 ...
              45 50 55 60 65 70 75 80 85 90];
CAL(1).mW  = [2.3  2.3 3.4 5.7 11.2 16.5 23 34.4 48.8 59.6 70.5 89.5 103 ...
              123 138 161 177 200 240 276 312 341 365 383 392 392 385 370];

% Re-measured 2026-08-05 across the full 0-85 % range.  This SUPERSEDES the
% earlier 1-50 % post-260723 table.  Unlike the pre_260723 table it starts at a
% true zero (0 mW at 0-1 %) rather than a 2.3 mW meter floor, and it covers the
% rollover: power peaks at 371 mW near 81 % and falls again by 85 %.
CAL(2).name      = 'post_260723';
CAL(2).validFrom = datetime(2026,7,23);
CAL(2).validTo   = datetime(2100,1,1);
CAL(2).pct = [0 0.5 1 2 3 4 5 6 8 10 12 15 18 20 22 25 27 30 32 35 37 40 ...
              42 45 47 50 55 60 65 70 75 80 81 85];
CAL(2).mW  = [0 0 0 0.1 0.5 1.2 2.2 3.5 6.8 11.4 17 27.3 39.5 48.9 59.5 ...
              76.5 88.5 108 122 143 158 180 195 216 231 253 286 315 337 ...
              354 365 370 371 368];

% --- explicit table name overrides the date --------------------------------
if (ischar(acqDate) || isstring(acqDate)) && any(strcmpi(acqDate, {CAL.name}))
    cal = CAL(strcmpi(acqDate, {CAL.name}));
    [mW, cal] = interp_clamped(pct, cal);
    return
end

% --- resolve the date ------------------------------------------------------
if ischar(acqDate) || isstring(acqDate)
    acqDate = datetime(acqDate);
elseif isnumeric(acqDate)
    acqDate = datetime(acqDate, 'ConvertFrom', 'datenum');
end
assert(isa(acqDate, 'datetime') && ~isnat(acqDate), ...
    'acqDate must be a datetime, datenum, or date string');

k = find(acqDate >= [CAL.validFrom] & acqDate < [CAL.validTo], 1);
assert(~isempty(k), 'No power calibration covers %s', datestr(acqDate));
cal = CAL(k);

[mW, cal] = interp_clamped(pct, cal);
end

function [mW, cal] = interp_clamped(pct, cal)
% pchip, not linear or spline: the curve is steeply convex at the low end so a
% linear chord overestimates (by up to 12.5 % near 1.7 %), and spline can
% overshoot where the point spacing changes.  pchip is shape-preserving and
% cannot make power fall as setpoint rises.  Above ~10 % the choice is worth
% <0.5 %; below 5 % the table spacing is the real limit either way.
lo = min(cal.pct);  hi = max(cal.pct);
cal.clamped = pct < lo | pct > hi;
if any(cal.clamped(:))
    warning('laser_power_calibration:clamped', ...
        ['%d of %d setpoint(s) fall outside the %s table (%g-%g %%) and were ' ...
         'CLAMPED to the endpoint. mW is not measured there -- treat those ' ...
         'points as unquantified, do not fit through them.'], ...
        nnz(cal.clamped), numel(pct), cal.name, lo, hi);
end
mW = interp1(cal.pct, cal.mW, min(max(pct, lo), hi), 'pchip');
end
