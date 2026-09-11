% laser_power_vs_z_260530.m
% -----------------------------------------------------------------------
%  Plot ScanImage laser power (%) vs z-depth for two map folders.
%  Reads P0/LC and the z-list from the first raw TIFF found in each folder
%  (all TIFFs in a session share the same depth-power schedule).
%  Model:  P(z) = P0 * exp((z - z0) / LC)
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =====================
folders = { ...
    'C:\Users\Admin\Desktop\live_251104_definitely_not_sst_tdtomato_fitc\Map'; ...
    'C:\Users\Admin\Desktop\260114_homo_cy5\map'};

doSave  = true;
outDir  = 'C:\Users\Admin\Desktop';
outName = 'laser_power_vs_z';

% --- calibration: ScanImage % -> actual mW (measured at sample) ---
cal_pct = [0.10 1 3 5 8 10 12 15 18 20 22 25 27 30 32 35 37 40 ...
           45 50 55 60 65 70 75 80 85 90];
cal_mW  = [2.3  2.3 3.4 5.7 11.2 16.5 23 34.4 48.8 59.6 70.5 89.5 103 ...
           123 138 161 177 200 240 276 312 341 365 383 392 392 385 370];
% =========================================================

set(0,'DefaultAxesFontName','Arial');
fig = figure('Color','w','Units','centimeters','Position',[2 2 18 13]); hold on;
cols = lines(numel(folders));
S = cell(numel(folders), 1);

for kk = 1:numel(folders)
    folderPath = folders{kk};
    [~, lbl] = fileparts(folderPath);
    tiffList = dir(fullfile(folderPath,'**','*.tif'));
    if isempty(tiffList)
        warning('No TIFF found under %s', folderPath); continue;
    end
    tp = fullfile(tiffList(1).folder, tiffList(1).name);
    fprintf('[%d] %s\n     reading meta from: %s\n', kk, lbl, tiffList(1).name);

    [z_motor, p_pct, info] = read_si_power(tp);
    if isempty(z_motor)
        warning('Could not extract power schedule from %s', tp); continue;
    end

    % depth (um from surface): assume z_motor(1) is surface or top of stack
    depth_um = z_motor - z_motor(1);
    % clamp into calibration range so interp doesn't extrapolate weirdly
    p_clamp = min(max(p_pct, min(cal_pct)), max(cal_pct));
    p_mW    = interp1(cal_pct, cal_mW, p_clamp, 'pchip');

    plot(depth_um, p_mW, '-o', 'Color', cols(kk,:), ...
         'MarkerFaceColor', cols(kk,:), 'MarkerSize', 4, ...
         'LineWidth', 1.5, 'DisplayName', ...
         sprintf('%s  (P_0=%.1f%%, L_C=%.0f \\mum)', lbl, info.P0, info.LC));

    S{kk} = struct('folder',folderPath,'tiff',tp, ...
                   'z_motor',z_motor,'depth_um',depth_um, ...
                   'p_pct',p_pct,'p_mW',p_mW, ...
                   'P0',info.P0,'P1',info.P1,'LC',info.LC);
end

xlabel('Depth from top of stack (\mum)');
ylabel('Laser power at sample (mW)');
title('ScanImage power vs depth (%-to-mW calibrated)');
grid on; box on;
legend('Location','northwest','Interpreter','tex');

if doSave
    base = fullfile(outDir, outName);
    exportgraphics(fig, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'S');
    fprintf('Saved %s.png/.pdf/.mat\n', base);
end

%% ===================== LOCAL FUNCTIONS =====================
function [z_motor, p_pct, info] = read_si_power(tiffPath)
% Returns motor-Z list and power(%) at each Z from ScanImage TIFF metadata.
    z_motor = []; p_pct = []; info = struct('P0',NaN,'P1',NaN,'LC',NaN);
    inf1 = imfinfo(tiffPath);
    meta = "";
    if isfield(inf1(1),'Software') && ~isempty(inf1(1).Software)
        meta = string(inf1(1).Software);
    elseif isfield(inf1(1),'ImageDescription') && ~isempty(inf1(1).ImageDescription)
        meta = string(inf1(1).ImageDescription);
    else
        t = Tiff(tiffPath,'r'); c = onCleanup(@() t.close());
        try, meta = string(t.getTag('ImageDescription'));
        catch, try, meta = string(t.getTag('Software')); catch, return; end
        end
    end
    meta = replace(meta, char(8629), newline);   % stray CR codes

    LC = local_get_num(meta, "SI.hBeams.lengthConstants");
    P0 = local_get_num(meta, "SI.hBeams.stackStartPower");
    P1 = local_get_num(meta, "SI.hBeams.stackEndPower");
    zs = local_get_vec(meta, "SI.hStackManager.zs");
    z0 = local_get_num(meta, "SI.hStackManager.stackZStartPos");
    z1 = local_get_num(meta, "SI.hStackManager.stackZEndPos");
    dz = local_get_num(meta, "SI.hStackManager.stackZStepSize");

    if ~isempty(zs)
        z_motor = zs(:);
    elseif all(isfinite([z0 z1 dz])) && dz~=0
        z_motor = (z0:dz:z1).';
    else
        return;
    end
    if ~isfinite(P0) || ~isfinite(LC) || LC == 0, return; end
    p_pct = P0 .* exp((z_motor - z_motor(1)) ./ LC);
    info.P0 = P0; info.P1 = P1; info.LC = LC;
end

function v = local_get_num(meta, key)
    pat = key + " = ";
    i = strfind(meta, pat);
    if isempty(i), v = NaN; return; end
    s = extractAfter(meta, i(1) + strlength(pat) - 1);
    line = extractBefore(s, newline);
    if isempty(line), line = s; end
    tok = regexp(line, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match', 'once');
    if isempty(tok), v = NaN; else, v = str2double(tok); end
end

function vec = local_get_vec(meta, key)
    vec = [];
    pat = key + " = ";
    i = strfind(meta, pat);
    if isempty(i), return; end
    s = extractAfter(meta, i(1) + strlength(pat) - 1);
    line = extractBefore(s, newline);
    if isempty(line), line = s; end
    m = regexp(line, '\[(.*)\]', 'tokens', 'once');
    if isempty(m), return; end
    nums = regexp(m{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    if isempty(nums), return; end
    vec = str2double(nums(:));
end
