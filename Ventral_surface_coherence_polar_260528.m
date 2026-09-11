% Ventral_surface_coherence_polar_260528.m
% -----------------------------------------------------------------------
%  All-ROI piecewise-phase x Ca-spike COHERENCE polar plot for the
%  Ventral_surface_summary dataset (ChAT + Vglut2 + Sst + Vgat).
%
%  PHASE REFERENCE (the input fed to Chronux):
%     phi(t) is piecewise-linear between detected breath events:
%        insp ONSET (foot, trough) -> phi = 0
%        insp PEAK                  -> phi = pi
%     reference signal = cos(phi(t))   (foot=+1, peak=-1)
%  Requires per-FOV  *breath_peak_data.mat  AND  *breath_insp_start_data.mat
%  FOVs missing the insp-start file are skipped.
%
%  EVERY ROI from EVERY FOV under <rootPath> on ONE polar axes. Dots are
%  colored by GROUP (cell type). For each FOV:
%     - detect breath-PSD peak -> FWHM-derived coherence band
%     - Chronux coherencyc(cos(phi), Ca-spike train) in-band
%     - r = mean(|C|), th = angle(mean(exp(-i*phi))), jackknife CI bar/arc
%  Style copied from chat_breath_coherence_polar_260526.m (Fig 1 polar):
%     theta-zero top, clockwise, R[0,1], filled colored markers + magnitude
%     CI radial bar + phase CI arc; dashed confC circle.
%
%  Layout assumed:
%     <rootPath>/<Group>/<Date>/<cell|IO>/<recording>/ca_spike_data.mat
%
%  Outputs (under outDir):
%     coherence_polar_all.png / .pdf      -- single combined polar
%     coherence_polar_data.mat            -- PP/labels/group colors
%
%  Dependencies: Chronux (coherencyc, mtspectrumc), detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath, 'coherence_polar_260528');

% group label -> color. IO is split out from ChAT: any folder under ChAT\
% whose cell-level dir name contains 'IO' counts as IO, the rest is true ChAT.
groups       = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst'};
group_colors = [0    0    0   ;     % IO      black
                0.85 0.10 0.10;     % ChAT    red
                0.10 0.65 0.20;     % Vglut2  green
                0.10 0.30 0.85;     % Vgat    blue
                0.55 0.20 0.75];    % Sst     purple

% which top-level rootPath folders to scan, and how each maps to a group
scan_dirs = {'ChAT', 'Vglut2', 'Vgat', 'Sst'};

nDrop           = 30;        % breath frames tossed up front (align to Ca)
fallback_fps    = 30;
minSpikes       = 2;         % include ROI if it has >= this many spikes
TW              = 4;         % multitaper TW for coherence
alpha_sig       = 0.001;     % primary significance level (jackknife err uses this)
alpha_sig2      = 0.05;      % secondary threshold (drawn as outer dashed circle)
ca_lag_sec      = 0.015;     % GCaMP rise compensation: shift sig ROIs earlier (s)

f_breath_search = [0.2 4];   % Hz, search band for breath PSD peak
fwhm_factor     = 0.6;       % coherence band = fwhm_factor x FWHM
min_bw          = 0.05;      % Hz, minimum coherence band width
fmin            = 0.05;      % Hz, PSD lower bound
fmax            = 15;        % Hz, PSD upper bound

doSave          = true;
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ---- single collector across ALL groups / FOVs / ROIs ----
PP = init_coll();
labels = {};
confC  = NaN;
% Chronux confC closed form: confC = sqrt(1 - alpha^(1/(K-1))), K = 2*TW-1
K_tap  = 2*TW - 1;
confC2 = sqrt(1 - alpha_sig2 ^ (1/(K_tap - 1)));
n_per_group = zeros(1, numel(groups));

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir)
        warning('Folder missing: %s', gdir); continue;
    end

    allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
    fprintf('\n=== [%s] %d recordings ===\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;
        recName    = folder_basename(folderPath);
        recDate    = date_from_path(folderPath, gdir);
        % assign group: any cell-level dir named "*IO*" -> IO (regardless of scan dir)
        gname = sname;
        if is_io_path(folderPath, gdir)
            gname = 'IO';
        end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), warning('No color for group "%s" -- skipping %s', gname, recName); continue; end
        try
            bp = dir(fullfile(folderPath, '*breath_peak_data.mat'));
            ip = dir(fullfile(folderPath, '*breath_insp_start_data.mat'));
            if isempty(bp), fprintf('  skip (no breath peaks): %s\n', recName); continue; end
            if isempty(ip), fprintf('  skip (no insp-start): %s\n', recName); continue; end

            fps = detect_session_fps(folderPath, fallback_fps);
            CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));
            nROI = numel(CA.roi_spikes);
            nCa  = numel(CA.roi_spikes(1).spike_train);

            BP = load(fullfile(bp(1).folder, bp(1).name));
            IP = load(fullfile(ip(1).folder, ip(1).name));
            nB = numel(BP.breath);
            bw = detrend(double(BP.breath(:)));
            bw(1:min(nDrop,numel(bw))) = [];
            bw = bw - mean(bw);

            % --- piecewise phase: foot=0, peak=pi  (post-toss frame index)
            peak_idx = round(BP.insp_onset_idx(:)) - nDrop;     % "insp_onset_idx" = breath PEAK in legacy naming
            foot_idx = round(IP.insp_start_idx(:)) - nDrop;     % insp-start (foot/trough)

            % --- per-session timing fix: Vglut2/1124 used rising-edge trigger
            %     instead of falling-edge -> breath samples sit 1 frame late.
            %     Shift breath events +1 frame to realign with Ca.
            if strcmpi(sname,'Vglut2') && strcmp(recDate,'1124')
                peak_idx = peak_idx + 1;
                foot_idx = foot_idx + 1;
                bw = [bw(1); bw(1:end-1)];   % delay breath waveform by 1 frame
            end

            T = min([numel(bw), nCa]);
            peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
            foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
            if numel(peak_idx) < 2 || numel(foot_idx) < 2
                fprintf('  skip (too few events): %s\n', recName); continue;
            end
            phi = piecewise_phase_local(peak_idx, foot_idx, T);
            ref = cos(phi); ref(isnan(ref)) = 0;
            ref = ref - mean(ref(:));
            bw  = bw(1:T);

            %% breath waveform PSD -> coherence band
            pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
            pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
            [Sb,fb] = mtspectrumc(bw, pB); Sb=Sb(:); fb=fb(:);
            m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
            [~,rl]=max(Sb(m)); ip=find(m,1)+rl-1; f_pk=fb(ip);
            h=Sb(ip)/2; lo=ip; while lo>1&&Sb(lo)>h, lo=lo-1; end
            hi=ip;            while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
            f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
            bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
            band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

            pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
            pc.fpass=band; pc.err=[2,alpha_sig];

            nInc = 0;
            for rid = 1:nROI
                st = double(CA.roi_spikes(rid).spike_train(:));
                st = st(1:min(T,numel(st)));
                if numel(st)<T, st(end+1:T)=0; end
                if sum(st) < minSpikes, continue; end
                nInc = nInc + 1;

                [PP, confC] = add_coh(PP, ref, st - mean(st), pc, band, gi, confC, f_pk);
                labels{end+1} = sprintf('%s/%s/%s/%d', gname, recDate, recName, rid); %#ok<SAGROW>
            end
            n_per_group(gi) = n_per_group(gi) + nInc;
            fprintf('  [%d] %-50s band [%.2f %.2f] Hz  %d/%d ROI\n', ...
                    kk, recName, band(1), band(2), nInc, nROI);
        catch ME
            warning('  ERROR %s: %s', recName, ME.message);
        end
    end
end

if isempty(PP.r), error('No ROIs collected.'); end

%% ====== GCaMP-rise compensation: shift sig ROIs earlier by ca_lag_sec ======
% Apply only to ROIs with r >= confC.  A time lead of dt at frequency f_pk
% maps to a phase advance of 2*pi*f_pk*dt.
sig_shift = PP.r >= confC;
PP.th_raw = PP.th;
PP.th(sig_shift) = wrapToPi(PP.th(sig_shift) - 2*pi*PP.f_pk(sig_shift)*ca_lag_sec);
fprintf('GCaMP comp: shifted %d significant ROIs by %.0f ms earlier.\n', ...
        sum(sig_shift), ca_lag_sec*1000);

%% ============================ POLAR =================================
fig = figure('Color','w','Name','Ventral surface coherence polar', ...
             'Units','centimeters','Position',[2 2 16 14]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
ax = polaraxes(fig,'Position',[0.08 0.08 0.78 0.80]);
plot_panel(ax, PP, group_colors, confC, confC2, ...
           'Ventral surface: cos(\phi)  x  Ca spikes   (foot=0, peak=\pi)');
add_group_legend(fig, ax, groups, group_colors, n_per_group);

sub = cellfun(@(g,n) sprintf('%s n=%d', g, n), groups, num2cell(n_per_group), ...
              'UniformOutput', false);
sgtitle(sprintf('%s   |   confC(\\alpha=%.3f)=%.2f (black), confC(\\alpha=%.3f)=%.2f (gray)   N=%d ROI total', ...
        strjoin(sub,'  '), alpha_sig, confC, alpha_sig2, confC2, numel(PP.r)));

%% ============ PRINT SIGNIFICANT ROIs (r >= confC) ====================
sig_mask = PP.r >= confC;
fprintf('\n==================== SIGNIFICANT ROIs (r >= %.3f) ====================\n', confC);
fprintf('  %d / %d ROIs significant\n', sum(sig_mask), numel(PP.r));
for gi = 1:numel(groups)
    idx = find(sig_mask & PP.colorIdx == gi);
    if isempty(idx), continue; end
    fprintf('\n-- %s (n=%d sig) --\n', groups{gi}, numel(idx));
    [~, ord] = sort(PP.r(idx), 'descend');
    for k = ord(:)'
        i = idx(k);
        fprintf('  r=%.3f  th=%+6.2f rad   %s\n', PP.r(i), PP.th(i), labels{i});
    end
end
fprintf('========================================================================\n\n');

% also drop a CSV so duplicate hunting is easy in Excel / pandas
if doSave
    sig_idx = find(sig_mask);
    fid = fopen(fullfile(outDir,'sig_rois.csv'),'w');
    fprintf(fid,'group,label,r,th_rad,th_deg,fov,roi\n');
    for ii = sig_idx(:)'
        lab = labels{ii};
        parts = regexp(lab,'#','split');
        fov_part = parts{1}; roi_part = parts{end};
        fprintf(fid,'%s,%s,%.6f,%.6f,%.3f,"%s",%s\n', ...
                groups{PP.colorIdx(ii)}, lab, PP.r(ii), PP.th(ii), ...
                rad2deg(PP.th(ii)), fov_part, roi_part);
    end
    fclose(fid);
    fprintf('Saved sig_rois.csv to %s\n', outDir);
end

if doSave
    exportgraphics(fig, fullfile(outDir,'coherence_polar_all.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'coherence_polar_all.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'coherence_polar_data.mat'), ...
         'PP','labels','confC','confC2','groups','group_colors','n_per_group', ...
         'TW','alpha_sig','alpha_sig2','f_breath_search','fwhm_factor','min_bw','minSpikes','nDrop');
    fprintf('\nSaved coherence_polar_all.png/.pdf + .mat to %s\n', outDir);
end
fprintf('Done. %d ROIs across %d groups.\n', numel(PP.r), sum(n_per_group>0));

%% ========================= LOCAL FUNCTIONS ==========================
function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[],'colorIdx',[],'f_pk',[]);
end

function [C, confC] = add_coh(C, x, y, pc, band, cidx, confC, f_pk)
% Chronux band-averaged magnitude + circular-mean phase + jackknife CI.
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f = f(:);
    mb = f>=band(1) & f<=band(2);
    if ~any(mb), mb = true(size(f)); end
    C.th(end+1,1)   = angle(mean(exp(1i*(-phi(mb)))));
    C.r(end+1,1)    = mean(Cxy(mb));
    C.rlo(end+1,1)  = max(0, mean(Cerr(1,mb)));
    C.rhi(end+1,1)  = min(1, mean(Cerr(2,mb)));
    C.dphi(end+1,1) = 1.96*mean(phistd(mb));
    C.colorIdx(end+1,1) = cidx;
    C.f_pk(end+1,1) = f_pk;
end

function plot_panel(ax, C, group_colors, confC, confC2, ttl)
% theta-zero top, clockwise, dashed confC circles (primary + secondary),
% marker per ROI colored by group with magnitude-CI bar + phase-CI arc.
    hold(ax,'on');
    thc = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(confC, 1,360), 'k--','LineWidth',1);
    polarplot(ax, thc, repmat(confC2,1,360), '--','Color',[0.5 0.5 0.5],'LineWidth',0.8);

    for k = 1:numel(C.r)
        ci = C.colorIdx(k);
        if isnan(C.th(k)) || isnan(C.r(k)) || ci<1 || ci>size(group_colors,1), continue; end
        col = group_colors(ci,:);
        sig = C.r(k) >= confC;                  % above confC dashed circle
        if sig
            polarplot(ax, [C.th(k) C.th(k)], [C.rlo(k) C.rhi(k)], ...
                      '-', 'Color',col, 'LineWidth',1.0);
            if ~isnan(C.dphi(k))
                arc = linspace(C.th(k)-C.dphi(k), C.th(k)+C.dphi(k), 30);
                polarplot(ax, arc, C.r(k)*ones(size(arc)), ...
                          '-', 'Color',col, 'LineWidth',1.0);
            end
        end
        if ci == 1     % IO: hollow black circle
            polarplot(ax, C.th(k), C.r(k), 'o', ...
                      'MarkerFaceColor','none', 'MarkerEdgeColor','k', ...
                      'MarkerSize',6, 'LineWidth',0.8);
        else
            polarplot(ax, C.th(k), C.r(k), 'o', ...
                      'MarkerFaceColor',col, 'MarkerEdgeColor','k', ...
                      'MarkerSize',6, 'LineWidth',0.4);
        end
    end
    ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
    title(ax, ttl, 'Interpreter','none');
end

function add_group_legend(fig, refAx, groups, group_colors, n_per_group)
% Color-key legend showing group -> hue + ROI count.
    pos = refAx.Position;
    legAx = axes(fig, 'Position', [pos(1)+pos(3)+0.01, pos(2)+0.30*pos(4), 0.12, 0.40*pos(4)]);
    hold(legAx,'on'); axis(legAx,'off');
    N = numel(groups);
    y = linspace(0.90, 0.10, max(N,2));
    for k = 1:N
        if k == 1     % IO swatch: hollow black, matches plot
            plot(legAx, 0.08, y(k), 'o', 'MarkerFaceColor','none', ...
                 'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',0.8);
        else
            plot(legAx, 0.08, y(k), 'o', 'MarkerFaceColor',group_colors(k,:), ...
                 'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',0.4);
        end
        text(legAx, 0.25, y(k), sprintf('%s (n=%d)', groups{k}, n_per_group(k)), ...
             'FontSize',8, 'Interpreter','none');
    end
    xlim(legAx,[0 1]); ylim(legAx,[0 1]);
end

function d = date_from_path(folderPath, groupRoot)
% Path component immediately under groupRoot (e.g. '0521' under .../ChAT/).
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    if isempty(parts), d=''; else, d = parts{1}; end
end

function tf = is_io_path(folderPath, groupRoot)
% True if the cell-level dir (the one right after the date dir under
% groupRoot) contains 'IO' (case-insensitive). Layout:
%   <groupRoot>\<date>\<cell-or-IO>\<recording>\
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; this rebuilds the full
% last path segment for folders with dots in the name.
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Piecewise-linear phase reference: FEET at 0/2pi/..., PEAKS at pi/3pi/...
% Linear ramps in time between consecutive events. NaN outside the range.
    phi = nan(T,1);
    events = [peak_idx(:); foot_idx(:)];
    types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
    [events, ord] = sort(events);
    types = types(ord);
    % strict alternation: drop adjacent same-type duplicates
    keep = true(size(events));
    for i = 2:numel(events)
        if types(i) == types(i-1), keep(i) = false; end
    end
    events = events(keep); types = types(keep);
    if numel(events) < 2, return; end
    phases  = nan(size(events));
    phi_cur = types(1) * pi;        % type=1 (peak)->pi; type=0 (foot)->0
    for i = 1:numel(events)
        phases(i) = phi_cur; phi_cur = phi_cur + pi;
    end
    for i = 1:numel(events)-1
        a = events(i); b = events(i+1);
        if a < 1 || b > T || b <= a, continue; end
        phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
    end
end
