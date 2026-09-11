% Ventral_surface_coherence_polar_phase_260529.m
% -----------------------------------------------------------------------
%  Re-do of Ventral_surface_coherence_polar_260528.m using the new
%  piecewise-linear phase reference (peak = 0, insp start = pi) instead of
%  the binary peak-onset train.
%
%  Walks <rootPath>\{ChAT, Vglut2, Vgat, Sst}\**\ca_spike_data.mat.  For
%  each FOV with BOTH peak data AND insp-start data, builds the phase
%  reference signal cos(phi(t)) and computes Chronux coherence against
%  every ROI's spike train.  Polar style and group colors mirror the
%  earlier script.
%
%  Exclusions:
%    - Vglut2/1124/IO  -> skipped entirely (per user request)
%    - other IOs (ChAT IO etc.) are kept and colored black hollow as before
%
%  Outputs land in <rootPath>\coherence_polar_phase_260529\
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath = 'D:\Ventral_surface_summary';
outDir   = fullfile(rootPath, 'coherence_polar_phase_260529');

% group label -> color. IO is split out from ChAT: any folder under ChAT\
% whose cell-level dir name contains 'IO' counts as IO, the rest is true ChAT.
groups       = {'IO', 'ChAT', 'Vglut2', 'Vgat', 'Sst'};
group_colors = [0    0    0   ;     % IO      black
                0.85 0.10 0.10;     % ChAT    red
                0.10 0.65 0.20;     % Vglut2  green
                0.10 0.30 0.85;     % Vgat    blue
                0.55 0.20 0.75];    % Sst     purple

scan_dirs = {'ChAT', 'Vglut2', 'Vgat', 'Sst'};

% --- exclude list: any recording whose folder contains any of these
%     substrings is skipped (case-insensitive). Vglut2 IO is intentionally
%     left out of this polar.
exclude_substrings = { fullfile('Vglut2','1124','IO') };

nDrop           = 30;
fallback_fps    = 30;
minSpikes       = 2;
TW              = 5;
alpha_sig       = 0.001;

f_breath_search = [0.2 4];
fwhm_factor     = 0.6;
min_bw          = 0.05;
fmin            = 0.05;
fmax            = 15;

doSave          = true;
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ---- single collector across ALL groups / FOVs / ROIs ----
PP = init_coll();
labels = {};
confC  = NaN;
n_per_group = zeros(1, numel(groups));

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), warning('Folder missing: %s', gdir); continue; end

    allMat = dir(fullfile(gdir, '**', 'ca_spike_data.mat'));
    fprintf('\n=== [%s] %d recordings ===\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;
        recName    = folder_basename(folderPath);

        % --- exclusion check ---
        if any(cellfun(@(p) ~isempty(strfind(lower(folderPath), lower(p))), exclude_substrings)) %#ok<STREMP>
            fprintf('  EXCLUDED: %s\n', folderPath); continue;
        end

        % assign group: within ChAT, split into IO vs ChAT by cell-level dir
        gname = sname;
        if strcmpi(sname,'ChAT') && is_io_path(folderPath, gdir)
            gname = 'IO';
        end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), warning('No color for group "%s"', gname); continue; end

        try
            % --- find peak + insp-start data (try stemmed and bare names) ---
            peak_file = find_data_file(folderPath, 'breath_peak_data.mat');
            ip_file   = find_data_file(folderPath, 'breath_insp_start_data.mat');
            if isempty(peak_file)
                fprintf('  skip (no peak): %s\n', recName); continue;
            end
            if isempty(ip_file)
                fprintf('  skip (no insp_start): %s\n', recName); continue;
            end

            fps = detect_session_fps(folderPath, fallback_fps);
            BP  = load(peak_file);
            IP  = load(ip_file);
            CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));

            % --- breath waveform + events (toss leading frames) ---
            breath = double(BP.breath(:));
            breath = breath(nDrop+1:end);
            T_b    = numel(breath);

            peak_idx = double(BP.insp_onset_idx(:)) - nDrop;
            peak_idx(peak_idx < 1 | peak_idx > T_b) = [];
            foot_idx = double(IP.insp_start_idx(:)) - nDrop;
            foot_idx(foot_idx < 1 | foot_idx > T_b) = [];

            nCa = numel(CA.roi_spikes(1).spike_train);
            T   = min(T_b, nCa);
            breath = breath(1:T);
            peak_idx(peak_idx > T) = [];
            foot_idx(foot_idx > T) = [];

            if numel(peak_idx) < 2 || numel(foot_idx) < 1
                fprintf('  skip (too few events): %s\n', recName); continue;
            end

            % --- piecewise-linear phase reference (peak=0, foot=pi) ---
            phi = piecewise_phase(peak_idx, foot_idx, T);
            if all(isnan(phi))
                fprintf('  skip (phase build failed): %s\n', recName); continue;
            end
            ref = cos(phi); ref(isnan(ref)) = 0;
            ref_m = ref - mean(ref);

            % --- breath PSD -> coherence band ---
            pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
            pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
            bw_ws = detrend(breath); bw_ws = bw_ws - mean(bw_ws);
            [Sb,fb] = mtspectrumc(bw_ws, pB); Sb=Sb(:); fb=fb(:);
            m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
            if ~any(m), fprintf('  skip (no PSD in band): %s\n', recName); continue; end
            [~,rl]=max(Sb(m)); ip2=find(m,1)+rl-1; f_pk=fb(ip2);
            h=Sb(ip2)/2; lo=ip2; while lo>1&&Sb(lo)>h, lo=lo-1; end
            hi=ip2;             while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
            f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
            bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
            band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

            pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
            pc.fpass=band; pc.err=[2,alpha_sig];

            % --- coherence per ROI ---
            nROI = numel(CA.roi_spikes);
            nInc = 0;
            for rid = 1:nROI
                st = double(CA.roi_spikes(rid).spike_train(:));
                st = st(1:min(T,numel(st)));
                if numel(st) < T, st(end+1:T) = 0; end
                if sum(st) < minSpikes, continue; end
                nInc = nInc + 1;

                [PP, confC] = add_coh(PP, ref_m, st - mean(st), pc, band, gi, confC);
                labels{end+1} = sprintf('%s/%s#%d', gname, recName, rid); %#ok<SAGROW>
            end
            n_per_group(gi) = n_per_group(gi) + nInc;
            fprintf('  [%d] %-50s band [%.2f %.2f] Hz  P=%d F=%d  %d/%d ROI  -> %s\n', ...
                    kk, recName, band(1), band(2), numel(peak_idx), numel(foot_idx), ...
                    nInc, nROI, gname);
        catch ME
            warning('  ERROR %s: %s', recName, ME.message);
        end
    end
end

if isempty(PP.r), error('No ROIs collected.'); end

%% ============================ POLAR =================================
fig = figure('Color','w','Name','Ventral surface phase polar (peak=0, foot=pi)', ...
             'Units','centimeters','Position',[2 2 16 14]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
ax = polaraxes(fig,'Position',[0.08 0.08 0.78 0.80]); hold(ax,'on');

thc = linspace(0,2*pi,360);
polarplot(ax, thc, repmat(confC,1,360), 'k--','LineWidth',1);

for k = 1:numel(PP.r)
    ci = PP.colorIdx(k);
    if isnan(PP.th(k)) || isnan(PP.r(k)) || ci<1 || ci>size(group_colors,1), continue; end
    col = group_colors(ci,:);
    sig = PP.r(k) >= confC;
    if sig
        polarplot(ax, [PP.th(k) PP.th(k)], [PP.rlo(k) PP.rhi(k)], '-', 'Color',col, 'LineWidth',1.0);
        if ~isnan(PP.dphi(k))
            arc = linspace(PP.th(k)-PP.dphi(k), PP.th(k)+PP.dphi(k), 30);
            polarplot(ax, arc, PP.r(k)*ones(size(arc)), '-', 'Color',col, 'LineWidth',1.0);
        end
    end
    if ci == 1     % IO hollow black
        polarplot(ax, PP.th(k), PP.r(k), 'o', 'MarkerFaceColor','none', ...
                  'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.8);
    else
        polarplot(ax, PP.th(k), PP.r(k), 'o', 'MarkerFaceColor',col, ...
                  'MarkerEdgeColor','k', 'MarkerSize',6, 'LineWidth',0.4);
    end
end
ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
ax.RAxisLocation=180; ax.FontSize=8;

% legend via a hidden cartesian axes (polaraxes won't host plot() handles)
legAx = axes(fig, 'Position',[0.88 0.10 0.10 0.78], 'Visible','off'); hold(legAx,'on');
hL = []; lL = {};
for gi = 1:numel(groups)
    nAll = n_per_group(gi); if nAll==0, continue; end
    col = group_colors(gi,:);
    if gi == 1
        hL(end+1) = plot(legAx, NaN, NaN, 'o', 'MarkerFaceColor','none', ...
                         'MarkerEdgeColor','k', 'MarkerSize',7, 'LineWidth',0.8); %#ok<SAGROW>
    else
        hL(end+1) = plot(legAx, NaN, NaN, 'o', 'MarkerFaceColor',col, ...
                         'MarkerEdgeColor','k', 'MarkerSize',7); %#ok<SAGROW>
    end
    lL{end+1} = sprintf('%s (n=%d)', groups{gi}, nAll); %#ok<SAGROW>
end
if ~isempty(hL), legend(legAx, hL, lL, 'Location','east','FontSize',8); end

title(ax, sprintf('Ventral surface phase polar  (insp start=0, peak=\\pi)\nconfC=%.2f, \\alpha=%.3f, N=%d ROI', ...
                  confC, alpha_sig, numel(PP.r)), 'Interpreter','tex');

if doSave
    exportgraphics(fig, fullfile(outDir,'coherence_polar_phase.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'coherence_polar_phase.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'coherence_polar_phase_data.mat'), ...
         'PP','labels','confC','groups','group_colors','n_per_group', ...
         'TW','alpha_sig','f_breath_search','fwhm_factor','min_bw','minSpikes','nDrop', ...
         'exclude_substrings');
    fprintf('\nSaved phase polar + .mat to %s\n', outDir);
end
fprintf('Done. %d ROIs across %d groups.\n', numel(PP.r), sum(n_per_group>0));

%% ========================= LOCAL FUNCTIONS ==========================
function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[],'colorIdx',[]);
end

function [C, confC] = add_coh(C, x, y, pc, band, cidx, confC)
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f = f(:); mb = f>=band(1) & f<=band(2); if ~any(mb), mb = true(size(f)); end
    C.th(end+1,1)   = angle(mean(exp(1i*(-phi(mb)))));
    C.r(end+1,1)    = mean(Cxy(mb));
    C.rlo(end+1,1)  = max(0, mean(Cerr(1,mb)));
    C.rhi(end+1,1)  = min(1, mean(Cerr(2,mb)));
    C.dphi(end+1,1) = 1.96*mean(phistd(mb));
    C.colorIdx(end+1,1) = cidx;
end

function phi = piecewise_phase(peak_idx, foot_idx, T)
% Piecewise-linear phase: FEET at 0/2pi/4pi..., PEAKS at pi/3pi/5pi...
% (inspiration start = 0, breath peak = pi).  NaN outside the event range.
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];   % 1=peak, 0=foot
[events, ord] = sort(events);
types = types(ord);
keep = true(size(events));
for i = 2:numel(events)
    if types(i) == types(i-1), keep(i) = false; end
end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases  = nan(size(events));
phi_cur = types(1) * pi;      % type=1 (peak) -> pi; type=0 (foot) -> 0
for i = 1:numel(events)
    phases(i) = phi_cur; phi_cur = phi_cur + pi;
end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end

function fp = find_data_file(folderPath, basename)
% Locate either [stem]_basename or bare basename in folderPath.
% basename example: 'breath_peak_data.mat' or 'breath_insp_start_data.mat'
hits = dir(fullfile(folderPath, ['*' basename]));
if ~isempty(hits)
    fp = fullfile(hits(1).folder, hits(1).name);
else
    bare = fullfile(folderPath, basename);
    if isfile(bare), fp = bare; else, fp = ''; end
end
end

function tf = is_io_path(folderPath, groupRoot)
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end
