% breath_amplitude_by_phase_260730.m
% -----------------------------------------------------------------------
%  What does the BREATH WAVEFORM look like in phase coordinates?
%
%  The phase used by Ventral_surface_polar_coh_vs_rayleigh_260729.m is defined
%  purely by LANDMARKS -- onset = 0, peak = pi, next onset = 2pi, linear in time
%  within each half. So the breath is guaranteed to peak "at pi" by construction.
%  What is NOT guaranteed, and is worth looking at, is the SHAPE: how sharp the
%  peak is, how symmetric the two halves are, and whether every recording agrees.
%  If a recording's waveform did not peak at pi here, its landmark detection is
%  wrong -- so this doubles as a QC on the peak/trough GUIs.
%
%  Breath trace = the PER-RECORDING SVD PC1 (BP.breath from breath_peak_pc1.mat),
%  Z-SCORED WITHIN EACH RECORDING. Deliberately not the acquisition-comparable
%  fixed-metric version: here we want the SHAPE of each recording's own cycle on
%  a common vertical scale, not amplitudes comparable between recordings.
%  CAVEAT: Vglut2/0728 sourced its trace from breath_fixedmetric.mat (50 s window)
%  rather than a per-video SVD, so those three are flagged in the printout.
%
%  Panels
%    1  mean z-scored breath vs phase, one line per recording + grand mean
%    2  the same in polar, to sit alongside the cell polar plot
%    3  OCCUPANCY: frames per phase bin -- the ~7x inspiration/expiration
%       imbalance that forces the occupancy correction in the Rayleigh
%    4  where the shuffle-confirmed cells sit, drawn on top of the mean waveform
%
%  Runqi Zhang / 2026-07-30
% -----------------------------------------------------------------------
clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(scriptDir);

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
outDir    = fullfile(rootPath, 'polar_coh_vs_rayleigh_260729');   % sits with the polar plot
cellData  = fullfile(outDir, 'polar_coh_vs_rayleigh_data.mat');   % for panel 4

scan_dirs = {'ChAT','Vglut2','Vgat','Sst','Sert'};
groups       = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};
group_colors = [0 0 0; 0.85 0.10 0.10; 0.10 0.65 0.20; ...
                0.10 0.30 0.85; 0.55 0.20 0.75; 0.90 0.45 0.10];

nDrop     = 30;          % breath frames tossed up front, to align with Ca
nPhaseBins= 72;          % 5 deg bins -- finer than the Rayleigh's 36, for shape
fallback_fps = 30;
doSave    = true;
excludeRecordings = coh_cfg_260727().excludeRecordings;
% =========================================================

edges = linspace(0, 2*pi, nPhaseBins+1);
ctrs  = edges(1:end-1) + diff(edges)/2;

% Field list and ORDER must match the struct() built in the loop exactly, or
% B(end+1) = struct(...) fails with "dissimilar structures".
B  = struct('rec',{},'gname',{},'date',{},'prof',{},'occ',{},'gi',{}, ...
            'fps',{},'cyc_s',{},'insp_frac',{},'flag',{},'nCyc',{},'fr_per_insp',{});

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), continue; end
    allMat = dir(fullfile(gdir,'**','ca_spike_data.mat'));
    for kk = 1:numel(allMat)
        fp = allMat(kk).folder;
        rec = folder_basename(fp);
        rel = regexprep(strrep(fp,gdir,''),'^[\\/]+','');
        pp  = regexp(rel,'[\\/]','split');
        if numel(pp) < 3, continue; end
        dateStr = pp{1};
        gname = sname;
        if ~isempty(regexpi(pp{2},'IO','once')), gname = 'IO'; end
        gi = find(strcmp(groups,gname),1);
        if isempty(gi) || any(strcmp(rec, excludeRecordings)), continue; end
        bpF = fullfile(fp,'breath_peak_pc1.mat');
        ipF = fullfile(fp,'breath_insp_start_pc1.mat');
        if ~isfile(bpF) || ~isfile(ipF), continue; end
        try
            fps = detect_session_fps(fp, fallback_fps);
            BP  = load(bpF);  IP = load(ipF);
            CA  = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
            nCa = numel(CA.roi_spikes(1).spike_train);

            bw = detrend(double(BP.breath(:)));
            bw(1:min(nDrop,numel(bw))) = [];
            T  = min(numel(bw), nCa);
            bw = bw(1:T);
            % Z-SCORE within this recording: shape on a common vertical scale.
            bz = (bw - mean(bw)) / max(std(bw), eps);

            peak_idx = round(BP.insp_onset_idx(:)) - nDrop;
            foot_idx = round(IP.insp_start_idx(:)) - nDrop;
            if strcmpi(sname,'Vglut2') && strcmp(dateStr,'1124')
                peak_idx = peak_idx + 1;  foot_idx = foot_idx + 1;
                bz = [bz(1); bz(1:end-1)];
            end
            peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
            foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
            if numel(peak_idx) < 2 || numel(foot_idx) < 2, continue; end

            phi  = piecewise_phase_local(peak_idx, foot_idx, T);
            phiW = mod(phi, 2*pi);
            v    = ~isnan(phiW);
            bidx = nan(T,1);  bidx(v) = discretize(phiW(v), edges);
            occ  = accumarray(bidx(v), 1, [nPhaseBins 1]);   % genuine frame density

            % ---- profile by PER-CYCLE INTERPOLATION, not by binning ----------
            % Binning frames into 5 deg bins aliases badly on the inspiratory
            % half: inspiration is ~21% of a ~1.5 s cycle, so ~9-10 frames span
            % 180 deg, i.e. ~18-20 deg between consecutive frames. Most 5 deg
            % bins are then EMPTY for any given recording, and which ones are
            % populated depends on that recording's frames-per-inspiration. The
            % grand mean at a sparse bin is therefore an average over a small,
            % biased subset of recordings -- that is what produced the spurious
            % "ditch" at 172 deg, computed from 5 of 36 recordings while its
            % neighbours had 15 and 36.
            %
            % Interpolating each CYCLE onto the common grid instead gives every
            % bin a value from every cycle: no empty bins, no subset bias, and no
            % dependence on how many frames happened to land in a cycle.
            fi = sort(foot_idx);
            prof_c = nan(nPhaseBins, numel(fi)-1);
            for q = 1:numel(fi)-1
                idx = fi(q):fi(q+1);
                % Use the UNWRAPPED phi minus this cycle's start, NOT phiW. phiW is
                % mod(phi,2pi), so the cycle's last frame wraps back to ~0 and the
                % first can land at ~2pi-eps from floating-point mod of 2*k*pi --
                % either breaks monotonicity and the cycle gets thrown away.
                % phi - phi(start) is exactly 0..2pi and monotone by construction.
                pq  = phi(idx) - phi(fi(q));
                bq  = bz(idx);
                good = ~isnan(pq);
                if nnz(good) < 4, continue; end
                pq = pq(good);  bq = bq(good);
                [pq, iu] = unique(pq, 'stable');  bq = bq(iu);
                if numel(pq) < 4 || any(diff(pq) <= 0), continue; end
                prof_c(:,q) = interp1(pq, bq, ctrs(:), 'linear', NaN);
            end
            prof = mean(prof_c, 2, 'omitnan');
            nCycUsed = nnz(any(~isnan(prof_c),1));

            flag = "";
            if isfield(BP,'src_file') && contains(string(BP.src_file),'fixedmetric')
                flag = "fixed-metric source, not per-video SVD";
            end
            if T < 0.5*nCa
                flag = flag + sprintf(" | breath covers only %.0f%% of the recording", 100*T/nCa);
            end

            B(end+1) = struct('rec',rec,'gname',gname,'date',dateStr, ...
                'prof',prof,'occ',occ,'gi',gi,'fps',fps, ...
                'cyc_s',median(diff(sort(foot_idx)))/fps, ...
                'insp_frac',sum(occ(1:nPhaseBins/2))/sum(occ), 'flag',flag, ...
                'nCyc',nCycUsed, ...
                'fr_per_insp', sum(occ(1:nPhaseBins/2))/max(numel(foot_idx)-1,1)); %#ok<SAGROW>
        catch ME
            warning('  %s: %s', rec, ME.message);
        end
    end
end
assert(~isempty(B), 'no recordings');

fprintf('\n%-40s %-7s %6s %7s %7s %6s %7s  %s\n', ...
        'recording','group','fps','cycle_s','insp_%','cycles','fr/insp','note');
for i=1:numel(B)
    fprintf('%-40s %-7s %6.1f %7.2f %6.0f%% %6d %7.1f  %s\n', B(i).rec(1:min(40,end)), ...
            B(i).gname, B(i).fps, B(i).cyc_s, 100*B(i).insp_frac, B(i).nCyc, ...
            B(i).fr_per_insp, B(i).flag);
end
fprintf(['\nframes per inspiration: median %.1f  -> ~%.0f deg between consecutive frames.\n' ...
         'With %d deg bins, binning would leave most inspiratory bins empty; the profile\n' ...
         'is built by per-cycle interpolation instead, so every bin has every cycle.\n'], ...
        median([B.fr_per_insp]), 180/median([B.fr_per_insp]), round(360/nPhaseBins));
P = cat(2, B.prof);                    % nBins x nRec
O = cat(2, B.occ);
gm = mean(P,2,'omitnan');  gs = std(P,0,2,'omitnan');
fprintf('\ninspiration occupies %.0f%% of the cycle on average (range %.0f-%.0f%%)\n', ...
        100*mean([B.insp_frac]), 100*min([B.insp_frac]), 100*max([B.insp_frac]));
[~,ip] = max(gm);  [~,it] = min(gm);
fprintf('grand-mean z-scored breath PEAKS at %.0f deg, TROUGHS at %.0f deg\n', ...
        rad2deg(ctrs(ip)), rad2deg(ctrs(it)));

%% ============================ FIGURE ============================
fig = figure('Color','w','Units','centimeters','Position',[1 1 30 20]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);

% --- 1  linear, per recording ---
ax1 = subplot(2,2,1); hold(ax1,'on');
for i=1:numel(B)
    plot(ax1, rad2deg(ctrs), B(i).prof, '-', 'Color',[group_colors(B(i).gi,:) 0.35], 'LineWidth',0.7);
end
fill(ax1, [rad2deg(ctrs) fliplr(rad2deg(ctrs))], [gm+gs; flipud(gm-gs)]', ...
     [0.2 0.2 0.2],'EdgeColor','none','FaceAlpha',0.18);
plot(ax1, rad2deg(ctrs), gm, 'k-','LineWidth',2);
xline(ax1,180,'k--'); yline(ax1,0,'k:');
xlim(ax1,[0 360]); xticks(ax1,0:90:360);
xlabel(ax1,'phase (deg)   0 = insp onset,  180 = insp peak');
ylabel(ax1,'breath (z, per recording)');
title(ax1, sprintf('breath waveform in phase coordinates  (%d recordings)', numel(B)));
grid(ax1,'on');

% --- 2  polar ---
ax2 = polaraxes(fig,'Position',[0.58 0.545 0.30 0.34]); hold(ax2,'on');
r0 = min(gm)-0.3;                       % shift so radius is positive
for i=1:numel(B)
    polarplot(ax2, [ctrs ctrs(1)], [B(i).prof; B(i).prof(1)]-r0, ...
              '-','Color',[group_colors(B(i).gi,:) 0.30],'LineWidth',0.7);
end
polarplot(ax2, [ctrs ctrs(1)], [gm; gm(1)]-r0, 'k-','LineWidth',2);
polarplot(ax2, [0 0], [0 max(gm)-r0], 'k--');
ax2.ThetaZeroLocation='right'; ax2.ThetaDir='counterclockwise'; ax2.RTickLabel={};
title(ax2,'same, polar   (radius offset for display)');

% --- 3  occupancy ---
ax3 = subplot(2,2,3); hold(ax3,'on');
Onorm = O ./ sum(O,1) * 100;
for i=1:numel(B)
    plot(ax3, rad2deg(ctrs), Onorm(:,i), '-','Color',[group_colors(B(i).gi,:) 0.35],'LineWidth',0.7);
end
plot(ax3, rad2deg(ctrs), mean(Onorm,2), 'k-','LineWidth',2);
yline(ax3, 100/nPhaseBins, 'r--','LineWidth',1);
text(ax3, 5, 100/nPhaseBins, ' uniform', 'Color','r','VerticalAlignment','bottom','FontSize',7);
xline(ax3,180,'k--'); xlim(ax3,[0 360]); xticks(ax3,0:90:360);
xlabel(ax3,'phase (deg)'); ylabel(ax3,'% of frames in bin');
title(ax3, sprintf('OCCUPANCY: inspiration is only %.0f%% of the cycle', 100*mean([B.insp_frac])));
grid(ax3,'on');

% --- 4  cells on top of the waveform ---
ax4 = subplot(2,2,4); hold(ax4,'on');
plot(ax4, rad2deg(ctrs), gm, 'k-','LineWidth',2);
xline(ax4,180,'k--'); yline(ax4,0,'k:');
xlabel(ax4,'phase (deg)'); ylabel(ax4,'breath (z)');
xlim(ax4,[0 360]); xticks(ax4,0:90:360); grid(ax4,'on');
if isfile(cellData)
    S = load(cellData);
    ok = (S.RY.p_sh < S.alpha_ray) & (S.RY.logZ >= S.logZ_crit);
    yyaxis(ax4,'right');
    ax4.YColor = [0.3 0.3 0.3];
    hb = [];  lb = {};
    for gi = [6 3 2 1]                        % Sert, Vglut2, ChAT, IO
        m = ok & S.RY.colorIdx==gi;
        if nnz(m) < 2, continue; end
        d = mod(rad2deg(S.RY.thMed(m)), 360);
        [n,e2] = histcounts(d, 0:20:360);
        hb(end+1) = stairs(ax4, e2(1:end-1), n, '-', ...
             'Color',group_colors(gi,:),'LineWidth',1.4); %#ok<SAGROW>
        lb{end+1} = sprintf('%s (n=%d)', groups{gi}, nnz(m)); %#ok<SAGROW>
    end
    ylabel(ax4,'shuffle-confirmed cells / 20\circ bin');
    if ~isempty(hb), legend(ax4, hb, lb, 'Location','northwest','FontSize',7,'Box','off'); end
    title(ax4,'where the cells fire, on the breath waveform');
    yyaxis(ax4,'left');
else
    title(ax4,'(run the polar script first for the cell overlay)');
end

sgtitle(sprintf(['breath amplitude across phase   |   per-recording SVD PC1, z-scored within recording' ...
                 '   |   %d bins of %.0f\\circ'], nPhaseBins, 360/nPhaseBins), 'FontSize',9);

if doSave
    if ~isfolder(outDir), mkdir(outDir); end
    exportgraphics(fig, fullfile(outDir,'breath_amplitude_by_phase.png'), ...
                   'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'breath_amplitude_by_phase.pdf'), ...
                   'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'breath_amplitude_by_phase.mat'), 'B','P','O','ctrs','gm','gs','nPhaseBins');
    fprintf('\nsaved breath_amplitude_by_phase.png/.pdf/.mat to %s\n', outDir);
end

%% ===================== local functions =====================
function name = folder_basename(p)
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);  name = [n e];
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
% Identical to the polar script's: feet at 0/2pi, peaks at pi, linear ramps in
% time between consecutive landmarks, strict alternation, NaN outside the range.
    phi = nan(T,1);
    events = [peak_idx(:); foot_idx(:)];
    types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
    [events, ord] = sort(events);  types = types(ord);
    keep = true(size(events));
    for i = 2:numel(events)
        if types(i) == types(i-1), keep(i) = false; end
    end
    events = events(keep); types = types(keep);
    if numel(events) < 2, return; end
    phases  = nan(size(events));
    phi_cur = types(1) * pi;
    for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
    for i = 1:numel(events)-1
        a = events(i); b = events(i+1);
        if a < 1 || b > T || b <= a, continue; end
        phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
    end
end
