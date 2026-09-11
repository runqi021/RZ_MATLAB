function cell_coherence_pooled_260727()
%% cell_coherence_pooled_260727  Breath x Ca coherence PER CELL, pooled across
%                                every recording in which that cell was matched.
% -----------------------------------------------------------------------
% STEP 3. Turns the per-observation container from cell_pool_260727.m into one
% coherence estimate per CELL, plus the per-recording estimates it was built
% from, so pooling can always be checked against its parts.
%
% ================= HOW THE POOLING WORKS (and why it is safe) =============
% The recordings of one cell differ in length, in frame rate, and in the animal's
% breathing rate. Three ways to pool them would each destroy something:
%   - concatenating the traces invents continuity across the seam between two
%     recordings taken minutes apart, and the seam lands in the breathing band;
%   - resampling to a common fps distorts spike timing, the very thing measured;
%   - forcing one shared frequency band ignores that breathing rate drifts
%     between recordings, so the band would sit off-peak for some of them.
% None of those happen here. Instead each recording is transformed on its OWN
% sampling grid, over ITS OWN breath band, and only the resulting cross-spectral
% quantities are added up:
%
%     A  = SUM over recordings, tapers, in-band frequencies of  Jx .* conj(Jy)
%     Bx = SUM of |Jx|^2          By = SUM of |Jy|^2
%     coherence = |A| / sqrt(Bx*By)        phase = angle(A)
%
% Jx, Jy are the multitaper transforms of the breath phase reference and of that
% cell's spike train in that recording. This is the standard multi-segment
% coherence estimator -- exactly what Chronux's own trialave does, generalized to
% segments of unequal length. Nothing is interpolated, resampled, or concatenated.
%
% WEIGHTING. Summing weights a recording by how much data it contributes. That is
% the default and is the statistically standard choice. The equal-per-recording
% alternative is ALSO computed and stored (r_pool_equal), so the two can be
% compared rather than argued about.
%
% SIGNIFICANCE. confC = sqrt(1 - alpha^(1/(dim-1))) with dim = K tapers x number
% of recordings -- so it is PER CELL, and a cell pooled from 4 recordings has a
% lower threshold than a cell seen once. That is the whole benefit of pooling and
% it means there is no single significance circle on the polar plot; each cell
% carries its own. Frequencies within the band add further degrees of freedom
% that this deliberately does NOT claim, so the threshold stays conservative.
% Cells seen in >=3 recordings additionally get a leave-one-recording-out
% jackknife CI on both magnitude and phase.
%
% COMPARABILITY. The foundation reports mean(coherence) over the band (a mean of
% ratios); this reports a ratio of sums. They are different estimators and do not
% have to agree to the third decimal. Both are computed per recording here
% (r_rec, r_rec_bandmean) so the pooled number can be traced back to the
% foundation's number for the same observation.
%
% Output (into cfg.cellDir):
%   cell_coherence_pooled.mat / .csv
%   cell_coherence_pooled.png / .pdf   polar + consistency panels
%
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
TW          = 4;        % multitaper TW  (match coherence_polar_general_260727.m)
alpha_sig   = 0.001;    % primary significance level
alpha_sig2  = 0.05;     % secondary level (outer dashed circle)
ca_lag_sec  = 0.1;      % GCaMP lead compensation: spikes shifted EARLIER this much
minSpikes   = 2;        % an observation contributes only with >= this many spikes
minObs      = 1;        % a cell is reported with >= this many contributing observations
doSave      = true;
% =====================================================================
K_tap = 2*TW - 1;

fprintf('\n=========== cell_coherence_pooled_260727 ===========\n');
assert(isfile(cfg.poolFile), ['cell_pool.mat not found:\n  %s\nRun cell_pool_260727.m first.'], cfg.poolFile);
P = load(cfg.poolFile, 'pool');  pool = P.pool;
rec = pool.rec;  obs = pool.obs;  cells = pool.cells;
fprintf('pool: %d cells, %d usable observations, %d usable recordings\n', ...
        numel(cells), pool.audit.n_obs_usable, pool.audit.n_rec_usable);

%% ---- per-cell pooled coherence ----
nCell = numel(cells);
Cres = struct('cell_id',{},'n_obs',{},'n_rec',{},'n_spikes',{}, ...
              'r_pool',{},'th_pool',{},'r_pool_equal',{},'confC',{},'confC2',{},'sig',{}, ...
              'r_jk_lo',{},'r_jk_hi',{},'th_jk_sd',{}, ...
              'r_rec',{},'th_rec',{},'r_rec_bandmean',{},'confC_rec',{},'n_sig_rec',{}, ...
              'th_circmean',{},'th_circstd',{},'th_R',{},'r_mean',{},'r_sd',{}, ...
              'f_pk',{},'band_lo',{},'band_hi',{},'rec_names',{},'roi_idx',{},'group',{});

fprintf('\ncomputing...\n');
for c = 1:nCell
    ii = cells{c};
    if isempty(ii), continue; end

    A = 0; Bx = 0; By = 0;                       % pooled accumulators
    rj = []; thj = []; rbm = []; Aj = []; Bxj = []; Byj = [];
    nsp = 0; recIdx = []; recNm = strings(0,1); roiIx = [];
    fpk = []; bandLo = []; bandHi = [];

    for q = 1:numel(ii)
        o = obs(ii(q));
        if ~o.usable, continue; end
        if o.n_spikes_used < minSpikes, continue; end
        k  = o.rec;  R = rec(k);
        Fs = R.fps;  T = R.T;

        x = R.ref(:);                                  % breath phase reference cos(phi)
        st = full(o.spikes(:));
        lag = round(ca_lag_sec*Fs);                    % GCaMP lead: shift spikes EARLIER
        if lag > 0 && lag < T
            st = [st(1+lag:end); zeros(lag,1)];
        end
        y = st - mean(st);
        if all(y == 0), continue; end

        [Jx, Jy] = seg_ffts(x, y, Fs, TW, R.band);
        if isempty(Jx), continue; end

        a_  = sum(sum(Jx .* conj(Jy)));                % sign convention matches the
        bx_ = sum(sum(abs(Jx).^2));                    %   foundation's th = -angle(C12)
        by_ = sum(sum(abs(Jy).^2));
        A = A + a_;  Bx = Bx + bx_;  By = By + by_;
        Aj(end+1,1) = a_;  Bxj(end+1,1) = bx_;  Byj(end+1,1) = by_; %#ok<AGROW>

        rj(end+1,1)  = abs(a_)/sqrt(bx_*by_); %#ok<AGROW>
        thj(end+1,1) = angle(a_); %#ok<AGROW>

        % the foundation's estimator on the same data, for traceability
        S12 = mean(conj(Jx).*Jy, 2);  S1 = mean(abs(Jx).^2, 2);  S2 = mean(abs(Jy).^2, 2);
        rbm(end+1,1) = mean(abs(S12)./sqrt(S1.*S2)); %#ok<AGROW>

        nsp = nsp + o.n_spikes_used;
        recIdx(end+1,1) = k; %#ok<AGROW>
        recNm(end+1,1)  = R.name; %#ok<AGROW>
        roiIx(end+1,1)  = o.roi_index; %#ok<AGROW>
        fpk(end+1,1)    = R.f_pk; %#ok<AGROW>
        bandLo(end+1,1) = R.band(1); bandHi(end+1,1) = R.band(2); %#ok<AGROW>
    end

    nSeg = numel(rj);
    if nSeg < minObs || nSeg == 0, continue; end

    r_pool  = abs(A)/sqrt(Bx*By);
    th_pool = angle(A);
    % equal-weight-per-recording alternative
    Aeq = sum(Aj ./ sqrt(Bxj.*Byj)) / nSeg;
    r_pool_equal = abs(Aeq);

    dim   = K_tap * nSeg;
    confC = sqrt(1 - alpha_sig ^(1/max(dim-1,1)));
    conf2 = sqrt(1 - alpha_sig2^(1/max(dim-1,1)));
    confC_rec = sqrt(1 - alpha_sig^(1/max(K_tap-1,1)));   % same threshold the foundation uses

    % leave-one-recording-out jackknife (needs >= 3 segments to mean anything)
    r_lo = NaN; r_hi = NaN; th_sd = NaN;
    if nSeg >= 3
        zAll = atanh(min(r_pool, 1-1e-12));
        zj = nan(nSeg,1);  tj = nan(nSeg,1);
        for j = 1:nSeg
            keep = true(nSeg,1); keep(j) = false;
            Ak = sum(Aj(keep)); Bxk = sum(Bxj(keep)); Byk = sum(Byj(keep));
            zj(j) = atanh(min(abs(Ak)/sqrt(Bxk*Byk), 1-1e-12));
            tj(j) = angle(Ak);
        end
        ps  = nSeg*zAll - (nSeg-1)*zj;
        se  = std(ps)/sqrt(nSeg);
        r_lo = tanh(mean(ps) - 1.96*se);   r_lo = max(r_lo, 0);
        r_hi = tanh(mean(ps) + 1.96*se);   r_hi = min(r_hi, 1);
        th_sd = circ_std_local(tj) * sqrt(nSeg-1);   % jackknife SE of the pooled phase
    end

    % consistency across the recordings of this cell
    [th_cm, th_cs, th_R] = circ_stats_local(thj);

    gname = cfg.genotype;
    if all(contains(recNm, 'IO', 'IgnoreCase', true)), gname = 'IO'; end

    Cres(end+1) = struct('cell_id',c, 'n_obs',nSeg, 'n_rec',numel(unique(recIdx)), ...
        'n_spikes',nsp, 'r_pool',r_pool, 'th_pool',th_pool, 'r_pool_equal',r_pool_equal, ...
        'confC',confC, 'confC2',conf2, 'sig',r_pool >= confC, ...
        'r_jk_lo',r_lo, 'r_jk_hi',r_hi, 'th_jk_sd',th_sd, ...
        'r_rec',{rj}, 'th_rec',{thj}, 'r_rec_bandmean',{rbm}, 'confC_rec',confC_rec, ...
        'n_sig_rec',nnz(rj >= confC_rec), ...
        'th_circmean',th_cm, 'th_circstd',th_cs, 'th_R',th_R, ...
        'r_mean',mean(rj), 'r_sd',std(rj), ...
        'f_pk',{fpk}, 'band_lo',{bandLo}, 'band_hi',{bandHi}, ...
        'rec_names',{recNm}, 'roi_idx',{roiIx}, 'group',string(gname)); %#ok<AGROW>
end

assert(~isempty(Cres), 'No cell produced a coherence estimate (check minSpikes / the pool audit).');
nObsV  = [Cres.n_obs]';
rPool  = [Cres.r_pool]';
thPool = [Cres.th_pool]';
sigV   = [Cres.sig]';
confV  = [Cres.confC]';

fprintf('\n---- results ----\n');
fprintf('  cells with an estimate : %d  (%d from >1 recording)\n', numel(Cres), nnz(nObsV>1));
fprintf('  significant (own confC): %d of %d  (%.1f%%)\n', nnz(sigV), numel(Cres), 100*mean(sigV));
fprintf('  confC ranges %.3f (%dx) .. %.3f (%dx) because it depends on cell size\n', ...
        min(confV), max(nObsV), max(confV), min(nObsV));
multi = nObsV > 1;
if any(multi)
    fprintf('  multi-recording cells: pooled r %.3f vs mean per-recording r %.3f\n', ...
            mean(rPool(multi)), mean([Cres(multi).r_mean]));
    fprintf('  phase agreement across recordings: median circular SD %.3f rad (%.1f deg)\n', ...
            median([Cres(multi).th_circstd]), rad2deg(median([Cres(multi).th_circstd])));
end

%% ---- figure ----
fig = figure('Color','w','Name','Pooled per-cell coherence', ...
             'Units','centimeters','Position',[2 2 32 13]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);

% --- panel 1: polar, one marker per cell, observations tethered to it ---
ax = polaraxes(fig,'Position',[0.03 0.10 0.30 0.78]);  hold(ax,'on');
thc = linspace(0,2*pi,360);
polarplot(ax, thc, repmat(median(confV),1,360), 'k--','LineWidth',1);
polarplot(ax, thc, repmat(min(confV),1,360), ':','Color',[.5 .5 .5],'LineWidth',0.8);
polarplot(ax, thc, repmat(max(confV),1,360), ':','Color',[.5 .5 .5],'LineWidth',0.8);
for k = 1:numel(Cres)
    col = cfg.genotype_color;  if Cres(k).group == "IO", col = [0 0 0]; end
    pale = col + (1-col)*0.65;                % lightened, not RGBA: polarplot's
    mid  = col + (1-col)*0.40;                % Color does not take a 4th alpha element
    if Cres(k).n_obs > 1                      % tether each observation to its cell
        for j = 1:Cres(k).n_obs
            polarplot(ax, [Cres(k).th_rec(j) Cres(k).th_pool], [Cres(k).r_rec(j) Cres(k).r_pool], ...
                      '-', 'Color',pale, 'LineWidth',0.4);
            polarplot(ax, Cres(k).th_rec(j), Cres(k).r_rec(j), '.', 'Color',mid, 'MarkerSize',4);
        end
    end
    ms = 4 + 2*min(Cres(k).n_obs, 5);
    if Cres(k).sig
        polarplot(ax, Cres(k).th_pool, Cres(k).r_pool, 'o', 'MarkerFaceColor',col, ...
                  'MarkerEdgeColor','k','MarkerSize',ms,'LineWidth',0.4);
    else
        polarplot(ax, Cres(k).th_pool, Cres(k).r_pool, 'o', 'MarkerFaceColor','none', ...
                  'MarkerEdgeColor',mid,'MarkerSize',ms,'LineWidth',0.5);
    end
end
ax.RLim = [0 1]; ax.ThetaZeroLocation = 'right'; ax.ThetaDir = 'counterclockwise';
ax.RAxisLocation = 180;
title(ax, {sprintf('%s: pooled per-cell coherence', cfg.genotype), ...
           'filled = significant vs its OWN confC; marker size = #recordings'}, 'FontSize',8);

% --- panel 2: pooled vs mean per-recording (does pooling change the answer?) ---
ax2 = axes(fig,'Position',[0.40 0.60 0.25 0.32]); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, [0 1],[0 1], 'k:');
rMeanV = [Cres.r_mean]';                      % column, to match rPool's orientation
scatter(ax2, rMeanV(~multi), rPool(~multi), 14, [.7 .7 .7], 'filled');
scatter(ax2, rMeanV(multi),  rPool(multi),  26, cfg.genotype_color, 'filled');
xlabel(ax2,'mean per-recording coherence'); ylabel(ax2,'pooled coherence');
title(ax2,'pooling vs its parts (grey = 1x cells, on the line by construction)','FontSize',7);
xlim(ax2,[0 1]); ylim(ax2,[0 1]);

% --- panel 3: phase agreement across recordings of the same cell ---
ax3 = axes(fig,'Position',[0.40 0.11 0.25 0.32]); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
thSdV = [Cres.th_circstd]';
if any(multi)
    scatter(ax3, nObsV(multi) + 0.12*(rand(nnz(multi),1)-0.5), thSdV(multi), ...
            26, cfg.genotype_color, 'filled', 'MarkerFaceAlpha',0.7);
end
xlabel(ax3,'# recordings the cell was matched in'); ylabel(ax3,'phase circular SD (rad)');
title(ax3,'phase agreement across recordings (low = stable locker)','FontSize',7);

% --- panel 4: coherence gain vs cell size, and the moving threshold ---
ax4 = axes(fig,'Position',[0.72 0.60 0.25 0.32]); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
scatter(ax4, nObsV + 0.12*(rand(numel(Cres),1)-0.5), rPool, 22, ...
        cfg.genotype_color, 'filled', 'MarkerFaceAlpha',0.6);
us = unique(nObsV);
plot(ax4, us, arrayfun(@(s) sqrt(1-alpha_sig^(1/max(K_tap*s-1,1))), us), 'k--','LineWidth',1);
xlabel(ax4,'# recordings pooled'); ylabel(ax4,'pooled coherence');
title(ax4,sprintf('dashed = confC(\\alpha=%.3f) for that many recordings', alpha_sig),'FontSize',7);

% --- panel 5: how many spikes went into each cell ---
ax5 = axes(fig,'Position',[0.72 0.11 0.25 0.32]); hold(ax5,'on'); box(ax5,'on'); grid(ax5,'on');
scatter(ax5, [Cres.n_spikes]', rPool, 22, cfg.genotype_color, 'filled', 'MarkerFaceAlpha',0.6);
set(ax5,'XScale','log');
xlabel(ax5,'# spikes pooled into the cell'); ylabel(ax5,'pooled coherence');
title(ax5,'estimate quality vs how much evidence it had','FontSize',7);

sgtitle(sprintf('%s  |  %d cells (%d multi-recording)  |  %d significant  |  TW=%d, %d tapers, %.0f ms GCaMP lead', ...
    cfg.genotype, numel(Cres), nnz(multi), nnz(sigV), TW, K_tap, ca_lag_sec*1000), 'FontSize',9);

%% ---- save ----
if doSave
    if ~isfolder(cfg.cellDir), mkdir(cfg.cellDir); end
    exportgraphics(fig, fullfile(cfg.cellDir,'cell_coherence_pooled.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(cfg.cellDir,'cell_coherence_pooled.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');

    T = table([Cres.cell_id]', [Cres.group]', nObsV, [Cres.n_rec]', [Cres.n_spikes]', ...
              rPool, thPool, rad2deg(thPool), [Cres.r_pool_equal]', confV, sigV, ...
              [Cres.r_jk_lo]', [Cres.r_jk_hi]', [Cres.th_jk_sd]', ...
              [Cres.r_mean]', [Cres.r_sd]', [Cres.th_circmean]', [Cres.th_circstd]', ...
              [Cres.th_R]', [Cres.confC_rec]', [Cres.n_sig_rec]', ...
              arrayfun(@(s) strjoin(cellstr(s.rec_names),'|'), Cres, 'UniformOutput',false)', ...
              arrayfun(@(s) strjoin(string(s.roi_idx),'|'), Cres, 'UniformOutput',false)', ...
        'VariableNames', {'cell_id','group','n_obs','n_rec','n_spikes','r_pooled','th_pooled_rad', ...
                          'th_pooled_deg','r_pooled_equalweight','confC','significant', ...
                          'r_jackknife_lo','r_jackknife_hi','th_jackknife_sd', ...
                          'r_perrec_mean','r_perrec_sd','th_perrec_circmean','th_perrec_circsd', ...
                          'th_perrec_R','confC_perrec','n_sig_perrec','recordings','roi_indices'});
    writetable(T, fullfile(cfg.cellDir,'cell_coherence_pooled.csv'));

    params = struct('TW',TW,'K_tap',K_tap,'alpha_sig',alpha_sig,'alpha_sig2',alpha_sig2, ...
                    'ca_lag_sec',ca_lag_sec,'minSpikes',minSpikes,'minObs',minObs); %#ok<NASGU>
    save(fullfile(cfg.cellDir,'cell_coherence_pooled.mat'), 'Cres','params','cfg','-v7.3');
    fprintf('\nSaved cell_coherence_pooled.mat/.csv/.png/.pdf to\n  %s\n', cfg.cellDir);
end
fprintf('Done.\n');
end

%% ========================= helpers =========================
function [Jx, Jy] = seg_ffts(x, y, Fs, TW, band)
% Multitaper transforms of one segment on its OWN grid, restricted to its OWN
% breath band. Same tapers/nfft/grid Chronux's coherencyc would use, so a 1x cell
% reproduces the foundation up to the estimator difference documented at the top.
Jx = []; Jy = [];
N = numel(x);
if N < 8 || any(~isfinite(band)) || band(2) <= band(1), return; end
tap  = dpsschk([TW, 2*TW-1], N, Fs);       % scaled by sqrt(Fs), as Chronux does
nfft = max(2^nextpow2(N), N);              % pad = 0
[f, findx] = getfgrid(Fs, nfft, band); %#ok<ASGLU>
if isempty(findx) || ~any(findx), return; end
Jx = mtfftc(x(:), tap, nfft, Fs);  Jx = Jx(findx,:);
Jy = mtfftc(y(:), tap, nfft, Fs);  Jy = Jy(findx,:);
end

function [m, s, R] = circ_stats_local(th)
% Circular mean, circular SD (Mardia) and resultant length.
z = mean(exp(1i*th(:)));
m = angle(z);
R = abs(z);
if R <= 0, s = Inf; else, s = sqrt(-2*log(R)); end
end

function s = circ_std_local(th)
[~, s] = circ_stats_local(th);
end
