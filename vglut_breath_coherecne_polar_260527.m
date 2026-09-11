% vglut_breath_coherecne_polar_260527.m
% -----------------------------------------------------------------------
%  Vglut2 analog of chat_breath_coherence_polar_260526.m. Plotting and
%  analysis are identical; only the dataset roots, the highlighted cells,
%  and a few labels change.
%
%  Produces FOUR figures (Fig 5 IO/reference dropped):
%    Fig 1  coherence_peakXpeak_polar           : breath-peak x Ca-spike
%                                                 coherence polar plot
%                                                 (only the 5 designated ROIs)
%    Fig 2  vglut5_dff_breath (+_full)          : dF/F + breath trace overlay,
%                                                 5 designated Vglut2 cells
%                                                 (windowed AND full versions)
%    Fig 3  vglut5_spectra_coherence_waveform   : WAVEFORM power spectra (L)
%                                                 + waveform coherence (R), 5x2
%    Fig 4  vglut5_triggered_avg_heatmap        : breath-sorted dF/F heatmap
%                                                 + triggered average (+/- 2 s)
% -----------------------------------------------------------------------
%  Only the 5 ROIs listed in chat_list are kept. The recording-discovery
%  loop still scans both rootPaths, but background ROIs are skipped via
%  cidx == 0 continue.
%
%  Method: Chronux coherencyc, TW tapers, band-averaged r=mean(C),
%  phase th=angle(mean(exp(-i*phi))), dashed confC circle, theta-zero=top,
%  clockwise. Band auto-set per recording from breath waveform PSD peak.
%
%  Alignment: breath cam 2P-triggered (fps = imaging); toss first nDrop
%  breath frames, truncate to common length, drop tail.
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
% Search across BOTH dataset roots (recursive ** discovery within each).
rootPaths = { ...
    'D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing', ...
    'D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing' };

outDir   = 'D:\batch_dffQC_test_260325\vglut_breath_coh_260527';   % all outputs land here

nDrop           = 30;        % breath frames to toss (match calcium)
fallback_fps    = 30;
minSpikes       = 2;         % include ROI if it has >= this many spikes (>1)
TW              = 4;         % multitaper TW for POLAR plot / significance (Fig 1)
TW_spec         = 6;         % multitaper TW for Fig 3 spectra (PSD + coherence)
alpha_sig       = 0.01;      % coherence error / significance level

f_breath_search = [0.2 4];   % Hz, search band for breath PSD peak
fwhm_factor     = 0.6;       % coherence band = fwhm_factor x FWHM
min_bw          = 0.05;      % Hz, minimum coherence band width
fmin            = 0.05;      % Hz, PSD lower bound
fmax            = 15;        % Hz, PSD / spectrum upper bound

% true Vglut2 neurons: {recording-folder substring, ROI index}
chat_list = { ...
    'roi5_1400-1230-0_x4.4_15lp_6000f_00001',     5; ...
    'left_pFN_roi1_z0_12x_00001',                 3; ...
    'pFN_roi4_z0_512x512_5x_6000f_00001',         7; ...
    'pFN_roi1_z0_512x512_6x_2000f_00001',         9; ...
    'pFN_roi2_z5_512x512_6x_6000f_00001',         6 ; ...
    'pFN_roi3_z20_512x512_3x_2000f_00001',        24};
chat_colors = [0.85 0.10 0.10;    % red
               0.10 0.45 0.85;    % blue
               0.10 0.65 0.20;    % green
               0.85 0.55 0.10;    % orange
               0.55 0.20 0.75;    % purple
               0.55 0.85 0.75];   

% Fig 2 time window (s) per Vglut2 cell, same row order as chat_list
chat_xlim  = [ 20  80;    % 1  roi5_1400-1230-0
               7  67;    % 2  left_pFN_roi1_z0_12x
               80 140;    % 3  pFN_roi4_z0_512x512_5x
               1  61;    % 4  pFN_roi1_z0_512x512_6x_2000f
               1  61    % 5  pFN_roi2_z5_512x512_6x_6000f
               2  62];  % 6  pFN_roi3_z20_512x512_3x_2000f

doSave          = true;
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end

%% ---- collector (one row per ROI) ----
PP = init_coll();
labels = {};
confC  = NaN;

%% ---- discover recordings across all rootPaths ----
allMat = [];
for r = 1:numel(rootPaths)
    h = dir(fullfile(rootPaths{r}, '**', 'ca_spike_data.mat'));
    allMat = [allMat; h]; %#ok<AGROW>
end
fprintf('Found %d recording(s) across %d root(s).\n', numel(allMat), numel(rootPaths));

for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    recName = folder_basename(folderPath);   % robust to dots in folder names
    try
        bp = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        if isempty(bp), fprintf('skip (no breath): %s\n', recName); continue; end

        fps = detect_session_fps(folderPath, fallback_fps);
        CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));
        nROI = numel(CA.roi_spikes);
        nCa  = numel(CA.roi_spikes(1).spike_train);

        BP = load(fullfile(bp(1).folder, bp(1).name));
        nB = numel(BP.breath);
        if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
            evt = double(BP.insp_onsets_train(:) ~= 0);
        else
            evt = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); evt(oi(oi>=1 & oi<=nB)) = 1;
        end
        evt(1:min(nDrop,numel(evt))) = [];
        bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);

        T = min([numel(evt), numel(bw), nCa]);
        evt = evt(1:T); bw = bw(1:T);
        if sum(evt) < 3, continue; end

        %% breath PSD (waveform) -> band
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0; pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(bw, pB); Sb=Sb(:); fb=fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        [~,rl]=max(Sb(m)); ip=find(m,1)+rl-1; f_pk=fb(ip);
        h=Sb(ip)/2; lo=ip; while lo>1&&Sb(lo)>h, lo=lo-1; end
        hi=ip; while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

        pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0; pc.fpass=band; pc.err=[2,alpha_sig];
        br_p = evt - mean(evt);

        nInc = 0;
        for rid = 1:nROI
            st = double(CA.roi_spikes(rid).spike_train(:));
            st = st(1:min(T,numel(st))); if numel(st)<T, st(end+1:T)=0; end
            if sum(st) < minSpikes, continue; end       % include >= minSpikes
            nInc = nInc + 1;

            cidx = 0;                                    % Vglut2 identity
            for c = 1:size(chat_list,1)
                if contains(recName, chat_list{c,1}) && rid == chat_list{c,2}, cidx = c; break; end
            end
            if cidx == 0, continue; end                  % keep only the 5 designated ROIs

            [PP, confC] = add_coh(PP, br_p, st - mean(st), pc, band, cidx, confC);
            labels{end+1} = sprintf('%s#%d', recName, rid); %#ok<SAGROW>
        end
        fprintf('[%d] %-40s band [%.2f %.2f] Hz  %d/%d ROI included\n', ...
                kk, recName, band(1), band(2), nInc, nROI);
    catch ME
        warning('  ERROR %s: %s', recName, ME.message);
    end
end

if isempty(PP.r), error('No ROIs collected.'); end

%% ============================ POLAR =================================
fig = figure('Color','w','Name','Breath peak x Ca peak coherence', ...
             'Units','centimeters','Position',[2 2 13 12]);
set(fig,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
ax = polaraxes(fig,'Position',[0.08 0.06 0.80 0.80]);
plot_panel(ax, PP, chat_colors, confC, 'breath peaks  x  Ca spikes');
sgtitle(sprintf('Vglut2 (color, w/ CI)   |   confC=%.2f, \\alpha=%.2f   N=%d ROI', ...
        confC, alpha_sig, numel(PP.r)));

if doSave
    exportgraphics(fig, fullfile(outDir,'coherence_peakXpeak_polar.png'), 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'coherence_peakXpeak_polar.pdf'), 'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'coherence_peakXpeak.mat'), 'PP','labels','confC', ...
         'chat_list','chat_colors','TW','alpha_sig','f_breath_search','fwhm_factor','min_bw','minSpikes');
    fprintf('\nSaved coherence_peakXpeak_polar.png/.pdf + .mat to %s\n', outDir);
end
fprintf('Done. %d ROIs (%d Vglut2).\n', numel(PP.r), sum(PP.colorIdx>0));

%% =================================================================== %%
%%   PER-Vglut2-CELL FIGURES (the 5 designated colored neurons)        %%
%%   Fig 2: dF/F + breath trace overlay (windowed + full)              %%
%%   Fig 3: power spectra (left) + breath x Ca coherence (right)       %%
%%   Fig 4: breath-sorted dF/F heatmap (top) + triggered average (bot) %%
%% =================================================================== %%
nC = size(chat_list,1);
clear S;
for c = 1:nC
    sub = chat_list{c,1}; roi = chat_list{c,2};
    hit = [];
    for r = 1:numel(rootPaths)
        h = dir(fullfile(rootPaths{r},'**',[sub '*'],'ca_spike_data.mat'));
        if ~isempty(h), hit = h; break; end
    end
    if isempty(hit), error('Could not find folder for %s', sub); end
    fp = hit(1).folder; recName = folder_basename(fp);

    df = dir(fullfile(fp, '*_ch1_dFF.mat'));
    bp = dir(fullfile(fp, '*breath_peak_data.mat'));
    assert(~isempty(df) && ~isempty(bp), 'missing dFF/breath in %s', fp);

    fps = detect_session_fps(fp, fallback_fps);
    D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
    BPc = load(fullfile(bp(1).folder, bp(1).name));
    CA  = load(fullfile(fp, 'ca_spike_data.mat'));

    bw = detrend(double(BPc.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
    if isfield(BPc,'insp_onsets_train') && numel(BPc.insp_onsets_train)==numel(BPc.breath)
        ev = double(BPc.insp_onsets_train(:) ~= 0);
    else
        ev = zeros(numel(BPc.breath),1); oi = round(BPc.insp_onset_idx(:));
        ev(oi(oi>=1 & oi<=numel(ev))) = 1;
    end
    ev(1:min(nDrop,numel(ev))) = [];

    dff_all = double(D.dFF);
    st_all  = double(CA.roi_spikes(roi).spike_train(:));
    T = min([size(dff_all,1), numel(st_all), numel(bw), numel(ev)]);
    dff = dff_all(1:T, roi); st = st_all(1:T); bw = bw(1:T); ev = ev(1:T);
    t   = (0:T-1)'/fps;

    pB2.Fs=fps; pB2.tapers=[TW_spec,2*TW_spec-1]; pB2.pad=0; pB2.fpass=[fmin,min(fmax,fps/2)]; pB2.err=[2,alpha_sig];

    % --- WAVEFORM spectra: breath waveform + dF/F derivative (dF/F per s) ---
    [Sbw,fbw,SbwErr] = mtspectrumc(bw, pB2);            Sbw=Sbw(:); fbw=fbw(:);
    [Sdd,fdd,SddErr] = mtspectrumc(diff(dff)*fps, pB2); Sdd=Sdd(:); fdd=fdd(:);

    % peak / FWHM / detection band from breath WAVEFORM PSD
    mm = fbw>=f_breath_search(1) & fbw<=f_breath_search(2);
    [~,rl2]=max(Sbw(mm)); ip2=find(mm,1)+rl2-1; f_pk=fbw(ip2);
    hh=Sbw(ip2)/2; lo2=ip2; while lo2>1&&Sbw(lo2)>hh, lo2=lo2-1; end
    hi2=ip2;               while hi2<numel(fbw)&&Sbw(hi2)>hh, hi2=hi2+1; end
    f_fwhm=[max(fbw(lo2),f_breath_search(1)), min(fbw(hi2),f_breath_search(2))];
    bwd2=max(diff(f_fwhm)*fwhm_factor, min_bw);
    band2=[max(f_pk-bwd2/2,fmin), min(f_pk+bwd2/2,fmax)];

    % --- PEAK/SPIKE train auto-spectra ---
    [Sbe,fbe] = mtspectrumc(ev-mean(ev), pB2); Sbe=Sbe(:); fbe=fbe(:);
    [Sse,fse] = mtspectrumc(st-mean(st), pB2); Sse=Sse(:); fse=fse(:);

    % --- coherences (both): waveform x waveform AND peak x spike ---
    pc2.Fs=fps; pc2.tapers=[TW_spec,2*TW_spec-1]; pc2.pad=0; pc2.fpass=[fmin,min(fmax,fps/2)]; pc2.err=[2,alpha_sig];
    [~,Cw,~,~,~,~,fcw,confCw,~,Cerrw] = coherencyc(bw,         dff-mean(dff), pc2);
    [~,Cp,~,~,~,~,fcp,confCp,~,Cerrp] = coherencyc(ev-mean(ev), st-mean(st),  pc2);

    S(c).recName=recName; S(c).fps=fps; S(c).t=t; S(c).dff=dff; S(c).bw=bw; S(c).ev=ev;
    S(c).f_pk=f_pk; S(c).band=band2;
    S(c).fbw=fbw; S(c).Sbw=Sbw; S(c).SbwErr=SbwErr; S(c).fdd=fdd; S(c).Sdd=Sdd; S(c).SddErr=SddErr;
    S(c).fcw=fcw(:); S(c).Cw=Cw(:); S(c).confCw=confCw; S(c).Cerrw=Cerrw;
    S(c).fbe=fbe; S(c).Sbe=Sbe; S(c).fse=fse; S(c).Sse=Sse;
    S(c).fcp=fcp(:); S(c).Cp=Cp(:); S(c).confCp=confCp; S(c).Cerrp=Cerrp;
    fprintf('  Vglut2[%d] %-42s fps=%.2f peak=%.2f Hz band[%.2f %.2f]\n', ...
            c, recName, fps, f_pk, band2(1), band2(2));
end

%% ---- Fig 2: dF/F + breath trace overlay (5 rows), windowed + full ----
f2     = dff_breath_fig('vglut5 dFF vs breath (window)', S, chat_colors, chat_xlim, true);
f2full = dff_breath_fig('vglut5 dFF vs breath (full)',   S, chat_colors, chat_xlim, false);

%% ---- Fig 3: WAVEFORM spectra + waveform coherence (5x2) ----
f3 = spec_coh_fig('vglut5 waveform spectra + coherence', S, chat_colors, fmin, fmax, f_breath_search, ...
        'fbw','Sbw','SbwErr','fdd','Sdd','SddErr','fcw','Cw','confCw','Cerrw', ...
        'breath waveform + dF/F'' PSD', 'breath x Ca (waveform)');

%% ---- Fig 4: sorted heatmap (top) + triggered average (bottom) ----
f4 = figure('Color','w','Name','vglut5 triggered avg + heatmap', ...
            'Units','normalized','Position',[0.06 0.06 0.9 0.82]);
for c = 1:nC
    fps = S(c).fps; dff = S(c).dff; ev = S(c).ev;
    win = round(2 * fps);                  % fixed +/- 2 s window
    tau = (-win:win)/fps;
    on  = find(ev > 0); on = on(on-win >= 1 & on+win <= numel(dff));
    E   = zeros(numel(on), 2*win+1);
    for k = 1:numel(on), E(k,:) = dff(on(k)-win : on(k)+win); end

    post = mean(E(:, tau>=0), 2);
    [~, si] = sort(post, 'descend'); Es = E(si, :);
    cl = prctile(Es(:), [5 99.5]);

    axT = subplot(2,nC,c);
    imagesc(axT, tau, 1:size(Es,1), Es); axis(axT,'tight'); xlim(axT,[-2 2]);
    colormap(axT, flipud(gray(256))); caxis(axT, cl);
    hold(axT,'on'); xline(axT, 0, 'r-', 'LineWidth',1.2); hold(axT,'off');
    set(axT,'YDir','reverse'); ylabel(axT,'breath # (sorted)');
    title(axT, sprintf('%s  (n=%d)', S(c).recName, numel(on)), 'Interpreter','none');
    cb = colorbar(axT); cb.Label.String = 'dF/F';

    axB = subplot(2,nC,nC+c); hold(axB,'on');
    mu = mean(E,1); sd = std(E,0,1);
    fill(axB, [tau fliplr(tau)], [mu+sd fliplr(mu-sd)], 0.30*chat_colors(c,:)+0.70, ...
         'EdgeColor','none');
    plot(axB, tau, mu, 'Color',chat_colors(c,:), 'LineWidth',1.5);
    xline(axB, 0, 'k--', 'LineWidth',0.8);
    xlim(axB,[-2 2]); ylim(axB,[-0.2 0.6]); xlabel(axB,'time from inspiration (s)');
    ylabel(axB,'dF/F'); grid(axB,'on');
end

if doSave
    exportgraphics(f2, fullfile(outDir,'vglut5_dff_breath.png'),                 'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f2, fullfile(outDir,'vglut5_dff_breath.pdf'));
    exportgraphics(f2full, fullfile(outDir,'vglut5_dff_breath_full.png'),        'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f2full, fullfile(outDir,'vglut5_dff_breath_full.pdf'));
    exportgraphics(f3, fullfile(outDir,'vglut5_spectra_coherence_waveform.png'), 'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f3, fullfile(outDir,'vglut5_spectra_coherence_waveform.pdf'));
    exportgraphics(f4, fullfile(outDir,'vglut5_triggered_avg_heatmap.png'),      'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f4, fullfile(outDir,'vglut5_triggered_avg_heatmap.pdf'));
    fprintf('Saved vglut5_dff_breath / spectra_coherence / triggered_avg .png to %s\n', outDir);
end

%% ========================= LOCAL FUNCTIONS ==========================
function name = folder_basename(p)
% Robust folder basename: fileparts treats "...dir.x" as filename + ".x"
% extension, so directories with dots (e.g. "...x4.4_15lp_6000f_00001")
% get truncated. This rebuilds the full last path segment.
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end

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

function plot_panel(ax, C, chat_colors, confC, ttl)
    hold(ax,'on');
    thc = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(confC,1,360), 'k--','LineWidth',1);
    for k = 1:numel(C.r)        % background first (under colored)
        if C.colorIdx(k)>0 || isnan(C.th(k)) || isnan(C.r(k)), continue; end
        polarplot(ax, C.th(k), C.r(k), 'o', 'MarkerEdgeColor','k', ...
                  'MarkerFaceColor','none', 'MarkerSize',4, 'LineWidth',0.5);
    end
    for k = 1:numel(C.r)        % Vglut2 on top, with CI
        ci = C.colorIdx(k);
        if ci==0 || isnan(C.th(k)) || isnan(C.r(k)), continue; end
        col = chat_colors(ci,:);
        polarplot(ax, [C.th(k) C.th(k)], [C.rlo(k) C.rhi(k)], '-', 'Color',col, 'LineWidth',1.5);
        if ~isnan(C.dphi(k))
            arc = linspace(C.th(k)-C.dphi(k), C.th(k)+C.dphi(k), 30);
            polarplot(ax, arc, C.r(k)*ones(size(arc)), '-', 'Color',col, 'LineWidth',1.5);
        end
        polarplot(ax, C.th(k), C.r(k), 'o', 'MarkerFaceColor',col, ...
                  'MarkerEdgeColor','k', 'MarkerSize',9);
    end
    ax.RLim=[0 1]; ax.ThetaZeroLocation='top'; ax.ThetaDir='clockwise';
    ax.RAxisLocation=180; ax.FontSize=8;
    title(ax, ttl, 'Interpreter','none');
end

function save_vector_pdf(fig, pdfPath)
    set(fig,'Units','inches'); p = get(fig,'Position');
    set(fig,'PaperUnits','inches','PaperSize',p(3:4),'PaperPosition',[0 0 p(3:4)]);
    print(fig, char(pdfPath), '-dpdf', '-painters');
    [d,n] = fileparts(char(pdfPath));
    print(fig, fullfile(d,[n '.svg']), '-dsvg', '-painters');
end

function fig = dff_breath_fig(figName, S, chat_colors, chat_xlim, useWindow)
% dF/F (color, real scale) + breath (gray, normalized) overlay, one row/cell.
% useWindow=true -> plot only chat_xlim demo window; false -> full recording.
    norm01 = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);
    nC  = numel(S);
    fig = figure('Color','w','Name',figName,'Units','normalized','Position',[0.04 0.04 0.92 0.92]);
    for c = 1:nC
        ax = subplot(nC,1,c);
        if useWindow, w = chat_xlim(c,:);
        else,         w = [S(c).t(1) S(c).t(end)]; end
        mw = S(c).t >= w(1) & S(c).t <= w(2);
        tw = S(c).t(mw);
        yyaxis(ax,'right');
        plot(ax, tw, norm01(S(c).bw(mw)), '-', 'Color',[0.55 0.55 0.55], 'LineWidth',0.8);
        ylim(ax,[-0.05 1.05]); ylabel(ax,'breath (norm)');
        yyaxis(ax,'left');
        plot(ax, tw, S(c).dff(mw), '-', 'Color',chat_colors(c,:), 'LineWidth',0.8);
        ylabel(ax,'dF/F');
        ax.YAxis(1).Color = chat_colors(c,:); ax.YAxis(2).Color = [0.55 0.55 0.55];
        xlim(ax, w);
        title(ax, sprintf('%s   |   breath (gray, norm) vs dF/F (color, real)   peak %.2f Hz', ...
              S(c).recName, S(c).f_pk), 'Interpreter','none');
        if c==nC, xlabel(ax,'Time (s)'); end
    end
end

function add_band(ax, band, f_pk)
    yl = ylim(ax);
    p = patch(ax, [band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.90 0.90 0.90], 'EdgeColor','none');
    xline(ax, f_pk, 'k--', 'LineWidth',1);
    uistack(p,'bottom'); ylim(ax, yl);
end

function fig = spec_coh_fig(figName, S, cols, fmin, fmax, cohXlim, ...
        fBf, SBf, SBerrf, fSf, SSf, SSerrf, fcf, Cff, conff, Cerrf, leftTtl, rightTtl)
% nCx2 figure: left = power spectra (gray breath, color Ca) in dB log-x with
% shaded jackknife CI (Serr); right = coherence spectrum with shaded CI
% (Cerr) + confC line. Breath peak (dashed) + band (shaded) on both.
    nC = numel(S);
    fig = figure('Color','w','Name',figName,'Units','normalized','Position',[0.06 0.04 0.78 0.92]);
    for c = 1:nC
        col = cols(c,:); bd = S(c).band;
        axL = subplot(nC,2,2*c-1); hold(axL,'on');
        fB = S(c).(fBf); fB=fB(:)'; eB = S(c).(SBerrf);
        fill(axL, [fB fliplr(fB)], 10*log10([eB(1,:) fliplr(eB(2,:))]), [0.85 0.85 0.85], ...
             'EdgeColor','none');
        fS = S(c).(fSf); fS=fS(:)'; eS = S(c).(SSerrf);
        fill(axL, [fS fliplr(fS)], 10*log10([eS(1,:) fliplr(eS(2,:))]), 0.30*col+0.70, ...
             'EdgeColor','none');
        plot(axL, fB, 10*log10(S(c).(SBf)), 'Color',[0.6 0.6 0.6], 'LineWidth',0.9);
        plot(axL, fS, 10*log10(S(c).(SSf)), 'Color',col,           'LineWidth',1.1);
        set(axL,'XScale','log'); xlim(axL,[fmin fmax]);
        xticks(axL,[0.1 0.3 1 3 10]); set(axL,'XMinorTick','off');
        add_band(axL, bd, S(c).f_pk); ylabel(axL,'power (dB)');
        pbaspect(axL,[1 1 1]);
        title(axL, sprintf('%s  %s', S(c).recName, leftTtl), 'Interpreter','none');
        if c==nC, xlabel(axL,'Frequency (Hz)'); end

        axR = subplot(nC,2,2*c); hold(axR,'on');
        fc = S(c).(fcf); fc = fc(:)';
        ce = S(c).(Cerrf);
        fill(axR, [fc fliplr(fc)], [ce(1,:) fliplr(ce(2,:))], 0.30*col+0.70, ...
             'EdgeColor','none');
        plot(axR, fc, S(c).(Cff), 'Color',col, 'LineWidth',1.0);
        yline(axR, S(c).(conff), 'k--', 'LineWidth',1);
        set(axR,'XScale','log'); xlim(axR,cohXlim); ylim(axR,[0 1]);
        xticks(axR, [0.25 0.5 1 2 4]); set(axR,'XMinorTick','off');
        add_band(axR, bd, S(c).f_pk); ylabel(axR,'coherence');
        pbaspect(axR,[1 1 1]);
        title(axR, sprintf('%s  confC=%.2f (shaded = jackknife CI)', ...
              rightTtl, S(c).(conff)), 'Interpreter','none');
        if c==nC, xlabel(axR,'Frequency (Hz)'); end
    end
end
