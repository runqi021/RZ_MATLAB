% chat_breath_coherence_polar_260526.m
% -----------------------------------------------------------------------
%  Produces FIVE figures for the ChAT_analysis dataset (Fig 5 = IO reference):
%    Fig 1  coherence_peakXpeak_polar           : pooled breath-peak x Ca-spike
%                                                 coherence polar plot (all ROIs)
%    Fig 2  chat3_dff_breath                    : dF/F + breath trace overlay,
%                                                 3 designated ChAT cells
%    Fig 3  chat3_spectra_coherence_waveform    : WAVEFORM power spectra (L,
%                                                 breath waveform + dF/F deriv)
%                                                 + waveform coherence (R), 3x2,
%                                                 breath peak (dashed) + band
%                                                 (shaded); log-x, dB power
%    Fig 4  chat3_triggered_avg_heatmap         : breath-sorted dF/F heatmap (top,
%                                                 inverted grayscale) + triggered
%                                                 average +/- SEM (bottom)
%    Fig 5  IO_dff_breath                       : IO reference recording (io_path):
%                                                 breath (top) + all-ROI stacked
%                                                 dF/F (bottom), single dataset
% -----------------------------------------------------------------------
%  Fig 1 detail: Breath PEAK train x Ca SPIKE train COHERENCE polar plot.
%
%  Pooled across all recordings. The 3 designated ChAT neurons are colored
%  and drawn with magnitude + phase 95% CI; all other ROIs (IO/non-ChAT) are
%  plain black circles with NO error bars.
%
%  Inclusion: every ROI with >= minSpikes detected calcium events (default 2,
%  i.e. more than 1 spike), NOT just the ifSpike-flagged cells.
%
%  Method = Chronux coherencyc, TW tapers, band-averaged r=mean(C),
%  phase th=angle(mean(exp(-i*phi))), dashed confC circle, theta-zero=top,
%  clockwise. Band auto-set per recording from breath waveform PSD peak.
%
%  ChAT_analysis specifics: breath camera is 2P-triggered (fps = imaging);
%  toss first nDrop breath frames, truncate to common length, drop tail.
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
rootPath        = 'C:\Users\Admin\Desktop\ChAT_analysis';

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

% true ChAT neurons: {recording-folder substring, ROI index}
chat_list = { 'roi5_7x_x-1200y200z-30_3000f_23lp', 1; ...
              'roi3_8x_x-1070y730z0_3000f_15lp',   1; ...
              'roi1_4x_x-900y700z-15_6000f_13lp',  1 };
chat_colors = [0.85 0.10 0.10;    % red
               0.10 0.45 0.85;    % blue
               0.10 0.65 0.20];   % green

% Fig 2 time window (s) per ChAT cell, same row order as chat_list
chat_xlim  = [ 35  95;    % roi5_7x (red)
               30  90;    % roi3_8x (blue)
              110 170];   % roi1_4x (green)

% IO reference recording for the breath + stacked-dF/F QC panel (Fig 5)
io_path = 'C:\Users\Admin\Desktop\ChAT_analysis\0124\IO\roi3_R_-1000_140_2x_34lp_512x256_00001';
io_xlim = [100 130];   % s, time window shown in Fig 5
% IO ROIs to plot in Fig 5 (QC keep-list); [] = all ROIs
io_sel  = [7, 8, 9, 12, 13, 14, 15, 17, 18, 19, 23, 22, 25, 26, 30, 29, 28, 31, 32, 34];
io_tbar_s = 5;         % Fig 5 time scale bar (s)
io_dffbar = 1;         % Fig 5 dF/F scale bar (dF/F units)

doSave          = true;
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% ---- collector (one row per ROI) ----
PP = init_coll();
labels = {};
confC  = NaN;

%% ---- discover recordings ----
allMat = dir(fullfile(rootPath, '**', 'ca_spike_data.mat'));
fprintf('Found %d recording(s).\n', numel(allMat));

for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    [~, recName] = fileparts(folderPath);
    try
        bp = dir(fullfile(folderPath, '*DLC*breath_peak_data.mat'));
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

            cidx = 0;                                    % ChAT identity
            for c = 1:size(chat_list,1)
                if contains(recName, chat_list{c,1}) && rid == chat_list{c,2}, cidx = c; break; end
            end

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
sgtitle(sprintf('ChAT (color, w/ CI) vs IO (black)   |   confC=%.2f, \\alpha=%.2f   N=%d ROI', ...
        confC, alpha_sig, numel(PP.r)));

if doSave
    exportgraphics(fig, fullfile(rootPath,'coherence_peakXpeak_polar.png'), 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(rootPath,'coherence_peakXpeak_polar.pdf'), 'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(rootPath,'coherence_peakXpeak.mat'), 'PP','labels','confC', ...
         'chat_list','chat_colors','TW','alpha_sig','f_breath_search','fwhm_factor','min_bw','minSpikes');
    fprintf('\nSaved coherence_peakXpeak_polar.png/.pdf + .mat to %s\n', rootPath);
end
fprintf('Done. %d ROIs (%d ChAT).\n', numel(PP.r), sum(PP.colorIdx>0));

%% =================================================================== %%
%%   PER-ChAT-CELL FIGURES (the 3 designated red/blue/green neurons)    %%
%%   Fig 2: dF/F + breath trace overlay                                 %%
%%   Fig 3: power spectra (left) + breath x Ca coherence (right)        %%
%%   Fig 4: breath-sorted dF/F heatmap (top) + triggered average (bot)  %%
%% =================================================================== %%
norm01 = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);

nC = size(chat_list,1);
clear S;
for c = 1:nC
    sub = chat_list{c,1}; roi = chat_list{c,2};
    hit = dir(fullfile(rootPath,'**',[sub '*'],'ca_spike_data.mat'));
    if isempty(hit), error('Could not find folder for %s', sub); end
    fp = hit(1).folder; [~, recName] = fileparts(fp);

    df = dir(fullfile(fp, '*_ch1_dFF.mat'));
    bp = dir(fullfile(fp, '*DLC*breath_peak_data.mat'));
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
    % 3rd output Serr = [2 x nf] jackknife lower/upper bounds (linear power).
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

    % --- PEAK/SPIKE train auto-spectra: breath onset train + Ca spike train ---
    [Sbe,fbe] = mtspectrumc(ev-mean(ev), pB2); Sbe=Sbe(:); fbe=fbe(:);
    [Sse,fse] = mtspectrumc(st-mean(st), pB2); Sse=Sse(:); fse=fse(:);

    % --- coherences (both): waveform x waveform AND peak x spike ---
    pc2.Fs=fps; pc2.tapers=[TW_spec,2*TW_spec-1]; pc2.pad=0; pc2.fpass=[fmin,min(fmax,fps/2)]; pc2.err=[2,alpha_sig];
    % coherencyc -> [C12(complex), C(magnitude), phi, S12,S1,S2, f, confC, phistd, Cerr]
    % Cerr = [2 x nf] jackknife lower/upper confidence bounds on C (err=[2,alpha]).
    [~,Cw,~,~,~,~,fcw,confCw,~,Cerrw] = coherencyc(bw,         dff-mean(dff), pc2);
    [~,Cp,~,~,~,~,fcp,confCp,~,Cerrp] = coherencyc(ev-mean(ev), st-mean(st),  pc2);

    S(c).recName=recName; S(c).fps=fps; S(c).t=t; S(c).dff=dff; S(c).bw=bw; S(c).ev=ev;
    S(c).f_pk=f_pk; S(c).band=band2;
    S(c).fbw=fbw; S(c).Sbw=Sbw; S(c).SbwErr=SbwErr; S(c).fdd=fdd; S(c).Sdd=Sdd; S(c).SddErr=SddErr;  % waveform
    S(c).fcw=fcw(:); S(c).Cw=Cw(:); S(c).confCw=confCw; S(c).Cerrw=Cerrw;
    S(c).fbe=fbe; S(c).Sbe=Sbe; S(c).fse=fse; S(c).Sse=Sse;   % peak/spike
    S(c).fcp=fcp(:); S(c).Cp=Cp(:); S(c).confCp=confCp; S(c).Cerrp=Cerrp;
    fprintf('  ChAT[%d] %-42s fps=%.2f peak=%.2f Hz band[%.2f %.2f]\n', ...
            c, recName, fps, f_pk, band2(1), band2(2));
end

%% ---- Fig 2: dF/F + breath trace overlay (3 rows) ----
%   Two versions: the SELECTED demo window (chat_xlim) and the FULL recording.
f2     = dff_breath_fig('chat3 dFF vs breath (window)', S, chat_colors, chat_xlim, true);
f2full = dff_breath_fig('chat3 dFF vs breath (full)',   S, chat_colors, chat_xlim, false);

%% ---- Fig 3: WAVEFORM spectra + waveform coherence (3x2) ----
f3 = spec_coh_fig('chat3 waveform spectra + coherence', S, chat_colors, fmin, fmax, f_breath_search, ...
        'fbw','Sbw','SbwErr','fdd','Sdd','SddErr','fcw','Cw','confCw','Cerrw', ...
        'breath waveform + dF/F'' PSD', 'breath x Ca (waveform)');

%% ---- Fig 4: sorted heatmap (top) + triggered average (bottom) ----
f4 = figure('Color','w','Name','chat3 triggered avg + heatmap', ...
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
    mu = mean(E,1); sd = std(E,0,1);                 % +/- 1 SD across breaths
    fill(axB, [tau fliplr(tau)], [mu+sd fliplr(mu-sd)], 0.30*chat_colors(c,:)+0.70, ...
         'EdgeColor','none');                          % solid (no alpha -> vector)
    plot(axB, tau, mu, 'Color',chat_colors(c,:), 'LineWidth',1.5);
    xline(axB, 0, 'k--', 'LineWidth',0.8);
    xlim(axB,[-2 2]); ylim(axB,[-0.2 0.3]); xlabel(axB,'time from inspiration (s)');   % hard-fixed window
    ylabel(axB,'dF/F'); grid(axB,'on');
end

%% ---- Fig 5: IO reference -- breath (top) + stacked dF/F (bottom) ----
f5 = [];
io_df = dir(fullfile(io_path, '*_ch1_dFF.mat'));
io_bp = dir(fullfile(io_path, '*DLC*breath_peak_data.mat'));
if ~isempty(io_df) && ~isempty(io_bp)
    [~, io_nm] = fileparts(io_path);
    Dio  = load(fullfile(io_df(1).folder, io_df(1).name), 'dFF');
    BPio = load(fullfile(io_bp(1).folder, io_bp(1).name));
    fps_io = detect_session_fps(io_path, fallback_fps);
    dffio  = double(Dio.dFF);

    bwio = detrend(double(BPio.breath(:))); bwio(1:min(nDrop,numel(bwio)))=[]; bwio = bwio - mean(bwio);
    if isfield(BPio,'insp_onsets_train') && numel(BPio.insp_onsets_train)==numel(BPio.breath)
        evio = double(BPio.insp_onsets_train(:) ~= 0);
    else
        evio = zeros(numel(BPio.breath),1); oi=round(BPio.insp_onset_idx(:)); evio(oi(oi>=1 & oi<=numel(evio)))=1;
    end
    evio(1:min(nDrop,numel(evio))) = [];
    Tio = min([size(dffio,1), numel(bwio), numel(evio)]);
    dffio = dffio(1:Tio,:); bwio = bwio(1:Tio); evio = evio(1:Tio);
    if isempty(io_sel), sel = 1:size(dffio,2);
    else,               sel = io_sel(io_sel>=1 & io_sel<=size(dffio,2)); end
    dffio = dffio(:, sel);
    tio = (0:Tio-1)'/fps_io; Nio = numel(sel);
    mio = tio >= io_xlim(1) & tio <= io_xlim(2);      % plot ONLY the window
    tio = tio(mio); bwio = bwio(mio); dffio = dffio(mio,:);

    rngio = max(dffio,[],1)-min(dffio,[],1); spio = max(prctile(rngio,80),0.3);
    f5 = figure('Color','w','Name',['IO ' io_nm],'Units','normalized','Position',[0.05 0.08 0.9 0.82]);
    ax1 = subplot(4,1,1); hold(ax1,'on');
    plot(ax1, tio, norm01(bwio), 'Color',[0 0.35 1], 'LineWidth',0.7);
    ylabel(ax1,'breath (norm)'); ylim(ax1,[-0.05 1.05]); xlim(ax1,io_xlim);
    title(ax1, sprintf('IO  %s  (fps=%.2f, %d ROIs)', io_nm, fps_io, Nio), 'Interpreter','none');
    ax2 = subplot(4,1,2:4); hold(ax2,'on');
    for r=1:Nio, plot(ax2, tio, dffio(:,r)+(r-1)*spio, 'LineWidth',0.5); end
    ylim(ax2,[-spio Nio*spio]); set(ax2,'YTick',(0:Nio-1)*spio,'YTickLabel',sel);
    xlabel(ax2,'Time (s)'); ylabel(ax2,'ROI dF/F (stacked)'); xlim(ax2,io_xlim);
    % scale bars (solid vector): time (horizontal) + dF/F (vertical), bottom-left
    xb = io_xlim(1) + 0.03*diff(io_xlim); yb = -0.7*spio;
    plot(ax2, [xb xb+io_tbar_s], [yb yb], 'k-', 'LineWidth',2);
    plot(ax2, [xb xb], [yb yb+io_dffbar], 'k-', 'LineWidth',2);
    text(ax2, xb+io_tbar_s/2, yb, sprintf('%g s', io_tbar_s), ...
         'Horizontal','center','Vertical','top', 'FontName','Arial','FontSize',9);
    text(ax2, xb, yb+io_dffbar/2, sprintf(' %g \\DeltaF/F', io_dffbar), ...
         'Horizontal','left','Vertical','middle', 'FontName','Arial','FontSize',9);
    linkaxes([ax1 ax2],'x');
else
    warning('IO panel (Fig 5) skipped: missing *_ch1_dFF.mat / *DLC*breath_peak_data.mat in %s', io_path);
end

if doSave
    exportgraphics(f2, fullfile(rootPath,'chat3_dff_breath.png'),                 'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f2, fullfile(rootPath,'chat3_dff_breath.pdf'));
    exportgraphics(f2full, fullfile(rootPath,'chat3_dff_breath_full.png'),        'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f2full, fullfile(rootPath,'chat3_dff_breath_full.pdf'));
    exportgraphics(f3, fullfile(rootPath,'chat3_spectra_coherence_waveform.png'), 'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f3, fullfile(rootPath,'chat3_spectra_coherence_waveform.pdf'));
    exportgraphics(f4, fullfile(rootPath,'chat3_triggered_avg_heatmap.png'),      'Resolution',150, 'BackgroundColor','white');
    save_vector_pdf(f4, fullfile(rootPath,'chat3_triggered_avg_heatmap.pdf'));
    if ~isempty(f5)
        exportgraphics(f5, fullfile(rootPath,'IO_dff_breath.png'),                'Resolution',150, 'BackgroundColor','white');
        save_vector_pdf(f5, fullfile(rootPath,'IO_dff_breath.pdf'));
    end
    fprintf('Saved chat3_dff_breath / chat3_spectra_coherence_waveform / chat3_triggered_avg_heatmap / IO_dff_breath .png to %s\n', rootPath);
end

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

function plot_panel(ax, C, chat_colors, confC, ttl)
    hold(ax,'on');
    thc = linspace(0,2*pi,360);
    polarplot(ax, thc, repmat(confC,1,360), 'k--','LineWidth',1);
    for k = 1:numel(C.r)        % IO first (under ChAT)
        if C.colorIdx(k)>0 || isnan(C.th(k)) || isnan(C.r(k)), continue; end
        polarplot(ax, C.th(k), C.r(k), 'o', 'MarkerEdgeColor','k', ...
                  'MarkerFaceColor','none', 'MarkerSize',4, 'LineWidth',0.5);
    end
    for k = 1:numel(C.r)        % ChAT on top, with CI
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
% True editable vector export for Illustrator. Uses the 'painters' renderer
% (NOT exportgraphics, which rasterizes complex axes into one embedded image).
% Also writes a sibling .svg, which Illustrator opens as fully editable paths.
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
    fig = figure('Color','w','Name',figName,'Units','normalized','Position',[0.04 0.08 0.92 0.84]);
    for c = 1:nC
        ax = subplot(nC,1,c);
        if useWindow, w = chat_xlim(c,:);
        else,         w = [S(c).t(1) S(c).t(end)]; end
        mw = S(c).t >= w(1) & S(c).t <= w(2);
        tw = S(c).t(mw);
        yyaxis(ax,'right');                              % breath, normalized in-window
        plot(ax, tw, norm01(S(c).bw(mw)), '-', 'Color',[0.55 0.55 0.55], 'LineWidth',0.8);
        ylim(ax,[-0.05 1.05]); ylabel(ax,'breath (norm)');
        yyaxis(ax,'left');                               % dF/F, real scale (on top)
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
% solid light-gray detection band + black dashed breath-peak line, sent to
% back. Solid (no FaceAlpha) so exportgraphics keeps it vector in PDF.
    yl = ylim(ax);
    p = patch(ax, [band(1) band(2) band(2) band(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.90 0.90 0.90], 'EdgeColor','none');
    xline(ax, f_pk, 'k--', 'LineWidth',1);
    uistack(p,'bottom'); ylim(ax, yl);
end

function fig = spec_coh_fig(figName, S, cols, fmin, fmax, cohXlim, ...
        fBf, SBf, SBerrf, fSf, SSf, SSerrf, fcf, Cff, conff, Cerrf, leftTtl, rightTtl)
% 3x2 figure: left = power spectra (gray breath, color Ca) in dB log-x with
% shaded jackknife CI (Serr); right = coherence spectrum with shaded CI
% (Cerr) + confC line.  Field names select waveform vs peak/spike data in S.
% Breath peak (dashed) + band (shaded) on both.
    nC = numel(S);
    fig = figure('Color','w','Name',figName,'Units','normalized','Position',[0.06 0.06 0.78 0.86]);
    for c = 1:nC
        col = cols(c,:); bd = S(c).band;
        axL = subplot(nC,2,2*c-1); hold(axL,'on');
        % solid lightened CI fills (no FaceAlpha -> stays vector in PDF)
        fB = S(c).(fBf); fB=fB(:)'; eB = S(c).(SBerrf);     % breath PSD CI
        fill(axL, [fB fliplr(fB)], 10*log10([eB(1,:) fliplr(eB(2,:))]), [0.85 0.85 0.85], ...
             'EdgeColor','none');
        fS = S(c).(fSf); fS=fS(:)'; eS = S(c).(SSerrf);     % Ca PSD CI
        fill(axL, [fS fliplr(fS)], 10*log10([eS(1,:) fliplr(eS(2,:))]), 0.30*col+0.70, ...
             'EdgeColor','none');
        plot(axL, fB, 10*log10(S(c).(SBf)), 'Color',[0.6 0.6 0.6], 'LineWidth',0.9);
        plot(axL, fS, 10*log10(S(c).(SSf)), 'Color',col,           'LineWidth',1.1);
        set(axL,'XScale','log'); xlim(axL,[fmin fmax]);
        xticks(axL,[0.1 0.3 1 3 10]); set(axL,'XMinorTick','off');  % ~x3 -> evenly spaced on log
        add_band(axL, bd, S(c).f_pk); ylabel(axL,'power (dB)');
        pbaspect(axL,[1 1 1]);                       % 1:1 plot box
        title(axL, sprintf('%s  %s', S(c).recName, leftTtl), 'Interpreter','none');
        if c==nC, xlabel(axL,'Frequency (Hz)'); end

        axR = subplot(nC,2,2*c); hold(axR,'on');
        fc = S(c).(fcf); fc = fc(:)';
        ce = S(c).(Cerrf);                          % [2 x nf] jackknife CI
        fill(axR, [fc fliplr(fc)], [ce(1,:) fliplr(ce(2,:))], 0.30*col+0.70, ...
             'EdgeColor','none');                    % solid (no alpha -> vector)
        plot(axR, fc, S(c).(Cff), 'Color',col, 'LineWidth',1.0);
        yline(axR, S(c).(conff), 'k--', 'LineWidth',1);
        set(axR,'XScale','log'); xlim(axR,cohXlim); ylim(axR,[0 1]);
        xticks(axR, [0.25 0.5 1 2 4]); set(axR,'XMinorTick','off');  % geometric -> evenly spaced on log
        add_band(axR, bd, S(c).f_pk); ylabel(axR,'coherence');
        pbaspect(axR,[1 1 1]);                       % 1:1 plot box
        title(axR, sprintf('%s  confC=%.2f (shaded = jackknife CI)', ...
              rightTtl, S(c).(conff)), 'Interpreter','none');
        if c==nC, xlabel(axR,'Frequency (Hz)'); end
    end
end
