% representative_traces_avgproj_260527.m
% -----------------------------------------------------------------------
%  Representative single-ROI panels across MANY dffQC FOVs.
%
%  For a hand-picked list of (FOV folder, ROI index) pairs, draws for EACH:
%     - AVG   : avg projection done like avg_proj_scalebar_260527.m
%               (crop crop_um, clip percentiles, gamma, burned scale bar),
%               with THAT ROI outlined in the trace's color.
%     - TRACE : that ROI's dF/F, FULL recording
%     - TRIG  : breath peak-triggered dF/F average +/- SD, +/-1 breath period,
%               square panel (1:1). Method copied from
%               chat_breath_coherence_polar_260526.m (Fig 4 bottom).
%  ...then stacks all of them into one combined summary figure.
%
%  Pairing (verified): dF/F came from the SINGLE-MC cpSAM, so
%     maskL  <- *_ch1_preproc_MC_cpSAM_output.mat
%     avg    <- *_ch1_preproc_MC_AVG_for_CP.tif   (single-MC mean; maskL space)
%  ROI index k == maskL label k == dFF column k.
%
%  Breath/dF/F alignment: the breath cam is 2P-triggered frame-per-frame on the
%  RAW movie, so breath is longer than dF/F by exactly the tossed frames. We
%  drop (len_breath - T_dff) frames off the FRONT of breath, then truncate.
%
%  Outputs:
%     <FOV>\repr_ROI##_trace_avgproj.pdf / .png      (one per FOV)
%     <outDir>\representative_summary.pdf / .png      (combined)
%
%  Standalone. Image Processing Toolbox (bwboundaries). Runqi Zhang / 2026.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath')); addpath(scriptDir);
addpath(genpath(fullfile(scriptDir,'chronux_2_12')));   % mtspectrumc (breath PSD peak)

%% ===================== USER-EDITABLE =====================
% {folder, ROI index, group label}
items = {
 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing\6x_-820-1070-z20_00001',          4, 'SST';
 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing\6x_-743-840-z30_00001',           2, 'SST';
 %'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\random\7x_-850-1032-z30_3000f_00001',             2, 'SST';
 'D:\batch_dffQC_test_260325\260330_sst_soma_g8s\phys\maybe_breathing\-1050-760-z5_5x_6000f_15lp_00001',             3, 'SST';
 'D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed\730-930-z25_9x_3000f_00001',                          3, 'Vgat';
 'D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed\760-750-z30_6x_3000f_00001',                          3, 'Vgat';
 'D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed\-500-990-z-30_8x_3000f_00001',                        1, 'Vgat';
};

% Per-trace plotting window (s) for the SELECTED-window figure, SAME row order
% as the (uncommented) items above. [NaN NaN] row -> use that trace's full range.
sel_xlim = [
   20  80;    % 1  6x_-820-1070-z20  ROI4  SST
   20  80;    % 2  6x_-743-840-z30   ROI2  SST
   20  80;    % 3  -1050-760-z5_5x   ROI3  SST
   35  95;    % 4  730-930-z25_9x    ROI3  Vgat
   20  80;    % 5  760-750-z30_6x    ROI3  Vgat
   28  88;    % 6  -500-990-z-30_8x  ROI1  Vgat
];

outDir   = 'D:\batch_dffQC_test_260325\representative_260527';

fallback_fps  = 30;
PixelSizeBase = 1.7778;     % um/px @ zoom 1 (fallback if no pixelSize_um)

% ---- avg-projection display (matches avg_proj_scalebar_260527.m) ----
clip_pct    = [0.5 99.9];   % intensity clip percentiles
gamma_val   = 0.6;          % display gamma (<1 brightens midtones)
crop_um     = 5;            % crop from EACH side (um); 0 = full FOV
scaleBar_um = 50;           % white scale bar on the projection (um)
outlineLW   = 1.4;          % ROI outline width

% ---- trace display ----
dffLW       = 0.8;          % dF/F line width

% ---- breath peak-triggered average (matches chat_breath_coherence_polar) ----
nDrop          = 30;          % breath frames tossed up front (align to dFF)
TW_breath      = 6;           % multitaper TW for breath-PSD peak
f_breath_search= [0.2 4];     % Hz, breath PSD peak search band
fmin_b         = 0.05;        % Hz, breath PSD lower bound
fmax_b         = 15;          % Hz, breath PSD upper bound
trigLW         = 1.2;         % triggered-average line width

doSave   = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if doSave && ~isfolder(outDir), mkdir(outDir); end
nItem  = size(items,1);
colors = roi_colormap(nItem);          % one hue per FOV/ROI (shared trace<->outline)

P = struct('nm',{},'grp',{},'roi',{},'t',{},'tr',{},'img',{},'B',{},'bar',{}, ...
           'tau',{},'mu',{},'sd',{},'non',{},'fpk',{});

for i = 1:nItem
    d=items{i,1}; roi=items{i,2}; grp=items{i,3};
    assert(isfolder(d),'Folder not found: %s', d);

    dff_hit = dir(fullfile(d,'*_ch1_dFF.mat'));
    sam_hit = dir(fullfile(d,'*_ch1_preproc_MC_cpSAM_output.mat'));
    avg_hit = dir(fullfile(d,'*_ch1_preproc_MC_AVG_for_CP.tif'));
    assert(~isempty(dff_hit),'No *_ch1_dFF.mat in %s', d);
    assert(~isempty(sam_hit),'No *_ch1_preproc_MC_cpSAM_output.mat in %s', d);
    assert(~isempty(avg_hit),'No *_ch1_preproc_MC_AVG_for_CP.tif in %s', d);
    stem = erase(dff_hit(1).name,'_ch1_dFF.mat');

    % --- fps + pixel size from meta ---
    fps=fallback_fps; px_um=NaN;
    mhit=dir(fullfile(d,[stem '_ch1_meta.mat'])); if isempty(mhit), mhit=dir(fullfile(d,[stem '_meta.mat'])); end
    if ~isempty(mhit)
        M=load(fullfile(mhit(1).folder,mhit(1).name));
        if isfield(M,'fps')&&isfinite(M.fps)&&M.fps>0, fps=M.fps; end
        if isfield(M,'pixelSize_um')&&isfinite(M.pixelSize_um)&&M.pixelSize_um>0, px_um=M.pixelSize_um;
        elseif isfield(M,'zoomFactor')&&isfinite(M.zoomFactor)&&M.zoomFactor>0, px_um=PixelSizeBase/M.zoomFactor; end
    end
    assert(isfinite(px_um),'No pixel size for %s', stem);

    % --- dF/F trace ---
    D=load(fullfile(dff_hit(1).folder,dff_hit(1).name),'dFF'); dff=double(D.dFF);
    assert(roi>=1&&roi<=size(dff,2),'ROI %d out of range (1..%d) in %s',roi,size(dff,2),stem);
    tr=dff(:,roi); T=numel(tr);
    t=(0:T-1)'/fps;

    % --- breath onsets + waveform, peak-triggered dF/F average ---
    % align: breath cam is 2P-triggered, longer by the tossed frames -> drop
    % nDrop off the FRONT of breath/onsets, then truncate to common length.
    tau=[]; mu=[]; sd=[]; non=0; fpk=NaN;
    bp_hit=dir(fullfile(d,'*breath_peak_data.mat'));
    if isempty(bp_hit)
        warning('No *breath_peak_data.mat in %s -- trig-avg panel blank.', d);
    else
        BP=load(fullfile(bp_hit(1).folder,bp_hit(1).name));
        bw=detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw)))=[]; bw=bw-mean(bw);
        if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==numel(BP.breath)
            ev=double(BP.insp_onsets_train(:)~=0);
        else
            ev=zeros(numel(BP.breath),1); oi=round(BP.insp_onset_idx(:)); ev(oi(oi>=1&oi<=numel(ev)))=1;
        end
        ev(1:min(nDrop,numel(ev)))=[];
        Tb=min([numel(ev),numel(bw),T]); ev=ev(1:Tb); bw=bw(1:Tb); trb=tr(1:Tb);
        % breath rate from waveform PSD peak -> +/-1 breath period window
        pB.Fs=fps; pB.tapers=[TW_breath,2*TW_breath-1]; pB.pad=0;
        pB.fpass=[fmin_b,min(fmax_b,fps/2)]; pB.err=0;
        [Sb,fb]=mtspectrumc(bw,pB); Sb=Sb(:); fb=fb(:);
        msk=fb>=f_breath_search(1)&fb<=f_breath_search(2);
        [~,rl]=max(Sb(msk)); ip=find(msk,1)+rl-1; fpk=fb(ip);
        win=max(1,round(2*fps)); tau=(-win:win)/fps;   % fixed +/-2 s window
        on=find(ev>0); on=on(on-win>=1 & on+win<=Tb);
        if isempty(on)
            warning('No valid breath onsets for %s -- trig-avg blank.', stem); tau=[];
        else
            E=zeros(numel(on),2*win+1);
            for k=1:numel(on), E(k,:)=trb(on(k)-win:on(k)+win); end
            mu=mean(E,1); sd=std(E,0,1); non=numel(on);
        end
    end

    % --- avg projection: crop -> clip -> gamma ---
    avg=double(imread(fullfile(avg_hit(1).folder,avg_hit(1).name)));
    cpx=round(crop_um/px_um);
    if cpx>0 && 2*cpx<min(size(avg)), avg=avg(cpx+1:end-cpx,cpx+1:end-cpx); else, cpx=0; end
    lo=prctile(avg(:),clip_pct(1)); hi=prctile(avg(:),clip_pct(2));
    img=min(max((avg-lo)/max(hi-lo,eps),0),1).^gamma_val;
    [H,W]=size(img);

    % --- ROI outline (shift for crop) ---
    S=load(fullfile(sam_hit(1).folder,sam_hit(1).name),'maskL');
    B=bwboundaries(S.maskL==roi,'noholes'); for k=1:numel(B), B{k}=B{k}-cpx; end
    if isempty(B), warning('ROI %d has no pixels in %s maskL.',roi,stem); end

    % --- scale-bar geometry (image coords) ---
    barLen=min(max(1,round(scaleBar_um/px_um)),W-2);
    margin=round(0.04*H); barThk=max(3,round(0.012*H));
    bar=[margin, H-margin-barThk, barLen, barThk];

    P(i)=struct('nm',stem,'grp',grp,'roi',roi,'t',t,'tr',tr,'img',img,'B',{B},'bar',bar, ...
                'tau',tau,'mu',mu,'sd',sd,'non',non,'fpk',fpk);
    fprintf('%d) %s ROI%d | T=%d @%.3g Hz | %.4f um/px | breath %.2f Hz n=%d\n', ...
            i,stem,roi,T,fps,px_um,fpk,non);

    % ---------------- per-FOV figure ----------------
    fA=figure('Color','w','Name',sprintf('%s ROI%d',stem,roi),'Units','normalized','Position',[0.06 0.25 0.86 0.42]);
    tl=tiledlayout(fA,1,6,'TileSpacing','compact','Padding','compact');
    ax1=nexttile(tl,1); draw_proj(ax1,img,bar); title(ax1,sprintf('%s | ROI%d',grp,roi));
    ax2=nexttile(tl,2,[1 4]); draw_trace(ax2,t,tr,colors(i,:),dffLW);
    title(ax2,stem,'Interpreter','none','FontSize',9);
    ax3=nexttile(tl,6); draw_trigavg(ax3,tau,mu,sd,colors(i,:),fpk,non,trigLW);
    if doSave
        exportgraphics(fA,fullfile(d,sprintf('repr_ROI%02d_trace_avgproj.pdf',roi)),'ContentType','vector','BackgroundColor','white');
        exportgraphics(fA,fullfile(d,sprintf('repr_ROI%02d_trace_avgproj.png',roi)),'Resolution',300,'BackgroundColor','white');
        % N per-input PNGs = AVG PROJECTION ONLY, native resolution, scalebar burned in
        ip=img; rr=bar(2)+1:bar(2)+bar(4); cc=bar(1)+1:bar(1)+bar(3); ip(rr,cc)=1;
        imwrite(uint16(round(ip*65535)), fullfile(outDir,sprintf('avgproj_%02d_%s_ROI%02d.png',i,stem,roi)));
    end
end

%% ===================== COMBINED SUMMARY =====================
%   Two versions: FULL trace and the SELECTED window (sel_xlim) per trace.
fFull = summary_fig('Representative summary (full)',   '\DeltaF/F full recording',  P, colors, dffLW, trigLW, sel_xlim, false);
fSel  = summary_fig('Representative summary (window)', '\DeltaF/F selected window',  P, colors, dffLW, trigLW, sel_xlim, true);

if doSave
    exportgraphics(fFull,fullfile(outDir,'representative_summary_full.pdf'),'ContentType','vector','BackgroundColor','white');
    exportgraphics(fFull,fullfile(outDir,'representative_summary_full.png'),'Resolution',300,'BackgroundColor','white');
    exportgraphics(fSel, fullfile(outDir,'representative_summary_window.pdf'),'ContentType','vector','BackgroundColor','white');
    exportgraphics(fSel, fullfile(outDir,'representative_summary_window.png'),'Resolution',300,'BackgroundColor','white');
    save(fullfile(outDir,'representative_traces.mat'),'items','sel_xlim','colors','clip_pct','gamma_val','crop_um','scaleBar_um', ...
         'P','nDrop','TW_breath','f_breath_search');
    fprintf('Saved full + window summaries + per-FOV panels to %s\n',outDir);
end

%% ========================= LOCAL FUNCTIONS =========================
function draw_proj(ax,img,bar)
    imshow(img,[0 1],'Parent',ax); colormap(ax,gray(256)); hold(ax,'on');
    rectangle(ax,'Position',bar,'FaceColor','w','EdgeColor','none'); hold(ax,'off');
end

function draw_trace(ax,t,tr,dffCol,dffLW,w)
% dF/F trace; w = [t0 t1] plot window (default / NaN -> full trace).
    if nargin<6 || isempty(w) || any(isnan(w)), w=[t(1) t(end)]; end
    m = t>=w(1) & t<=w(2);                          % clip to window (no off-axis paths in PDF)
    plot(ax,t(m),tr(m),'-','Color',dffCol,'LineWidth',dffLW);
    set(ax,'YColor','k'); ylabel(ax,'\DeltaF/F');
    xlim(ax,w); xlabel(ax,'Time (s)'); box(ax,'off');
end

function fig = summary_fig(figName, ttl, P, colors, dffLW, trigLW, sel_xlim, useWindow)
% Combined nItem x 6 summary: avg proj | dF/F trace (full or windowed) | trig-avg.
    nItem = numel(P);
    fig = figure('Color','w','Name',figName,'Units','normalized','Position',[0.02 0.04 0.94 0.92]);
    tl  = tiledlayout(fig,nItem,6,'TileSpacing','compact','Padding','compact');
    title(tl, ['Representative ROIs   |   avg proj + ' ttl ' + breath-triggered avg'],'FontWeight','bold');
    for i=1:nItem
        p=P(i);
        if useWindow && i<=size(sel_xlim,1), w=sel_xlim(i,:); else, w=[p.t(1) p.t(end)]; end
        ax1=nexttile(tl,(i-1)*6+1); draw_proj(ax1,p.img,p.bar);
        ax2=nexttile(tl,(i-1)*6+2,[1 4]); draw_trace(ax2,p.t,p.tr,colors(i,:),dffLW,w);
        if i<nItem, xlabel(ax2,''); set(ax2,'XTickLabel',[]); end
        text(ax2,0,1,p.nm,'Units','normalized','VerticalAlignment','bottom','FontSize',7,'Interpreter','none','Color',[.3 .3 .3]);
        ax3=nexttile(tl,(i-1)*6+6); draw_trigavg(ax3,p.tau,p.mu,p.sd,colors(i,:),p.fpk,p.non,trigLW);
        if i<nItem, xlabel(ax3,''); end
    end
end

function draw_trigavg(ax, tau, mu, sd, col, fpk, non, lw)
% breath peak-triggered dF/F average +/- SD, square (1:1) panel.
% Method: chat_breath_coherence_polar_260526.m Fig 4 bottom.
    if isempty(tau)
        text(ax,0.5,0.5,'no breath','Units','normalized','HorizontalAlignment','center', ...
             'Color',[.6 .6 .6],'FontSize',8); axis(ax,'off'); return;
    end
    mu=mu(:)'; sd=sd(:)'; tau=tau(:)';
    hold(ax,'on');
    fill(ax,[tau fliplr(tau)],[mu+sd fliplr(mu-sd)],0.30*col+0.70,'EdgeColor','none');  % solid -> vector
    plot(ax,tau,mu,'-','Color',col,'LineWidth',lw);
    xline(ax,0,'k--','LineWidth',0.8);
    xlim(ax,[-2 2]); ylim(ax,[-0.2 0.3]);   % hard-fixed window + dF/F range
    box(ax,'off'); pbaspect(ax,[1 1 1]);  % 1:1 height/width
    xlabel(ax,'t - insp (s)'); ylabel(ax,'\DeltaF/F');
    title(ax,sprintf('%.2f Hz, n=%d',fpk,non),'FontSize',8,'FontWeight','normal');
end

function cmap = roi_colormap(N)
    if N<=0, cmap=zeros(0,3); return; end
    gr=0.618033988749895; h=mod((0:N-1)*gr,1);
    cmap=hsv2rgb([h(:) 0.85*ones(N,1) 0.90*ones(N,1)]);
end
