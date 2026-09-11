%% Run_SVD_mode_psd_260811.m
% Companion to Run_SVD_recon_260811.m.
%
% Reads the saved SVD_result.mat of each variant and plots the power
% spectrum of the leading temporal modes v_k(t). This is what tells you
% WHAT the top modes are: a dipole eigenimage whose v_k peaks near the
% heart rate (~8-12 Hz) is pulsation; near the breath rate (~3-8 Hz) is
% respiration; a peak pinned at Nyquist is a scan/alias artifact.
%
% No re-run of the SVD is needed — this only touches the saved MATs.

clear; clc;

%% USER PARAMS
outRoot   = "D:\Ventral_surface_summary\Vgat\0730\cell1\roi1_3.2x_x1100y1050_z6_6000f_12lp_00001\svd_check_260811";
labels    = {'raw','mcmc'};
nShowMode = 6;        % leading modes to spectrally analyse
winSec    = 17;       % pwelch window (s) -> ~0.06 Hz resolution at 30 fps

outRoot = char(outRoot);
nV = numel(labels);

f = figure('Color','w','Visible','off','Position',[60 60 1500 260*nShowMode]);
tl = tiledlayout(f, nShowMode, 1, 'TileSpacing','compact', 'Padding','compact');
cols = lines(nV);

pk = nan(nShowMode, nV);

for k = 1:nShowMode
    ax = nexttile(tl); hold(ax,'on');
    for iv = 1:nV
        S = load(fullfile(outRoot, labels{iv}, 'SVD_result.mat'));
        R = S.SVD_result;
        fps = R.fps;

        v   = double(R.V(:,k));
        nw  = min(round(winSec*fps), numel(v));
        [pxx, fx] = pwelch(v - mean(v), hamming(nw), round(nw/2), 4096, fps);

        plot(ax, fx, 10*log10(pxx), '-', 'LineWidth',1.3, 'Color',cols(iv,:), ...
             'DisplayName', sprintf('%s (%.2f%% var)', labels{iv}, 100*R.varExp(k)));

        keep = fx > 0.3;                       % ignore the DC/drift shoulder
        [~, im] = max(pxx(keep));
        fk = fx(keep);
        pk(k, iv) = fk(im);
    end
    xlim(ax, [0 fps/2]);
    grid(ax,'on');
    ylabel(ax, sprintf('v_{%d}  PSD (dB)', k));
    legend(ax, 'Location','northeast');
    title(ax, sprintf('mode %d   peak: %s', k, ...
          strjoin(arrayfun(@(x) sprintf('%.2f Hz', x), pk(k,:), ...
                           'UniformOutput', false), ' / ')), 'FontSize',9);
    if k == nShowMode; xlabel(ax, 'frequency (Hz)'); end
end

title(tl, 'temporal modes v_k(t) — power spectra', 'FontWeight','bold');
outPNG = fullfile(outRoot, 'mode_psd.png');
exportgraphics(f, outPNG, 'Resolution',180);
close(f);

fprintf('\n[PSD] peak frequency (Hz), rows = mode, cols = %s\n', strjoin(labels,' / '));
disp(pk);
fprintf('[PSD] Saved -> %s\n', outPNG);
