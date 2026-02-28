% Batch
folderPath = 'D:\RUNQI\phys\processed\sound';
list = dir(folderPath);

ct = 0;
for k = 1:numel(list)
    dlc_csv = []; SAM =[];
    if strcmp(list(k).name,'.') || strcmp(list(k).name,'..'); continue;
    else
        currentPath = fullfile(folderPath, list(k).name);
        sublist = dir(currentPath);
        for i = 1:numel(sublist)
             if endsWith(sublist(i).name, 'cpSAM_output.mat')
                SAM = fullfile(currentPath, sublist(i).name); 
            end
        end

        if isempty(SAM); continue; end
               
        %% Run
        load(SAM); ct = ct + 1;

        n_pulse = 10;
        T_base = 10;
        T_pulse = .25;
        T_cycle = 10;
        
        fs = 30 % Hz
        toss = 5; % Sec
        F(1:fs*toss, :) = [];
        F(end-98:end, :) = [];
        
        [T, N] = size(F);
        t = (0:T-1) / fs;
        
        dFF = dFF_RZ(F).dFF;
        
        t_pulse = linspace(T_base,T_base+T_cycle-1, n_pulse);
        
        %% stackDFF and stim overlay
        stackDFF(dFF);
        
        for i = 1:numel(t_pulse)
            x1 = t_pulse(i);
            x2 = x1 + T_pulse;
            yl = ylim;
            h = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1 1 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7);
            uistack(h,'bottom');
        end
        
        outFig = fullfile(currentPath, 'stackdFF_Stim.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        outFig = fullfile(currentPath, 'stackDFF_Stim.fig')
        savefig(gcf, outFig);
        
        close all;
        %fprintf('Saved:\n  %s\n', outFig);
        
        %%
        hmapDFF(dFF);
        for i = 1:numel(t_pulse)
            xline(t_pulse(i), 'LineWidth',3);
            % x1 = t_pulse(i);
            % x2 = x1 + T_pulse;
            % yl = ylim;
            % h = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1 1 0], ...
            % 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            % %uistack(h,'bottom');
        end
        
        outFig = fullfile(currentPath, 'heatMap_Stim.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        outFig = fullfile(currentPath, 'heatMap_Stim.fig')
        savefig(gcf, outFig);
        close all;
        
        %% full trace average
        dFF_avg = mean(dFF, 2);
        dFF_std = std(dFF,0,2);
        plot(t, dFF_avg, 'Color','k'); hold on;
        y1 = dFF_avg - dFF_std; y2 = dFF_avg + dFF_std;
        y1 = y1'; y2 = y2';
        X = [t, fliplr(t)];
        Y = [y1, fliplr(y2)];
        fill(X, Y, [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        for i = 1:numel(t_pulse)
            x1 = t_pulse(i);
            x2 = x1 + T_pulse;
            yl = ylim;
            h = patch([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [1 1 0], ...
            'EdgeColor', 'none', 'FaceAlpha', 0.7);
            uistack(h,'bottom');
        end
        
        outFig = fullfile(currentPath, 'timeSerires_avg.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        outFig = fullfile(currentPath, 'timeSeries_avg.fig')
        savefig(gcf, outFig)
        
        close all;
        
        %% triggered average
        pulse_onsets = t_pulse * fs;
        
        win = 1*fs;
        nEv = numel(pulse_onsets);
        
        fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);
        
        dff_seg_all = cell(N,1);
        
        for i = 1:N
            roi_trace = dFF(:, i);          % RAW ΔF/F
        
            seg = zeros(nEv, win+1);
            for e = 1:nEv
                c = pulse_onsets(e);
                seg(e,:) = roi_trace(c : c+win);
            end
        
            dff_seg_all{i} = seg;
        end
        
        t = (0:win)/fs;
        ct=0;
        figure('Color', 'white');
        for i = 1:N
            seg = dff_seg_all{i};
            avg = mean(seg,1);
            plot(t, avg, 'color', [0.8, 0.8, 0.8]); hold on;
            ct=ct+1
            seg_avg(i,:) = avg;
        end
        
        pop_avg = mean(seg_avg,1); hold on;
        pop_std = std(seg_avg, 0, 1); hold on;
        plot(t, pop_avg, 'Color', 'k');
        y1 = pop_avg - pop_std; y2 = pop_avg + pop_std;
        
        X = [t, fliplr(t)];
        Y = [y1, fliplr(y2)];
        
        fill(X, Y, [0.7 0.7 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        xlabel('Time (s)');
        ylabel('ΔF/F');
        title('Event-Triggered Average with Standard Deviation');
        axis square
        
        outFig = fullfile(currentPath, 'triggered_avg.png')
        exportgraphics(gcf, outFig, 'ContentType','vector');
        outFig = fullfile(currentPath, 'triggered_avg.fig')
        savefig(gcf, outFig)
        
        close all;
        %fprintf("%d ROI average overlay \n", ct);
       
%%
    end
end