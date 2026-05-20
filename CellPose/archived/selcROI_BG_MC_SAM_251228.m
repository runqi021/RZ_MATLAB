clear; clc; delete(gcp('nocreate'));

filePath = "C:\Users\zhang\Desktop\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001_ch1\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001_ch1.tif";
F_raw = loadtiff(filePath);
F = F_raw;
%%
trimL = 10;
trimR = 10;
%% channel
F(:,1:trimL,:) = []; F(:,end-trimR+1:end,:) = [];
%F1 = F(:, :, 1:2:end); 
%F3 = F(:, :, 2:2:end);
F1 = F;

darkOut1 = estimate_dark_from_movie(F1, ...
    'MaxFrames', 500, ...   % sample up to 200 frames
    'Subsample', 4);        % every 4th pixel in x,y

dark_mean1 = darkOut1.dark_mean;   % dark current
dark_std1  = darkOut1.dark_std;    % dark noise

%darkOut3 = estimate_dark_from_movie(F3, ...
%    'MaxFrames', 500, ...   % sample up to 200 frames
%    'Subsample', 4);        % every 4th pixel in x,y
%
%dark_mean3 = darkOut3.dark_mean;   % dark current
%dark_std3  = darkOut3.dark_std;    % dark noise

%% Dark subtraction
F1 = F1 - uint16(dark_mean1); %F1 = F1 - uint16(dark_mean1);
%F3 = F3 - int16(dark_mean3);

F1 = max(F1, 0); % negative mask
%F3 = max(F3, 0);

%% per pixel

thr1 = 3 * dark_std1;
meanF1 = mean(F1, 3);
mask1 = meanF1 > thr1;
fracKeep = nnz(mask1) / numel(mask1) * 100;
fprintf('[darkMask] Keeping %.2f %% of pixels (mean > %.2f)\n', fracKeep, thr1);

F1(~mask1(:,:,ones(1,size(F1,3)))) = 0;

%%
sumProj = sum(F1, 3);

figure('Color','w');
imagesc(sumProj);
axis image off;
colormap(gray);
hold on;

badMask = ~mask1;
[r,c] = find(badMask);
scatter(c, r, 2, 'r', 'filled');

title(sprintf('Sum projection with masked pixels (%.2f%%)', ...
    100*nnz(badMask)/numel(badMask)));
hold off;

%%
F1_sat = 65535 - dark_mean1;  % unit16: 65535,  int16: 32767
F1_plot = F1;
F1_plot(1:30,:)=[];
Fhist(F1_plot, sat = F1_sat);
%% save trimmed video
[folderPath, baseName, ~] = fileparts(char(filePath));
outPath1 = fullfile(folderPath, [baseName '_ch1_minusDark.tif']);  % pure chars
%outPath3 = fullfile(folderPath, [baseName '_ch3_minusDark.tif']);  % pure chars

%%
opts = struct('overwrite', true, 'message', true);
saveastiff(F1, outPath1, opts);
%saveastiff(F3, outPath3, opts);

fprintf('Saved dark-subtracted TIFF in :\n%s\n', folderPath);

%% now run MC on the dark-corrected movie
mcOut1 = run_rigid_mc(outPath1);
%mcOut3 = run_rigid_mc(outPath3);

%%
sh1 = mcOut1.shifts;

T = numel(sh1);
sx1 = zeros(T,1);
sy1 = zeros(T,1);

for t = 1:T
    s1 = sh1(t).shifts;
    sx1(t) = s1(1);
    sy1(t) = s1(2);
end

r1 = hypot(sx1, sy1);

% sh3 = mcOut3.shifts;
% 
% T = numel(sh3);
% sx3 = zeros(T,1);
% sy3 = zeros(T,1);
% 
% for t = 1:T
%     s3 = sh3(t).shifts;
%     sx3(t) = s3(1);
%     sy3(t) = s3(2);
% end
% 
% r3 = sqrt(sx3.^2 + sy3.^2);

%%
figure; plot(r1,'LineWidth',1.5); ylabel('|shift| (pixel)'); 
xlabel('Frame'); grid on

%%
outPath1 = fullfile(folderPath, [baseName '_ch1_minusDark_MC.tif']);  % pure chars

out_MC = cp_sam_extract_F_cli(outPath1, ...
    'FPS',30, ...
    'Diameter',30, ...
    'PythonExe','C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe', ...
    'UseGPU',true, ...
    'FlowThreshold',0.50, ...
    'CellprobThreshold',0);

