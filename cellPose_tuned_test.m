addpath(genpath('C:\Users\Admin\Desktop\RZ_MATLAB'));                                                                                                  
                                                                                                                                                   
modelPath = 'D:\cellpose_training\models\cellpose_1775197773.3387594';
baseDir   = 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed';
saveDir   = fullfile(baseDir, 'cellpose_finetuned_test');
if ~isfolder(saveDir), mkdir(saveDir); end

d = dir(fullfile(baseDir, '**', '*_preproc_MC_MC.tif'));
fprintf('Found %d FOVs, saving to:\n  %s\n', numel(d), saveDir);

for i = 1:numel(d)
  tiffPath = fullfile(d(i).folder, d(i).name);
  fprintf('\n[%d/%d] %s\n', i, numel(d), d(i).name);
  try
      out = cp_sam_extract_F_cli_finetuned(tiffPath, ...
          'ModelPath', modelPath, ...
          'DiameterUm', 16, ...
          'Savedir', saveDir);
  catch ME
      fprintf(2, '  >> FAILED: %s\n', ME.message);
  end
end
fprintf('\nDone.\n');