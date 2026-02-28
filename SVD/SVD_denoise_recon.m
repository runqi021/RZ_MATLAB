k = 20;
kList = 1:k;

U_F = SVD_result.U; S_F = SVD_result.S; V_F = SVD_result.V;

Y_recon = U_F(:,kList) * S_F(kList,kList) * V_F(:,kList)'; 
Y_recon = Y_recon + SVD_result.mu;
%lims = [prctile(Y_recon(:), 0.5), prctile(Y_recon(:), 99.5)];
%%
Hc = SVD_result.Hc; Wc = SVD_result.Wc; T = SVD_result.T; 

Y_recon_stack = zeros(Hc, Wc, T);

for i = 1:T
    Y_recon_stack(:,:,i) = reshape(Y_recon(:,i), [Hc, Wc]);
end

Y_uint16 = uint16(Y_recon_stack);

srcTif = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\raw\fov2_-200_L_30lp_1100nm_5x_00003.tif";
[folderPath, baseName, ~] = fileparts(srcTif);
outPath = fullfile(folderPath, baseName + "_SVDk20_denoised.tif");

opts = struct('overwrite', true, 'message', true);

saveastiff(Y_uint16, char(outPath), opts);   % <-- char() fixes saveastiff delete bug
