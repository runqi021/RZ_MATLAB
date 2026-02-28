clear
gcp;  % optional parallel pool

% ===== USER PATH =====
imgDir = "E:\RZ\data\shi\260115_homo_cal590\test";  % your folder with TIFFs
% ======================

files = [dir(fullfile(imgDir, '*.tif')); dir(fullfile(imgDir, '*.tiff'))];

for k = 1:numel(files)
    inFile = fullfile(imgDir, files(k).name);
    [~, base, ext] = fileparts(inFile);
    outFile = fullfile(imgDir, [base '_MC' ext]);

    if exist(outFile, 'file')
        fprintf('Skipping %s (already motion-corrected)\n', files(k).name);
        continue
    end

    fprintf('\n=== Processing %s ===\n', files(k).name);

    % --- read image stack ---
    tic;
    Y = read_file(inFile);
    Y = single(Y);
    Y = Y - min(Y(:));   % ensure nonnegative
    T = size(Y, ndims(Y));
    fprintf('Loaded %d frames in %.2f s\n', T, toc);

    % --- set rigid NoRMCorre parameters ---
    options_rigid = NoRMCorreSetParms( ...
        'd1', size(Y,1), ...
        'd2', size(Y,2), ...
        'bin_width', 200, ...
        'max_shift', 15, ...
        'us_fac', 50, ...
        'init_batch', 200);

    % --- perform motion correction ---
    tic;
    [M1, ~, ~, ~] = normcorre(Y, options_rigid);
    fprintf('Motion correction done in %.2f s\n', toc);

    % --- save as 16-bit TIFF ---
    mc = M1;
    mc(mc < 0) = 0;
    mc(mc > 65535) = 65535;
    mc = uint16(mc);

    fprintf('Saving %s ...\n', outFile);
    T = size(mc,3);
    for t = 1:T
        if t == 1
            imwrite(mc(:,:,t), outFile, 'tif', 'Compression', 'none');
        else
            imwrite(mc(:,:,t), outFile, 'tif', 'WriteMode','append','Compression','none');
        end
    end
    fprintf('Saved motion-corrected file: %s\n', outFile);
end

fprintf('\nAll done.\n');
