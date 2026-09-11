function breath_cache_diffmap(nFr)
% breath_cache_diffmap  Compute the full-resolution inspiration-baseline diff map
% (top nFr highest-signal frames - bottom nFr lowest) + mean camera frame from
% each breath_crop.mat cube, and append them into the sibling breath_pc1.mat as
% diffImg / camImg. The peak GUI then displays the cached full-res diffImg next
% to the eigenimage (no 178 MB cube load at GUI time).
if nargin<1 || isempty(nFr), nFr=100; end
roots={'D:\Ventral_surface_summary\Vglut2', ...
       'D:\Ventral_surface_summary\ChAT'};
hits=[];
for ri=1:numel(roots), hits=[hits; dir(fullfile(roots{ri},'**','breath_pc1.mat'))]; end %#ok<AGROW>
for j=1:numel(hits)
    folder=hits(j).folder;
    cropf=fullfile(folder,'breath_crop.mat');
    if ~isfile(cropf), fprintf('skip (no cube): %s\n', folder); continue; end
    C=load(cropf,'mov'); mov=C.mov;
    P=load(fullfile(folder,'breath_pc1.mat'),'breathTrace'); bt=double(P.breathTrace(:));
    [~,~,T]=size(mov); T2=min(T,numel(bt)); mov=mov(:,:,1:T2); bt=bt(1:T2);
    z=(bt-mean(bt))/std(bt); if mean(z.^3)<0, bt=-bt; end      % inspiration = high side
    nf=min(nFr,floor(T2/2));
    [~,ord]=sort(bt,'descend'); top=ord(1:nf); bot=ord(end-nf+1:end);
    camImg  = mean(single(mov),3);                             %#ok<NASGU>
    diffImg = mean(single(mov(:,:,top)),3) - mean(single(mov(:,:,bot)),3); %#ok<NASGU>
    save(fullfile(folder,'breath_pc1.mat'),'diffImg','camImg','-append');
    [~,leaf]=fileparts(folder); fprintf('cached diffImg -> %s\n', leaf);
end
fprintf('Done.\n');
end
