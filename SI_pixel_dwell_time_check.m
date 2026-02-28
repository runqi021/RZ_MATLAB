fn = "C:\Users\Admin\Downloads\shi_250115\test\preBotZ_1200_400_z50-700_15-55lp_00001.tif";
info = imfinfo(fn);

S = string({info.Software});
S = replace(S, "↵", newline);

% Quick check: are all headers identical?
sameAll = all(S == S(1));
fprintf("All pages have identical Software header: %d\n", sameAll);

% If not identical, show which pages differ
if ~sameAll
    idx = find(S ~= S(1));
    fprintf("First differing pages (up to 10): %s\n", mat2str(idx(1:min(10,end))));
end


%%
info = imfinfo(fn);
hdr  = string(info(1).Software);
hdr  = replace(hdr,"↵",newline);

m = regexp(hdr, "SI\.hScan2D\.mask\s*=\s*\[(.*?)\]", "tokens", "once");
body = m{1};
body = strrep(body, ";", " ");
body = strrep(body, ",", " ");
body = strrep(body, newline, " ");

mask = sscanf(body, "%f");

fprintf("mask length = %d\n", numel(mask));


%%
info = imfinfo(fn);
hdr = replace(string(info(1).Software),"↵",newline);

% parse mask
m = regexp(hdr,"SI\.hScan2D\.mask\s*=\s*\[(.*?)\]","tokens","once");
mask = sscanf(strrep(m{1},";"," "),"%f");
mask = mask(:)';

% compute mean projection quickly from a subset
H = info(1).Height; W = info(1).Width;
idx = round(linspace(1,numel(info),200));  % 200 frames sample
acc = zeros(H,W,'double');
for k = idx
    acc = acc + double(imread(fn,k,"Info",info));
end
proj = acc/numel(idx);

col = mean(proj,1);
col = col/mean(col);
maskn = mask/mean(mask);

figure; plot(col,'k'); hold on;
plot(maskn,'r'); plot(1./maskn,'b');
legend('raw column profile','mask (normalized)','1/mask (normalized)');
title('Which one matches raw?');
