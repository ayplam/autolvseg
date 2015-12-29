% Simple LV Localization on SSFP CINE images
% In the case that the cine image is not the highest intensity, you may
% have to choose a different kmeans group.
k = 3;
totImgMask = kmeans2(sum(cine,3),k);
stdImg = std(cine,[],3);

% The blood pool is the brightest and will be the highest k-group

[B,L] = bwboundaries(totImgMask == k);
L2 = (totImgMask == k) .* stdImg;
Lstd = zeros(length(B),1);

for i = 1:length(B)    
    blob = L == i;
    Lstd(i) = sum(sum(blob .* L2));
end


figure;
set(gcf,'Position',[450 450 1010 460]);
[~,ind] = max(Lstd);
subplot(1,3,1)
imagesc(L2)
title('Blobs highlighting standard deviations across time')
axis equal
subplot(1,3,2)
plot(Lstd)
title('Plot of sum of std blobs')
subplot(1,3,3)
imagesc(L == ind);
title('Actual LV location')
axis equal

LVpos = (L==ind);

%% After you get the LV do some basic LV localization across the rest of 
% the CINE image
k = 3
sz = size(cine);
LVlocal = zeros(size(cine));
se = strel('disk',2)
for i = 1:sz(3)
    tmpK = kmeans2(cine(:,:,i),k);
    [~,L] = bwboundaries(tmpK == k);
    inds = find(LVpos);
    bdr = L == mode(L(inds));
    
    % Erode to ensure you remain inside the LV cavity.
%     [B,L] = bwboundaries(imerode(bdr,se));
    [B,L] = bwboundaries(bdr);
    dt = delaunayTriangulation(B{1}(:,1),B{1}(:,2));
    ch = convexHull(dt);
    LVlocal(:,:,i) = roipoly(LVlocal(:,:,1), ...
        round(dt.Points(ch,2)),round(dt.Points(ch,1)));
end
    
    


%% Possible to do some simple localization using 
cine = dcmread('Z:\Images\NeoStem AMR\AMR-002_132019 - Sorted\001-002\Screening - 1.2.840.113970.3.57.1.25749326.20120119.1140621\DCM9021_TRUFI_CINE_SA_stack');
test = zeros(size(cine));
for i = 1:size(cine,3)
    test(:,:,i) = kmeans2(cine(:,:,i),2);
end