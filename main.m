%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Automated LV segmentation based on Lynch & Whelan 2004
% -This script is meant to serve as a base for endocardial segmentation of
% short-axis cardiac images to allow calculation of ejection fraction.
% There are two main steps to the segmentation: (1) Smoothing and (2)
% Clustering the pixel intensities.
%
% Written by: Adrian Lam
% ayplam@gmail.com
%
%



%% Read in dicom images
close('all')
dcmdir = 'DCM0028_trufi_cine_shortaxis';
files = dir(dcmdir);
files(1:2) = [];
fullfilenames = cellfun(@(x)fullfile(dcmdir,x),{files.name},'uniformoutput',0);
tmp = dicomread(fullfilenames{1});
cine = zeros([size(tmp) length(fullfilenames)]);
for i = 1:length(fullfilenames)
    cine(:,:,i) = dicomread(fullfilenames{i});
end
sz = size(cine);

% Images to display as an example
d = round(linspace(1,30,3));

% Automatically identify the location of the LV
LVlocal = cineLVLocalize(cine);

opengl software
f1 = figure(1);
for i = 1:length(d)
    subplot(1,3,i)
    imagesc(cine(:,:,d(i)));
    axis image
    hold on
    tmp = zeros(sz(1),sz(2),3);
    tmp(:,:,1) = LVlocal(:,:,d(i));
    overlay = imagesc(tmp);
    set(overlay,'AlphaData',LVlocal(:,:,d(i)) * 0.7);
end
colormap('gray')
set(gcf,'Position', [560 609 954 339]);

hb = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off', ...
    'Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.975,'\bf Original DICOM Images with Simple LV Overlay', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment', 'top', ...
    'FontSize',13);

% Step 1
% Perform smoothing on the cine images for improved image clustering
cineSmoothed = adaptiveSmoothing(cine);

% Step 2
% Use a k-means to start clustering like-segments together
kInfo = getKClusters(cineSmoothed,22);

opengl software
f2 = figure(2);
for i = 1:length(d)
    subplot(2,3,i)
    imagesc(cineSmoothed(:,:,d(i)));
    subplot(2,3,i+3)
    imagesc(kInfo.mask(:,:,d(i)));
end

colormap('jet')
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off', ...
    'Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.975,'Smoothed Images', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment', 'top', ...
    'FontSize',12);
text(0.5, 0.525,'KMeansClustered Images', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment', 'top', ...
    'FontSize',12);

% Step 2b
% Combine all related clusters
combinedClusters = kMeansClusterCombine(kInfo);

% Finalize the cluster combining by extracting which cluster is actually 
% the LV and doing some convex rounding
LVseg = finalClusterCombine(combinedClusters,LVlocal);

opengl software
f3 = figure(3);
im = imagesc(cine(:,:,1));

hold on
axis image

tmp_ov = zeros(sz(1),sz(2),3);
ov = imagesc(tmp_ov);

% Review all segmented images
for i = 1:sz(3)
    
    set(im,'CData',cine(:,:,i));
    tmp_ov(:,:,1) = LVseg(:,:,i);
    set(ov,'AlphaData',LVseg(:,:,d(i)) * 0.7,'CData',tmp_ov);    
    
    colormap('gray')
    title(['Segmented Short-Axis Image: #',numstr(i)])
    
    drawnow;
end
