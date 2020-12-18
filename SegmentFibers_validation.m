% SegmentFibers_validation.m
% This script runs the 'segmentfibers' function on an example image with
% known fiber lengths. We compare the output of the segmentation method
% with the fiber length data.
%
% Thien-Khoi N. Phung (December 17, 2020)

%% Load in data
% Fiber Image
% We use a synthetic stress fiber image created using 
%       Generate_FiberImage.m
% SFIMG  = imread('localalign_example/DUUAD_gauss.tif');
SFIMG  = imread('randomalign_example/SBXIQ_gauss.tif');

% Fiber Data
% From Generate_FiberImage.m, we load in the known fiber geometry
% information for later comparison to segmentation method.
% load('localalign_example/DUUAD_fiberinfo.mat')
load('randomalign_example/SBXIQ_fiberinfo.mat')

%% Setup parameters for segmentation
% Image resolution
params.pixres = 1; % um/pixel (side length)

% Expected fiber width range
params.fdiam = [1 4]; % min and max (um)
  
%% Run segmentation
% This step takes a while...
[fiberpx,fiberd,fiberlab] = segmentfibers(SFIMG,params,true);

%% Display results
figure('WindowStyle','docked','NumberTitle','off','name','Fiber Image')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
imshow(SFIMG,'Parent',axt)

figure('WindowStyle','docked','NumberTitle','off','name','Segmentation')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
imshow(fiberlab,'Parent',axt)

figure('WindowStyle','docked','NumberTitle','off','name','Segmentation Overlay')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
skoverlay = labeloverlay(SFIMG,...
                         imdilate(fiberlab>0,strel('disk',1)),...
                         'Colormap',[1 0 1; 0 0 0],...
                         'Transparency',0);
imshow(skoverlay,[],'Parent',axt)

% imwrite(skoverlay,'localalign_example/DUUAD_segmentfibers.tif','tif')

%% Compare segmentation to data
figure
t = tiledlayout(2,1,'TileSpacing','Compact');

nexttile
histogram(fiberinfo.lengths,linspace(0,80,100))
title('Data- Fiber Lengths')
ylabel('Number of Fibers')
xlabel('Fiber Lengths (pixels)')
xlim([0 80])

nexttile
histogram(fiberd,linspace(0,80,100))
title('Segmentation Results')
ylabel('Number of Fibers')
xlabel('Fiber Lengths (pixels)')
xlim([0 80])