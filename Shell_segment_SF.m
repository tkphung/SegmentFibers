% Shell_segment_SF.m
% This script takes in a *.tif and *.mat file to segments the stress fibers.
%
% Thien-Khoi N. Phung (August 15, 2022)

%% Load in dataset
% Stress Fiber Image
sffilename = 'test_images/U13_72_C2_4_SF.tif';
SFIMG      = imread(sffilename);
SFinfo     = imfinfo(sffilename);
pr         = 1/SFinfo.XResolution;

% Basal Cell Segmentation File
load('test_images/U13_72_C2_4_basalcell.mat'); % basalbody



%% Setup parameters for segmentation
% Image resolution
params.pixres = pr; % um/pixel (side length)

% Expected fiber width range
params.fdiam = [0.5 2]; % min and max (um)

% (OPTIONAL) include basal cell bodies to only search for stress fiber
% within basal cells. Usually only do this for CONTROL samples.
params.basalbody = basalbody;


%% Run segmentation
% This step takes a while...
[fiberpx,fiberd,fiberlab] = segmentfibers(SFIMG,params,true);


%% Display results

figure('WindowStyle','docked','NumberTitle','off','name','Fiber Image')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
imshow(SFIMG,[],'Parent',axt)

figure('WindowStyle','docked','NumberTitle','off','name','Segmentation')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
imshow(fiberlab,'Parent',axt)

figure('WindowStyle','docked','NumberTitle','off','name','Segmentation Overlay')
axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
skoverlay = labeloverlay(imadjust(SFIMG),...
                         imdilate(fiberlab>0,strel('disk',1)),...
                         'Colormap',[1 0 1; 0 0 0],...
                         'Transparency',0);
imshow(skoverlay,[],'Parent',axt)
