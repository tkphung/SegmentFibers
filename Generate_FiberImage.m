% Generate_FiberImage.m
% This script generates a synthetic image of fibers for testing a fiber
% segmentation algorithm.
% The images are meant to mimic epithelial basal cell stress fibers.
%
% Thien-Khoi N. Phung (December 17, 2020)

%% Image size parameters
% Set the image size
x = 1920; % width
y = 1216; % height

% Image resolution
pixres = 0.293; % um/pixel (side length of pixel)

% Create local regions of alignment?
localalignment_flag = false;

%% Fiber parameters
% Number of fibers
nfibs = 20000;

% Fiber width (aka fiber diameter)
fdm  = 0.6; % mean (um)
fdsd = 0.2; % std (um)
    % Metrics for log normal distribution
    fdmu = log((fdm^2)/sqrt(fdsd^2 + fdm^2));
    fdsigma = sqrt(log((fdsd^2)/(fdm^2)+1));
    
% Fiber lenghths
flm  = 6; % mean (um)
flsd = 3; % std (um)
    % Metrics for log normal distribution
    flmu = log((flm^2)/sqrt(flsd^2 + flm^2));
    flsigma = sqrt(log((flsd^2)/(flm^2)+1));

%% Fiber orientation distribution
if localalignment_flag
    % Create interpolation points to define fiber angles
    % Number of points
    nangs = 100;
    
    % Coordinates of points
    xangs = rand(nangs,1).*x;
    yangs = rand(nangs,1).*y;
    % Angles at points
    angs = rand(nangs,1).*180; % (degrees)
    
    % Use STDev to determine heterogeneity
    alignstd = 10; % degrees
    angsdev = ones(nangs,1).*alignstd; % high alignemnt 
    
    % Make half the points have random alignment by increased STDev
    randstd =  45; % degrees
    angsdev(round(nangs/2):end) = randstd; % random alignment

    % Create an interpolation functions for angles and stdev
    angint = scatteredInterpolant(xangs,yangs,angs,'linear','linear');
    stdint = scatteredInterpolant(xangs,yangs,angsdev,'linear','linear');
end


%% Generate image
% Position of fiber (centroid)
fx = rand(nfibs,1)*x;
fy = rand(nfibs,1)*y;

% Fiber length (pixels)
fl = lognrnd(flmu,flsigma,nfibs,1)/pixres;

% Fiber width (pixels)
fd = lognrnd(fdmu,fdsigma,nfibs,1)/pixres;

% Fiber angle distribution metrics
if localalignment_flag
    famean = angint(fx,fy);
    fastd  = stdint(fx,fy);
    % Cap the std between    
    fastd(fastd<alignstd) = alignstd;
    fastd(fastd>randstd)  = randstd;
end

fa = [];
fimg = zeros(y,x);
for ff = 1:nfibs
    disp(['Generating Fiber ' num2str(ff)])
    
    % Fiber orientation
    if localalignment_flag
        % Draw random fibers based on Gaussian distribution
        fa(ff) = normrnd(famean(ff),fastd(ff));
    else
        % Uniform distribution of angles
        fa(ff) = rand(1)*180;
    end
    
    % Create fiber skeleton
    se    = strel('line',fl(ff),fa(ff));
    fibse = double(se.Neighborhood);
    fibse = padarray(fibse,[5 5]);

    % Generate fiber width
    fibse = imdilate(fibse,strel('disk',round(fd(ff)/2)));

    % Add fiber to image
    fibtemp = zeros(y,x);
    fibtemp(ceil(fy(ff)),ceil(fx(ff))) = 1;
    fibtemp = conv2(fibtemp,fibse,'same');
    fimg    = fimg+fibtemp;   
end

%% Image post-processing
% Visualize Image
figure('WindowStyle','docked','NumberTitle','off','name','Gen. Fibers')
imagesc(fimg)
axis equal tight

% Binary fiber image
binfiber = fimg>0;
figure('WindowStyle','docked','NumberTitle','off','name','Binary')
imshow(binfiber,[])
axis equal tight

% Increase contrast between background and fiber
confiber   = fimg + (fimg>0).*5;
% Convert to grayscale
grayfiber  = mat2gray(confiber,[0 max(confiber(:))]);
% Blur fiber image
gaussfiber = imgaussfilt(grayfiber);
figure('WindowStyle','docked','NumberTitle','off','name','Gauss Blur')
imshow(gaussfiber,[])
axis equal tight

%% Save image and information
% Create a random name for the fiber data
alpha = 'A':'Z';
imtag = alpha(randi(numel(alpha),1,5)); 

% Save binary image
imwrite(binfiber,[imtag '_binary.tif'],'tif')

% Save blurred image
imwrite(gaussfiber,[imtag '_gauss.tif'],'tif')

% Save information about fibers
fiberinfo.n         = nfibs; % number of fibers
fiberinfo.lengths   = fl;    % lengths (pixels)
fiberinfo.diameters = fd;    % diameters (pixels)
fiberinfo.angles    = fa;    % angles (degrees)
fiberinfo.rawimg    = fimg;  % raw fiber image data
fiberinfo.tag       = imtag; % filenamecode
save([imtag '_fiberinfo.mat'],'fiberinfo')