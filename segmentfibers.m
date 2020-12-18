function [fiberpx,fiberd,fiberlab] = segmentfibers(IMG,params,progress_flag)
%segmentfibers: this script takes in an image containing fibers and
%segments the fibers.
% [fiberpx,fiberd,fiberlab] = segmentfibers(IMG,params)
% INPUTS:
%       IMG- image file (2D matrix)
%       params- struct containing:
%             .pixres- pixel resolution (um/pixel)
%             .fdiam- expected min and max fiber diameter (um)
%       progress_flag- logical of whether to display progress of
%                      segmentation
% OUTPUTS:
%       fiberpx- cell array containing coordinates for all segmented fibers
%       fiberd- vector containing lengths of all fibers (in pixels)
%       fiberlab- 2D matrix image of the segmented fibers
%
% The fiber segmentation algorithm is based on two papers:
% 1) We implement a fiber tracing method developed by
%    Rogge, H., Artelt, N., Endlich, N. & Endlich, K. Automated segmentation
%    and quantification of actin stress fibres undergoing experimentally
%    induced changes. J. Microsc. 268, 129-140 (2017).
% 2) We reconstruct the fiber pieces using a pairing algorithm developed by
%    Zhang, Z., Xia, S. & Kanchanawong, P. An integrated enhancement and
%    reconstruction strategy for the quantitative extraction of actin stress
%    fibers from fluorescence micrographs. BMC Bioinformatics (2017).
%
% Thien-Khoi N. Phung (December 17, 2020)

control_flag = false; % TODO add in control condition flag

%% STEP 1. Parameters
    % Pixel resolution
    pixres = params.pixres; % um/pixel (side length)

    % Expected fiber width range
    Wmin = round(params.fdiam(1)/pixres); % pixels
    Wmax = round(params.fdiam(2)/pixres); % pixels
    
%% STEP 2. Illumination correction
    % Tophat filter for illumination correction
    % Define size of structuring element based on pixel resolution
    tophatSE  = strel('disk',round(7/pixres));
    OriginImg = imtophat(im2uint8(IMG),tophatSE);

    % Adjust image intensity (saturate 1% of low and high intensities)
    % This is essentially Global Thresholding to define "background"
    OriginImg = imadjust(OriginImg);

%% Step 3. Image pre-processing
    % Gauss Filter (gaussimg)
    gaussimg = imgaussfilt(OriginImg,0.75*Wmin);


    % Path Opening (poimg)
    % Create a set of line structuring elements (SE) with length (lpo) and 
    % angles distributed over 180 degrees
    lpo   = 2*(Wmin + (Wmax-Wmin)/2);
    abins = 8; % number of angles
    angs  = linspace(0,180,abins+1); % degrees
    
    % For each SE, erode & dilate the image
    poimgs = zeros([size(gaussimg) abins]);
    for by = 1:abins
        lse = strel('line',lpo,angs(by)); % line SE
        % Erosion
        eroimg = imerode(gaussimg,lse);
        % Dilation
        dilimg = imdilate(eroimg,lse);
        % Store PO image
        poimgs(:,:,by) = dilimg;
    end    
    % Take the maximum intensity projection of the path opened images
    poimg = max(poimgs,[],3);

    
    % Directional Top-Hat (thimg)
    % Using same set of line SE angles as Path Opening, but change length
    lth = Wmax;
    % For each SE, erode & dilate the image
    thimgs = zeros(size(poimgs));
    for by = 1:abins
        lse = strel('line',lth,angs(by)); % line SE
        % Erosion
        eroimg = imerode(poimg,lse);
        % Dilation
        dilimg = imdilate(eroimg,lse);
        % Store Directional Top-Hat image
        thimgs(:,:,by) = poimg - dilimg;
    end
    % Sum together the top-hat images
    thimg = sum(thimgs,3);

    
    % Hysteresis Threshold (create mask for fibers) (hthimg)
    % Object size minimum
    Amin = 3*lth;
    
    % Determine range of non-zero pixels pixel intensities to set threshold
    lowpi  = min(thimg(thimg>0));
    highpi = prctile(thimg(thimg>0),90); % nth prctile
    range  = lowpi:highpi;
    
    % Use a nested for loop with 'break' to search thresholds
    hthcount = nan(numel(lowpi:highpi));
    flag_break = 0;
    for hh = lowpi:highpi         % set high threshold
        if flag_break
            break
        end
        
        % Above high threshold
        [rr,cc] = find(thimg>hh);
        
        for ll = (hh-1):-1:lowpi  % set low threshold
            % Above low threshold
            lowimg  = thimg>ll;
            
            % Obtain all connected regions in lowimg that includes points
            % from [rr,cc]
            hthimg = bwselect(lowimg,cc,rr,8);

            % Count objects larger than Amin
            hthcount(hh,ll) = sum(cell2mat(struct2cell(...
                                  regionprops(hthimg,'Area')))>Amin);
            
            % If ll==(hh-1) the first low threshold for this new high thr.
            %    AND the count is less than the previous maximum
            % Stop the loops
            if ll==(hh-1) && hthcount(hh,ll)<max(hthcount(:)) && max(hthcount(:))>100
                flag_break = 1;
                break
            end
            
            % If count is 0
            %    OR lower threshold changes number of count objects 
            % move to the next high threshold
            if hthcount(hh,ll)==0 || hthcount(hh,ll)<hthcount(hh,hh-1)
                break
            end
        end
    end
    % Set high threshold based to find a maximum number of objects of size
    % larger than Amin
    % Decrease low threshold as long as number of connected objects is
    % unchanged (to prevent merging of objects)
    [hh,ll]    = find(hthcount==max(hthcount(:)));
    lowthresh  = range(ll);
    highthresh = range(hh);
    
    % Perform final Hysteresis Threshold
    % Above low threshold
    lowimg  = thimg>lowthresh;

    % Above high threshold
    [rr,cc] = find(thimg>highthresh);

    % Obtain all connected regions in lowimg that includes points in
    % [rr,cc]
    hthimg = bwselect(lowimg,cc,rr,8);

    % Remove any objects smaller than Amin
    hthimg = bwareaopen(hthimg,Amin);

% For CONTROL cond, remove any regions that are not within cell boundaries
    if control_flag
        % Only search for fibers in the Basal Cell bodies
        hthimg = hthimg & BCBMASK;
    end
    
    
%% Step 4. Trace algorithm
    % Sobel filter to calculate grayscale gradients
    [Gx,Gy] = imgradientxy(OriginImg,'Sobel');

    % Calculate Orientation field    
        % Components
        gsx = Gx.^2 - Gy.^2;
        gsy = 2.*Gx.*Gy;
        % Smooth
        bo  = ones(Wmax); % average over block the same size as Wmax
        gsx = conv2(gsx,bo,'same');
        gsy = conv2(gsy,bo,'same');
        ofield = 1/2.*atan2(gsy,gsx) + (pi/2);

    % Angular coherence Map
        % Window size 
        bco = round(3*Wmax); 

        % Pad ofield
        padofield = nan(size(ofield)+bco-1);
        padofield((1:size(ofield,1))+floor(bco/2),...
                  (1:size(ofield,2))+floor(bco/2)) = ofield;

        % Performed the equation:
        %       abs(cos(ofield(x0,y0) - ofield(x,y))) for x,y in bco square
        ofieldstack = zeros([size(ofield) bco^2]);
        lyr = 0;
        for rr = 1:bco
            for cc = 1:bco
                lyr = lyr+1;
                ofieldstack(:,:,lyr) = abs(cos(ofield -...
                                       padofield((1:size(ofield,1))+(rr-1),...
                                                 (1:size(ofield,2))+(cc-1))));
            end
        end
        angcoh = nansum(ofieldstack,3);

    % Parameter- Trace Step Size
    tracestep = 1.5*(Wmin + (Wmax-Wmin)/2); % pixels
    curvethr  = (22/Wmax)*pi/180; % radians


    % Create an image containing relevant angular coherance belonging to fibers
    % Mask the non-fiber regions of the angcohfib with nan
    angcohfib          = angcoh;
    angcohfib(~hthimg) = nan;

    % Create scatteredInterpolant for angles- Interpolate vector components
    [frow,fcol] = find(hthimg);
    fvecx       = scatteredInterpolant(frow,fcol,cos(2.*ofield(hthimg)),'linear');
    fvecy       = scatteredInterpolant(frow,fcol,sin(2.*ofield(hthimg)),'linear');

    flag_segfiber = true;
    count = 0;

    fmask = zeros(size(angcohfib));
    while flag_segfiber
        % Find max ang coh to intiate tracing (for plotting, y-rows, x-columns)
        [ro,co] = find(angcohfib==max(angcohfib(:)),1);
        
        if progress_flag
            disp(['Traced Count ' num2str(count) '; Checking (' ...
              num2str(co) ',' num2str(ro) ')-- Search Space Remaining '...
              num2str(sum(~isnan(angcohfib(:)))) ' Pixels'])
        end
        
        % Isolate the fiber mask for segmentation
        fibermask = bwselect(~isnan(angcohfib),co,ro);
            % Check if fibermask is small
            if sum(fibermask(:))<Wmax && sum(~isnan(angcohfib(:)))>0
                % Remove from angcohfib
                 angcohfib(fibermask) = nan;

                % Don't trace, continue to next iteration
                continue
            elseif sum(~isnan(angcohfib(:)))==0
                flag_segfiber = false;
                continue
            end
        fxlog = [];
        fylog = [];
        falog = [];
        flag_trace   = true; % Are we still tracing the same fiber?
        flag_forward = true; % Have we finished checking one direction of fiber?
        % flag_curvy   = false;
        while flag_trace
        % Find borders and correct point
            % Find angle at the point
            ao = 1/2*atan2(fvecy(ro,co),fvecx(ro,co)); % radians
                % Flip angle to be within -90 to +90 deg of the prior angle
                % This ensures the sides of the borders do not swap
                if ~isempty(falog) && abs(ao-falog(end))>(pi/2)
                    if ao<0 % ao in quadrants I or II
                        ao = ao + pi;
                    else    % ao in quadrants III or IV
                        ao = ao - pi;
                    end
                    % ao should still be within -180 to 180 deg
                end
            % Calculate orthogonal angle
            aortho = ao + pi/2;

            % Grab profile of edgesimg to find border
            % Border 1
            [b1x,b1y,b1i] = improfile(fibermask,[co co+tracestep*cos(aortho)],...
                                                [ro ro+tracestep*sin(aortho)]);
            b1i(isnan(b1i)) = 0;
            b1id            = find([b1i;0]==0,1)-1;
            % Border 2
            [b2x,b2y,b2i] = improfile(fibermask,[co co-tracestep*cos(aortho)],...
                                                [ro ro-tracestep*sin(aortho)]);
            b2i(isnan(b2i)) = 0;
            b2id            = find([b2i;0]==0,1)-1;
            % [Border1 CorrectPoint Border2]
            co = round((b1x(b1id)+b2x(b2id))/2);
            ro = round((b1y(b1id)+b2y(b2id))/2);
            fx = [b1x(b1id)   co   b2x(b2id)];
            fy = [b1y(b1id)   ro   b2y(b2id)];


        % Corrected fiber angle
        ao = 1/2*atan2(fvecy(ro,co),fvecx(ro,co)); % radians
            % Flip angle to be within -90 to +90 deg of the prior angle
            % This ensures the sides of the borders do not swap
            if ~isempty(falog)
                if abs(ao-lastao)>(pi/2)
                    if ao<0 % ao in quadrants I or II
                        ao = ao + pi;
                    else    % ao in quadrants III or IV
                        ao = ao - pi;
                    end
                    % ao should still be within -180 to 180 deg
                end
            end

        % Record data and Step along fiber direction 
        if flag_forward
            % Record fiber information
            fxlog = [fxlog; fx];
            fylog = [fylog; fy];
            falog = [falog; ao];
            lastao = ao;
            % Step Forward
            co = round(co - tracestep*cos(ao));
            ro = round(ro - tracestep*sin(ao));
        else % Reverse direction
            % Record fiber information
            fxlog = [fx; fxlog];
            fylog = [fy; fylog];
            falog = [ao; falog];
            lastao = ao;
            co = round(co + tracestep*cos(ao));
            ro = round(ro + tracestep*sin(ao));
        end

        % Check if predicted point is in the fiber body?
        [~,~,io] = improfile(fibermask,co,ro); % sample the point
            if isnan(io)
                io = 0;
            end
        % Check if predicted point is unique (not already in fiber)
            if ismember([co ro],[fxlog(:,2) fylog(:,2)],'rows')
                io = 0; % Point is already traced
            end
        if ~logical(io) % NOT in body
            if flag_forward % Forward
                % Reverse direction
                flag_forward = false;

                % Reset to the other side of the fiber
                co = round(fxlog(1,2) + tracestep*cos(ao));
                ro = round(fylog(1,2) + tracestep*sin(ao));

                % Check if this point is in the fiber body?
                [~,~,io] = improfile(fibermask,co,ro); % sample the point
                    if isnan(io)
                        io = 0;
                    end
                if ~logical(io) % Reverse point is not in body
                    flag_trace = false; % Stop tracing
                end
            else % Already reverse
                flag_trace = false; % Stop tracing
            end
        end
        end

        % Create mask to remove this segmentation iteration
            if size(fxlog,1)>1 % Fiber was segmented
                % Create Mask for segmented fiber
                byefiber = poly2mask([fxlog(:,1); flipud(fxlog(:,3))],...
                                 [fylog(:,1); flipud(fylog(:,3))],...
                                 size(angcohfib,1),size(angcohfib,2));
                % If the mask is empty (the fiber and borders all coincide
                if mean(byefiber(:))==0
                    byefiber(round(fylog(:)),round(fxlog(:))) = 1;
                end

                % Count fiber and store information
                count = count+1;
                fibers{count}.x = fxlog;
                fibers{count}.y = fylog;
                fix = fxlog(:,2);
                fiy = fylog(:,2);
                % If the pixels are out of frame, round to nearest frame pixel
                fiy(fiy<1) = 1;
                fiy(fiy>size(fmask,1)) = size(fmask,1);
                fix(fix<1) = 1;
                fix(fix>size(fmask,2)) = size(fmask,2);
                % Fill in pixels to form a continuous fiber
                [fillx,filly,~] = improfile(fmask,fix,fiy);
                fmask(sub2ind(size(fmask),round(filly),round(fillx))) = 1;

            else % Fiber was not segmented, remove starting point
                % Create mask of the searched point
                byefiber = zeros(size(angcohfib));
                [roo,coo] = find(angcohfib==max(angcohfib(:)),1);
                byefiber(roo,coo) = 1;
            end
        % Remove this search region from angcohfib
            byefiber = imdilate(byefiber,strel('disk',round((Wmin+Wmax)/2)));
            angcohfib(logical(byefiber)) = nan;

        % Clean up search space- remove small objects
            keeparea = bwareaopen(~isnan(angcohfib),Amin);
            angcohfib(~keeparea) = nan;
    end % flag_segfiber
    
    
%% Step 5. Calculate fiber fragment endpoint orientations
    % Remove single pixel fiber fragments
    linfrags = bwmorph(fmask,'clean');

    % Connection parameters
        % Set search distance for fiber tips
        dmax   = 3*Wmax; % pixels
        % Maximum allowed fiber kinking
        kmax   = 25/(3*Wmax)*pi/180; % radians

    % Find endpoints of line fragments
    linends       = bwmorph(linfrags,'endpoints');
    [lerow,lecol] = find(linends); % row, col locations of endpoints

    % Associate endpoints with fragments
    linlabel  = bwlabel(linfrags);
    lelab     = linlabel(sub2ind(size(linlabel),lerow,lecol)); % frag ID for each endpoint

    % Pixel locations for all fragments
    pxlist    = struct2cell(regionprops(linfrags,'pixellist'));
    % Pad all pixel lists to be the same size
    padnum    = max(cellfun(@(x) size(x,1),pxlist));
    pixellist = cellfun(@(x) padarray(x,padnum-size(x,1),inf,'post'),...
                        pxlist,'UniformOutput',false);
    % Reshape pixellist into px and py (pixel x y coordinate matrices)
    % each row is a fragment (the row index corresponds to lelab)
    % each col lists a coordinate for a pixel of the frag
    pixelxy   = cell2mat(pixellist');
    px        = reshape(pixelxy(:,1),padnum,max(linlabel(:)))';
    py        = reshape(pixelxy(:,2),padnum,max(linlabel(:)))';

    % For each end point, grab its corresponding fragment coordinates
    % each row here is an end point
    % each col lists the coordinates for its associated fragment
    lepx      = px(lelab,:);
    lepy      = py(lelab,:);
    % Calculate distance from endpoints to its fragment pixels
    lepdist   = ((lepx - repmat(lecol,1,padnum)).^2 + ...
                 (lepy - repmat(lerow,1,padnum)).^2).^(1/2);

    % Filter fragment pixels below distance threshold
    keepp        = lepdist<dmax;  % Distance threshold
    lepx(~keepp) = nan;           % Assign nan to all other pixels
    lepy(~keepp) = nan;

    % Calculate centroid of remaining fragment pixels for each endpoint
    centx     = nanmean(lepx,2);
    centy     = nanmean(lepy,2);

    % Calculate vector between centroid -> endpoint
    levecx    = lecol - centx;
    levecy    = lerow - centy; % note +Y is downward in axis IJ

    % Calculate orientation angle using full -180 to 180 degrees
    % Use -levecy to define +angle in counterclockwise direction in IJ coord
    leang    = atan(-levecy./levecx).*(180/pi); % [-90 to 90 deg]
    % Correct angles up to +/-180 degrees
    leang(levecx<0 & (-levecy)>=0) = leang(levecx<0 & (-levecy)>=0) + 180;
    leang(levecx<0 & (-levecy)<0)  = leang(levecx<0 & (-levecy)<0) - 180;


%% STEP 6. Pair tips (aka endpoints)
    % Keep track of hmidx & nbridx
    % These are the indices for the potential home-neighbor tip pairings

    % PROXIMITY: Filter potential hm-nbr tip pairs by distance
        threshhndist    = dmax; % radius from home tip to search (pixels)

        % Calculate all distances between home and neighbors
        % each row refers to the home tip
        % each col refers to the nbr tip
        % This returns an upper triangular matrix of distances; this is to avoid
        % pairing a-b with b-a; we only look at nbr tips with and index greater
        % than the home tip.
        hndist            = triu(squareform(pdist([lecol,lerow])),1);
        hndist(hndist==0) = nan; % set lower triangle to nan

        % Keep only nbr tips with distance <= threshhndist
        hndist(hndist>threshhndist) = nan;

        % This matrix should be very sparse- convert to (home, nbr) pairs
        [hmidx,nbridx] = find(~isnan(hndist));

        % Remove any hm-nbr pairs from the same fragment (two ends of the same
        % fragment)
        samefrag = lelab(hmidx) == lelab(nbridx);
        hmidx(samefrag)  = [];
        nbridx(samefrag) = [];

        % PROXIMITY metric- distance between hm and nbr
        % This was already used as the first filter step; however, we calculate
        % update the list here after FAN SEARCH.
        proxdist = hndist(sub2ind(size(hndist),hmidx,nbridx));


    % CONTINUITY: Filter potential hm-nbr tip pairs by fan angle from hm tip
        threshfanangle  = 20; % angle around home tip orientation to search (deg)

        % Calculate vector between home -> nbr tips
        hnvecx = lecol(nbridx) - lecol(hmidx);
        hnvecy = lerow(nbridx) - lerow(hmidx); % note +Y is downward in axis IJ

        % Calculate angle of home-neighbor vector using full -180 to 180 degrees
        % Use -hnvecy to define +angle in counterclockwise direction in IJ coord
        hnang  = atan(-hnvecy./hnvecx).*(180/pi); % [-90 to 90 deg]
        % Correct angles up to +/-180 degrees
        hnang(hnvecx<0 & (-hnvecy)>=0) = hnang(hnvecx<0 & hnvecy>=0) + 180;
        hnang(hnvecx<0 & (-hnvecy)<0)  = hnang(hnvecx<0 & hnvecy<0) - 180;

        % CONTINUITY metric- angle difference between hm angle and hm->nbr angle
        contang = abs(hnang - leang(hmidx));
        % If angle > 180, take 360 - angle
        contang(contang>180) = 360 - contang(contang>180);

        % Remove hm-nbr pairs with angles outside fan angle window
        outfan = contang>(threshfanangle/2);
        % Update all hn vectors
        hmidx(outfan)    = [];
        nbridx(outfan)   = [];
        contang(outfan)  = [];
        proxdist(outfan) = []; % Update proximity metric

    % SIMILARITY metric- (supplementary) angle difference between hm angle and nbr angle
        % leang- hm angles for all tips
        % Difference between hm and nbr tip angle
        simang = abs(leang(nbridx) - leang(hmidx));
        % If angle > 180, take 360 - angle
        simang(simang>180) = 360 - simang(simang>180);
        % Angles closer to 180 deg mean more similar
        % Flip metric by taking supplementary angle
        % simang closer to 0 means more similar
        simang = 180-simang;

        % Filter by SIMILARITY metric
        % This is making sure the hm and nbr tips are pointing (somewhat)
        % towards each other
        % threshold for SIMILARITY angle
        threshsimang = 45; % degrees (aka 2*threshsimang deg fan region)

        % Remove hn pairs based on threshold
        outsim = simang>threshsimang;
        % Update all hn vectors
        hmidx(outsim)    = []; % home indices
        nbridx(outsim)   = []; % nbr indices
        proxdist(outsim) = []; % proximity metric
        contang(outsim)  = []; % continuity metric
        simang(outsim)   = []; % similarity metric


    % PEAR tips based on PROXIMITY, CONTINUITY, and SIMILARITY
        % So far the hn lists only contain h-n pairs where hmidx<nbridx
        % Now we need to consider the reverse n-h pairs as well when pairing
        % Here we assume no fibers will branch
        % Duplicate the hn vectors to contain the nh combination as well
        hnnhidx  = [hmidx nbridx; nbridx hmidx];
        proxdist = [proxdist; proxdist];
        contang  = [contang; contang];
        simang   = [simang; simang];

        % Sort hnnhidx
        [~,hsort] = sort(hnnhidx(:,1));
        hnnhidx   = hnnhidx(hsort,:);
        proxdist  = proxdist(hsort);
        contang   = contang(hsort);
        simang    = simang(hsort);

        % Count number of tips per home
        % Unique hm indices
        [uhmidx,~,ic] = unique(hnnhidx(:,1));
        % Count how many times each unique home index occurs
        hmcount = accumarray(ic,1);

        % Divide hn vectors into cells
        % Each cell is the nbr metrics for a single hm listed in uhmidx
        gnbridx   = mat2cell(hnnhidx(:,2),hmcount);
        gproxdist = mat2cell(proxdist,hmcount);
        gcontang  = mat2cell(contang,hmcount);
        gsimang   = mat2cell(simang,hmcount);

        % Calculate metric score
        % Weights for scores
        proxwt = 1;
        contwt = 1;
        simwt  = 1;
        % Normalize PROXIMITY, CONTINUITY, and SIMILARITY 
        % *If angle max is 0, add 1 to avoid divide by zero error.
        ngproxdist = cellfun(@(x) x./max(x),gproxdist,'UniformOutput',false);
        ngcontang  = cellfun(@(x) x./(max(x) + 1*(max(x)==0)),gcontang,'UniformOutput',false);
        ngsimang   = cellfun(@(x) x./(max(x) + 1*(max(x)==0)),gsimang,'UniformOutput',false);

        % Add scores
        pcsscore  = cellfun(@(p,c,s) proxwt.*p + contwt.*c + simwt.*s,...
                              ngproxdist,ngcontang,ngsimang,'UniformOutput',false);

        % Find min scoring pair for each hm tip
        %   In case of a tie in scores- just take the first one
        minscoreid = cellfun(@(score) find(score==min(score),1),pcsscore,'UniformOutput',false);
        % Pull out the information about hte min score tip
        npear = cellfun(@(minid,nbr,p,c,s) ...
                        [nbr(minid)  p(minid)  c(minid)  s(minid)],...
                        minscoreid,gnbridx,gproxdist,gcontang,gsimang,'UniformOutput',false);

        % Combine home and neighbor pear information
        % This is the most likely pair for each tip
        % [HtipID  NtipID Proximity Continuity Similarity]
        hnnhpear = [uhmidx cell2mat(npear)];

        % Assume that each tip can only be paired with one other tip
        % To confirm they are a good match, the combination A-B and B-A should
        % both be present in hnnhpear(:,1:2)
        % Sort the columns so that lowest index is the home
        hnnhpear(:,1:2) = sort(hnnhpear(:,1:2),2);

        % Unique pears
        [upear,~,ic] = unique(hnnhpear,'rows');
        % Count how many times each pear is in hnnhpear
        upearcount   = accumarray(ic,1);

        % Keep pears with upearcount==2
        upear(upearcount~=2,:) = []; % Remove all other upear
        
%% Step 7. Check bridge for each pair
    % Determine if bridge between each paired tips has a potential fiber by
    % looking at the pixel intensities
    % Take in upear, and remove any pairs that do not contain a valid bridge

    % Create an image with pixel intensities
    % Determine the pixel [Min Mean Max] intensity of the linfrags
    fragi = [cell2mat(struct2cell(regionprops(linfrags,OriginImg,'MinIntensity')))' ...
             cell2mat(struct2cell(regionprops(linfrags,OriginImg,'MeanIntensity')))'... 
             cell2mat(struct2cell(regionprops(linfrags,OriginImg,'MaxIntensity')))'];

    % Determine the mean pixel intensity between endpoints
        % X Y coordinates for each tip pair
        xx = lecol(upear(:,1:2)); % Each row is a pair, each col is coord
        yy = lerow(upear(:,1:2));

        bridgei = zeros(size(xx,1),3);
        for pairjz = 1:size(xx,1)
            % Sample intensities across the bridge
            [~,~,ii] = improfile(OriginImg,xx(pairjz,:),yy(pairjz,:));

            % Calculate [Min Mean Max] intensity
            bridgei(pairjz,:) = [min(ii) mean(ii) max(ii)];
        end

    % Check bridge pixel intensity compared to fragments
        % Compile each bridge's respective fragment pixel intensity information
        % Each row is a bridge
        pairedfragi = [fragi(lelab(upear(:,1)),:) fragi(lelab(upear(:,2)),:)];

        % Is mean bridge pixel intensity brighter than the min PI of the two
        % fiber fragments?
        yesbridge = min(pairedfragi,[],2)<bridgei(:,2);

    % Remove dim bridges
    upear(~yesbridge,:) = [];
    
%% STEP 8. Connect fragments based on tip pairing
    % LINE FRAGMENT VARIABLES
    %   linfrags- original fiber fragments mask
    %   pxlist-   each cell contains the pixels (X,Y) for a fragment
    % TIP PAIR VARIABLES
    %   upear-    unique pairs of home-neighbor tips
    %             [HtipID   NtipID   Proximity   Continuity   Similarity]
    %             HtipID and NtipID are the indices used for lerow/col/lab
    %   lerow-    tip row index (y coordinate)
    %   lecol-    tip col index (x coordinate)
    %   lelab-    associated fragment ID for each tip
    %             fragment IDs are the indices for cells of pxlist
    % 
    % Note that connecting tips may cross fibers, so the fibers are stored as
    % pixel coordinates. A mask of the fiber fragments cannot distinguish
    % individual fibers if they cross (they are considered one connected
    % component).

    % Group together the new fragments
        % Replace the Tip ID with the Fragment ID
        % These pairs of fragments are considered a fiber.
        fragpear = [lelab(upear(:,1)) lelab(upear(:,2))];
        % Sort the fragpears so that the lowest ID appears first
        fragpear = sort(fragpear,2);

        % Vector of the Fragment IDs
        fragid   = 1:numel(pxlist);

        % Replace fragid entries to be the lowest Frag ID in the Grouped Fiber
        % By the end of the replacements, the IDs in fragpear(:,2) should no
        % longer be present in fragid vector
        while sum(ismember(fragid,fragpear(:,2)))>0
            % Find the fragid entries that need to be paired
            [replaceme,withwhom] = ismember(fragid,fragpear(:,2));

            % Replace those entries with its lower ID pair
            fragid(replaceme) = fragpear(withwhom(replaceme),1);
        end
        % fragid(jz) now tells you which fragment ID=jz belongs to

    % Generate the pixels for the connecting line between each pair of tips
        % X Y coordinates for each tip pair
        xx = lecol(upear(:,1:2)); % Each row is a pair, each col is coord
        yy = lerow(upear(:,1:2));

        % Distance between points (aka Proximity)
        peardist = upear(:,3);    
        % Use city-block distance to get number of pixels for each bridge
        cityblockdx   = ceil(abs(xx(:,1) - xx(:,2)));
        cityblockdy   = ceil(abs(yy(:,1) - yy(:,2)));
        cityblocks    = max([cityblockdx cityblockdy],[],2);
        % Stepsize between city blocks along the 
        citystep      = peardist./cityblocks;

        % Define the bridging pixels for each tip pair
        % (this excludes the tip pairs)
        % each cell includs the [X Y] coordinate pixels 
        for jz = 1:size(upear,1)
            bridgexy{jz} = round(interp1([0 peardist(jz)],[xx(jz,:)' yy(jz,:)'],...
                       citystep(jz):citystep(jz):(peardist(jz)-citystep(jz))));
        end

        % Identify which fragment group each bridge belongs to
        bridgefragid = fragid(fragpear(:,1));

    % Combine fragments and bridges for each fiber
        % Make a mask for the fibers
        fiberlab = zeros(size(linfrags));

        % Start a new numbering for the fibers
        fibercount = 0;
        fiberpx    = [];
        fiberd     = [];
        % Cycle through each fiber
        for jz = unique(fragid)
            if progress_flag
                disp(['Fiber ' num2str(jz) ' of ' num2str(max(unique(fragid)))])
            end
            
            % Identify all the fragments belonging to this fiber
            groupedfrag = fragid==jz;

            % Identify all the bridges belonging to this fiber
            groupedbridge = bridgefragid==jz;

            % Combine all of the pixels
            fiberpts = [cell2mat(pxlist(groupedfrag)');      % frags
                        cell2mat(bridgexy(groupedbridge)')]; % bridges     

            % Create new mask of all fiber points
            groupmask = zeros(size(linfrags));
            for ptid = 1:size(fiberpts,1)
                groupmask(fiberpts(ptid,1),fiberpts(ptid,2)) = ptid;
            end

            % Find endpoints of fragment
            groupepmask = bwmorph(groupmask>0,'endpoints');
            endpts      = groupmask(groupepmask);
                % If no endpoints- likely a ring instead of fiber
                if numel(endpts)<2
                    continue
                end

            % Connect the pixels
            % Distance between all points
            nbrmat   = squareform(pdist(fiberpts));
            % Make adjacency matrix for pts (numbered by row in fiberpts)
            % Direct neighbor length is 1 or sqrt(2) because unit grid   
            nbredges = ismember(nbrmat,[1 sqrt(2)]);
            nbrmat(~nbredges) = 0;
            % Find shortestpath between the endpoints
            % Note that this may remove some unecessary pixels to create the
            % fiber (this is still a valid fiber with the removed pixels)
            [ptorder,d] = shortestpath(graph(nbrmat),endpts(1),endpts(2));

            % Increase fiber count
            fibercount = fibercount + 1;

            % Store ordered points & distance
            fiberpx{fibercount} = fiberpts(ptorder,:);
            fiberd(fibercount)  = d; % Store distance for optional filtering


            % Add fiber to the mask
            fiberlab(sub2ind(size(fiberlab),fiberpx{fibercount}(:,2),fiberpx{fibercount}(:,1))) = fibercount;
        end

