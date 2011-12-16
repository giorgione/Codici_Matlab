% function [Tracks, CTracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, detdir, KLT, ext, fstep, framerate, iZ, sparams, firstframe, lastframe)
% Tracking from 1 video Sequence:
%  PARAMETERS
% - imgdir:  directory of images extracted from the video sequence.
% - detdir:  directory of detection results for each frame of the video sequence
% - KLT:     directory containing detections of KLF features for the GROUND
% - ext:     filename NOT USED
% - fstep: 
% - framerate, 
% - iZ : Initial Camera Model
% - sparams:               All Parameters in Model
% - firstframe
% - lastframe 
function [Tracks, CTracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, detdir, KLT, ext, fstep, framerate, iZ, sparams, firstframe, lastframe)
%Detection Threshold
detth = sparams.detth;

Tracks = [];
hTracks = [];
CTracks = [];
hCTracks = [];
camTrack = [];

imfiles = dir([imgdir '*.png']);
if length(imfiles) == 0
    imfiles = dir([imgdir '*.jpg']);
    if length(imfiles) == 0
        imfiles = dir([imgdir '*.bmp']);
    end
end
speed = [];
stang = [];
yaw = [];
latacc = [];
longacc = [];
% 
% addpath('C:\Users\wgchoi\VisionLabWork\FordProject\Dataset\Ford_U of M');
% [imfiles, speed, stang, yaw, latacc, longacc] = readData(imgdir, 'C:\Users\wgchoi\VisionLabWork\FordProject\Dataset\Ford_U of M\Test17.CSV');

NFRAMES = length(imfiles);
if nargin < 8
    firstframe = 1;
    lastframe = length(imfiles);
elseif nargin < 9
    lastframe = length(imfiles);
end
    disp(['Dir ' imgdir ': ' num2str(NFRAMES) ' video frames']);

Im = imread([imgdir imfiles(1).name]);
isz = size(Im);

%% Track init/term parameters
debug = 0; % 1; %1;
nTrackTerm = 4 + 1;
nInitTrack1 = 5;
nInitTrack2 = 3;

%% cparams structure
% should depend on the fstep size
cparams.appth = -log(.55); % 0.7
cparams.ovth = -log(.45); %-log();
cparams.covth = -log(.3); %-log();

cparams.alpha = 1;
cparams.beta = 1;
cparams.adf = 0.6;

%Size related to MeanShift Tracker
mstsize = 1.3;

%% mstrack parameters
%  GENERATE nkernel KERNELs for meanshift TRACKING  composed by increasing
%  number of samples  in [-1 1] x [-1 1] 
nBit = 6;
nkernel = 61;
szkernel = floor(logspace(-log(2)/log(10), log(4)/log(10), nkernel) * 64) * 2;
for s = 1:length(szkernel)
  %                            npoints x width-height
  kernels(s) = buildKernel( szkernel(s)/2, szkernel(s));
end
simth = 0.9;
%Draw some Kernels
figure(4)
for kidx=1:9
    kidx1=kidx*6,
    subplot(3,3,kidx)
    plot3(kernels(kidx1).xs,kernels(kidx1).ys,kernels(kidx1).K,'ob')
end

cparams.nperv = 5;
% cparams.ncar = 5;
cparams.ncamv = 5;
cparams.notrial = 10;

% intial number of people
targetcnt = 1;
trackcnt = 1;
% initial number of cars
targetccnt = 1;
trackccnt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init Z containing all VARIABLES in the MODEL
Z = iZ;
% clear Theta, sparams_cam;
truecnt = []; TP = []; FP = []; FN = []; F = [];
% imfiles(1:150) = [];

faildetect = 0;
icam = [];
% 
% %%%%%%%%%% test remove the first n frames
imfiles([1:(firstframe-1), (lastframe+1):NFRAMES]) = [];
if(size(KLT.x, 2) >= NFRAMES)
    KLT.x(:, [1:(firstframe-1), (lastframe+1):NFRAMES]) = [];
    KLT.y(:, [1:(firstframe-1), (lastframe+1):NFRAMES]) = [];
end
% speed([1:(firstframe-1), (lastframe+1):NFRAMES]) = [];
%%%%%%%%%%%
figure(8); %Previous Frame
figure(9); %Current Frame

%
% PROCESSING FRAMES by FRAMES
for i = 1:fstep:length(imfiles)
    disp(['Processing ' num2str(i) 'th frame']);
    if i > 1
        %Show previous image
        PrevIm=Im; 
        %Save Previous State
        PrevZ=Z;
       
        figure(8);
        imshow(PrevIm);hold on
        plot(KLT.x(PrevZ.gfidx,i-fstep),KLT.y(PrevZ.gfidx,i-fstep),'r*') % RETAINED KLT for frame 1
        
    end
    Im = imread([imgdir imfiles(i).name]);
    %Show Current image
    figure(9);
    imshow(Im);hold on
    plot(KLT.x(:,i),KLT.y(:,i),'r*') %all KLT
    
    
    fileNameNoExt=imfiles(i).name(1:find(imfiles(i).name=='.',1,'last')-1 );
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % READ DETECTIONS in Current Frame of PEOPLE AND CAR
    % X=  Detected PEOPLE
    % Xc= Detected CAR
    [X, Xc, det] = getDets(imgdir, detdir, i-1, detth, isz);
    %Draw detections
    RectangleFaceAlpha(X,9);
    RectangleFaceAlpha(Xc,9);
    
    if exist([imgdir fileNameNoExt '_ant.mat'])
        load([imgdir fileNameNoExt '_ant.mat']);
        truecnt(i) = length(oneFrameAnnotation);
    else
        truecnt(i) = 0; % truecnt(i - fstep);
        oneFrameAnnotation = {};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %Get The TRACKS  from the previous frame FRAMES of i
    [Y] = getMStracks(Im, Z, KLT, i, sparams, kernels, szkernel, mstsize, nBit, simth);
    %%%%%%%%%%%%%% MS track end  correspondence
    %Compute Corrispondences between DETECTED PERSON (X) and PERSON
    %VARIABLES stored in  the model
    [corres] = getCorrespondence2(Z, X, Y, cparams, sparams);
    zcorres = -1 * ones(1, Z.nTargets);
    for j = 1:Z.nTargets
        temp = find(corres == j);
        if ~isempty(temp)
            zcorres(j) = temp;
        end
    end
    
    %Compute Corrispondences between DETECTED CAR (Xc) and CAR VARIABLES 
    %stored in the model
    [corres_car] = getCorrespondenceCars(Z, Xc, cparams, sparams);
    % take MCMC samples! 
    zcorres_car = -1 * ones(1, Z.nCarTargets);
    for j = 1:Z.nCarTargets
        temp = find(corres_car == j);
        if ~isempty(temp)
            zcorres_car(j) = temp;
        end
    end
    %%%%%%%%%%%%%%%%%%% correspondence
    %% %%%%%%%%% get KLT match
    %
    % 
    [Z, tKLT, gfmatch] = KLTmatch(Z, KLT, i, sparams, isz);
    KLTused(i).tKLT = tKLT;
    if i > 1
    % View KLT matching features
    MergedIm=[PrevIm Im];
        figure(10)
        
        imshow(MergedIm);hold on
        plot(KLT.x(Z.gfidx,i-fstep),KLT.y(Z.gfidx,i-fstep),'r*') % RETAINED KLT for frame 1
        plot(2040+KLTused(i).tKLT.x ,KLTused(i).tKLT.y,'g*') % RETAINED KLT for frame 1
        %draw line between corrispondence
        for klti=1:length(Z.gfidx) 
            Xlinea=[2040+KLTused(i).tKLT.x(klti) ;KLT.x(Z.gfidx(klti),i-fstep)];
            Ylinea=[KLTused(i).tKLT.y(klti);KLT.y(Z.gfidx(klti),i-fstep)];
            line(Xlinea,Ylinea)
        end
   end  
    
    %%%%%%%%%%% done
    if i == 1,   sparams.dt = 0;
    else         sparams.dt = fstep/framerate;
    end
    
    %% %%%%%% test
    if 1 % i < fstep * 10 % get better initialization.
        sparams.burinin = 1500;  sparams.thinning = 50;
    else
        sparams.burinin = 2000;  sparams.thinning = 200;
    end
    %% %%%%%%%%%%%%% MCMC part
    tic;
    X.speed = 0; %speed(i);
 
    [Z, zcorres, zcorres_car] = MCMCSamplesJointStatesParametrerizationWCar...
        (Z, X, Xc, Y, tKLT, zcorres, zcorres_car, sparams, i);
 
    toc;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% display results!!!
    disp(['mean camera : ' num2str(mean(Z.cam, 2)')]);
    disp(['mean car : ' num2str(mean(Z.car, 2)')]);
    disp(['Camera Variance' num2str(diag(Z.V{1,1})')])
    disp(['avg feature conf :' num2str(sum(Z.gfeat(3:3:end, :), 2)')]);
    disp(['horizon : ' num2str(Z.cam(4)) ', variance : ' num2str(Z.V{1,1}(4,4))])
    heights = mean(Z.per(5:6:end, :), 2)'
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Get new features and remove invalid features
    [Z, KLT, KLTused] = postProcessFeatures(Z, KLT, KLTused, i, sparams, isz, gfmatch);
    %Draw ground features retained for current frame
    figure(9)       
    imshow(Im)
    plot(KLT.x(Z.gfidx,i),KLT.y(Z.gfidx,i),'b*') % RETAINED KLT for frame 1
    
    
    
    %% Update Model 
    [Z, targetcnt] = updateMStracker(Im, Z, X, targetcnt, corres, kernels, szkernel, nBit, cparams);
    
    [Z, targetccnt] = updateCarTracks(Im, Z, X, targetccnt, corres_car);
    
    %% Manage tracks! filter out obvious miss tracks
    [Z, hTracks, Tracks, hCTracks, CTracks, zcorres, zcorres_car, trackcnt, trackccnt] = ManageTracks(Z, hTracks, Tracks, hCTracks, CTracks, zcorres, zcorres_car, trackcnt, trackccnt, i, sparams, fstep, nInitTrack1, nInitTrack2, nTrackTerm);
    % add function for cars.
    
    %% %%%%%%%%%%% occlusion reasoning
    [Tracks] = inferOcclusion(Tracks, Z, i, sparams);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%% draw trajectories!
    if 1
%         try
            [TP(i), FP(i), FN(i)] = showOneFrame(Im, i, Z, Tracks, CTracks, oneFrameAnnotation, sparams, KLT, tKLT, det, Xc, Y, 0);
            close(2);
            drawnow;
%         catch
%         end    
        recall = sum(TP(1:4:end)) / sum((TP(1:4:end) + FN(1:4:end)))
        FPPI = sum(FP(1:4:end)) / length(FP(1:4:end))
        disp(['vel diff = ' num2str(X.speed - Z.cam(6)) ', (obs : ' num2str(X.speed) ', est : ' num2str(Z.cam(6)) ')'])
    end
    Zs(i)=Z;
end
end

function [Z, targetcnt] = updateMStracker(Im, Z, X, targetcnt, corres, kernels, szkernel, nBit, cparams)
for j = 1:length(corres)
    if corres(j) == 0, continue; end

    [dummy, kidx] = min(abs( szkernel - X.obs(4, j)));
    patch = getImgPatch(Im, X.obs(:, j)');
    patch = uint8(imresize(patch, [szkernel(kidx) szkernel(kidx)/2]));
    
    Qc = bitshift(reshape(patch, [], 3), nBit-8);

    Z.model(corres(j)).timg = patch;
    Z.model(corres(j)).qS = (1 - cparams.adf) * buildHist( Qc, kernels(kidx), nBit ) + cparams.adf * Z.model(corres(j)).qS;
    % update this by Z estimation..
end

newtracks = find(corres == 0);
for j = newtracks
    [dummy, kidx] = min(abs(szkernel - X.obs(4, j)));

    patch = getImgPatch(Im, X.obs(:, j)');
    patch = uint8(imresize(patch, [szkernel(kidx) szkernel(kidx)/2]));

    Qc = bitshift(reshape(patch, [], 3), nBit-8);

    model.timg = patch;
    model.qS = buildHist( Qc, kernels(kidx), nBit );
    model.ppos = [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];        
    Z.model = [Z.model, model];

    Z.peridx = [Z.peridx, targetcnt];
    Z.nTargets = Z.nTargets + 1;
    targetcnt = targetcnt + 1;
end

end
%function [X, Xc, det, detc] = getDets(imgdir, detdir, idx, detth, isz)
%
% Read Detections of PEOPLE and CARS in Frame idx and return DETECTION data
% filtered and unfiltered
% PARAMETERS INPUT:
%
%  - imgdir: directory containing the frames
%  - detdir: directory containing the outputs of Dector for each frame
%  - idx: frame index
%  - detth:  Detection Threshold
%  - isz: size of the frame
%
% PARAMETERS OUTPUT:
% 
%  - X:  struct data contaning Detections set to []
%  - Xc: struct data contaning filtered Detections
%           idx:  id of Detection
%           obs:  Detection Data (uo,vo,w,h)
%           pobj: Probability of Detection for the scale
%
%  - det:  [ ]
%  - detc: filtered Detections (uo,vo,w,h,Scale)

function [X, Xc, det, detc] = getDets(imgdir, detdir, idx, detth, isz)
det = [];
detc = [];

% Xc DETECTION
%
if exist([detdir '/conf_' num2str(idx, '%08d') '.conf'])
    %                                      (down-left)-loc.     upper-right loc.                                  
    %                                                    |       | 
    %Detection Vector (num. detection in frame) x 5: [ Uo , Vo , U1 , V1 , Scale]
    detc = load([detdir '/conf_' num2str(idx, '%08d') '.conf']);

%     detc(1, :) = [];
    
    %% Detection Filtering
    %
    % 1) Remove Detection under Threshold
    detc(detc(:, 5) < detth, :) = [];
    %Calculate the W and H of the box
    detc(:, 3:4) = detc(:, 3:4) - detc(:, 1:2) + 1;

    % 2) Reduce BBOX-width
    for j = 1:size(detc, 1)
        %If width of bbox > img width - 15 reduce the bounding box width
        if detc(j, 1) + detc(j, 3) > isz(2) - 15
            detc(j, 3) = isz(2) - 15 - detc(j, 1);
        end
    end
    
    Xc.obs = [];
    Xc.idx = [];
    % 3) Calculate for each detection couple (j,k) by getOverlap(Dj,Dk) in getCorrispondence.m 
    %    the OVERLAP betwee bounding bbox
    % 
    omat = zeros(size(detc, 1), size(detc, 1));
    for j = 1:size(detc, 1)
        for k = j + 1:size(detc, 1)
            omat(j, k) = getOverlap(detc(j, 1:4), detc(k, 1:4));
        end
    end
     % 4) Estraction of Overlapping Detections > 0.4 
    [rid, cid]=find(omat > 0.4);

    filterout = [];
    
    % 5) Scale Analysis: Select Overlapping Detection with the Minimum
    %    Scale value
    for j = 1:length(rid)
        if detc(rid(j), 5) > detc(cid(j), 5)
            % set to be deleted
            filterout(end + 1) = cid(j);
        else
            filterout(end + 1) = rid(j);
        end
    end
    %delete selected Detection
    detc(filterout, :) = [];

    %Xc struct:
    % 
    %  idx:  id of Detection
    %  obs:  Detection Data (uo,vo,w,h)
    %  pobj: Probability of Detection for the scale
    Xc.idx = 1:size(detc, 1);
    for j = 1:size(detc, 1)
        Xc.obs(:, j) = detc(j, 1:4)';
        Xc.pobj(j) = 1/(1+exp(-detc(j, 5)));
    end
end

% X DETECTION
%
% if exist([imgdir fileNameNoExt ext])
%     det = load([imgdir fileNameNoExt ext]);
%     det(1, :) = [];
%     det(det(:, 6) < detth, :) = [];
% 
%     X.obs = [];
%     X.idx = [];
% 
%     omat = zeros(size(det, 1), size(det, 1));
%     for j = 1:size(det, 1)
%         for k = j + 1:size(det, 1)
%             omat(j, k) = getOverlap(det(j, 1:4), det(k, 1:4));
%         end
%     end
%     [rid, cid]=find(omat > 0.4);
% 
%     filterout = [];
%     for j = 1:length(rid)
%         if det(rid(j), 6) > det(cid(j), 6)
%             filterout(end + 1) = cid(j);
%         else
%             filterout(end + 1) = rid(j);
%         end
%     end
%     det(filterout, :) = [];
% 
%     X.idx = 1:size(det, 1);
%     for j = 1:size(det, 1)
%         X.obs(:, j) = det(j, 1:4)';
%         X.pobj(j) = 1/(1+exp(-det(j, 6)));
%     end
% else
%     error;
% end

X.idx = [];
X.obs = [];
X.pobj = [];

end
%function [Y] = getMStracks(Im, Z, KLT, frameidx, sparams, kernels, szkernel, mstsize, nBit, simth)
% Get the TRACKS hypotesis
% 
% PARAMETERS INPUT
% - Im: Frame 
% - Z: current model
% - KLT: klt features
% - frameidx: frame index 
% - kernels: KERNELs of different numbers of point
% - szkernel: Number of points in the i-Kernel --> (szkernel/2,szkernel)
%
%  PARAMETERS OUTPUT:
% - Y: Vector 3x(Number of Target) containing TARGET LOCATION (u,v) and Kernel size
%      describing each TARGET TRACKS
function [Y] = getMStracks(Im, Z, KLT, frameidx, sparams, kernels, szkernel, mstsize, nBit, simth)
%%%%% get MS tracker
Y = zeros(3, length(Z.model));
%calculate pan angle estimate
tpan = getRoughCamPan(Z, KLT, frameidx, sparams);

%Only after processing First frame
% For each TARGET j do Meanshift Tracker
for j = 1:length(Z.model)
 
    %Get BOUNDING BOX of TARGETs given the current VARIABLE STATES
    [bb] = getImageProjections(Z, j, sparams, 0, tpan);
    %Get a mean value for the BBOX
    bb = mean(bb, 2);
    % ppos feature vector (central point of bbox) x,y,w,h
    ppos = [bb(1:2) + bb(3:4)/2; bb(3:4)];
    
    %% Get candidate scales values for start Tracking 
    % Minimum scale value: ppos(4) / mstsize --> bbox height/mstsize
    % Maximum scale value: ppos(4) * mstsize --> bbox height*mstsize
    cands = find(szkernel >= (ppos(4) / mstsize) & szkernel <= (ppos(4) * mstsize));
    bestSim = 0;        
    %Tracking with different kernels and save the best RESULT 
    for s = cands
        %sim: similiarity measure between TARGET MODEL and TARGET CANDIDATE
        [p, pos, Ic, sim] = kernelTrack(Im, Z.model(j).qS, ppos(1:2)', kernels(s), nBit);
        
        if sim > bestSim, best={p, pos, Ic, s}; bestSim = sim; end
    end
    %If the Similiarity measure is up to the threshoold save in Y: 
    % 1) Location of the Target --> best{2}
    % 2) Kernel size associated to best response
    if bestSim > simth
        Y(1:2, j) = best{2};
        Y(3, j) = szkernel(best{4});
    end
end
%%%%%%%%%%%%%% MS track end
end

%function [Z, tKLT, gfmatch] = KLTmatch(Z, KLT, frameidx, sparams, isz)
% MATCHING between KLT Features.
% The features are already matched in the DATASET, so we can serach
% by index the CORRISPONDENCE between frames
%
% PARAMETERS OUT
%
% - Z: Model Data Updated with MATCHED KLT features
% 
% - tKLT: MATCHED KLT FEATURES in a struct format (x,y,idx)
%
% - gfmatch: INDEX of MATCHED KLT FEATURES in 
function [Z, tKLT, gfmatch] = KLTmatch(Z, KLT, frameidx, sparams, isz)
%Filtering of KLT FEATURES with imposed spatial costraint
%KLTidx: index of VALID KLT FEATURES
KLTidx = getValidKLT2(KLT, frameidx, Z, isz(2), sparams)';
removelist = [];
gfmatch = [];

%MATCHING PROCEDURE: for the first frame Z.gfdix=[]
for j = 1:length(Z.gfidx)
    %search corrispondence betwen features in Z and KLT
    temp = find(KLT.vidx(KLTidx) == Z.gfidx(j));
   
    if ~isempty(temp)
        %store matched features
        gfmatch(end + 1) = temp;
    else
        % store features to be removed
        removelist(end+1) = j;
    end
end

% temp KLT: features matched
tKLT.x = KLT.x(KLTidx(gfmatch), frameidx);
tKLT.y = KLT.y(KLTidx(gfmatch), frameidx);
tKLT.idx = KLT.vidx(KLTidx(gfmatch));

% remove out non-detected features!
% First frame, removelist is []
Z = FilterOutGFeats(Z, removelist, sparams);
end


%function [Z, KLT, KLTused] = postProcessFeatures(Z, KLT, KLTused, frameidx, sparams, isz, gfmatch)
% Post Processing of KLT FEATURES: save the KLT features with high
% confidence level
%
% 
function [Z, KLT, KLTused] = postProcessFeatures(Z, KLT, KLTused, frameidx, sparams, isz, gfmatch)

% filter out features!
removelist = find(sum(Z.gfeat(3:3:end, :), 2) < 0.2 * Z.nSamples);
% filterout from KLTs
removeKLT = [];
KLTused(frameidx).invalid = Z.gfidx(removelist);

for j = removelist(:)'
    removeKLT(end + 1) = find(KLT.vidx == Z.gfidx(j));
end
KLT.x(removeKLT, :) = [];
KLT.y(removeKLT, :) = [];
KLT.val(removeKLT, :) = [];
KLT.vidx(removeKLT) = [];
Z = FilterOutGFeats(Z, removelist(:)', sparams);

%KLT FEATURES INIT (when processing first frame)
if Z.nFeats < sparams.nfeatuse
    %Filter KLT features 
    KLTidx = getValidKLT2(KLT, frameidx, Z, isz(2), sparams, 0)';
    %check for some error
    if sum(KLT.x(KLTidx, frameidx) == 0) > 0
        keyboard;
    end

    % Sort KLT values response in Ascending Order
    [dummy, rp] = sort(KLT.val(KLTidx));
    rp = rp(end:-1:1); %Index of KLT Values response in Descending Order
%         rp = KLTidx(rp);

    %Conservo le KLT FEATURES che hanno risposta MAGGIORE
    cnt = Z.nFeats; %counter
    j = 0;

    while(cnt < sparams.nfeatuse)
        j = j + 1;
        if j > length(rp) %exit from while
            break;
        end
        %Get index of Features with MAXIMUM RESPONSE
        idx1 = rp(j); 
        if ~isempty(find(gfmatch == idx1))
            continue;
        end
        %check error
        if sum([KLT.x(KLTidx(idx1), frameidx); KLT.y(KLTidx(idx1), frameidx)] <= 0) > 0
            keyboard;
        end
        %Analyze all CAMERA STATES :
        %
        for k = 1:Z.nSamples
            %Get k-Camera Sample Status:  each camera configuration contain
            % size(tempcam, 2) samples of Camera Status
            tempcam = Z.cams{k};
            %Get temp Ground Feauture 3 x 30
            tempgf = zeros(3, size(tempcam, 2));
            
            %For each Status in tempgf :
            %  BackProject the KLT feature in the current frame to the 
            %  3D World and generate GROUND FEATURE on Vector tempgf containing
            %  different kind of information data:
            
            for l = 1:size(tempcam, 2)
                %                                       ground Im. point                   not used    cam. status   prob. confidence
                %                                           |                                  |             |           |
                tempgf(:, l) = [getGIProj([KLT.x(KLTidx(idx1), frameidx); KLT.y(KLTidx(idx1), frameidx)], tempcam(:, l)); 0.9];
            end
            
            %GROUND FEATURE MEAN VALUE for k-Sample --> params.ngfeat = 3;
            %Append by row Z.gefeat the new mean value of tempgf
            Z.gfeat((Z.nFeats*sparams.ngfeat + 1):((Z.nFeats+1)*sparams.ngfeat), k) = mean(tempgf, 2);
        
            %generate the GROUND FEATURE COVARIANCE only for (x,z) data
            Z.gfV{Z.nFeats+1, k} = cov(tempgf(1:2, :)') + (0.15)^2*eye(2); % safe guard...
        end
        %Save the index of KLT Features for the current STATE
        Z.gfidx(end + 1) = KLT.vidx(KLTidx(idx1));
        Z.gfcnt(end + 1) = 1;
        %counter
        cnt = cnt + 1;
        
        Z.nFeats = Z.nFeats + 1;
    end
end
end


function [Z, targetccnt] = updateCarTracks(Im, Z, X, targetccnt, corres_car)

newtracks = find(corres_car == 0);
for j = newtracks
    Z.caridx = [Z.caridx, targetccnt];
    Z.nCarTargets = Z.nCarTargets  + 1;
    targetccnt = targetccnt + 1;
end

end

function [Z, hTracks, Tracks, hCTracks, CTracks, zcorres, zcorres_car, trackcnt, trackccnt] = ManageTracks(Z, hTracks, Tracks, hCTracks, CTracks, zcorres, zcorres_car, trackcnt, trackccnt, frameidx, sparams, fstep, nInitTrack1, nInitTrack2, nTrackTerm)
%% Filter out tracks
filterout = [];
for j = 1:Z.nTargets
    if sum(Z.per((sparams.nperv + 1) * j, :)) < Z.nSamples * 0.2
        filterout = [filterout, j];
    end
end
[Z, hTracks, Tracks] = RemoveTracks(Z, hTracks, Tracks, filterout, sparams);
zcorres(filterout) = [];
%% Filter out cars
filterout = [];
for j = 1:Z.nCarTargets
    if sum(Z.car((sparams.ncarv + 1) * j, :)) < Z.nSamples * 0.2
        filterout = [filterout, j];
    end
end

[Z, hCTracks, CTracks] = RemoveCarTracks(Z, hCTracks, CTracks, filterout, sparams);
zcorres_car(filterout) = [];
%% initiate and terminate tracks
targetidx = Z.peridx;
targetidx = targetidx(:)';
detidx = setdiff(targetidx, targetidx(zcorres == -1));

for j = detidx(:)'
    matched = 0;
    for k = 1:length(hTracks)
        if hTracks(k).tid == j
            matched = 1;
            hTracks(k).det(end + 1) = frameidx;
            break;
        end
    end
    for k = 1:length(Tracks)
        if Tracks(k).tid == j
            matched = 1;
            Tracks(k).det(end + 1) = frameidx;
            break;
        end
    end
    if matched == 0
        track.id = -1;
        track.tid = j;
        track.det = frameidx;
        track.term = -1;
        track.occlusion = [];
        hTracks = [hTracks, track];
    end
end

removelist = [];
for j = 1:length(hTracks)
%         if length(hTracks(j).det) >= 5
    if sum(hTracks(j).det > frameidx-(nInitTrack1*fstep)) >= nInitTrack2
        hTracks(j).id = trackcnt;
        Tracks = [Tracks, hTracks(j)];
        % remove from hypothesis
        removelist = [removelist, j];        
        trackcnt = trackcnt + 1;
    end
%         end
    if (hTracks(j).det(end) < frameidx - (nTrackTerm*fstep))
        zidx = find(Z.peridx == hTracks(j).tid);
        targetidx = unique(Z.peridx);

        Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), :) = [];
        mid = find(targetidx == hTracks(j).tid);
        Z.model(mid) = [];
        Z.peridx(mid) = [];
        Z.tcnt(mid) = [];
        Z.nTargets = Z.nTargets - 1;
        %%%%% remove interaction too!!!!!
        Z.beta(:, mid, :, :) = [];
        Z.beta(mid, :, :, :) = [];

        removelist = [removelist, j];
    end
end

hTracks(removelist) = [];
for j = 1:length(Tracks)
     if (Tracks(j).term == -1) && (Tracks(j).det(end) < frameidx - (nTrackTerm*fstep))
        zidx = find(Z.peridx == Tracks(j).tid);
        targetidx = unique(Z.peridx);

        Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), :) = [];
        mid = find(targetidx == Tracks(j).tid);
        Z.model(mid) = [];
        Z.peridx(mid) = [];
        Z.tcnt(mid) = [];
        Z.nTargets = Z.nTargets - 1;
        %%%%% remove interaction too!!!!!
        Z.beta(:, mid, :, :) = [];
        Z.beta(mid, :, :, :) = [];

        Tracks(j).term = frameidx;
    end
end
%% initiate and terminate car tracks
caridx = Z.caridx;
caridx = caridx(:)';
detidx = setdiff(caridx, caridx(zcorres_car == -1));

for j = detidx(:)'
    matched = 0;
    for k = 1:length(hCTracks)
        if hCTracks(k).tid == j
            matched = 1;
            hCTracks(k).det(end + 1) = frameidx;
            break;
        end
    end
    
    for k = 1:length(CTracks)
        if CTracks(k).tid == j
            matched = 1;
            CTracks(k).det(end + 1) = frameidx;
            break;
        end
    end
    
    if matched == 0
        track.id = -1;
        track.tid = j;
        track.det = frameidx;
        track.term = -1;
        track.occlusion = [];
        hCTracks = [hCTracks, track];
    end
end

removelist = [];
for j = 1:length(hCTracks)
%         if length(hTracks(j).det) >= 5
    if sum(hCTracks(j).det > frameidx-(nInitTrack1*fstep)) >= nInitTrack2
        hCTracks(j).id = trackccnt;
        CTracks = [CTracks, hCTracks(j)];
        % remove from hypothesis
        removelist = [removelist, j];        
        trackccnt = trackccnt + 1;
    end
%         end
    if (hCTracks(j).det(end) < frameidx - (nTrackTerm*fstep))
        zidx = find(Z.caridx == hCTracks(j).tid);
        targetidx = unique(Z.caridx);

        Z.car((1 + (zidx-1) * (sparams.ncarv + 1)):(zidx * (sparams.ncarv + 1)), :) = [];
        
        mid = find(targetidx == hCTracks(j).tid);
        Z.caridx(mid) = [];
        Z.tccnt(mid) = [];
        Z.nCarTargets = Z.nCarTargets - 1;
        
        removelist = [removelist, j];
    end
end
hCTracks(removelist) = [];

for j = 1:length(CTracks)
     if (CTracks(j).term == -1) && (CTracks(j).det(end) < frameidx - (nTrackTerm*fstep))
        zidx = find(Z.caridx == CTracks(j).tid);
        targetidx = unique(Z.caridx);

        Z.car((1 + (zidx-1) * (sparams.ncarv + 1)):(zidx * (sparams.ncarv + 1)), :) = [];
        mid = find(targetidx == CTracks(j).tid);
        Z.caridx(mid) = [];
        Z.tccnt(mid) = [];
        Z.nCarTargets = Z.nCarTargets - 1;

        Tracks(j).term = frameidx;
    end
end
end
