clear;
addpath(genpath('common'))
addpath(genpath('mexFunctions'))

% imgdir = '../test4/';
% imfiles = dir([imgdir '*.jpg']);
% ext = '_det.txt';
% fstep = 2;
% framerate = 30;
% Z.nSamples = 1;
% Z.cam = [500 1.7 360 180 0 0 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-50 50 1e-50 .2 1e-10 1e-10]).^2;
% Z.W = 1;

imgdir = '../Dataset/ETH/seq02/';
imfiles = dir([imgdir '*.png']);
ext = '_det06.txt';
fstep = 1;
framerate = 14;
Z.nSamples = 1;
Z.cam = [500 0.94 320 220 0 1.3 0 0]';
Z.V{1, 1} = diag([1 0.01 1e-50 20 1e-6 .4 1e-6 1e-6]).^2;
Z.W = 1;
% 
% imgdir = '../Dataset/ETH/seq03/';
% imfiles = dir([imgdir '*.png']);
% ext = '_det06.txt';
% fstep = 1;
% framerate = 14;
% Z.nSamples = 1;
% Z.cam = [500 1.1 320 220 0 .5 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-6 5 1e-6 .2 1e-6 1e-6]).^2;
% Z.W = 1;
% 
% imgdir = '../Dataset/ETH/seq04/';
% imfiles = dir([imgdir '*.png']);
% ext = '_det06.txt';
% fstep = 1;
% framerate = 15;
% Z.nSamples = 1;
% Z.cam = [500 0.94 320 235 0 1.2 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-6 10 1e-6 .1 1e-6 1e-6]).^2;
% Z.W = 1;

Im = imread([imgdir imfiles(1).name]);
isz = size(Im);
[KLT.x, KLT.y, KLT.val] = klt_read_featuretable([imgdir, 'features.txt']);
KLT = preprocessKLT(KLT);

debug = 0; %1;
nTrackTerm = 4 + 1;
nInitTrack1 = 5;
nInitTrack2 = 3;

detth = -0.5;
% should depend on the fstep size
cparams.appth = -log(.5); % 0.7
cparams.ovth = -log(.4); %-log();
cparams.alpha = 2;
cparams.beta = 1;
cparams.adf = 0.5;
mstsize = 1.3;

sparams = InitParams(3, fstep, framerate);

% sparams_cam.asr = 100;
%%%% mstrack parameters
nBit = 4;
nkernel = 61;
szkernel = floor(logspace(-log(2)/log(10), log(4)/log(10), nkernel) * 64) * 2;
for s = 1:length(szkernel)
  kernels(s) = buildKernel( szkernel(s)/2, szkernel(s));
end
simth = 0.94;

cparams.nperv = 5;
cparams.ncamv = 5;
cparams.notrial = 10;

Tracks = [];
hTracks = [];

targetcnt = 1;
trackcnt = 1;
% 
% clear Z;
% 
% sparams_cam.A = eye(5);
% sparams_cam.R = diag([25, 0.0001, 4, 4, 0]);
% sparams_cam.lenState = 5;
% sparams_cam.numPerTarget = 200; %sparams.nretain;
% Theta = SampleIndependentParticles(Theta, sparams_cam);

Z.per = zeros(0, Z.nSamples);
Z.gfeat = zeros(0, Z.nSamples);
Z.gfidx = [];
Z.peridx = [];
Z.tcnt = [];
Z.gfcnt = [];
Z.model = [];
Z.nTargets = 0;
Z.nFeats = 0;
Z.beta = [];
% clear Theta, sparams_cam;
truecnt = []; TP = []; FP = []; FN = [];

% imfiles(1:150) = [];

dbgklt = [];
camswt = 1;
camvacc = 0;
for i = 1:fstep:length(imfiles)
    Im = imread([imgdir imfiles(i).name]);
    fileNameNoExt=imfiles(i).name(1:find(imfiles(i).name=='.',1,'last')-1 );
    if exist([imgdir fileNameNoExt ext])
        det = load([imgdir fileNameNoExt ext]);
        det(1, :) = [];
        det(det(:, 6) < detth, :) = [];

        X.obs = [];
        X.idx = [];
        
        omat = zeros(size(det, 1), size(det, 1));
        for j = 1:size(det, 1)
            for k = j + 1:size(det, 1)
                omat(j, k) = getOverlap(det(j, 1:4), det(k, 1:4));
            end
        end
        [rid, cid]=find(omat > 0.55);
        
        filterout = [];
        for j = 1:length(rid)
            if det(rid(j), 6) > det(cid(j), 6)
                filterout(end + 1) = cid(j);
            else
                filterout(end + 1) = rid(j);
            end
        end
        det(filterout, :) = [];
        
        X.idx = 1:size(det, 1);
        for j = 1:size(det, 1)
            X.obs(:, j) = det(j, 1:4)';
            X.pobj(j) = 1/(1+exp(-det(j, 6)));
        end
    else
        error;
    end
    if exist([imgdir fileNameNoExt '_ant.mat'])
        load([imgdir fileNameNoExt '_ant.mat']);
        truecnt(i) = length(oneFrameAnnotation);
    else
        truecnt(i) = 0; % truecnt(i - fstep);
        oneFrameAnnotation = {};
    end

    Y = zeros(3, length(Z.model));
    for j = 1:length(Z.model)
        [bb] = getImageProjections(Z, j, sparams, 1);
        bb = mean(bb, 2);
        ppos = [bb(1:2) + bb(3:4)/2; bb(3:4)];
        % get candidate scales        
%         cands = find(szkernel >= (Z.model(j).ppos(4) / mstsize) & szkernel <= (Z.model(j).ppos(4) * mstsize));
        cands = find(szkernel >= (ppos(4) / mstsize) & szkernel <= (ppos(4) * mstsize));
        bestSim = 0;        
        for s = cands
            [p, pos, Ic, sim] = kernelTrack(Im, Z.model(j).qS, ppos(1:2)', kernels(s), nBit);
%             [p, pos, Ic, sim] = kernelTrack(Im, Z.model(j).qS, Z.model(j).ppos(1:2)', kernels(s), nBit);
            if sim > bestSim, best={p, pos, Ic, s}; bestSim = sim; end
        end
%         szkernel(best{4})
        % super care needed!
        if bestSim > simth
            Y(1:2, j) = best{2};
            Y(3, j) = szkernel(best{4});
%             Z.model(j).ppos = [Y(1, j) - Y(3, j)/4, Y(2, j) - Y(3, j)/2, Y(3, j)/2, Y(3, j)]';% [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];
        end
    end
    [corres] = getCorrespondence2(Z, X, Y, cparams, sparams);

    %%%%% show corrspondences....debug;
    if debug == 1
        figure(3); 
        clf;
        K = 10;        
        for j = 1:min(K, length(Z.model))
            subplot(2, K, j);
            imshow(Z.model(j).timg);
            
            if sum(corres == j) > 0
                title('matched')
            else
                title('unmatched')
            end
        end
        
        cnt = sum(corres ~= 0) + 1;
        for j = 1:min(K, size(X.obs, 2))
            patch = getImgPatch(Im, X.obs(:, j)');
            if corres(j) == 0
                if cnt > K
                    continue;
                end
                subplot(2, K, K + cnt);
                imshow(patch);
                cnt = cnt + 1;
                title('unmatched')
            else
                if corres(j) > K
                    continue;
                end
                subplot(2, K, K + corres(j));
                imshow(patch);
                title('matched')
            end
        end
        pause(1);
    end
    
    % take MCMC samples! 
    zcorres = -1 * ones(1, Z.nTargets);
    for j = 1:Z.nTargets
        temp = find(corres == j);
        if ~isempty(temp)
            zcorres(j) = temp;
        end
    end
    
    KLTidx = getValidKLT2(KLT, i, Z, isz(2), sparams)';
    removelist = [];
    gfmatch = [];
    for j = 1:length(Z.gfidx)
        temp = find(KLT.vidx(KLTidx) == Z.gfidx(j));
        if ~isempty(temp)
            gfmatch(end + 1) = temp;
        else
            removelist(end+1) = j;
        end
    end
    
    % temp
    tKLT.x = KLT.x(KLTidx(gfmatch), i);
    tKLT.y = KLT.y(KLTidx(gfmatch), i);
    tKLT.idx = KLT.vidx(KLTidx(gfmatch));
    % remove out non-detected features!
    Z = FilterOutGFeats(Z, removelist, sparams);
%     if Z.nFeats > 0
%         removeidx = [];
%         for j = removelist
%             removeidx = [removeidx, [(1:sparams.ngfeat) + (j-1)*sparams.ngfeat]];
%         end
%         Z.gfeat(removeidx, :) = [];
%         Z.gfV(removelist, :) = [];
%         Z.gfidx(removelist) = [];
%         Z.gfcnt(removelist) = [];
%         Z.nFeats = length(gfmatch);
%     end    
    % need to deal with indexing (correspondence)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %%%%%% temp
%     corres
%     zcorres
    if i == 1
        sparams.dt = 0;
%         sparams.mPertCam = 36;
%         sparams.mPertPer = 4;
    else
        sparams.dt = fstep/framerate;
%         sparams.mPertCam = 64;
%         sparams.mPertPer = 2;
    end
    
%     if i < 3
%         sparams.burnin = 1000;
%         params.thinning = 100;
%     else
%         params.burnin = 1000;
%         params.thinning = 100;
%     end
%     
    if sum(tKLT.y < 200) > 0
        a=1;
    end
    tic;
    [Z, zcorres] = MCMCSamplesJointStatesParametrerization(Z, X, Y, tKLT, zcorres, sparams);
    toc;
    disp(['avg feature conf :' num2str(sum(Z.gfeat(3:3:end, :), 2)')]);
    
%     sqrt(diag(Z.V{2, 1}(3:4, 3:4))') * framerate/fstep
    % filter out features!
    removelist = find(sum(Z.gfeat(3:3:end, :), 2) < 0.2 * Z.nSamples);
    % filterout from KLTs
    removeKLT = [];
    for j = removelist(:)'
        removeKLT(end + 1) = find(KLT.vidx == Z.gfidx(j));
    end
    KLT.x(removeKLT, :) = [];
    KLT.y(removeKLT, :) = [];
    KLT.val(removeKLT, :) = [];
    KLT.vidx(removeKLT) = [];
    Z = FilterOutGFeats(Z, removelist(:)', sparams);
    
    % initialize features!
    if Z.nFeats < sparams.nfeatuse
        KLTidx = getValidKLT2(KLT, i, Z, isz(2), sparams, 0)';
        
        if sum(KLT.x(KLTidx, i) == 0) > 0
            keyboard;
        end
        
        rp = randperm(length(KLTidx));
        cnt = Z.nFeats; j = 0;
        
        while(cnt < sparams.nfeatuse)
            j = j + 1;
            if j > length(rp)
                break;
            end
            idx1 = rp(j);
            if ~isempty(find(gfmatch == idx1))
                continue;
            end

            if sum([KLT.x(KLTidx(idx1), i); KLT.y(KLTidx(idx1), i)] <= 0) > 0
                keyboard;
            end
            for k = 1:Z.nSamples
                tempcam = Z.cams{k};
                tempgf = zeros(3, size(tempcam, 2));
                for l = 1:size(tempcam, 2)
                    tempgf(:, l) = [getGIProj([KLT.x(KLTidx(idx1), i); KLT.y(KLTidx(idx1), i)], tempcam(:, l)); 0.97];
                end
                Z.gfeat((Z.nFeats*sparams.ngfeat + 1):((Z.nFeats+1)*sparams.ngfeat), k) = mean(tempgf, 2);
                Z.gfV{Z.nFeats+1, k} = cov(tempgf(1:2, :)') + (0.15)^2*eye(2); % safe guard...
            end

            Z.gfidx(end + 1) = KLT.vidx(KLTidx(idx1));
            Z.gfcnt(end + 1) = 1;
            cnt = cnt + 1;
            Z.nFeats = Z.nFeats + 1;
        end
    end
    sqrt(diag(Z.V{1,1})')
    Z.per(6:6:end, :)
    mean(Z.per(5:6:end, :), 2)'
    mean(Z.cam, 2)'
    %%% update model ?
    for j = 1:length(corres)
        if corres(j) == 0, continue; end
        
        [dummy, kidx] = min(abs( szkernel - X.obs(4, j)));
        patch = getImgPatch(Im, X.obs(:, j)');
        patch = uint8(imresize(patch, [szkernel(kidx) szkernel(kidx)/2]));
        
        Qc = bitshift(reshape(patch, [], 3), nBit-8);
        
        Z.model(corres(j)).timg = patch;
        Z.model(corres(j)).qS = (1 - cparams.adf) * buildHist( Qc, kernels(kidx), nBit ) + cparams.adf * Z.model(corres(j)).qS;
        % update this by Z estimation..
%         temp = getImageProjections(Z, corres(j), sparams);
%         temp = mean(temp, 2);
%         Z.model(corres(j)).ppos = [temp(1:2) + temp(3:4) / 2; temp(3:4)];
%         Z.model(corres(j)).ppos = [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];
    end

    newtracks = find(corres == 0);
    for j = newtracks
        [dummy, kidx] = min(abs( szkernel - X.obs(4, j)));

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
    
    % filter out obvious miss tracks
    filterout = [];
    for j = 1:Z.nTargets
        if sum(Z.per((sparams.nperv + 1) * j, :)) < Z.nSamples * 0.35
            filterout = [filterout, j];
        end
    end
    removelist = [];
    for k = 1:length(hTracks)
        idx = find(hTracks(k).tid == Z.peridx(filterout));
        if ~isempty(idx)
            removelist = [removelist, k]; 
        end
    end
    hTracks(removelist) = [];
    for k = 1:length(Tracks)
        idx = find(Tracks(k).tid == Z.peridx(filterout));
        if ~isempty(idx)
            Tracks(k).term = i;
        end
    end
    
    idx = [];
    for j = filterout
        idx = [idx, ((sparams.nperv + 1) * (j - 1) + 1):((sparams.nperv + 1) * j)];
    end
    Z.per(idx, :) = [];
    Z.model(filterout) = [];
    Z.peridx(filterout) = [];
    Z.tcnt(filterout) = [];
    Z.nTargets = Z.nTargets - length(filterout);
    %%%%% remove interaction too!!!!!
    Z.beta(:, filterout, :, :) = [];
    Z.beta(filterout, :, :, :) = [];
    
    zcorres(filterout) = [];
    
    % initiate and terminate tracks
    targetidx = Z.peridx;
    targetidx = targetidx(:)';
    detidx = setdiff(targetidx, targetidx(zcorres == -1));
%     newdets = sum(corres == 0);
%     detidx = [detidx, targetidx((end - newdets + 1):end)];
    for j = detidx(:)'
        matched = 0;
        for k = 1:length(hTracks)
            if hTracks(k).tid == j
                matched = 1;
                hTracks(k).det(end + 1) = i;
                break;
            end
        end
        for k = 1:length(Tracks)
            if Tracks(k).tid == j
                matched = 1;
                Tracks(k).det(end + 1) = i;
                break;
            end
        end
        if matched == 0
            track.id = -1;
            track.tid = j;
            track.det = i;
            track.term = -1;
            hTracks = [hTracks, track];
        end
    end
    
    removelist = [];
    for j = 1:length(hTracks)
%         if length(hTracks(j).det) >= 5
        if sum(hTracks(j).det > i-(nInitTrack1*fstep)) >= nInitTrack2
            hTracks(j).id = trackcnt;
            Tracks = [Tracks, hTracks(j)];
            % remove from hypothesis
            removelist = [removelist, j];        
            trackcnt = trackcnt + 1;
        end
%         end
        if (hTracks(j).det(end) < i - (nTrackTerm*fstep))
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
         if (Tracks(j).term == -1) && (Tracks(j).det(end) < i - (nTrackTerm*fstep))
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
            
            Tracks(j).term = i;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% draw trajectories!
    try
        [TP(i), FP(i), FN(i)] = showOneFrame(Im, i, Z, Tracks, oneFrameAnnotation, sparams, KLT, tKLT, det, Y, 0);
%         tKLT.y'
%         Z.gfidx'
        
        figure(1); F(i) = getframe;
        figure(2); tF2 = getframe;

        stf1 = size(F(i).cdata);
        tF2.cdata = imresize(tF2.cdata, [stf1(1) stf1(1)]);
        stf2 = size(tF2.cdata);
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 1) = [tF2.cdata(:,:,1); zeros(stf1(1) - stf2(1), stf2(2))];
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 2) = [tF2.cdata(:,:,2); zeros(stf1(1) - stf2(1), stf2(2))];
        F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 3) = [tF2.cdata(:,:,3); zeros(stf1(1) - stf2(1), stf2(2))];
    catch
    end
    
    Zs(i) = Z;
    
    recall = sum(TP) / sum((TP + FN))
    FPPI = sum(FP) / length(TP)
%     
%     if length(Z.peridx) < 2
%         camvacc = camvacc + (25*fstep/framerate)^2;
%         sparams.Qcam(4,4) = 0.1; % turn off horizon...
% %         keyboard
%     else
%         sparams.Qcam(4,4) = camvacc; % (25*fstep/framerate)^2;
%         camvacc = (25*fstep/framerate)^2;
%     end
end

close all;
% movie2avi(F, 'test_mstrack_100s.avi');

for i = 1:length(F)
    ttt(i) = 0;
    for j = 1:length(Tracks) 
        if Tracks(j).det(1) <= i && Tracks(j).term > i, ttt(i) = ttt(i) + 1;  end; 
    end
end
figure(1); plot(truecnt);
hold on; plot(ttt, ':'); hold off
title('Number of people in the scene'); legend({'Annotation', 'Tracked'})

figure(2); plot(TP ./(TP+FN)); title('Recall rate in each frame'); figure(3); plot(TP ./(TP+FP)); title('Precision in each frame');

recall = sum(TP ./(TP + FN)) / length(TP)
FPPI = sum(FP) / length(TP)