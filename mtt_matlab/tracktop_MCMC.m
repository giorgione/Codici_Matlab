clear;
addpath(genpath('common'))

imgdir = '../test4/';


imfiles = dir([imgdir '*.jpg']);
ext = '_det.txt';

Im = imread([imgdir imfiles(1).name]);
isz = size(Im);
[KLT.x, KLT.y, KLT.val] = klt_read_featuretable([imgdir, 'features.txt']);
% load testinit_uniform.mat

% Theta.state = [Theta.state; zeros(1, size(Theta.state, 2))]; % add 0 panning

debug = 0; %1;

nTrackTerm = 7;

fstep = 1;
dt = fstep/15;
% 
% sparams.nretry = 5;
% sparams.nretain = 50;
% sparams.nsamples = 10;
% sparams.burnin = 200;
% sparams.thinning = 10;

sparams.nretry = 5;
sparams.nsamples = 10;
sparams.burnin = 400;
sparams.thinning = 10;

sparams.Qdet = diag([16 32 64]) * (1/4)^2;
sparams.Qktrack = diag([9 18 36]) * (1/3)^2;

sparams.nperv = 5;
sparams.ncamv = 5;

% perturbation covariance
sparams.Aper = eye(5) + [0,0,dt,0,0; 0,0,0,dt,0; zeros(3, 5)];
sparams.Qper1 = diag([1e-3, 1e-3, 0.04, 0.04, 0.0001]);
sparams.Qper2 = diag([1e-6, 1e-6, 0.01, 0.01, 0.000001]);
sparams.Acam = eye(5);
% sparams.Qcam = diag([1e-10, 0.0001, 1e-10, 1, (pi/180*1/6)^2]);
sparams.Qcam = diag([1e-10, 0.0001, 1e-10, 4, (pi/180*1)^2]);

sparams.camPert = sparams.Qcam / 8; %diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;

sparams.perPert0 = diag([0.01, 0.01, 0, 0, 0.01/4]);
sparams.perPert1 = sparams.Qper1 / 16;
sparams.perPert2 = sparams.Qper2 / 16; % diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]) / 16;

% sparams_cam.asr = 100;
%%%% mstrack parameters
nBit = 4;
szkernel = logspace(-log(2)/log(10), log(2)/log(10), 21) * 128;
for s = 1:length(szkernel)
  kernels(s) = buildKernel( szkernel(s)/2, szkernel(s));
end
simth = 0.9;
%%%%
% weightmask = fspecial('gaussian',[64 64], 21);
% weightmask = fspecial('gaussian',[32 64], 10);
% weightmask = imresize(weightmask, [128, 64]);
% weightmask = weightmask(:) * 1000; % prevent error due to small value max(w) * 1000 = 0.4...
% edges = 32 * (0:8);

detth = 0;
cparams.appth = -log(.55); % 0.7
cparams.ovth = 10; %-log();

cparams.nperv = 5;
cparams.ncamv = 5;

cparams.alpha = 1;
cparams.beta = 2;
cparams.adf = 0.5;
cparams.notrial = 30;

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

Z.nSamples = 200;
Z.cam = diag([50 0.1 1e-5 15 1e-8]) * randn(5, Z.nSamples) + repmat([550 1.7 360 200 0]', 1, Z.nSamples);
Z.per = zeros(0, Z.nSamples);
Z.peridx = [];
Z.tcnt = [];
Z.model = [];
Z.nTargets = 0;

% clear Theta, sparams_cam;

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
        [rid, cid]=find(omat > 0.5);
        
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

    Y = zeros(3, length(Z.model));
    for j = 1:length(Z.model)
        % get candidate scales
        cands = find(szkernel >= (Z.model(j).ppos(4) / 1.2) & szkernel <= (Z.model(j).ppos(4) * 1.2));
        bestSim = 0;        
        for s = cands
            [p, pos, Ic, sim] = kernelTrack(Im, Z.model(j).qS, Z.model(j).ppos(1:2)', kernels(s), nBit);
            if sim > bestSim, best={p, pos, Ic, s}; bestSim = sim; end
        end
        % super care needed!
        if bestSim > simth
            Y(1:2, j) = best{2};
            Y(3, j) = szkernel(best{4});
%             Z.model(j).ppos = [Y(1, j) - Y(3, j)/4, Y(2, j) - Y(3, j)/2, Y(3, j)/2, Y(3, j)]';% [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];
        end
    end
    [corres] = getCorrespondence2(Z, X, Y, cparams);

    %%%%% show corrspondences....debug;
    if debug == 1
        figure(3)
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
    
    tic;
%     [Z, zcorres] = MCMCSamplesJointStates(Z, X, Y, zcorres, sparams);
    if i == 1
        tKLT = [];
    else
        KLTidx = getValidKLT(KLT, i, Z, isz(2), sparams.nperv);
        
%         KLTidx = KLTidx;
        tKLT.x = KLT.x(KLTidx, (i-1):i);
        tKLT.y = KLT.y(KLTidx, (i-1):i);
    end
    [Z, zcorres] = MCMCSamplesJointStates2(Z, X, Y, tKLT, zcorres, sparams);
%     sum(Z.cam, 2) / Z.nSamples'
    toc;
    %%% update model ?
    for j = 1:length(corres)
        if corres(j) == 0, continue; end
        
        patch = getImgPatch(Im, X.obs(:, j)');
        patch = uint8(imresize(patch, [128 64]));
        Qc = bitshift(reshape(patch, [], 3), nBit-8);
        
        Z.model(corres(j)).timg = patch;
        Z.model(corres(j)).qS = (1 - cparams.adf) * buildHist( Qc, kernels(11), nBit ) + cparams.adf * Z.model(corres(j)).qS;
        % update this by Z estimation..
        Z.model(corres(j)).ppos = [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];
    end

    newtracks = find(corres == 0);
    for j = newtracks
        patch = getImgPatch(Im, X.obs(:, j)');
        patch = uint8(imresize(patch, [128 64]));
        Qc = bitshift(reshape(patch, [], 3), nBit-8);
        model.timg = patch;
        model.qS = buildHist( Qc, kernels(11), nBit );
        model.ppos = [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];        
        Z.model = [Z.model, model];
        
        Z.peridx = [Z.peridx, targetcnt];
        Z.nTargets = Z.nTargets + 1;
        targetcnt = targetcnt + 1;
    end
    
    % filter out obvious miss tracks
    filterout = [];
    for j = 1:Z.nTargets
        if sum(Z.per((sparams.nperv + 1) * j, :)) < Z.nSamples * 0.05
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
        if sum(hTracks(j).det > i - 5) >= 3
            hTracks(j).id = trackcnt;
            Tracks = [Tracks, hTracks(j)];
            % remove from hypothesis
            removelist = [removelist, j];        
            trackcnt = trackcnt + 1;
        end
%         end
        if (hTracks(j).det(end) < i - nTrackTerm)
            zidx = find(Z.peridx == hTracks(j).tid);
            targetidx = unique(Z.peridx);
            
            Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), :) = [];
            mid = find(targetidx == hTracks(j).tid);
            Z.model(mid) = [];
            Z.peridx(mid) = [];
            Z.tcnt(mid) = [];
            Z.nTargets = Z.nTargets - 1;
                        
            removelist = [removelist, j];
        end
    end
    hTracks(removelist) = [];
    for j = 1:length(Tracks)
         if (Tracks(j).term == -1) && (Tracks(j).det(end) < i - nTrackTerm)
            zidx = find(Z.peridx == Tracks(j).tid);
            targetidx = unique(Z.peridx);
            
            Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), :) = [];
            mid = find(targetidx == Tracks(j).tid);
            Z.model(mid) = [];
            Z.peridx(mid) = [];
            Z.tcnt(mid) = [];
            Z.nTargets = Z.nTargets - 1;
            
            Tracks(j).term = i;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    targetidx = unique(Z.peridx);
    figure(1);
    imshow(Im);
    figure(2);
    clf;
    xbound = 15;
    ccc = 'rgbck';
    lst = {'-', '-', ':', '-.'};
    
    for j = 1:5
        sidx = ceil(rand * Z.nSamples);

        figure(1);
        rectangle('Position', [1, Z.cam(4, sidx), isz(2), 1], 'LineStyle', ':', 'EdgeColor', 'r');
        figure(2);
        hold on;
        plot([0 xbound], [0 2 * Z.cam(1, sidx) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        plot([0 -xbound], [0 2 * Z.cam(1, sidx) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        hold off;

        for k = targetidx(:)'
            trackid = 0;
            for l = 1:length(Tracks)
                if (Tracks(l).term == -1) && (k == Tracks(l).tid)
                    trackid = l;
                    break;
                end
            end
            if trackid == 0, continue; end

            zidx = find(Z.peridx == k);
            if Z.per((sparams.nperv + 1) * zidx, sidx) == 0
                continue;
            end
            
            person = Z.per((1 + (zidx-1) * (sparams.nperv + 1)):(zidx * (sparams.nperv + 1)), sidx);
            [dummy, bb] = getProjection(person, Z.cam(:, sidx));
            figure(1);
            rectangle('Position', bb, 'LineWidth', 1, 'EdgeColor', ccc(mod(trackid, 5)+1), 'LineStyle',  lst{mod(ceil(trackid / 5), 4) + 1});

            hold on;
            quiver(bb(1) + bb(3)/2, bb(2) + bb(4)/2, person(3) * 30, -person(4) * 30, [ccc(mod(trackid, 5)+1)]);
            hold off;

            figure(2);
            hold on
            scatter(person(1), person(2), [ccc(mod(trackid, 5)+1), '.']);
            quiver(person(1), person(2), person(3), person(4), ccc(mod(trackid, 5) + 1));
            hold off
        end
    end
    
    figure(1)
    if i > 1
        hold on
        scatter(KLT.x(KLTidx, i), KLT.y(KLTidx, i), 'b.');
        hold off
    end
    
    for j = 1:size(det, 1)
        rectangle('Position', det(j, 1:4), 'LineStyle', ':', 'LineWidth', 1, 'EdgeColor', 'k');
    end
    title(['frame' num2str(i)]);
    drawnow;
    F(i) = getframe;
    figure(2)
    axis([-xbound xbound 0 2*xbound*3/4]);
%     axis equal
    grid on;
    drawnow;

    tF2 = getframe;
    stf1 = size(F(i).cdata);
    tF2.cdata = imresize(tF2.cdata, [stf1(1) stf1(1)]);
    stf2 = size(tF2.cdata);
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 1) = [tF2.cdata(:,:,1); zeros(stf1(1) - stf2(1), stf2(2))];
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 2) = [tF2.cdata(:,:,2); zeros(stf1(1) - stf2(1), stf2(2))];
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)+stf2(2), 3) = [tF2.cdata(:,:,3); zeros(stf1(1) - stf2(1), stf2(2))];
end

close all;

% movie2avi(F, 'test_mstrack_100s.avi');