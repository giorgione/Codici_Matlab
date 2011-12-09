clear;
addpath(genpath('common'))

load testinit_uniform.mat

debug = 0; %1;

nTrackTerm = 10;

initN = 100;
sparams_per.numPerTarget = initN;
sparams_cam.numPerTarget = initN;

fstep = 1;
dt = fstep/30;
sparams_per.A = eye(5) + [0,0,dt,0,0;0,0,0,dt,0;zeros(3, 5)];
sparams_per.R = diag([0, 0, 0.01, 0.01, 0.000001]);
sparams_per.asr = 20;
sparams_per.rKernel = 1;

sparams_cam.A = eye(5);
sparams_cam.R = diag([25, 0.0001, 4, 4, (pi/180*5/6)^2]);
Theta.state = [Theta.state; zeros(1, size(Theta.state, 2))]; % add 0 panning
sparams_cam.lenState = 5;
% sparams_cam.asr = 100;

%%%% mstrack parameters
nBit = 4;
szkernel = logspace(-log(2)/log(10), log(2)/log(10), 21) * 128;
for s = 1:length(szkernel)
  kernels(s) = buildKernel( szkernel(s)/2, szkernel(s));
end
simth = 0.9;
%%%%

Z.idx = [];
Z.w = [];
Z.state = [];
Z.pest = [];
Z.model = [];

% weightmask = fspecial('gaussian',[64 64], 21);
% weightmask = fspecial('gaussian',[32 64], 10);
% weightmask = imresize(weightmask, [128, 64]);
% weightmask = weightmask(:) * 1000; % prevent error due to small value max(w) * 1000 = 0.4...
% edges = 32 * (0:8);

detth = 0;
wparams_per.Qdet = diag([16 64 256]);
wparams_per.Qktrack = diag([9 36 144]);

cparams.appth = -log(.55); % 0.7
cparams.ovth = 10; %-log();

cparams.alpha = 1;
cparams.beta = 2;
cparams.adf = 0.5;
cparams.notrial = 30;

Tracks = [];
hTracks = [];

targetcnt = 1;
trackcnt = 1;
% Theta.idx = ones(1, initN);
% Theta.w = ones(1, initN) / initN;
% % Theta.state = [randn(1, initN) * 50 + 400; randn(1, initN) * 0.1 + 1.5; isz(2) / 2 * ones(1, initN); isz(1)/2 + 30 * randn(1, initN)];
% Theta.state = [randn(1, initN) * 25 + 900; randn(1, initN) * 0.0025 + 1.45; isz(2) / 2 * ones(1, initN); 210 + 2 * randn(1, initN)];
for i = 1:fstep:length(imfiles)
    Im = imread([imgdir imfiles(i).name]);
    fileNameNoExt=imfiles(i).name(1:find(imfiles(i).name=='.',1,'last')-1 );
    if exist([imgdir fileNameNoExt ext])
        det = load([imgdir fileNameNoExt ext]);
        det(1, :) = [];
        det(det(:, 6) < detth, :) = [];

        a = cumsum(Theta.w);
        filterout = [];
        for j = 1:size(det, 1)
            ImObs = det(j, :);
            ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
            b = [];
            for k = 1:100
                samp = sum(a < rand * a(end)) + 1;
                b = [b, getIProjection(ImObs + mvnrnd([0;0;0], wparams_per.Q)', Theta.state(:, samp))];
            end

            % filter out errornous dets.
            if sum(b(5, :) > 2.1 | b(5, :) < 1.3) > 70
                filterout = [filterout, j];
            end
        end
        det(filterout, :) = [];
        
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
        end
    else
        error;
    end
    unique(Z.idx)
    Z = SampleIndependentParticles(Z, sparams_per);
    Theta = SampleIndependentParticles(Theta, sparams_cam);
    unique(Z.idx)
    %%%%%%
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
    [corres] = getCorrespondence2(Z, Theta, X, Y, Im, cparams);

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
        
%         Z.model(corres(j)) = model;
    end
    %%% get correspondence
%     corres = [0 0 0 0];
    %%%%%%
    [ Z, Theta ] = ComputeWeight(Z, Theta, X, Y, corres, wparams_per);
    newtracks = find(corres == 0);
    for j = newtracks
        ImObs = X.obs(:, j);
        ImObs = [ImObs(1) + ImObs(3) / 2; ImObs(2) + ImObs(4); ImObs(4)];
        
        % sample noise / cam parameters
        % get inverse projections
        a = cumsum(Theta.w);
        for k = 1:initN
            samp = sum(a < rand * a(end)) + 1;
            Z.state = [Z.state, getIProjection(ImObs + mvnrnd([0;0;0], wparams_per.Q)', Theta.state(:, samp)) + [0;0;mvnrnd([0;0], diag([.1, .1]))';0]];
        end
        
        Z.state(5, (end - initN + 1):end) = 1.7 + 0.1 * randn(1, initN);
        Z.idx = [Z.idx, targetcnt * ones(1, initN)];
        Z.w = [Z.w, ones(1, initN) / initN];
        
        patch = getImgPatch(Im, X.obs(:, j)');
        patch = uint8(imresize(patch, [128 64]));
        Qc = bitshift(reshape(patch, [], 3), nBit-8);
        model.timg = patch;
        model.qS = buildHist( Qc, kernels(11), nBit );
        model.ppos = [X.obs(1:2, j) + X.obs(3:4, j) / 2; X.obs(3:4, j)];        
        Z.model = [Z.model, model];
        targetcnt = targetcnt + 1;
    end
    % initiate and terminate tracks
    targetidx = unique(Z.idx);
    targetidx = targetidx(:)';
    detidx = targetidx(corres(corres > 0));
    newdets = sum(corres == 0);
    detidx = [detidx, targetidx((end - newdets + 1):end)];
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
            zidx = find(Z.idx == hTracks(j).tid);
            targetidx = unique(Z.idx);
            
            Z.w(zidx) = [];
            Z.idx(zidx) = [];
            Z.state(:, zidx) = [];
%             Z.model(:, zidx(end) / sparams_per.numPerTarget) = [];
            mid = find(targetidx == hTracks(j).tid);
            Z.model(mid) = [];
            
            removelist = [removelist, j];
        end
    end
    hTracks(removelist) = [];
    for j = 1:length(Tracks)
         if (Tracks(j).term == -1) && (Tracks(j).det(end) < i - nTrackTerm)
            zidx = find(Z.idx == Tracks(j).tid);
            targetidx = unique(Z.idx);
            
            Z.w(zidx) = [];
            Z.idx(zidx) = [];
            Z.state(:, zidx) = [];
            
%             Z.model(:, zidx(end) / sparams_per.numPerTarget) = [];
            mid = find(targetidx == Tracks(j).tid);
            Z.model(mid) = [];
            
            Tracks(j).term = i;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    targetidx = unique(Z.idx);
    figure(1)
    imshow(Im);
    figure(2)
    clf;
    a = cumsum(Theta.w);
    
    xbound = 15;
    for j = 1:3
        samp = sum(a < rand * a(end)) + 1;
        figure(1);
        rectangle('Position', [1, Theta.state(4, samp), isz(2), 1], 'LineStyle', ':', 'EdgeColor', 'r');
        figure(2);
        hold on
        plot([0 xbound], [0 2 * Theta.state(1, samp) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        plot([0 -xbound], [0 2 * Theta.state(1, samp) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        hold off;
    end
    
    ccc = 'rgbck';
    lst = {'-', '-', ':', '-.'};
    for j = targetidx(:)'
        trackid = 0;
        for k = 1:length(Tracks)
            if (Tracks(k).term == -1) && (j == Tracks(k).tid)
                trackid = k;
                break;
            end
        end
        
        if trackid == 0, continue; end
        
        zidx = find(Z.idx == j);
        
        b = cumsum(Z.w(zidx));
        for k = 1:5
            csamp = sum(a < rand * a(end)) + 1;
            zsamp = sum(b < rand * b(end)) + 1;
            
            [dummy, bb] = getProjection(Z.state(:, zidx(zsamp)), Theta.state(:, csamp));
            figure(1);
            rectangle('Position', bb, 'LineWidth', 1, 'EdgeColor', ccc(mod(trackid, 5)+1), 'LineStyle',  lst{mod(ceil(trackid / 5), 4) + 1});

            hold on
            quiver(bb(1) + bb(3)/2, bb(2) + bb(4)/2, Z.state(3, zidx(zsamp)) * 30, -Z.state(4, zidx(zsamp)) * 30, [ccc(mod(trackid, 5)+1)]);
            hold off
            
            figure(2);
            hold on
            scatter(Z.state(1, zidx(zsamp)), Z.state(2, zidx(zsamp)), [ccc(mod(trackid, 5)+1), '.']);
            quiver(Z.state(1, zidx(zsamp)), Z.state(2, zidx(zsamp)), Z.state(3, zidx(zsamp)), Z.state(4, zidx(zsamp)), ccc(mod(trackid, 5) + 1));
            hold off
        end
    end
    figure(1)
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
