clear;
addpath(genpath('common'))

load testinit_uniform.mat

initN = 500;
sparams_per.numPerTarget = 500;
sparams_cam.numPerTarget = 500;

Z.idx = [];
Z.w = [];
Z.state = [];
Z.pest = [];
Z.model = [];

% weightmask = fspecial('gaussian',[64 64], 21);
weightmask = fspecial('gaussian',[32 64], 10);
weightmask = imresize(weightmask, [128, 64]);
weightmask = weightmask(:) * 1000; % prevent error due to small value max(w) * 1000 = 0.4...
edges = 32 * (0:8);

detth = 0;
wparams_per.Q = diag([9 36 144]);

dt = 1/30;
sparams_per.A = eye(5) + [0,0,dt,0,0;0,0,0,dt,0;zeros(3, 5)];
sparams_per.R = diag([0, 0, 0.01, 0.01, 0.000001]);

sparams_cam.R = diag([25, 0.0001, 4, 4]);

cparams.appth = 0.8;
cparams.ovth = 1.5;

cparams.alpha = 1;
cparams.beta = 1;
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


for i = 1:length(imfiles)
    Im = imread([imgdir imfiles(i).name]);
    fileNameNoExt=imfiles(i).name(1:find(imfiles(i).name=='.',1,'last')-1 );
    
    if exist([imgdir fileNameNoExt ext])
        det = load([imgdir fileNameNoExt ext]);
        det(1, :) = [];
        det(det(:, 6) < detth, :) = [];
        
%         det(2:end, :) = [];
        % temp
%         det = det([1 2 4], :);
        %%%%%%%%%%%%%%%
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

            % filter out error dets.
            if sum(b(5, :) > 2.1 | b(5, :) < 1.3) > 70
                filterout = [filterout, j];
            end
        end
        det(filterout, :) = [];
        
        X.obs = [];
        X.idx = [];
        
        X.idx = 1:size(det, 1);
        for j = 1:size(det, 1)
            X.obs(:, j) = det(j, 1:4)';
        end
    else
        error;
    end
    
    Z = SampleIndependentParticles(Z, sparams_per);
    Theta = SampleIndependentParticles(Theta, sparams_cam);
    % get estimation %%%%%%%%%%
%     [msid, msm] = meanshift(Theta.state', 15);
%     Theta.pest.loc = msm(1, :);
%     
%     dcnt = 1;
%     Z.pest = [];

%     for j = targetidx(:)'
%         zidx = find(Z.idx == j);
%         
%         [msid, msm] = meanshift(Z.state(:, zidx)', 0.3);
%         Z.pest(dcnt).loc = msm(1, :);
%         
%         dcnt = dcnt + 1;
%     end
%     if i == 65
%         keyboard
%     end
    
    tic;
    [corres] = getCorrespondence(Z, Theta, X, Im, cparams);
    corres
%     if i == 60
%         keyboard;
%     end
    
    [ Z, Theta ] = ComputeWeight(Z, Theta, X, corres, wparams_per);
    
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
        
        % filter out error dets.
%         if sum(Z.state(5, (end - initN + 1):end) > 2.1 | Z.state(5, (end - initN + 1):end) < 1.3) > initN / 2
%             Z.state(:, (end - initN + 1):end) = [];
%             corres(j) = -1;
%             continue;
%         else
            Z.state(5, (end - initN + 1):end) = 1.7 + 0.1 * randn(1, initN);
%         end
        
        Z.idx = [Z.idx, targetcnt * ones(1, initN)];
        Z.w = [Z.w, ones(1, initN) / initN];
%         Z.state = [Z.state, [20 * rand(1, initN) - 10; 50 * rand(1, initN); mvnrnd([0;0], [1, 0; 0, 1], initN)'; 1.7 + 0.1 * randn(1, initN)]];
        
%         patch = Im(X.obs(2, j):(X.obs(2, j) + X.obs(4, j) - 1), X.obs(1, j):(X.obs(1, j) + X.obs(3, j) - 1), :);
        patch = getImgPatch(Im, X.obs(:, j)');
        
        patch = imresize(patch, [128 64]);
        patch = patch(:); patch = [patch(1:128*64), patch((128*64+1):(128*64*2)), patch((128*64*2+1):(128*64*3))];
        patch(patch < 0 | patch > 255) = 0;
        h = histc_nD(patch, edges, weightmask); h=h(:);
        Z.model(:, end + 1) = h;
        
        targetcnt = targetcnt + 1;
    end
%     Theta.pest.loc = sum(Theta.state, 2) / size(Theta.state, 2);
%     targetidx = unique(Z.idx);
%     dcnt = 1;
%     for j = targetidx(:)'
%         zidx = find(Z.idx == j);
%         Z.pest(dcnt).loc = sum(Z.state(:, zidx), 2) / length(zidx);
%         dcnt = dcnt + 1;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if i == 9
%         keyboard;
%     end
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
        if length(hTracks(j).det) >= 5
            if sum(hTracks(j).det > i - 5) >= 3
                hTracks(j).id = trackcnt;
                Tracks = [Tracks, hTracks(j)];
                % remove from hypothesis
                removelist = [removelist, j];        
                trackcnt = trackcnt + 1;
            end
        end
        
        if hTracks(j).det(end) < i - 4
            zidx = find(Z.idx == hTracks(j).tid);
            
            Z.w(zidx) = [];
            Z.idx(zidx) = [];
            Z.state(:, zidx) = [];
            Z.model(:, zidx(end) / sparams_per.numPerTarget) = [];
            
            removelist = [removelist, j];
        end
    end
    hTracks(removelist) = [];
    
    for j = 1:length(Tracks)
         if Tracks(j).term == -1 & Tracks(j).det(end) < i - 4
            zidx = find(Z.idx == Tracks(j).tid);
            
            Z.w(zidx) = [];
            Z.idx(zidx) = [];
            Z.state(:, zidx) = [];
            Z.model(:, zidx(end) / sparams_per.numPerTarget) = [];
            
            Tracks(j).term = i;
        end
    end
   
%     Tracks
%     hTracks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    targetidx = unique(Z.idx);
    figure(1)
    imshow(Im);
    figure(2)
    clf;
    a = cumsum(Theta.w);
    
    xbound = 7;
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
        isvalid = 0;
        for k = 1:length(Tracks)
            if (Tracks(k).term == -1) && (j == Tracks(k).tid)
                isvalid = 1;
                break;
            end
        end
        
        if isvalid == 0, continue; end
        
        zidx = find(Z.idx == j);
        
        b = cumsum(Z.w(zidx));
        for k = 1:5
            csamp = sum(a < rand * a(end)) + 1;
            zsamp = sum(b < rand * b(end)) + 1;
            
            [dummy, bb] = getProjection(Z.state(:, zidx(zsamp)), Theta.state(:, csamp));
            figure(1);
            rectangle('Position', bb, 'LineWidth', 1, 'EdgeColor', ccc(mod(j, 5)+1), 'LineStyle',  lst{mod(ceil(j / 5), 4) + 1});

            hold on
            quiver(bb(1) + bb(3)/2, bb(2) + bb(4)/2, Z.state(3, zidx(zsamp)) * 30, -Z.state(4, zidx(zsamp)) * 30, [ccc(mod(j, 5)+1)]);
            hold off
            
            figure(2);
            hold on
            scatter(Z.state(1, zidx(zsamp)), Z.state(2, zidx(zsamp)), [ccc(mod(j, 5)+1), '.']);
            quiver(Z.state(1, zidx(zsamp)), Z.state(2, zidx(zsamp)), Z.state(3, zidx(zsamp)), Z.state(4, zidx(zsamp)), ccc(mod(j, 5) + 1));
            hold off
        end
        
%         [msid, msm] = meanshift(bb', 10);
%         rectangle('Position', msm(1, :), 'LineWidth', 1, 'EdgeColor', ccc(j));
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
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 1) = [tF2.cdata(:,:,1), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 2) = [tF2.cdata(:,:,2), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
    F(i).cdata(1:stf1(1), stf1(2)+1:stf1(2)*2, 3) = [tF2.cdata(:,:,3), zeros(stf2(1), stf1(2) - stf2(2)); zeros(stf1(1) - stf2(1), stf1(2))];
%     
%     figure(3)
%     scatter3(Theta.state(1, :), Theta.state(2, :), Theta.state(4, :), '.'); drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


