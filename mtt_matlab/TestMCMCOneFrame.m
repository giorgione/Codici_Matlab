clear;
addpath(genpath('common'))

load testinit_uniform.mat

Theta.state = [Theta.state; zeros(1, size(Theta.state, 2))]; % add 0 panning

debug = 0; %1;

nTrackTerm = 10;

fstep = 1;
dt = fstep/15;

sparams.nretry = 10;
sparams.nretain = 50;
sparams.nsamples = 5;
sparams.burnin = 100;
sparams.thinning = 10;

sparams.Qdet = diag([16 32 64]) * (1/4)^2;
sparams.Qktrack = diag([9 18 36]) * (1/3)^2;

sparams.nperv = 5;
sparams.ncamv = 5;

% perturbation covariance
sparams.Aper = eye(5) + [0,0,dt,0,0; 0,0,0,dt,0; zeros(3, 5)];
sparams.Qper = diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]);
sparams.Acam = eye(5);
sparams.Qcam = diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]);

sparams.camPert = diag([4, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;
sparams.perPert = diag([0.01, 0.01, 0, 0, 0.01/4]); % diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]) / 16;

detth =  -1.0;

targetcnt = 1;
trackcnt = 1;

clear Z;

sparams_cam.A = eye(5);
sparams_cam.R = diag([25, 0.0001, 4, 4, (pi/180*5/6)^2]);
sparams_cam.lenState = 5;

sparams_cam.numPerTarget = 200; % sparams.nretain;
Theta = SampleIndependentParticles(Theta, sparams_cam);

Z.nSamples = size(Theta.state, 2);
Z.cam = Theta.state;
Z.per = zeros(0, Z.nSamples);
Z.peridx = [];
Z.model = [];
Z.nTargets = 0;

clear Theta, sparams_cam;
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

    Z.per = zeros(0, Z.nSamples);
    tic;
    [Z] = MCMCSamplesOneFrame(Z, X, sparams);
    toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     targetidx = unique(Z.peridx);
    figure(1)
    imshow(Im);
    figure(2)
    clf;
    xbound = 15;
    ccc = 'rgbck';
    lst = {'-', '-', ':', '-.'};
    
    for j = 1:5
        sidx = ceil(rand * Z.nSamples);
        
        figure(1);
        rectangle('Position', [1, Z.cam(4, sidx), isz(2), 1], 'LineStyle', ':', 'EdgeColor', 'r');
        figure(2);
        hold on
        plot([0 xbound], [0 2 * Z.cam(1, sidx) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        plot([0 -xbound], [0 2 * Z.cam(1, sidx) / isz(2) * xbound], 'r:', 'LineWidth', 2);
        hold off;
        
        for k = 1:size(X.obs, 2)
            if Z.per((sparams.nperv + 1) * k, sidx) == 0
                continue;
            end
            person = Z.per((1 + (k - 1) * (sparams.nperv + 1)):(k * (sparams.nperv + 1)), sidx);
            [dummy, bb] = getProjection(person, Z.cam(:, sidx));
            figure(1);
            rectangle('Position', bb, 'LineWidth', 1, 'EdgeColor', ccc(mod(k, 5)+1), 'LineStyle',  lst{mod(ceil(k / 5), 4) + 1});

            hold on
            quiver(bb(1) + bb(3)/2, bb(2) + bb(4)/2, person(3) * 30, -person(4) * 30, [ccc(mod(k, 5)+1)]);
            hold off

            figure(2);
            hold on
            scatter(person(1), person(2), [ccc(mod(k, 5)+1), '.']);
            quiver(person(1), person(2), person(3), person(4), ccc(mod(k, 5) + 1));
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
    
    
%     retainidx = ceil(rand(1, sparams.nretain) * Z.nSamples);
%     Z.nSamples = sparams.nretain;
%     Z.cam = Z.cam(:, retainidx);
%     Z.per = Z.per(:, retainidx);
%     tic;
end

close all;

% movie2avi(F, 'test_mstrack_100s.avi');