clear
addpath(genpath('common'))
addpath(genpath('mexFunctions'))

% imgdir = '../Dataset/test4/';
% ext = '_det06th.txt';
% fstep = 2;
% framerate = 30;
% Z.nSamples = 1;
% Z.cam = [500 1.5 360 230 0 0 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-50 20 1e-50 .2 1e-10 1e-10]).^2;
% Z.W = 1;
% sparams = InitParams(1, fstep, framerate);
% params.mpcam = 230;
% params.vpcam = 10^2;

% imgdir = '../test17/';
% ext = '_det06th.txt';
% fstep = 2;
% framerate = 30;
% Z.nSamples = 1;
% Z.cam = [500 1.5 360 200 0 0 0 0]';
% % Z.cam = [500 1.6 360 210 0 0 0 0]';
% Z.V{1, 1} = diag([10 0.2 1e-50 20 1e-50 .2 1e-10 1e-10]).^2;
% Z.W = 1;
% sparams = InitParams(1, fstep, framerate);
% 
imgdir = '../Dataset/seq02/';
ext = '_det06th.txt';
fstep = 1;
framerate = 14;
Z.nSamples = 1;
Z.cam = [500 0.94 320 220 0 1.3 0 0]';
Z.V{1, 1} = diag([1 0.01 1e-6 10 1e-6 .4 1e-6 1e-6]).^2;
Z.W = 1;
sparams = InitParams(3, fstep, framerate);

% 
% imgdir = '../Dataset/seq03/';
% imfiles = dir([imgdir '*.png']);
% ext = '_det06th.txt';
% fstep = 1;
% framerate = 14;
% Z.nSamples = 1;
% Z.cam = [500 1.1 320 220 0 .5 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-6 5 1e-6 .2 1e-6 1e-6]).^2;
% Z.W = 1;
% sparams = InitParams(3, fstep, framerate);

% imgdir = '../Dataset/ETH/seq04/';
% imfiles = dir([imgdir '*.png']);
% ext = '_det06.txt';
% fstep = 1;
% framerate = 15;
% Z.nSamples = 1;
% Z.cam = [500 0.94 320 235 0 1.2 0 0]';
% Z.V{1, 1} = diag([10 0.1 1e-6 10 1e-6 .1 1e-6 1e-6]).^2;
% Z.W = 1;

Z.per = zeros(0, Z.nSamples); Z.gfeat = zeros(0, Z.nSamples); Z.gfidx = [];
Z.peridx = []; Z.tcnt = []; Z.gfcnt = []; Z.model = []; Z.nTargets = 0;
Z.nFeats = 0; Z.beta = [];

sparams.camscnt = 6;
sparams.perscnt = 4;
sparams.featscnt = 2.5;
% sparams.nodraw = 1;

[KLT.x, KLT.y, KLT.val] = klt_read_featuretable([imgdir, 'features.txt']);
KLT = preprocessKLT(KLT);

%% 
for th = -.5
    sparams.detth = th;
    [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
    for i = 1:length(Zs)
        Zs(i).model = [];
        Zs(i).cams = [];
    end
    save(['seq02_res_th_' num2str(th, '%.02f') '.mat'], 'Zs', 'Tracks', 'TP', 'FP', 'FN', 'sparams',  'imgdir',  'framerate', 'KLTused');
end
% 
% sparams.detth = -0.1;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n1.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused
% 
% sparams.detth = -0.2;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n2.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused
% 
% 
% sparams.detth = -0.3;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n3.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused
% 
% sparams.detth = -0.4;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n4.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused
% 
% sparams.detth = -0.5;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n5.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused
% 
% sparams.detth = -0.6;
% [Tracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, KLT, ext, fstep, framerate, Z, sparams);
% for i = 1:length(Zs)
%     Zs(i).model = [];
%     Zs(i).cams = [];
% end
% save seq02_res_th_n65.mat Zs Tracks TP FP FN sparams imgdir framerate KLTused




