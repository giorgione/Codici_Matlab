% function TrackTop(imgdir, detdir, option, ihorizon, prefix, thlist)
% Main function for tracking a video.
%
% - imgdir: directory of images extracted from the video sequence.
%
% - detdir: directory of detection results for each frame of the video
%          sequence
%
% - option: options
% 
% - ihorizon: Initial Orizon Line estimate
%
% - prefix: 
% 
% - thlist: THRESHOLD LIST --> list of possible value of the threshold for
%                              the DETECTION PROCESS
function TrackTop(imgdir, detdir, option, ihorizon, prefix, thlist)
% ALREADY ADDED in the path
%addpath(genpath('common'))
%addpath(genpath('mexFunctions'))

%image resolution
imsize = [2040 2040];
% ihorizon = 1052;

%%%%%%%%%%%%%%%%%%%%%%%%
%number of frames to skip
fstep = 2;
framerate = 14.9;
dt = fstep / framerate;

% use only 3 for now: Initialize Parameter for the MODEL
% sparams: struct containing MODEL PARAMETERS
sparams = InitParams(3, fstep, framerate);

% motion uncertainty in camera parameters: MATRIX 8x8
sparams.Qcam = diag((dt * [1e-6, 1e-6, 1e-6, 5/dt, 2*pi/180/dt, 10, 1e-6, 1e-6]).^2);
% perturbation matrix
sparams.camPert = sparams.Qcam / 256; %diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;
sparams.nfeatuse = 15;

% prior on camera's height
sparams.mpcam = ihorizon; % mean
sparams.vpcam = 3^2; % variance

ext = '_det06.txt';

sparams.useInteraction = 0;

% dimension of one car state. (x, z, vx, vz, h, w) + 1 (class)
sparams.ncarv = 6; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialize samples.
sparams.camscnt = 20;
sparams.perscnt = 10;
sparams.carscnt = 10;
sparams.featscnt = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STARTING PARAMETERS SETTING for the WHOLE MODEL

%Initial numbers of samples for the Camera Model
Z.nSamples = 1;

%% CAMERA initial paramters:  Variable 8x1 parameters
% focal length, camera height(meter), x camera center(pixel), horizon
% (pixel), initial angle, initial speed, camera's x location, camera's z location
Z.cam = [1200 1.0 imsize(1)/2 ihorizon 0 11 0 0]';
%Camera MATRIX PROCESS covariance
Z.V{1, 1} = diag([1 0.01 1e-6 10 1e-6 1 1e-6 1e-6]).^2;
% Initial weight
Z.W = 1;

%% PERSON DATA 
Z.per = zeros(0, Z.nSamples); 
Z.peridx = []; 

%% GROUND DATA
Z.gfeat = zeros(0, Z.nSamples); % ground features
Z.gfidx = [];                   % ground features IDX
Z.gfcnt = [];                   % ground features count

%% 
Z.model = []; 
%% 
Z.nTargets = 0;
Z.nFeats = 0; 
Z.nCarTargets = 0;

Z.beta = [];
Z.tccnt = []; 

%% CAR DATA
Z.car = zeros(0, Z.nSamples);
Z.caridx = []; 

% sparams.nodraw = 1;
max_frames = 500;

% KLT = sift_read_featuretable(imgdir, max_frames);
%load KLT features for the video
load('KLT500.mat');

for th = thlist
    %Fix the Treshold for the Detection
    sparams.detth = th;
    disp(['Begin Processing ' imgdir ' with threshold = ' num2str(th)]);
    
    [Tracks, CTracks, camTrack, Zs, TP, FP, FN, F, KLTused] = TrackOneVideo(imgdir, detdir, KLT, ext, fstep, framerate, Z, sparams, 1, max_frames);
    
    % don't need these later. Erase these..
    for i = 1:length(Zs)
        Zs(i).model = [];
        Zs(i).cams = [];
    end
    
    %save([prefix num2str(th, '%.02f') '_' num2str(option) '.mat'], 'Zs', 'Tracks', 'CTracks', 'TP', 'FP', 'FN', 'sparams',  'imgdir',  'framerate', 'KLTused');
end

