function params = InitParams(mode, fstep, framerate)
dt = fstep / framerate;
if mode == 1
    % MCMC sampling parameters
    params.nretry = 1;
    params.nsamples = 400;
    params.burnin = 1000;
    params.thinning = 200;
%     
%     params.mpcam = 220;
%     params.vpcam = 5^2;
    
    % MCMC proposal distribution.
    params.mPertCam = 64;
    params.mPertPer = 2;
    params.mPertGF = 1;

    params.mPertCam = 32;
    params.mPertPer = 4;
    params.mPertGF = 4;
    
    % observation noise.
    params.Qdet = diag([16 32 64]) * (2/4)^2;
    params.Qktrack = diag([9 18 36]) * (2/3)^2;

    % dimension of the states
    params.nperv = 5;
    params.ncamv = 8;
    params.ngfeat = 3;

    % perturbation covariance
    params.Aper = eye(5) + [0,0,dt,0,0; 0,0,0,dt,0; zeros(3, 5)];
    params.Qper1 = diag([1e-3, 1e-3, 3, 3, 0.01].^2);

    %params.Qper2 = diag((dt * [1e-1, 1e-1, .5, .5, 0.001/dt]).^2);
    params.Qper2 = diag((dt * [1e-1, 1e-1, 1, 1, 0.001/dt]).^2);
    
    params.Acam = eye(5);
    params.Qcam = diag((dt * [2, 0.0001, 1e-2, 100, 30*pi/180 1e-3 1e-3 1e-3]).^2);

    params.camPert = params.Qcam / 4; %diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;
    params.perPert0 = diag([0.01, 0.01, 0, 0, 0.01/4]);
    params.perPert1 = params.Qper1 / 4;
    params.perPert2 = params.Qper2 / 4; % diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]) / 16;
    
    params.mh = 1.7; % 1.8;
    params.sh = 0.1;
    params.normh = 1/sqrt(2 * pi * params.sh^2);
    
    params.sv = 0.5;
    params.normv = 1/sqrt(2 * pi * params.sv^2);
    
    params.noHprob = 0.1; % normpdf(2.6 * params.sh, 0, params.sh); % *normpdf(3 * params.sv, 0, params.sv);
    
%     params.hdet = 0.99;
%     params.nhdet = 0.00001; 
    params.hdet = 0.6;
    params.nhdet = 0.001; 

    params.flipObj = 0.1;
    
    % KLT parameters
    kltstd = 2;
    params.kltPrec = eye(2) / kltstd^2; % sx = 1/2, sy = 1/2; 
    params.lkltNorm = log(1/sqrt(2 * pi * kltstd^2));
    params.lkltuprob = log(normpdf(2*kltstd, 0, kltstd)^2); % at 2 sigma
    params.flipFeat = 0.2;
    %Minimum Number of GROUND FEATURES to use
    params.nfeatuse = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    params.useCamAnnealing = 1;
    params.AnnealingConst = 0.9999;
    params.KLTmargin = 45;
elseif mode == 2
    %%%%%%% moving camera
    % MCMC sampling parameters
    params.nretry = 5;
    params.nsamples = 200;
    params.burnin = 2000;
    params.thinning = 20;

    % MCMC proposal distribution.
    params.mPertCam = 49;
    params.mPertPer = 9;

    % observation noise.
    params.Qdet = diag([16 32 64]) * (1/4)^2;
    params.Qktrack = diag([9 18 36]) * (1/3)^2;

    % dimension of the states
    params.nperv = 5;
    params.ncamv = 8; % [f, yc, vc, uc, phi, vel, xc, zc]
    params.ngfeat = 3;
    params.dt = dt;
    
    % perturbation covariance
    params.Aper = eye(5) + [0,0,dt,0,0; 0,0,0,dt,0; zeros(3, 5)];
%     params.Qper1 = diag([1e-3, 1e-3, 1, 1, 0.01].^2);
    params.Qper2 = diag((dt * [1e-3, 1e-3, .5, .5, 0.01]).^2);
    params.Qper1 = diag([1e-3, 1e-3, 1, 1, 0.1].^2);
%     params.Qper2 = diag((dt * [1e-3, 1e-3, .5, .5, 0.1]).^2);
    
    params.Qcam = diag((dt * [2, 1e-3, 1e-4, 10, 6*pi/180, .5, 1e-3, 1e-3]).^2);
    params.camPert = params.Qcam / 16; %diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;
    params.perPert0 = diag([0.01, 0.01, 0, 0, 0.01/4]);
    params.perPert1 = params.Qper1 / 16;
    params.perPert2 = params.Qper2 / 16; % diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]) / 16;
    
    params.mh = 1.7;
    params.sh = 0.1;
    params.normh = 1/sqrt(2 * pi * params.sh^2);
    
    params.sv = 0.5;
    params.normv = 1/sqrt(2 * pi * params.sv^2);
    
    params.noHprob = normpdf(3 * params.sh, 0, params.sh)*normpdf(3 * params.sv, 0, params.sv);
    
    params.hdet = 0.6;
    params.nhdet = 0.1; 

    params.flipObj = 0.1;
    
    % KLT parameters
    params.kltPrec = eye(2) * 2; % sx = 1/2, sy = 1/2; 
%     params.kltuprob = 0.1353; % at 2 sigma : normpdf(2) / normpdf(0)
    params.klttrunc = 3^2; % outside of 3 sigma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    params.useInteraction = 0;
    %%%%%% interaction parameters
    params.c1_pair = 0.3;
    params.c2_pair = 0.3;

    params.i1_pair = 3;   % slope of the soft step
    params.i2_pair = 1.5; % threshold of the soft step

    params.lr1_pair = 1; % 
    params.lr2_pair = 4; % 

    params.lg1_pair = 3;  %
    params.lg2_pair = 32; %
    
    params.betatrns = 0.985 * eye(3) + 0.005 * ones(3,3); % transition of interactions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Margin in width for the KLT FEATURES extracted from FRAMES:
    %    sparams.KLTmargin < x < imwidth - sparams.KLTmargin
    %    y is > mcam(4) <--> horizon line + 30 pixel
    params.KLTmargin = 20;
elseif mode == 3
    % temp
    %%%%%%% moving camera
    % MCMC sampling parameters
    
    %Number of Iteration of Sampling 
    params.nretry = 10;
    %Numbers of generated samples
    params.nsamples = 150;
    %Burning Period
    params.burnin = 5000; % 
    params.thinning = 100;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Z.cam = [500 0.94 320 220 0 1.3 0 0]';
%     Z.V{1, 1} = diag([1 0.01 1e-50 20 1e-6 .4 1e-6 1e-6]).^2;
    params.mpcam = 220;
    params.vpcam = 5^2;
    
    % MCMC proposal distribution.
    params.mPertCam = 100;    % Perturbation CAMERA
    params.mPertPer = 16;     % Perturbation PERSON
    params.mPertCar = 25*4;   % Perturbation CAR
    params.mPertGF = 4;       % Perturbation GROUND

    % Not Usefull
    params.useCamAnnealing = 0;%1;
    params.AnnealingConst = 0.995;
    
    % observation noise: MATRIX 3x3
    params.Qdet = diag([16 32 64]) * (2/4)^2;
    %params.Qcardet = diag([6^2 4^2 10^2 4^2]) * (1/4)^2;
    %params.Qcardet = diag([4^2 4^2 6^2 4^2]) * (1/4)^2;
   
    %Car detector noise: MATRIX 4x4
    params.Qcardet = diag([4^2 4^2 4^2 4^2]) * (1/4)^2;
    %Track noise: MATRIX 3x3
    params.Qktrack = diag([9 18 36]) * (2/3)^2;

    %% STATES VARIABLE DIMENSION:
    %Vector 5D Describing the PERSON - TARGET VARIABLE
    params.nperv = 5;
    %Vector 6D Describing the CAR - CAR VARIABLE
    params.ncarv = 6;
    %Vector 8D Describing the CAMERA:[f, yc, vc, uc, phi, vel, xc, zc]
    params.ncamv = 8; 
    %Vector 3D Describing the PERSON
    params.ngfeat = 3;
    params.dt = dt;
    
    %% COVARIANCE for Random Process
    % perturbation covariance for PERSON: 5x5 Matrix
    params.Aper = eye(5) + [0,0,dt,0,0;
                            0,0,0,dt,0; 
                            zeros(3, 5)];

    % perturbation covariance for PERSON: 6x6 Matrix   
    params.Qper2 = diag((dt * [1e-3, 1e-3, .55, .55, 0.04/dt]).^2);
    %params.Qper2 = diag((dt * [1e-3, 1e-3, .2, .2, 0.01]).^2);
    %params.Qper1 = diag([1e-3, 1e-3, 1.2, 1.2, 0.1].^2);
    %params.Qper1 = diag([1e-3, 1e-3, 1.0, 1.0, 0.15].^2); 
    
    % perturbation covariance for CAR: 6x6 Matrix
    params.Acar = eye(6) + [0,0,dt,0,0,0; 0,0,0,dt,0,0; zeros(4, 6)];
    params.Qcar2 = diag((dt * [1e-3, 1e-3, 1, 3, 0.03/dt, 0.005/dt]).^2);
    
    %Pertubatio
    params.Qcam = diag((dt * [1e-6, 1e-6, 1e-6, 5/dt, .5*pi/180/dt, .7, 1e-6, 1e-6]).^2);
    
    params.camPert = params.Qcam / 16; %diag([1e-10, 0.0001, 1e-10, 1, (pi/180*5/6)^2]) / 4;
    params.perPert0 = diag([0.01, 0.01, 0, 0, 0.01/4]);
    params.carPert0 = diag([0.01, 0.01, 0, 0, 0.16, 0.01/4]);
    
%     params.perPert1 = params.Qper1 / 16;
%     params.perPert2 = params.Qper2 / 16; % diag([1e-10, 1e-10, 0.01, 0.01, 0.0001]) / 16;
    
    params.mh = 1.7;
    params.sh = 0.1;
    params.normh = 1/sqrt(2 * pi * params.sh^2);

    params.mch = 1.3;
    params.sch = 0.1;
    params.normch = 1/sqrt(2 * pi * params.sch^2);
    params.noCprob = normpdf(2.8 * params.sch, 0, params.sch);
    
    params.sv = 1;
    params.normv = 1/sqrt(2 * pi * params.sv^2);
    params.noHprob = normpdf(2.8 * params.sh, 0, params.sh); % * normpdf(3 * params.sv, 0, params.sv);
    
    params.hdet = 0.6; % probability to detect wehn human
    params.nhdet = 0.01; % probability to detect when non human - false alarm.

    params.flipObj = 0.1;
    
    %% KLT parameters
    kltstd = 2;
    params.kltPrec = eye(2) / kltstd^2; % sx = 1/2, sy = 1/2; 
    params.lkltNorm = log(1/sqrt(2 * pi * kltstd^2));
    params.lkltuprob = log(normpdf(2*kltstd, 0, kltstd)^2); % at 2 sigma
    params.flipFeat = 0.1;
    params.nfeatuse = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.KLTmargin = 7;
end
params.dt = dt;

params.useInteraction = 1;
%% INTERACTIONS PARAMETERS
params.c1_pair = 0.9;
params.i1_pair = 3;   % slope of the soft step
params.i2_pair = 1.5; % threshold of the soft step

params.lr2_pair = 6; % 
params.lg1_pair = 5; %
params.lg2_pair = 1/2 * 4^2; % 0.125...
params.betatrns = 0.9 * eye(2) + 0.05 * ones(2,2); % transition of interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Precompute normalization parameters!!!! for faster computation.
params.Prec1 = inv(params.Qdet);     % observation noise
params.Prec2 = inv(params.Qktrack);  % ms track noise.
params.CPrec = inv(params.Qcardet);  % observation noise

params.PrecCam = inv(params.Qcam); 
params.lnormDet = log(1/ sqrt((2 * pi)^3 * det(params.Qdet)));
params.lnormKtrack = log(1/ sqrt((2 * pi)^3 * det(params.Qktrack)));
params.lnormCDet = log(1/ sqrt((2 * pi)^4 * det(params.Qcardet)));
end