% function campan = getRoughCamPan(Z, KLT, frame, param)
%  Calculate an Rough estimate of the Cam Pan angle by MATCHING KLT FEATURES
%  extracted from the GROUND PLANE
%  
%  - Z:     VARIABLE STATUS in Model
%  - KLT:   KLT features for the video loaded from FILE
%  - frame: frame index
%  - param:  DATA of the MODEL 
%
%  - campan: Estimate of CAMERA PAN ENGLE
function campan = getRoughCamPan(Z, KLT, frame, param)

if frame < -1 % 50
    campan = 0;
    return;
end

%GROUND FEATURES INDEX and MATCH INDEX:
%at first frame they are empty --> []
gfidx = find(Z.gfcnt > 1); 

%MATCH INDEX MATRIX :
%              0  -->  if the FEATURE is matched with current
% gfmatch(j) = 
%              FEATURE INDEX  --> if the FEATURE is matched with current
%                                 FRAME
gfmatch = zeros(1, length(gfidx));

%For each ground KLT features search MATCHING with Current Frame
for j = 1:length(gfidx)
    %search for matching index betwen:
    %
    % KLT: data containing all Features
    %
    % Z.gfidx: FEATURES index in the Current FRAME
    temp = find(KLT.vidx == Z.gfidx(gfidx(j)));
    %No matching: frame=0
    if length(temp) == 0
        continue;
    end
    if ((KLT.x(temp, frame) > 0) && (KLT.y(temp, frame) > 0))
        %save the Matching index in KLT
        gfmatch(j) = temp;
    end
end

%Compute GROUND FEATURES mean values
mfeat = mean(Z.gfeat, 2);

%3D Point
gflist = reshape(mfeat, 3, length(mfeat) / 3);

gflist = gflist(:, gfidx);

% Get New STATE ESTIMATE of CAMERA's LOCATION from CAMERA DYNAMIC MODEL
cam = getNextCamera(Z.cam(:, 1), param);
% Compute Distance betwen camera Position
best = computeDist(gflist, gfmatch, cam, KLT, frame);

% Analyze 30 possible Variiation for The PAN ANGLE given the current
% camera status and compute the total distance between KLT features point
campan = 0; %pan angle
cnt = 0;    %count variable
step = pi/180/5;
direction = 1;


while(cnt < 30)
    %Temp camera pan configuration
    tpan = campan + direction * step;
    %Temp camera configuration 
    tcam = cam;
    tcam(5) = tcam(5) + campan + tpan;
    %Compute the TOTAL distance on the ground plane between matched points
    tdist = computeDist(gflist, gfmatch, tcam, KLT, frame);
    %Pick the minimum value of the distance
    if tdist < best
        %THE BEST CAMPPAN ANGLE ESTIMATE
        campan = campan + tpan;
        best = tdist;
    else
        %Aalyze Panning Angle in campan+[-pi/2 pi/2]
        direction = direction * -1;
    end
    cnt = cnt + 1;
end

end
%function dist = computeDist(gflist, gfmatch, cam, KLT, frame)
% gflist:  
% gfmatch: Vector containing index of MATCHING features in Z
% cam
% KLT
% frame
function dist = computeDist(gflist, gfmatch, cam, KLT, frame)
dist = 0;
%Matching GROUND FEATURES
for i = 1:size(gflist, 2)
    if gfmatch(i) > 0
        %GET the Projection
        temp = getGProj(gflist(:, i), cam);
        %Extract (x,y) data of Current FEATURE
        xx = [KLT.x(gfmatch(i), frame), KLT.y(gfmatch(i), frame)];
        %
        dist = dist + sum((temp - xx) .^ 2);
    end
end

end