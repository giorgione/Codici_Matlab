% function campan = getRoughCamPan(Z, KLT, frame, param)
%  Calculate an Rough estimate of the Cam Pan angle through KLT FEATURES
%  extracted from the GROUND PLANE
%  
%  - Z:     Camera Model
%  - KLT:   KLT features for the video loaded from FILE
%  - frame: frame index
%  - param:
function campan = getRoughCamPan(Z, KLT, frame, param)

if frame < -1 % 50
    campan = 0;
    return;
end

%ground features index and match value
gfidx = find(Z.gfcnt > 1);
gfmatch = zeros(1, length(gfidx));

%For each ground features
for j = 1:length(gfidx)
    %search for matching index
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
gflist = reshape(mfeat, 3, length(mfeat) / 3);
gflist = gflist(:, gfidx);

campan = 0;
cnt = 0;
step = pi/180/5;
direction = 1;

% Get New STATE ESTIMATE of CAMERA's LOCATION from CAMERA DYNAMIC MODEL
cam = getNextCamera(Z.cam(:, 1), param);
% Compute 
best = computeDist(gflist, gfmatch, cam, KLT, frame);

while(cnt < 30)
    tpan = campan + direction * step;
    tcam = cam;
    tcam(5) = tcam(5) + campan + tpan;
    tdist = computeDist(gflist, gfmatch, tcam, KLT, frame);
    if tdist < best
        campan = campan + tpan;
        best = tdist;
    else
        direction = direction * -1;
    end
    cnt = cnt + 1;
end

end

function dist = computeDist(gflist, gfmatch, cam, KLT, frame)
dist = 0;

for i = 1:size(gflist, 2)
    if gfmatch(i) > 0
        temp = getGProj(gflist(:, i), cam);
        xx = [KLT.x(gfmatch(i), frame), KLT.y(gfmatch(i), frame)];
        dist = dist + sum((temp - xx) .^ 2);
    end
end

end