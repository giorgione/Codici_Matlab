function campan = getRoughCamPan(Z, KLT, frame, param)

if frame < -1 % 50
    campan = 0;
    return;
end

gfidx = find(Z.gfcnt > 1);
gfmatch = zeros(1, length(gfidx));

for j = 1:length(gfidx)
    temp = find(KLT.vidx == Z.gfidx(gfidx(j)));
    if length(temp) == 0
        continue;
    end
    if ((KLT.x(temp, frame) > 0) && (KLT.y(temp, frame) > 0))
        gfmatch(j) = temp;
    end
end

mfeat = mean(Z.gfeat, 2);
gflist = reshape(mfeat, 3, length(mfeat) / 3);
gflist = gflist(:, gfidx);

campan = 0;
cnt = 0;
step = pi/180/5;
direction = 1;
cam = getNextCamera(Z.cam(:, 1), param);
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