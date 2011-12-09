function [Tracks] = inferOcclusion(Tracks, Z, frame, sparams)

bbs = [];

for i = 1:length(Z.peridx)
    bbs(:, i) = getImageProjections(Z, i, sparams, 0, 0);
end

occlusion = zeros(size(bbs, 2), size(bbs, 2));
for i = 1:size(bbs, 2)
    for j = 1:size(bbs, 2)
        validTrack = 0;
        for k = 1:length(Tracks)
            if Tracks(k).tid == Z.peridx(i)
                validTrack = 1;
            end
        end
        if validTrack == 1
            occlusion(i, j) = getOcclusion(bbs(:, i), bbs(:, j));
        end
    end
end
% occlusion
occlusion = sum(occlusion, 1) > 0;
for i = 1:length(Tracks)
    idx = find(Z.peridx == Tracks(i).tid);
    if length(idx) == 0
        continue;
    end
    
    if occlusion(idx) == 1
        Tracks(i).occlusion(end + 1) = frame;
    end
end

end


function occ = getOcclusion(bb1, bb2)
occ = 0;

bb1(1) = bb1(1) + bb1(3) / 8; bb1(3) = bb1(3) * 3 / 4;
bb2(1) = bb2(1) + bb2(3) / 8; bb2(3) = bb2(3) * 3 / 4;

A1 = bb1' + [0, 0, bb1(1:2)' - 1];
B1 = bb2' + [0, 0, bb2(1:2)' - 1];

obox = [max(A1(1), B1(1)), max(A1(2), B1(2)), min(A1(3), B1(3)), min(A1(4), B1(4))];

oa = max(0, obox(3)-obox(1)+1) * max(0, obox(4)-obox(2)+1);
ol = oa / (bb2(3) * bb2(4));

if (ol > 0.55) && (A1(4) > B1(4))
    occ = 1;
end

end