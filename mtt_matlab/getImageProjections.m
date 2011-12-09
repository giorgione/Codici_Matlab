function [bb] = getImageProjections(Z, id, param, pred, campan)
if nargin < 4
    pred = 0;
    campan = 0;
elseif nargin < 5
    campan = 0;
end

for sidx = 1:Z.nSamples
    person = Z.per((1 + (id-1) * (param.nperv + 1)):(id * (param.nperv + 1)-1), sidx);
    
    if pred == 1
        person = param.Aper * person;
        cam = getNextCamera(Z.cam(:, sidx), param);
        cam(5) = cam(5) + campan;
    else
        cam = Z.cam(:, sidx);
    end
    
    [dummy, temp] = getProjection(person, cam);
    bb(:, sidx) = temp';
end
end