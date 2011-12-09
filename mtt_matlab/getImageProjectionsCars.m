function [bb] = getImageProjectionsCars(Z, id, param, pred, campan)
if nargin < 4
    pred = 0;
    campan = 0;
elseif nargin < 5
    campan = 0;
end

for sidx = 1:Z.nSamples
    car = Z.car((1 + (id-1) * (param.ncarv + 1)):(id * (param.ncarv + 1)-1), sidx);
    
    if pred == 1
        car = param.Acar * car;
        cam = getNextCamera(Z.cam(:, sidx), param);
        cam(5) = cam(5) + campan;
    else
        cam = Z.cam(:, sidx);
    end
    
    [dummy, temp] = getProjectionCar(car, cam);
    bb(:, sidx) = temp';
end
end