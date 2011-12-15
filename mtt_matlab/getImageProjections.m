% function [bb] = getImageProjections(Z, id, param, pred, campan)
% 
function [bb] = getImageProjections(Z, id, param, pred, campan)
if nargin < 4
    %First frame processing
    pred = 0;
    campan = 0;
else if nargin < 5
    %Other frames and camera pannning angle not specified
     campan = 0;
end
%Consider each GENERATED SAMPLE for the PERSON VARIABLE
% sidx: sample index
for sidx = 1:Z.nSamples
    ind1=(1 + (id-1) * (param.nperv + 1));
    ind2=id * (param.nperv + 1)-1; 
    
    %Extract the person 
    person = Z.per((ind1:ind2), sidx);
    
    %Generate the STATE VALUES for CAMERA VARIABLE cam
    
    if pred == 1
        % First frame of the video sequence processing
        %Evolution of Person
        person = param.Aper * person;
        %Evolution of CAMERA
        cam = getNextCamera(Z.cam(:, sidx), param);
        %Update CAMERA PAN
        cam(5) = cam(5) + campan;
    else
        %Other frames: 
        cam = Z.cam(:, sidx);
    end
    
    %Project person on the image plane
    [dummy, temp] = getProjection(person, cam);
    %Person BBOX 
    bb(:, sidx) = temp';
end
end