%function cam = getNextCamera(cam, params)
% CAMERA DYNAMIC MODEL encoding the temporal relationship between camera
% parameters
% 
% INPUT PARAMETERS:
%
% cam: focal length, camera height(meter), x camera center(pixel), horizon
% (pixel), initial angle, initial speed, camera's X location, camera's Z location
%
% params:  Other variables in the model (only params.dt is needed)'
%
% OUTPUT PARAMETERS:
%
% cam: the new state for (X,Z) camera's location
function cam = getNextCamera(cam, params)
    if params.ncamv == 8
        cam(7) = cam(7) - cam(6) * sin(cam(5)) * params.dt; % + ? - ?
        cam(8) = cam(8) + cam(6) * cos(cam(5)) * params.dt;
    end
end