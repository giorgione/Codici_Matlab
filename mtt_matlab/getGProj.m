%function [x] = getGProj(z, cam)
%
% Calculates the GROUND projection x of the 3D point z on the GROUND
% in image Plane according to the camera model.
%
%  EXTRINSIC CAMERA MODEL:
%
%  Zw: Z in World coordinate system
%  Zc: Z in Camera coordinate system
%
%          Rcw                            Rwc
%           |                              |   -1
%  Zw = R(theta)*Zc + T     <-->  Zc = R(theta)    (Zw-T)
%
%  Rwc: WORLD COORDINATE to CAMERA COORDINATE Rotation Matrix  
%            -1                                    -1
%  Rwc =  Rcw    =   [ cos(cam(5)) ,  -sin(cam(5) ]   
%                    [ sin(cam(5)) ,   cos(cam(5) ] 
%
%  PROJECTION MODEL:
% 
% u = f*x/z + uo
% v = f*y/z + vo   cam(2) = height = y of the classical pinhole model
function [x] = getGProj(z, cam)
%         f                 uo       vo     tetha
% cam = [focal, height, xcenter, yhorizon, panning, Xth, Zth];
% x = [x, y, imhiehgt]              --> Image Data
% z = [x, z, vx, vz, heigth]        --> 3D Point in WORLD COORDINATE SYSTEM

%Camera Alignment with World cordinate System
%

Rwc = [cos(cam(5)), sin(cam(5));
    -sin(cam(5)), cos(cam(5))];
if(length(cam) == 8)
    % Inverse Translation
    z(1:2) = z(1:2) - cam(7:8);
end

% CAMERA COORDINATE SYSTEM : WORLD COORDINATE --> CAMERA COORDINATE
z(1:2) = Rwc * z(1:2);
%z(3) is hz

% PROJECTIVE TRANSFORM to get pixel (u,v)
% u = f*x/z + uo
% v = f*y/z + vo   cam(2) = height = y of the classical pinhole model
x(1) = cam(1) / z(2) * z(1) + cam(3);
x(2) = cam(1) / z(2) * cam(2) + cam(4);

end