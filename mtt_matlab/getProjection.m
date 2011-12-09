%function [x] = getProjection(z, cam)
%
% Calculates the projection x of the 3D point z 
% in image Plane according to the camera model.
%
% INPUT:
%                     
%   - cam = [focal, height, xcenter, yhorizon, panning, absolute velocity,];
% 
%   - z = [x, z, vx, vz, heigth] --> 3D Point
%
% OUTPUT:
% 
%   - x = [x, y, imhiehgt] --> Image Point
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
function [x, bb] = getProjection(z, cam)


R = [cos(cam(5)), sin(cam(5)); 
    -sin(cam(5)), cos(cam(5))];
if(length(cam) == 8)
    z(1:2) = z(1:2) - cam(7:8);
end
%Camera Alignment with World cordinate System
z(1:2) = R * z(1:2);
%z(3) is hz

%Projective Transform to get pixel (u,v)
% u = f*x/z + uo
% v = f*y/z + vo   cam(2) = height = y of the classical pinhole model

x(1) = cam(1) / z(2) * z(1) + cam(3);
x(2) = cam(1) / z(2) * cam(2) + cam(4);
x(3) = cam(1) / z(2) * z(5);

if nargout > 1
    bb(1) = x(1) - x(3) / 4;
    bb(2) = x(2) - x(3);
    bb(3) = x(3) / 2;
    bb(4) = x(3);
end

end