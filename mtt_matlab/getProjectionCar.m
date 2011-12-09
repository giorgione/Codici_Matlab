function [x, bb] = getProjectionCar(z, cam)

% cam = [focal, height, xcenter, yhorizon, panning];
% x = [x, y, imhiehgt]
% z = [x, z, vx, vz, heihgt]

R = [cos(cam(5)), sin(cam(5)); -sin(cam(5)), cos(cam(5))];
if(length(cam) == 8)
    z(1:2) = z(1:2) - cam(7:8);
end
z(1:2) = R * z(1:2);

x(1) = cam(1) / z(2) * z(1) + cam(3);
x(2) = cam(1) / z(2) * cam(2) + cam(4);
x(3) = cam(1) / z(2) * z(5);
x(4) = cam(1) / z(2) * z(6);

if nargout > 1
    bb(1) = x(1) - x(3) / 2;
    bb(2) = x(2) - x(4);
    bb(3) = x(3);
    bb(4) = x(4);
end

end