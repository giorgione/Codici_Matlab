%Tracking MCMC
%imgdir='/home/gemiggio/Documents/MATLAB/ToyotaData_20nly/left_images/';
%detdir='/home/gemiggio/Documents/MATLAB/ToyotaData_20nly/detections/';

imgdir='/home/giorgio/Documents/Codici_Matlab/ToyotaData_20nly/left_images/';
detdir='/home/giorgio/Documents/Codici_Matlab/ToyotaData_20nly/detections/';

option=[];
ihorizon=230;
prefix=[];
thlist=-0.5;
TrackTop(imgdir, detdir, option, ihorizon, prefix, thlist)

%DRAwing
hold on;
plot(KLT.x(:,1),KLT.y(:,1),'r*') %all KLT
plot(KLT.x(Z.gfidx,1),KLT.y(Z.gfidx,1),'r*') % RETAINED KLT for frame 1
