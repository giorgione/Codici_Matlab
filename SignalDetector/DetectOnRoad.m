close all;clc
FrameI=imread('strada2.jpg');
FrameI=imresize(FrameI,0.5);
[mF nF cF]=size(FrameI);
FrameI=im2double(FrameI);
Npoints=mF *nF;
figure(1)
imshow(FrameI);
R=FrameI(:,:,1);
G=FrameI(:,:,2);
B=FrameI(:,:,3);
R=reshape(R,1,Npoints);
G=reshape(G,1,Npoints);
B=reshape(B,1,Npoints);


Ihsv=rgb2hsv(FrameI);
figure(6)
H=Ihsv(:,:,1);
imshow(H);


FrameI=[R;G;B;];
NumFig=6 %Numero finestre attuali
for i=1:NModes
    %ImTmp=mahal(FrameI.',Cluster(i).Data);
    ImTmp=mvnpdf(FrameI.',Cluster(i).Mean,Cluster(i).Cov);
    
    figure(5)
    hist(ImTmp,100)
    ImTmp=ImTmp./max(ImTmp);
   
    %
    figure(2)
    imshow(reshape(ImTmp,mF,nF))
   
    figure(3)
    idx=(ImTmp > 0.13);
    ImF=zeros(mF*nF,1);
    ImF(idx)=1;
    ImF=reshape(ImF,mF,nF);
    imshow(ImF)
    
    figure(4)
    ImColor=zeros(10,10,3);
    ImColor(:,:,1)=Cluster(i).Mean(1);
    ImColor(:,:,2)=Cluster(i).Mean(2);
    ImColor(:,:,3)=Cluster(i).Mean(3);
    imshow(ImColor )

    pause;
    close(2);close(3);close(4)
    
end