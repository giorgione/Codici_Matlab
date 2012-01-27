close all;clc
BWF=zeros(m,n);
NumFig=6 %Numero finestre attuali
for i=1:NModes
    %ImTmp=mahal(FrameI.',Cluster(i).Data);
    ImTmp=mvnpdf(FrameI.',Cluster(i).Mean,Cluster(i).Cov)
    %
    figure(3)
    subplot(1,NumFig,2)
    Max=max(ImTmp);
    %ImTmp=ImTmp./Max;
    imshow(Ic)
    subplot(1,NumFig,1)
    imshow(reshape(ImTmp,m,n))
   
    %
    subplot(1,NumFig,3)
    ImColor=zeros(m,n,3);
    ImColor(:,:,1)=Cluster(i).Mean(1);
    ImColor(:,:,2)=Cluster(i).Mean(2);
    ImColor(:,:,3)=Cluster(i).Mean(3);
    imshow(ImColor) 
    
    %figure(5)
    subplot(1,NumFig,4)
    ImColor=zeros(m,n,3);
    ImColor(:,:,1)=Cluster(i).Mode(1);
    ImColor(:,:,2)=Cluster(i).Mode(2);
    ImColor(:,:,3)=Cluster(i).Mode(3);
    imshow(ImColor) 
    
    figure(4)
    hist(ImTmp,linspace(0,1,100));hold on;
    plot(mvnpdf( Cluster(i).Mean, Cluster(i).Mean,Cluster(i).Cov)  *ones(1,3),linspace(0,9000,3),'r')
    plot(linspace(0,1,3),ones(1,3)*Cluster(i).Npixel  ,'m')
    
    Y= hist(ImTmp,linspace(0,1,100));
    Ycum=cumsum(Y);
    idx=find(Ycum > Cluster(i).Npixel);
    plot(.8*ones(1,3),linspace(0,9000,3),'g')
    
    figure(3)
    subplot(1,NumFig,5)
    idx=(ImTmp > 0.8);
    ImF=zeros(m*n,1);
    ImF(idx)=1;
    ImF=reshape(ImF,m,n)
    imshow(ImF)
    
    
    subplot(1,NumFig,6)
    BW = edge(ImF,'canny');
    imshow(BW);
    BWF=BWF+BW;
    pause;
    clf(4);
end

