close all;clc

%Per ogni immagine faccio una segmentazione
% sul colore
%
% Genero sui colori il contorno
%
% Sommo i contorni

 
for i = 1:length(files)
	filename = ['./' files(i).name];
   
    IrgbTmp=imread(filename);
    
    [mF,nF,cF]=size(IrgbTmp);
    disp(['Image' filename ': ' num2str(m) ' x ' num2str(n)])
    Npoints=mF*nF;
    
     Irgb=im2double(IrgbTmp);
    

     R=Irgb(:,:,1);
     G=Irgb(:,:,2);
     B=Irgb(:,:,3);
     R=reshape(R,1,Npoints);
     G=reshape(G,1,Npoints);
     B=reshape(B,1,Npoints);

     FrameI=[R;G;B];

     
 
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
 

end