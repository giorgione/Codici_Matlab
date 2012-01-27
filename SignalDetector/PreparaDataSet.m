%Script that reads a content files and displays annotation info
%/larsson@isy.liu.se 2011
clc;clear;close all
%Set paths to the images and to the annotation file
base = './';
imagePath = fullfile(base,'Set2Part0');
annotationFn = fullfile(imagePath,'annotations.txt');

 %Get content struct from annotation file
content = parseSignAnnotations(annotationFn);      
content1=EstraiSegnali(content,'70_SIGN');
%For all images with annotation, display content
N = length(content1);


%X: Punti di Colore estratti
X=[];
NTrain=0;
for I = 1:N
	
    %title(content1(i).name);
	
    %Estraggo la patch
    %Numero segnali al tempo t
    nrSigns=length(content1(I).signs);
   
    for i=1:nrSigns
        if(i==1)
           fn = fullfile(imagePath,content1(I).name);
           img = im2double(imread(fn));    
           %imagesc(img);axis image;
           %title(content1(i).name);
        end
       %plot bounding box: ordine cui sono memorizzati i segnali
       % 1   4
       % 
       % 2   3
       % 
		bb = round(content1(I).signs(i).signBB);
		x = [bb(1) bb(1) bb(3) bb(3) bb(1)] ;
		y = [bb(2) bb(4) bb(4) bb(2) bb(2)] ;   
%         hold on;
% 		for j=1:4
%         
%             plot(x(j),y(j),'or')%,'lineWidth',lw);
%             pause
%         end
%         hold off;
        
        Irgb=img(bb(2):bb(4),bb(1):bb(3),:);
        [m,n,c]=size(Irgb);
        Npoints=m*n;

        R=Irgb(:,:,1);
        G=Irgb(:,:,2);
        B=Irgb(:,:,3);
        R=reshape(R,1,Npoints);
        G=reshape(G,1,Npoints);
        B=reshape(B,1,Npoints);
        
        X=cat(2,X,[R;G;B]);
         
        
        NTrain=NTrain+1;
        
        
%         hold off;
%         drawnow
%         pause
    end
    if(NTrain >15)
            break;
        end
    
	 
end

[m,Npoints]=size(X);
   
  figure; hold on
  for i=1:Npoints
    plot3(X(1,i),X(2,i),X(3,i),'o','MarkerFaceColor',X(:,i), ...
        'MarkerEdgeColor',X(:,i)); 
        grid on; box on; 
  end
  
  
 Cluster=ColorMeanShift(X)
 lenght()
