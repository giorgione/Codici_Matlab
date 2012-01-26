% Cerca i contorni nel immagine


clc;clear;close all;
Irgb_orig=imread('./limit30.png');
Irgb=Irgb_orig;

Igray=rgb2gray(Irgb);
figure(1)
BW = edge(Igray,'canny');
imshow(BW);
figure(2);
hold on;
s=size(BW);

k=0;
X=zeros(s);

for i = 1:s(1)
   
  for j=2:s(2)
    
      if(BW(i,j)==1 & X(i,j)==0)

           contour = bwtraceboundary(BW, [i, j], 'W',8);
           
           if(~isempty(contour) )
              
              % elimina contorni duplicati
              if(k>0)
                  [r c]=find(X==0);
                  idx1=r*s(2)+j; %non visitati
                  idx2=contour(:,1)*s(2)+contour(:,2); %pixel di contorno trovati
                  
                  S=setdiff(idx2,idx1);
                  %se ho piu di 40 pixel non visitati allora poso inserire
                  %un nuovo contorno
                  if(length(S)<40)
                       continue;
                  end
              end
              %Marca i pixel gia processati
              np=size(contour(:,2));
              for l=1:np
                X(contour(l,1),contour(l,2))=1;
              end
              %aggiorna il numero di contorni
              k=k+1;
              disp(['contur ' num2str(k) ': ' num2str(size(contour,1)) ]);
              contorno(k).x=contour(:,2);
              contorno(k).y=contour(:,1);
              plot(contour(:,2),contour(:,1),'g','LineWidth',2);
              plot(contour(1,2),contour(1,1),'or','MarkerFaceColor','r','MarkerEdgeColor','r');
              pause
              plot(contour(end,2),contour(end,1),'ob','MarkerFaceColor','b','MarkerEdgeColor','b');
              pause  
              
           end
      end
      
   end
end