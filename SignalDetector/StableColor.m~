%Esempio del Mean Shift CLUSTERING applicato alle Immagini
%
%Genero dei punti da una Mistura di gaussiane 2D
clc;clear;close all;
Irgb_orig=imread('./limit30.png');
Irgb=Irgb_orig;
% Irgb=Irgb_orig(1:4:end,1:4:end,:);
% figure(1);
% imshow(Irgb);
% R=Irgb(:,:,1);
% G=Irgb(:,:,2);
% B=Irgb(:,:,3);
% 
% [m,n,c]=size(Irgb);
% Npoints=m*n;
% 
% R=reshape(R,1,Npoints);
% G=reshape(G,1,Npoints);
% B=reshape(B,1,Npoints);
% Xi=double([R;G;B])/255.0;
% 
% %Immagine Clusterizzata
% Rc=zeros(1,Npoints);
% Gc=zeros(1,Npoints);
% Bc=zeros(1,Npoints);
% 
% %% Visualizzo lo Spazio Colore
% % figure(2);hold on;
% % for i=1:Npoints
% %         %Normalizzo
% %         color=Xi(:,i);
% %         %'MarkerSize',Sizes(i)
% %         plot3(Xi(1,i),Xi(2,i),Xi(3,i),'o','MarkerFaceColor',color, ...
% %         'MarkerEdgeColor',color); 
% %         grid on; box on;
% %         
% % end
% 
% %Etichette da assegnare ai Punti
% Labels=zeros(1,Npoints);
% 
% h=0.3;
% %Eseguo T iterazioni per vedere dove va il Mean shift
% T=50;
% t=2;
%  
% Xp=zeros(2,T);
% %Salvo la lista dei Punti Visitati
% NVPoints=1:Npoints;
% 
% Soglia=0.25;
% 
% %Numero Mode
% NModes=0;
% Modes=[];
% while isempty(NVPoints)==0
%     
%     %% Mean Shift sul singolo punto X
%     
%     %Ottieni il punto iniziale dall'insieme dei punti non ancora visitati.
%     Xindex=NVPoints(1);
%     X=Xi(:,Xindex);
%     %plot(X(1),X(2),'og');hold on
% 
%     % Valore Iniziale per il criterio  di Arresto
%     E=10;
%     
%     %Punti Visitati
%     VPoints=[];
%     
%     %t: Numero di Iterazioni
%     t=2;
%     Xp=X;
%     
%     while t< T &&  E > Soglia
%        %Calcolo il Mean Shift Vector
%        %Circ= (Z-X).'*(Z-X)-h^2;
%        %ezplot(Circ,[-10^2 10^2])
%        
%        %Cerco i punti che cadono nell'IperSfera Centrata in X
%        Sx=[];
%        Nx=0; 
% 
%        for i=1:Npoints
%             if (X-Xi(:,i)).'*(X-Xi(:,i))< h^2
%                 %Insieme dei Pixel che cadono nell'intorno
%                 Sx=[Sx Xi(:,i)];
%                 Nx=Nx+1;
%                 %indici dei punti visitati
%                 VPoints=union(VPoints,i);
%                 
%                 
%             end        
%        end
% 
%        X=repmat(X,1,Nx);
%        %Il mean shift è un Media Locale
%        MeanShift(:,t)=sum(Sx-X,2)/Nx;
%        %Fhu=Nx/(N*h^2*pi);
% 
%        %Nuovo Punto
%        Xp(:,t)=MeanShift(:,t)+Xp(:,t-1);
%        E=norm(Xp(:,t)-Xp(:,t-1),2);
%        
%        %Punto Iniziale Shiftato: al termine del ciclo è la MODA
%        X=Xp(:,t);
%        
%        %Traiettoria generata durante la Procedura di Shifting verso la Moda
%        %plot(Xp(1,1:t),Xp(2,1:t),'og-');
%        t=t+1;
%     end
%     
%     %Aggiorna l' insieme dei Punti NON VISITATI
%     NVPoints=setdiff(NVPoints,VPoints);
%     %plot(Xi(1,VPoints),Xi(2,VPoints),'+m');
%      
%     %Definizione delle Mode
%     if isempty(Modes)
%         Modes=X;
%         NModes=NModes+1;
%         LabelIndex=1;
%     else
%         %Cerco la Moda più vicina al volore cui ho avuto Convergenza
%         for i=1:NModes
%             Res(i)=norm(X-Modes(:,i),2);
%         end
%         %calcolo la moda + vicna e l'assegno
%         [val,j]=min(Res);
%         
%         if(val<=h)
%             %Moda già trovata
%             LabelIndex=j;
%         else
%             %Nuova Moda
%             NModes=NModes+1;
%             %Aggiungo la Nuova Moda alla lista delle Mode
%             Modes=[Modes X];
%             %Label della Nuova Moda
%             LabelIndex=NModes;
%         end
%     end
%     
%     %Setta tutti i punti visitati col valore della Moda trovata
%     Label(VPoints)=LabelIndex;
%   
%     
% end
% 
% 
% %Genero l'immagine clusterizzata
% [m,n,c]=size(Irgb);
% Npoints=m*n;
% R1=reshape(Rc,m,n);
% G1=reshape(Gc,m,n);
% B1=reshape(Bc,m,n);
% 
% Ic=zeros(m,n,3);
% Ic(:,:,1)=R1;
% Ic(:,:,2)=G1;
% Ic(:,:,3)=B1; 
% figure(2)
% %imshow(Ic)
% 
% Gray=round(linspace(1,255,NModes));
% Label1=Label
% for i=1:NModes
%     idx=(Label==i);
%     Label1(idx)=Gray(i);
% end
% Igray1=reshape(Label1,m,n)
% 
% BW1 = edge(Igray1,'canny')
% figure(2)
% imshow(BW1);

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



