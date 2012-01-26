%Esempio del Mean Shift CLUSTERING applicato alle Immagini
%
% Meanshift to cluster image in Color
%
%
%
clc;clear;close all;

files = dir([ './*.png']);
%Tutte le immagini Concatenate
Irgb=[];
%Pixel Concatenati
Xi=[];
%par
for i = 2:length(files)
	filename = ['./' files(i).name];
   
    Irgb_orig=imread(filename);
   % imshow(Irgb_orig);
   % pause;
    IrgbTmp=imresize(Irgb_orig,[64 64]);
    [m,n,c]=size(IrgbTmp);
    disp(['Image' filename ': ' num2str(m) ' x ' num2str(n)])
    Npoints=m*n;
    
    %Normalize
    IrgbTmp=im2double(IrgbTmp);
    Irgb=cat(1,Irgb,IrgbTmp);

end
figure(1);
imshow(Irgb);
drawnow
pause
[m,n,c]=size(Irgb);
Npoints=m*n;

R=Irgb(:,:,1);
G=Irgb(:,:,2);
B=Irgb(:,:,3);
R=reshape(R,1,Npoints);
G=reshape(G,1,Npoints);
B=reshape(B,1,Npoints);

Xi=double([R;G;B]);%/norm([255.0 ;255 ;255]);
    


%Frame dove fare detection: stored as 3 x Numpixels
FrameI=Xi;

%Immagine Clusterizzata
Rc=zeros(1,Npoints);
Gc=zeros(1,Npoints);
Bc=zeros(1,Npoints);

%% Visualizzo lo Spazio Colore
% figure(2);hold on;
% for i=1:Npoints
%         %Normalizzo
%         color=Xi(:,i);
%         %'MarkerSize',Sizes(i)
%         plot3(Xi(1,i),Xi(2,i),Xi(3,i),'o','MarkerFaceColor',color, ...
%         'MarkerEdgeColor',color); 
%         grid on; box on;
%         
% end

%Etichette da assegnare ai Punti
Labels=zeros(1,Npoints);


%Eseguo T iterazioni per vedere dove va il Mean shift
T=50;
t=2;
 
Xp=zeros(2,T);
%Salvo la lista dei Punti Visitati
NVPoints=1:Npoints;

h=0.55;
Soglia=0.25;

%Numero Mode
NModes=0;
Modes=[];
while isempty(NVPoints)==0
    
    %% Mean Shift sul singolo punto X
    
    %Ottieni il punto iniziale dall'insieme dei punti non ancora visitati.
    Xindex=NVPoints(1);
    X=Xi(:,Xindex);
    %plot(X(1),X(2),'og');hold on

    % Valore Iniziale per il criterio  di Arresto
    E=10;
    
    %Punti Visitati
    VPoints=[];
    
    %t: Numero di Iterazioni
    t=2;
    Xp=X;
    
    while t< T &&  E > Soglia
       %Calcolo il Mean Shift Vector
       %Circ= (Z-X).'*(Z-X)-h^2;
       %ezplot(Circ,[-10^2 10^2])
       
       %Cerco i punti che cadono nell'IperSfera Centrata in X
       Sx=[];
       Nx=0; 

       for i=1:Npoints
            if (X-Xi(:,i)).'*(X-Xi(:,i))< h^2
                %Insieme dei Pixel che cadono nell'intorno
                Sx=[Sx Xi(:,i)];
                Nx=Nx+1;
                %indici dei punti visitati
                VPoints=union(VPoints,i);
                
                
            end        
       end

       X=repmat(X,1,Nx);
       %Il mean shift è un Media Locale
       MeanShift(:,t)=sum(Sx-X,2)/Nx;
       %Fhu=Nx/(N*h^2*pi);

       %Nuovo Punto
       Xp(:,t)=MeanShift(:,t)+Xp(:,t-1);
       E=norm(Xp(:,t)-Xp(:,t-1),2);
       
       %Punto Iniziale Shiftato: al termine del ciclo è la MODA
       X=Xp(:,t);
       
       %Traiettoria generata durante la Procedura di Shifting verso la Moda
       %plot(Xp(1,1:t),Xp(2,1:t),'og-');
       t=t+1;
    end
    
    %Aggiorna l' insieme dei Punti NON VISITATI
    NVPoints=setdiff(NVPoints,VPoints);
    %plot(Xi(1,VPoints),Xi(2,VPoints),'+m');
     
    %Definizione delle Mode
    if isempty(Modes)
        Modes=X;
        NModes=NModes+1;
        LabelIndex=1;
    else
        %Cerco la Moda più vicina al volore cui ho avuto Convergenza
        for i=1:NModes
            Res(i)=norm(X-Modes(:,i),2);
        end
        %calcolo la moda + vicna e l'assegno
        [val,j]=min(Res);
        
        if(val<=h)
            %Moda già trovata
            LabelIndex=j;
        else
            %Nuova Moda
            NModes=NModes+1
            %Aggiungo la Nuova Moda alla lista delle Mode
            Modes=[Modes X];
            %Label della Nuova Moda
            LabelIndex=NModes;
        end
    end
     
    %Setta tutti i punti visitati col valore della Moda trovata
    Label(VPoints)=LabelIndex;
  
    
end

figure(6);hold on; grid on
for i=1:NModes
    [j]=find(Label==i);
    colore=Modes(:,i);
    X=zeros(length(j),3);
    plot3(Xi(1,j),Xi(2,j),Xi(3,j),'o','MarkerFaceColor',colore, ...
        'MarkerEdgeColor',colore); 
        grid on; box on; 
    
    %Costruisco una distribuzione di colore sulla Moda
    X(:,1)=R(j).';
    X(:,2)=G(j).';
    X(:,3)=B(j).';
    
    
    Cluster(i).Mean=mean(X);
    Cluster(i).Cov=cov(X);
    Cluster(i).Data=X;
    Cluster(i).Mode=colore.';
    Cluster(i).Dmax=norm(Cluster(i).Mean-Cluster(i).Mode);
    Cluster(i).Npixel=size(X,1);
    
    disp(['Cluster ' num2str(i) ' : '  num2str(Cluster(i).Npixel)])
    disp('Mean:')
    disp(Cluster(i).Mean  )
    disp('Moda:')
    disp(Cluster(i).Mode  )
    disp('Moda:')
    disp(Cluster(i).Dmax  )
    disp('Cov:')
    disp(Cluster(i).Cov )
    
    Rc(j)=Modes(1,i);
    Gc(j)=Modes(2,i);
    Bc(j)=Modes(3,i);
end

%Genero l'immagine clusterizzata
[m,n,c]=size(Irgb);
Npoints=m*n;
R1=reshape(Rc,m,n);
G1=reshape(Gc,m,n);
B1=reshape(Bc,m,n);

Ic=zeros(m,n,3);
Ic(:,:,1)=R1;
Ic(:,:,2)=G1;
Ic(:,:,3)=B1; 
figure(2)
imshow(Ic)
drawnow
figure(3)

%Filtragggio su Colore
% 1) Seleziona i colori che discostano da una certa soglia dal valore medio
%    di ogni cluster rilevato in precedenza
for i=1:NModes
    ImTmp=mahal(FrameI.',Cluster(i).Data);

    %
    figure(3)
    subplot(1,4,1)
    Max=max(ImTmp);
    ImTmp=ImTmp./Max;
    imshow(Ic)
    %imshow(reshape(1-ImTmp,m,n))
   
    %
    subplot(1,4,2)
    ImColor=zeros(m,n,3);
    ImColor(:,:,1)=Cluster(i).Mean(1);
    ImColor(:,:,2)=Cluster(i).Mean(2);
    ImColor(:,:,3)=Cluster(i).Mean(3);
    imshow(ImColor) 
    
    %figure(5)
    subplot(1,4,3)
    ImColor=zeros(m,n,3);
    ImColor(:,:,1)=Cluster(i).Mode(1);
    ImColor(:,:,2)=Cluster(i).Mode(2);
    ImColor(:,:,3)=Cluster(i).Mode(3);
    imshow(ImColor) 
    
    figure(4)
    hist(ImTmp,linspace(0,1,100));hold on;
    plot(Cluster(i).Dmax *ones(1,3),linspace(0,9000,3),'r')
    plot(linspace(0,1,3),ones(1,3)*Cluster(i).Npixel  ,'m')
    
    Y= hist(ImTmp,linspace(0,1,100));
    Ycum=cumsum(Y);
    idx=find(Ycum > Cluster(i).Npixel);
    plot(ImTmp(idx(1)) *ones(1,3),linspace(0,9000,3),'g')
    
    figure(3)
    subplot(1,4,4)
    idx=(ImTmp < ImTmp(idx(1) ));
    ImF=zeros(m*n,1);
    ImF(idx)=1;
    imshow(reshape(ImF,m,n))
    pause;
    clf(4);
end



