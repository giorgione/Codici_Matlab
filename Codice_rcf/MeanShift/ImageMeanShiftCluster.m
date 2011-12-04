%Esempio del Mean Shift CLUSTERING applicato alle Immagini
%
%Genero dei punti da una Mistura di gaussiane 2D
clc;clear;close all;
Irgb=imread('lena_256.jpg');
Irgb=Irgb(1:4:end,1:4:end,:);
figure(1);
imshow(Irgb);
R=Irgb(:,:,1);
G=Irgb(:,:,2);
B=Irgb(:,:,3);

[m,n,c]=size(Irgb);
Npoints=m*n;

R=reshape(R,1,Npoints);
G=reshape(G,1,Npoints);
B=reshape(B,1,Npoints);
Xi=double([R;G;B])/255.0;

%Immagine Clusterizzata
Rc=zeros(1,Npoints);
Gc=zeros(1,Npoints);
Bc=zeros(1,Npoints);

%%Visualizzo lo Spazio Colore
figure(2);hold on;
for i=1:Npoints
        %Normalizzo
        color=Xi(:,i);
        %'MarkerSize',Sizes(i)
        plot3(Xi(1,i),Xi(2,i),Xi(3,i),'o','MarkerFaceColor',color, ...
        'MarkerEdgeColor',color); 
        grid on; box on;
        
end

%Etichette da assegnare ai Punti
Labels=zeros(1,Npoints);

x = sym('x','real');
y = sym('y','real');
Z=[x;y]

h=0.1;

%Eseguo T iterazioni per vedere dove va il Mean shift
T=50;
t=2;
 
Xp=zeros(2,T);
%Salvo la lista dei Punti Visitati
NVPoints=1:Npoints;

Soglia=10^-2;

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
            NModes=NModes+1;
            %Aggiungo la Nuova Moda alla lista delle Mode
            Modes=[Modes X];
            %Label della Nuova Moda
            LabelIndex=NModes;
        end
    end
    
    %Setta tutti i punti visitati col valore della Moda trovata
    Label(VPoints)=LabelIndex;
  
    
end

for i=1:NModes
    [j]=find(Label==i);
    colore=Modes(:,i);
    %'MarkerSize',Sizes(i)
        plot3(Xi(1,j),Xi(2,j),Xi(3,j),'o','MarkerFaceColor',colore, ...
        'MarkerEdgeColor',colore); 
        grid on; box on;
        
        Rc(j)=Modes(1,i);
        Gc(j)=Modes(2,i);
        Bc(j)=Modes(3,i);
end

%Genero l'immagine clusterizzata
[m,n,c]=size(Irgb);
Npoints=m*n;
R=reshape(Rc,m,n);
G=reshape(Gc,m,n);
B=reshape(Bc,m,n);

Ic=zeros(m,n,3);
Ic(:,:,1)=R;
Ic(:,:,2)=G;
Ic(:,:,3)=B; 
figure(4)
imshow(Ic)

%Pixel Totali Elaborati
TotPixel=0;
%Numero totale di regioni (CONNESSE) clusterizzate con MeanShift.
TotRegioni=0;
%Numero Regioni Accettate
AcceptedRegion=0;

for Mi=1:NModes
    [j]=find(Label==Mi);
    Rc=zeros(1,Npoints);
    Gc=zeros(1,Npoints);
    Bc=zeros(1,Npoints);
      
    Rc(j)=Modes(1,Mi);
    Gc(j)=Modes(2,Mi);
    Bc(j)=Modes(3,Mi);

    R=reshape(Rc,m,n);
    G=reshape(Gc,m,n);
    B=reshape(Bc,m,n);

    Ic=zeros(m,n,3);
    Ic(:,:,1)=R;
    Ic(:,:,2)=G;
    Ic(:,:,3)=B;
    figure(5)
    imshow(Ic)
    
    %Informazioni di Spazio
    Pi=logical(zeros(1,Npoints));    
    %Genero i Pattern Spaziali
    Pi(j)=1;
    
    Pi=reshape(Pi,m,n);
    %Cerco le componenti connnesse
    [L ,numR] = bwlabel(Pi,8);
    
    %Numero totale di regioni (CONNESSE) clusterizzate con MeanShift.
    TotRegioni=TotRegioni+numR;
    
    disp(['Moda ' num2str(Mi) ' : ' num2str(numR) ' Regioni'])
    figure(6);hold on
    
    for i=1:numR
        Ri=(L==i);
        [r, c] = find(Ri);
        rc = [c r];
        npx=length(r);
        disp(['Regione ' num2str(i) ': ' num2str(npx) ' Pixel']) 
        TotPixel=TotPixel+npx;
        if npx>5
            %Genero l'ellisse attorno ai dati
            Media=mean(rc).';
            S=cov(rc);
            [V l]=eig(S);
            
            Ellissi=@(theta)V*3.0*sqrt(l)*[cos(theta); -sin(theta)]+repmat(Media,1,length(theta));
            theta=linspace(0,2*pi,50);
            E=Ellissi(theta);
            
            %%Disegno il layer
            figure(7)            
            imshow(Ri)
    
            % Disegno l'ellisse
            figure(6);
            plot(E(1,:),E(2,:),'Color',Modes(:,Mi));            
            plot(Media(1),Media(2),'+','Color',Modes(:,Mi))
            p=patch(E(1,:),E(2,:),ones(1,50));
            set(p,'FaceColor',Modes(:,Mi))
            AcceptedRegion=AcceptedRegion+1;
        end
    end
end

disp(['Regioni Totali: ' num2str(TotRegioni)])
disp(['Regioni Trattenute: ' num2str(AcceptedRegion)])
