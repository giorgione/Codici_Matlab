%Esempio del Mean Shift CLUSTERING
%
%Genero dei punti da una Mistura di gaussiane 2D
clc;clear;close all;
mu1 = [20 20];
mu2 = [ 30 10];
SIGMA1 = [10  0.5;
          0.5 4];

SIGMA2 = [10 -0.6;
          -0.6 8 ];
      
Npoints=100;
Xi =[ mvnrnd(mu1,SIGMA1,Npoints);
      mvnrnd(mu2,SIGMA2,Npoints)];
Xi=Xi.';
Npoints=2*Npoints;

%Etichette da assegnare ai Punti
Labels=zeros(1,Npoints);

%Disegno  il punto iniziale
plot(Xi(1,:),Xi(2,:),'+');hold on;

%Disegno le Mode della distribuzione
plot(mu1(1),mu1(2),'*r');
plot(mu2(1),mu2(2),'*r');

syms x y;
Z=[x;y]
%f=exp(-x^2/2);
%df=diff(f,'x')

g=@(x)0.5*(exp(-(x^2)));
h=6;

%Eseguo T iterazioni per vedere dove va il Mean shift
T=50;
t=2;
 
Xp=zeros(2,T);
%Salvo la lista dei Punti Visitati

NVPoints=1:Npoints;


Soglia=10^-2;
I=1;

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
       Circ= (Z-X).'*(Z-X)-h^2;
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
    colore=rand(3,1);
    colore=colore/sum(colore);
    plot(Xi(1,j),Xi(2,j),'+','color',colore)
end