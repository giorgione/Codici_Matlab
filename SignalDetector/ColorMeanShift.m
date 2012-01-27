function Cluster=ColorMeanShift(Xi)
[m,Npoints]=size(Xi);


%Etichette da assegnare ai Punti
Labels=zeros(1,Npoints);


%Eseguo T iterazioni per vedere dove va il Mean shift
T=15;
t=2;
 
Xp=zeros(3,T);
%Salvo la lista dei Punti Visitati
NVPoints=[];

h=0.55;
Soglia=0.25;

%Numero Mode
NModes=0;
Modes=[];
Pindex=1:Npoints;
for Xindex=1:Npoints
        
    %% Mean Shift sul singolo punto X
    
    %Ottieni il punto iniziale dall'insieme dei punti non ancora visitati.
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
       
       % Ottimizzato per MATLAB: Molto + veloce
       %Cerco i punti che cadono nell'IperSfera Centrata in X
       Sx=[];
       Sindex=[];
       Nx=0; 
       Xrep=repmat(X,1,Npoints);
       Xtmp=Xi-Xrep;
       Xtmp=Xtmp.^2;
       Xtmp=sum(Xtmp)-h^2;
       j=Xtmp<0;
       Sindex=Pindex(j);      %Indice punti nell'intorno
       Sx=Xi(:,Sindex);       %Colore dei punti
       Nx=length(Sindex);     %Numero di Punti nell'intorno
       
% VERSIONE NON OTTIMIZZATA       
%        for i=1:Npoints
%             if (X-Xi(:,i)).'*(X-Xi(:,i))< h^2
%                 %Insieme dei Pixel il cui COLORE cade nell'intorno del
%                 %COLORE corrente X
%                 Sx=[Sx Xi(:,i)];
%                 Sindex=[Sindex i];
%                 Nx=Nx+1;
%                 %indici dei punti visitati
%                 %VPoints=union(VPoints,i);
%                 
%                 
%             end        
%        end
       
       %Per ogni punto che cade vado a vedere se:
       %
       % 1) e' stato gia visitato
       %
       % 2) qual'e' la moda + ricorrente
       %
       % Assegno la Moda + ricorrente nell'intorno al pixel
       H=hist(Labels(Sindex),0:NModes);
       H=H/sum(H);
       
       [MaxMode j]=max(H);
       if j>1 | sum(H(2:end)) > 0.6
           [MaxMode j]=max(H(2:end));
           Labels(Xindex)=j;
           break;
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
    %NVPoints=setdiff(NVPoints,VPoints);
    NVPoints=[NVPoints Xindex];
    
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
            ModaCur=Modes(:,j);
        else
            %Nuova Moda
            NModes=NModes+1
            %Aggiungo la Nuova Moda (X) alla lista delle Mode
            ModaCur=X;
            Modes=[Modes ModaCur];
            %Label della Nuova Moda
            LabelIndex=NModes;
        end
        
         
        
        
    end
     
    %Setta tutti i punti visitati col valore della Moda trovata
    Labels(Xindex)=LabelIndex;
  
    disp([ num2str(Xindex) ' : ' num2str(LabelIndex) ])
end

for i=1:NModes
    [j]=find(Labels==i);
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
    Cluster(i).Dmax=mvnpdf(Cluster(i).Mode,Cluster(i).Mean,Cluster(i).Cov)
    norm(Cluster(i).Mean-Cluster(i).Mode);
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
    
    Rc(j)=Cluster(i).Mean(1);
    Gc(j)=Cluster(i).Mean(2);
    Bc(j)=Cluster(i).Mean(3);
end


