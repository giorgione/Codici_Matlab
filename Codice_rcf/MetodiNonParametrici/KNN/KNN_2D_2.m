% KNN_2D_2.m
%
% Simulazione Algoritmo Knn con Kn=1 (Metodod del Punto + Vicino)per la 
% Stima di una Distribuzione Gaussiana 2-D N(0,1)
%
% Idea:
% 
%       1) Estraggo n=64 campioni da una Distribuzione Normale N(0,1)
%
%       2) Eseguo l' algoritmo K-NN con kn=1 per stimare N
%
%       3) Visualizzo i risultati.
%
% 
clc;clear;close all

    n=64;
    kn=1;

    %Estraggo i campioni da una distribuzione Normale N(0,1)
    TrainingSet=randn(2,n)
    
    %Genero la Griglia in vado a Stimare la Normale
    [x y]=meshgrid(linspace(-3,3,10));
    [r c]=size(x);
    
    figure;
    plot(TrainingSet(1,:),TrainingSet(2,:),'ob','MarkerSize',2,'MarkerFaceColor','b');
    hold on;    
    plot(x(:),y(:),'or','MarkerSize',2,'MarkerFaceColor','r')
    legend('Training Set','Punti x in cui vado a stimare P')
    title(['Numero di Campioni: ' num2str(n)])
    
    %Per ogni (x,y) in TestSet  calcolo la distanza dai vettori del training
    %Set (n)
    Distanza=zeros(r,c,n);
    for i=1:r
        for j=1:c
            %calcolo la distanza del punto (x(i);y(j)) dagli n vettori del training Set
            for k=1:n
                Distanza(i,j,k)=norm([x(i,j);y(i,j)]-TrainingSet(:,k),2);
            end
       end
    end
    
    %Ordina le righe di Distanza in modo crescente: 
    %   - I contiene gli indici
    %   - D contiene la Matrice delle Distanze Oridnate
    [D,I]=sort(Distanza,3);

    %Per ogni (x,y) in TestSet calcolo il Volume contente k-vicini 
    V=zeros(r,c);
    for i=1:r
        for j=1:c
             % D(i,j,kn) specifica il Raggio della Circonferenza di cui calcolo l' Area          
             V(i,j)=(pi)*D(i,j,kn)^2; % calcola il volume in funzione dei k vicini
             
             %Disegno i cerchi
             DisegnaCirconferenza([x(i,j);y(i,j)],D(i,j,kn),20,1);

        end  
            end

    axis equal;

    figure;
    P=zeros(r,c);
    P=kn./(n.*V);
    surf(x,y,P)
    hold on;
    
    pause;
    close all;