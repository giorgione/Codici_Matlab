% KNN_1D_2.m
%
% Simulazione Algoritmo Knn per la Stima di una Distribuzione Gaussiana 1-D.
% N(0,1)
%
% Idea:
%
% for I=2:2:12
% 
%       1) Estraggo n=2^I campioni da una Distribuzione Normale N(0,1)
%
%       2) Eseguo l' algoritmo K-NN con kn=sqrt(n) per stimare N
%
%       3) Visualizzo i risultati.
%
% end
%
% per n-> +inf 0sservo che : 
% 
% - la Distribuzione di Probabilità Stimata Pn converge a P(x)
%
% - l' integrale di Pn -> 1
%
clc;clear;close all

for I=2:2:14
    %numero Campioni del Training Set
    n=2^I;
    kn=sqrt(n);

    %Estraggo i campioni da una distribuzione Normale N(0,1)
    TrainingSet=randn(1,n);
    TrainingSet=sort(TrainingSet);
    
    x=linspace(-3,3,40);
    [r c]=size(x);
    
    figure;
    subplot(1,2,1);
    plot(TrainingSet,zeros(1,n),'ob','MarkerSize',2,'MarkerFaceColor','b');
    hold on;    
    plot(x,zeros(1,c),'or','MarkerSize',2,'MarkerFaceColor','r')
    legend('Training Set','Punti x in cui vado a stimare P')
    title(['Numero di Campioni: ' num2str(n)])
    
    %Per ogni x  calcolo la distanza dai vettori del training Set
    Distanza=zeros(c,n);
    for i=1:c
        %calcolo la distanza dai vettori del training Set
        for j=1:n
            Distanza(i,j)=abs(x(i)-TrainingSet(j));
        end
    end
    %Ordina le righe di Distanza in modo crescente: 
    %   - I contiene gli indici
    %   - D contiene la Matrice delle Distanze Oridnate
    [D,I]=sort(Distanza,2);

    %Per ogni x in TestSet calcolo il Volume contente k-vicini 
    V=zeros(c,1);
    for i=1:c
        V(i)=2*D(i,kn); % calcola il volume in funzione dei k vicini
        DisegnaCirconferenza([x(i);0],V(i),20,1);
        %pause;
    end 
    axis equal;

    subplot(1,2,2)
    P=zeros(1,c);
    P=kn./(n.*V);
    plot(x,P)
    hold on;
    %Disegno la Gaussiana che sto stimando mediante KNN
    Media=0;
    Var=1;
    fun=@(x) exp( (-.5)*( ( (x-Media)/Var ).^2 ))./(sqrt(2*pi)*Var);
    t=linspace(-3,3,50);
    plot(t,fun(t),'r')
    
    %Verifico che la Funzione Stimata è una distribuzione di Probabilità:
    %Il suo integrale deve essere 1
    Fun=@(x) Distribuzione_KNN(x,TrainingSet);
    Integrale=quad(Fun,-10,10);
    
    legend(['Distribuzione P(x) Stimata Int=' num2str(Integrale)],'Distribuzione Originale')

    pause;
end
    
pause;
close all;