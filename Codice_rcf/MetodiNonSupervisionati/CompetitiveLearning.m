%Esempio di Competitive Learning su un Set di dati 2-D ed C=4 clusters
%
clc;clear;close all;
x1=-2+randn(1,10);
y1=-2+randn(1,10);

x2=2+randn(1,10);
y2=2+randn(1,10);

%Insericsco i dati in un unico Vettore
Dati=[x1 x2;y1 y2];
n=2*length(x1);

plot(x1,y1,'or','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(x2,y2,'or','MarkerSize',5,'MarkerFaceColor','r')
title('Campioni nello Spazio Originale')

%Spazio Augmented
Dati=[ones(1,n);Dati];
%Normalizzo i Dati
for i=1:n
    Dati(:,i)=Dati(:,i)/norm(Dati(:,i));
end

figure;
plot3(Dati(1,:),Dati(2,:),Dati(3,:),'or','MarkerSize',5,'MarkerFaceColor','r');
title('Campioni nello Spazio Augmented')

%Inizializzo i Pesi della Neurone 1
w1=randn(1,3);
w1=w1/norm(w1);
w1=[0 0 0;w1];
%Inizializzo i Pesi della Neurone 2
w2=randn(1,3);
w2=w2/norm(w2);
w2=[0 0 0;w2];

W=[w1(2,:) ;w2(2,:)];
K=2;

alpha=@(x) .4*x;
tol=1e-3;
maxiter=200;
while(norm(w1(K,:)-w1(K-1,:)) >tol || norm(w2(K,:)-w2(K-1,:)) >tol && K<maxiter)
    IndiciRand=randperm(n);
    for i=1:n
        I=IndiciRand(i);
        X=Dati(:,I);

        Z=W*X;
        %Selezione il Neourone avente massima risposta
        [Val,J]=max(Z);
        %Memorizzo in Cluster il risultato della Classificazione
        ClusterClusterOut(K-1,I)=J;
        
        W(J,:)=W(J,:)+alpha(K)*X.';
        W(J,:)=W(J,:)/norm(W(J,:));

        w1(K+1,:)=W(1,:);
        w2(K+1,:)=W(2,:);
    end
    K=K+1;
end

figure;
DisegnaClusters3D(Dati(1,:),Dati(2,:),Dati(3,:),ClusterClusterOut(end,:))
title('Campioni Clusterizzati nello spazio Augmented')
plot3(w1(end,1),w1(end,2),w1(end,3),'or','MarkerSize',5,'MarkerFaceColor','b');
plot3(w2(end,1),w2(end,2),w2(end,3),'or','MarkerSize',5,'MarkerFaceColor','g');

figure;
DatiNew=[x1 x2;y1 y2];
DisegnaClusters(Dati(1,:),Dati(2,:),ClusterClusterOut(end,:))
title('Campioni nello Spazio Originale')



