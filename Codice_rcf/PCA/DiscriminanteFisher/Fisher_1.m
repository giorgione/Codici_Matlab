clc;clear ;close all
% carattere di newline
s=sprintf('\n');

%Genero un un Set di N dati 3-Dimensionali 
N=20;

a=10;
b=25;
x1=a+(b-a)*randn(1,N);
x2=a+(b-a)*randn(1,N);

X=[x1;x2];

D1=X(:,1:10);
D2=X(:,11:end);

Media_D1=mean(D1,2);
Media_D2=mean(D2,2);


D1=X(:,1:10);
D2=X(:,11:end);

%Disegno i dati di D1
plot(D1(1,:),D1(2,:),'or','MarkerFaceColor','r');hold on;

%Disegno la Media di D1, intesa come rappresentazione sintetica
plot(Media_D1(1),Media_D1(2),'og','MarkerFaceColor','g');
grid on;

%Disegno i dati di D2
plot(D2(1,:),D2(2,:),'ob','MarkerFaceColor','b');hold on;

%Disegno la Media di D2, intesa come rappresentazione sintetica
plot(Media_D2(1),Media_D2(2),'og','MarkerFaceColor','g');
grid on;

%Disegno la linea congiungente le Medie
line([Media_D1(1) Media_D2(1)],[Media_D1(2) Media_D2(2)])
pause
%Il discriminante di Fisher Lineare ricerca la retta di direttori w che meglio 
%separa le proiezioni y dei Dati x:
%      t
% y = w  * x
%
%

%Calcolo la Matrici di Scatter per D1 e D2
S1=MatriceScatter(D1,Media_D1);
S2=MatriceScatter(D2,Media_D2);

%Calcolo la Matrice Whitin-Scatter
Sw=S1+S2

%Calcolo la Matrice Between-Scatter
Sb=(Media_D1-Media_D2)*(Media_D1-Media_D2).';

%Calcolo il Direttore w della retta che separa meglio le Proiezioni
w=inv(Sw)*(Media_D1-Media_D2);
%normalizzo il vettore w
w=w/norm(w,2);

Retta=[w*50  -w*50];
%Disegno la Retta
plot(Retta(1,:),Retta(2,:),'LineWidth',2);


%Calcolo le proiezioni
y1=w.'*D1;
y2=w.'*D2;

for i=1:10
    Proj_y1(:,i)=y1(i)*w;
    Proj_y2(:,i)=y2(i)*w;
end


%Disegno le proizioni dei punti sulla Retta
%figure;
plot(Proj_y1(1,:),Proj_y1(2,:),'or','MarkerFaceColor','r','MarkerSize',5);hold on
plot(Proj_y2(1,:),Proj_y2(2,:),'ob','MarkerFaceColor','b','MarkerSize',5);
%Vettore3D(w(1),w(2),0,[1 0 0])

%Disegno le rette congiungenti i punti e le proiezioni sul vettore w
Proj=[];
for i=1:N/2
    
    Punti=[D1(:,i) Proj_y1(:,i)];
    l=line(Punti(1,:),Punti(2,:));
    set(l,'Color','r','LineWidth',1)
    
    Punti=[D2(:,i) Proj_y2(:,i)];
    l=line(Punti(1,:),Punti(2,:));
    set(l,'Color','b','LineWidth',1)
end

w.'*Proj_y1
w.'*Proj_y2
