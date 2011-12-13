clc;close all
%Assi Mondo
Aw=100*eye(3);

figure(1)
plot3(0,0,0,'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
coloreM=[1 0 0];
colorX=[1 0 0];
colorY=[0 1 0];
colorZ=[0 0 1];
%Sistema di Riferimento Mondo centrato nell'origine
Vettore3D(Aw(1,1),Aw(2,1),Aw(3,1),colorX,1) %X in red
Vettore3D(Aw(1,2),Aw(2,2),Aw(3,2),colorY,1) %Z in green
Vettore3D(Aw(1,3),Aw(2,3),Aw(3,3),colorZ,1) %Y in blue
plot3(0,0,0,'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);

xlabel('X')
ylabel('Y')
zlabel('Z') 
title('WORLD REFERENCE SYSTEM')
%Fisso la distanza focale a 10
f=1;

%%    ROTOTRASLAZIONE che lega sistema di riferimento CAMERA e MONDO
%
% Zc= R*Zm + T

% Posizione del Centro Ottico C0
Xth=80;
Zth=80;  %profondita
Hth=80;  %altezza della camera
T=[Xth Zth Hth].';
CO=T;
plot3(T(1),T(2),T(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);

%Matrici di Rotazione tra Mondo e Camera e determinata dall' Orientazione
%della Camera.
%Considero un Punto Zm in Coordinate Mondo ed assumo che la camera sia
%puntata verso di Esso

%Disegno un punto in Coordinate Mondo
Zw=[30;20;20]
plot3(Zw(1),Zw(2),Zw(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
%Lina passante per il CO
Prj=[CO Zw];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')
Prj=[[0;0;0] CO];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')

%punto Zw sul piano di terra
Zwo=[Zw-[0;0;Zw(3)]];
plot3(Zwo(1),Zwo(2),Zwo(3),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',4);
%Disegno le proiezioni sugli assi
Zx=[Zw(1); 0 ;0];
Zy=[0; Zw(2); 0];
Zz=[0 ;0; Zw(3)];
plot3(Zx(1),Zx(2),Zx(3),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',4);
plot3(Zy(1),Zy(2),Zy(3),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',4);
plot3(Zz(1),Zz(2),Zz(3),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',4);

%Disegno le proiezioni sugli assi
Prj=[Zwo Zx];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','g','LineStyle','-')
Prj=[Zwo Zy];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','g','LineStyle','-')
Prj=[Zwo Zw]; 
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','g','LineStyle','-')


%Traslo l'asse Ottico nell'origine in modo da calcolare il Vettore Normale 
%al piano Focale
N=-T;
%Vettore3D(N(1),N(2),N(3),colorX,1)
%Calcolo il Piano Ortogonale ad N: attraverso i 2 generatori del piano ed
%il punto per il quale esso passa (T)
% Calcolo i generatori
N=N/norm(N);
Ac=null(N.');
%Matrice di Rotazone Per allineare Mondo a Camera: i vettori colonna sono
%normalizzati
Rw2c=[Ac(:,1) -Ac(:,2) N];
% Vettore3D(R(1,1),R(2,1),R(3,1),colorX,1) %X in red
% Vettore3D(R(1,2),R(2,2),R(3,2),colorY,1) %Y in verde  
% Vettore3D(R(1,3),R(2,3),R(3,3),colorZ,1) %Z in blue

%Vettori generatori del Piano Immagine in Coordinate MONDO
Ac=10*Rw2c;
% Sistema di RIFERIMENTO CAMERA centrato in T
% Ac--> Assi Camera di lunghezza f: sono VETTORI APPLICATI nell'origine
%       Sommando T ottengo i punti sugli ASSI a distanza f dal Centro
%       Ottico

%Disegno i Vettori Ac applicati in T
Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX)
Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY)
Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),colorZ)

%Considero i Punti in CORDINATE MONDO sugli assi del Piano Focale
Ac1=Ac+repmat(T,1,3);
plot3(Ac1(1,1),Ac1(2,1),Ac1(3,1),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
plot3(Ac1(1,2),Ac1(2,2),Ac1(3,2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
plot3(Ac1(1,3),Ac1(2,3),Ac1(3,3),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);

%Disegno il Piano Immagine-Camera
w=10;
h=8;
%Centro della Camera: concide sul punto a distanza f da CO sul raggio
%OTTICO
CCam=Ac1(:,3);

%Assi del Piano immagine
AsseU=w*Ac(:,1)./norm(Ac(:,1));
AsseV=h*Ac(:,2)./norm(Ac(:,2));
DisegnaPiano(AsseU,AsseV,CCam);
%Disegno la box della camera
Prj=[CO CCam+AsseU+AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[CO CCam+AsseU-AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[CO CCam-AsseU+AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[CO CCam-AsseU-AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')

%Recupero la Matrice di Rotazione che Lega i Sistemi di Riferimento
% imponendo che:
% - l' asse Z sia concidente con N
% - l' asse X sia l'asse Orizzontale del piano immagine
% - l' asse Y sia l'asse Verticale del piano immagine
%                    -1
%   Zc=R*Zw+T   --> R  (Zc-T)=Zw
%
%ROTAZIONE CAMERA TO MONDO
Rc2w=inv(Rw2c);
ku=10;
kv=10;
%MATRICE INTRINSECI
K=[f*ku  0      w
   0    f*kv    h
   0     0      1];

%Proietto nel Piano Immagine
P=K*(Rw2c*Zw+T)
%Divido per z
P2D=P/P(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DISEGNO IL MODELLO IN COORDINATE CAMERA (camera centrata nello 0)
%
% Per fare questo devo applicare la trasformazione INVERSA alla CAMERA
% ovvero allineare la camera a l MONDO
%
%     1)   Sottrarre il Vettore T per centrare il Centro Ottico nell'origine
%          lavorando sui vettori del piano non mi serve
%    
%     2)    Ruotare gli assi della Camera in Modo che siano allineati
%           secondo il MODELLO Pinhole:
%               - l' asse Z sia concidente con N
%               - l' asse X sia l'asse Orizzontale del piano immagine
%                - l' asse Y sia l'asse Verticale del piano immagine
%            x   y  z
Assi =     [ 1   0  0;
             0   0  1;
             0  -1  0];
%Matrice di Rotazione che mi consente di disegnare la camera centrata
%nell'origine
R=Assi;
%Verifico che R ruota gli assi della camera e li dispone come richiesto
% dal modello pinhole
% figure(3)
% Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,3) %X in red
% Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,3) %Y in verde  
% Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,3) %Z in blue
% Ac=R*Ac;
% Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,3) %X in red
% Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,3) %Y in verde  
% Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,3) %Z in blue

    
%Assi Camera allineati con il sistema di Riferimento
Ac2=5*Assi;
    
figure(2);
hold on
title('CAMERA REFERENCE SYSTEM')
%L'origine del sistema di Riferimento Mondo si trovera ad una distanza
%norm(T) dal Centro Ottico lungo l'asse ottico
dis=norm(T)
Ow=[0;dis;0]

Ow1=R*(Rc2w*([0;0;0]-T));
plot3(Ow(1),Ow(2),Ow(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);
%Centro Ottico
Co1=[0;0;0];
Co2=R*(Rc2w*(T-T));
plot3(Co1(1),Co1(2),Co1(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);
%Centro Immagine
CI1=Co1+[0;f;0];
plot3(CI1(1),CI1(2),CI1(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);
Prj=[Ow Co1];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')

%Assi Sistema di Riferimento CAMERA:
%X : asse Orizzontale
%Y : Asse Verticale
%Z : Profondita'
Vettore3D(Ac2(1,1),Ac2(2,1),Ac2(3,1),colorX,2) %X in red
Vettore3D(Ac2(1,2),Ac2(2,2),Ac2(3,2),colorY,2) %Y in verde  
Vettore3D(Ac2(1,3),Ac2(2,3),Ac2(3,3),colorZ,2) %Z in blue
 
AsseU=w*Ac2(:,1)./norm(Ac2(:,1));
AsseV=h*Ac2(:,2)./norm(Ac2(:,2));
%Piano Immagine
DisegnaPiano(AsseU,AsseV,CI1);
%Disegno la box della camera
Prj=[Co1 CI1+AsseU+AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[Co1 CI1+AsseU-AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[Co1 CI1-AsseU+AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
Prj=[Co1 CI1-AsseU-AsseV];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')

%Considero gli Assi del sistema MONDO trasformati nel nuovo sistema di
%Riferimento
Ac1=R*Rc2w*(Aw-repmat(T,1,3));

plot3(Ac1(1,1),Ac1(2,1),Ac1(3,1),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
plot3(Ac1(1,2),Ac1(2,2),Ac1(3,2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
plot3(Ac1(1,3),Ac1(2,3),Ac1(3,3),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
%Disegno la box della camera
Prj=[Ow Ac1(:,1)];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')
Prj=[Ow Ac1(:,2)];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','g','LineStyle','-')
Prj=[Ow Ac1(:,3)];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
 
Zw1=R*Rc2w*(Zw-T);
%Disegno il punto Zw
plot3(Zw1(1),Zw1(2),Zw1(3),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
title('CAMERA REFERENCE SYSTEM')
Prj=[Co1 Zw1];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')

figure(3);hold on
plot(P2D(1),P2D(2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
[x,y,z] = sphere(30);
figure(2)
[m,n]=size(x);
x=reshape(x,1,m*n);
y=reshape(y,1,m*n);
z=reshape(z,1,m*n);
P=10*[x;y;z]+10
P1=R*Rc2w*(P-repmat(T,1,m*n));
x=reshape(P1(1,:),m,n);
y=reshape(P1(2,:),m,n);
z=reshape(P1(3,:),m,n);
surf(x,y,z)  % sphere centered at origin

figure(3)
P=K*(Rw2c*P+repmat(T,1,m*n));
u=P(1,:)./P(3,:);
v=P(2,:)./P(3,:);
plot(u,v,'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
axis([0 w*ku 0 h*kv])
axis ij