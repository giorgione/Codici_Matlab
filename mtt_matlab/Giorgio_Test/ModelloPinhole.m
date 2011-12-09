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

%%    ROTOTRASLAZIONE che lega sistema di riferimento CAMERA e MONDO
%
% Zc= R*Zm + T

% Posizione del Centro Ottico
Xth=80;
Zth=80;  %profondita
Hth=80;  %altezza della camera
T=[Xth Zth Hth].';
plot3(T(1),T(2),T(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);

%Matrici di Rotazione tra Mondo e Camera e determinata dall' Orientazione
%della Camera.
%Considero un Punto Zm in Coordinate Mondo ed assumo che la camera sia
%puntata verso di Esso

%Disegno un punto in Coordinate Mondo
Zw=[30;20;20]
plot3(Zw(1),Zw(2),Zw(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
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

%%Disegno la la linea tra 0w e C: Asse OTTico
%Prj=[[0;0;0] T];
%line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','m','LineStyle','-')

%Traslo l'asse Ottico nell'origine in modo da calcolare il Vettore Normale 
%al piano Focale
N=T;
Vettore3D(N(1),N(2),N(3),colorX,1)
%Calcolo il Piano Ortogonale ad N: attraverso i 2 generatori del piano ed
%il punto per il quale esso passa (T)
% Calcolo i generatori
Ac=null(N.');
Ac=[Ac -N/norm(N)];
% Disegno il piano passante per T di generatori V
% DisegnaPiano(Ac(:,1),Ac(:,2),T)
Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX)
Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY)
Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),colorZ)
% Verifico che la Rotazione Rw2c ruota gli assi nel mondo allineandoli con
% quelli della Camera
%Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,1)
%Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,1)
%Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,1)
%Considero i Punti in CORDINATE MONDO sugli assi del Piano Focale
Ac1=Ac+repmat(T,1,3);

plot3(Ac1(1,1),Ac1(2,1),Ac1(3,1),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
plot3(Ac1(1,2),Ac1(2,2),Ac1(3,2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
plot3(Ac1(1,3),Ac1(2,3),Ac1(3,3),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
%Recupero la Matrice di Rotazione che Lega i Sistemi di Riferimento
% imponendo che:
% - l' asse Z sia concidente con N
% - l' asse X sia l'asse Orizzontale del piano immagine
% - l' asse Y sia l'asse Verticale del piano immagine
%                    -1
%   Zc=R*Zw+T   --> R  (Zc-T)=Zw
%
Rw2c=Ac;

%ROTAZIONE CAMERA TO MONDO
Rc2w=inv(Rw2c);

%Disegno gli assi nel SISTEMA di riferimento CAMERA trasformando i
%corrispettivi nel SISTEMA MONDO
figure(2);hold on
title('CAMERA REFERENCE SYSTEM')

%Considero gli Assi del sistema di RIFERIMENTO
Ac1=Rw2c+repmat(T,1,3);
plot3(Ac1(1,1),Ac1(2,1),Ac1(3,1),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);
plot3(Ac1(1,2),Ac1(2,2),Ac1(3,2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',5);
plot3(Ac1(1,3),Ac1(2,3),Ac1(3,3),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);

%Vettori applicati all'origine senza Traslare
Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,2)
Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,2)
Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,2)
%Assi del Sistema di riferimento Camera
Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX)
Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY)
Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),colorZ)

%Disegno l'origine
plot3(T(1),T(2),T(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);

%Disegno il punto Zw nel sistema di riferimento CAMERA
Zc=Rw2c*Zw+T;
plot3(Zc(1),Zc(2),Zc(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);

title('CAMERA REFERENCE SYSTEM')


