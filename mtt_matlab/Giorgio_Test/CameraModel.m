%Modello Camera adottato nel PAPER

clc;close all; clear
 

K=@(f,uc,vc)[f 0 uc
             0 f vc
             0 0  1];



Rc2w =@(tx) [cos(tx)  sin(tx)  0; 
            -sin(tx)  cos(tx)  0;
                0          0   1];
    
Rw2c=@(tx)inv(Rc2w(tx));

%World to Camera
C2W=@(tx,Zc,T) Rc2w(tx)*(Zc-T);

%Camera to World
W2C=@(tx,Zw,T) Rw2c(tx)*Zw+T;

ProJ=@(Zw,f,uc,vc,tx,T) K(f,uc,vc )*W2C(tx,Zw,T) ;
%Fisso i PARAMETRI INTRINSECI
fo=1;
Uo=10;
Vo=10;
w=20; ku=1;
h=20;  kv=1;
%PARAMETRI ESTRINSECI
tx=1;
 
% Posizione del Centro Ottico C0
Xt=10;
Zt=10; %profondita
Ht=0;  %altezza della camera
T=[Xt Zt Ht].';
CO=T;

subplot(1,2,1); hold on
plot3(CO(1),CO(2),CO(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);

%Assi Mondo
Aw=100*eye(3);
plot3(0,0,0,'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
coloreM=[1 0 0];
colorX=[1 0 0];
colorY=[0 1 0];
colorZ=[0 0 1];
%Sistema di Riferimento Mondo centrato nell'origine
Vettore3D(Aw(1,1),Aw(2,1),Aw(3,1),colorX,1) %X in red
Vettore3D(Aw(1,2),Aw(2,2),Aw(3,2),colorY,1) %Z in green
Vettore3D(Aw(1,3),Aw(2,3),Aw(3,3),colorZ,1) %Y in blue
plot3(T(1),T(2),T(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);
xlabel('X')
ylabel('Y')
zlabel('Z') 
title('WORLD REFERENCE SYSTEM')

%Fissato un angolo di Rotazione vado a vedere la relazione che intercorre
%trai 2 sistemi di riferimento
%Verifico che RW2C e' una ROTAZIONE sull' asse Z
Ac=W2C(tx,Aw,repmat([0;0;0],1,3));
%Sistema di Riferimento Mondo centrato nell'origine
% Vettore3D(Ac(1,1),Ac(2,1),Ac(3,1),colorX,1) %X in red
% Vettore3D(Ac(1,2),Ac(2,2),Ac(3,2),colorY,1) %Y in green
% Vettore3D(Ac(1,3),Ac(2,3),Ac(3,3),colorZ,1) %Z in blue

%Disegno i Vettori applicati in C0
%La camera e allineata al piano di terra e l'asse focale e l'asse Y
AsseU=w*Ac(:,1)./norm(Ac(:,1));
AsseV=h*Ac(:,3)./norm(Ac(:,3));
%Centro Camera
CC=CO+fo*Ac(:,2)/norm(Ac(:,2));
plot3(CC(1),CC(2),CC(3),'oy','MarkerFaceColor','y','MarkerEdgeColor','y','MarkerSize',5);


DisegnaPiano(AsseU,AsseV,CC,1,colorZ);
%Normalizzo gli assi Camera
Ac(:,1)=Ac(:,1)/norm(Ac(:,1));
Ac(:,2)=Ac(:,2)/norm(Ac(:,2));
Ac(:,3)=Ac(:,3)/norm(Ac(:,3));
Ac=10*Ac;
Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX)
Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY)
Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),colorZ)
drawnow
axis equal

%Matrici di Rotazione tra Mondo e Camera e determinata dall' Orientazione
%della Camera.
%Considero un Punto Zm in Coordinate Mondo ed assumo che la camera sia
%puntata verso di Esso

%Disegno un punto in Coordinate Mondo
Zw=[60;230;30]
plot3(Zw(1),Zw(2),Zw(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
%Lina passante per il CO
Prj=[CO Zw];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')
Prj=[[0;0;0] CO];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','r','LineStyle','-')
drawnow

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

Pu=ProJ(Zw,fo,Uo,Vo,tx,T);
Pn=Pu./Pu(3);
subplot(1,2,2)
plot(Pn(1),Pn(2),'og','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',4); hold on

[x,y,z] = sphere(30);
[m,n]=size(x);
x=reshape(x,1,m*n);
y=reshape(y,1,m*n);
z=reshape(z,1,m*n)+10;
Zw=10*[x;y;z]+repmat(Zw,1,m*n);
 
x=reshape(Zw(1,:),m,n);
y=reshape(Zw(2,:),m,n);
z=reshape(Zw(3,:),m,n);
subplot(1,2,1)
surf(x,y,z)  % sphere centered at origin
axis equal

subplot(1,2,2)
Pu=ProJ(Zw,fo,Uo,Vo,tx,repmat(T,1,m*n));
Pn=[Pu(1,:)./Pu(3,:);Pu(2,:)./Pu(3,:)]
plot(Pn(1,:),Pn(2,:),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5)
axis([0 w*ku 0 h*kv])
axis ij
axis equal