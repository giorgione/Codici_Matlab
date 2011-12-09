% syms r1 r2 r3 r4 r5 r6 r7 r8 r9 t1 t2 t3;
% R=[ r1 r2 r3; r4 r5 r6; r7 r8 r9];
% T=[t1  t2 t3].';
% 
% Z=[z1 z2 z3].';
% Rt= [R zeros(3,1)];
% Rt=[Rt ; [0 0 0 1]];
% Zn=Rt*Z+T

clc;close all
syms f uc vc x y z yc tx

X=[x y z 1].';

K=[f 0 uc
   0 f vc
   0 0  1]./z;

R=[1      0    0     ;   
   0 cos(tx) -sin(tx);
   0 sin(tx)  cos(tx)];

T=[0 yc 0].';

Estr=[R T];

P=K*Estr*X;
u=P(1);
v=P(2);

pretty(simplify(u))
pretty(simplify(v))
pretty(simplify(P(3)))


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
xlabel('X')
ylabel('Y')
zlabel('Z') 
title('WORLD REFERENCE SYSTEM')

%%   Distanza Camera - Origine Mondo
%
% la camera  e':
%
% TRASLATA di T rispetto al sistem a mondo
%
% RUOTATA di theta mediante R(theta) rispetto al sistema mondo (theta= pi/2)
Xth=80;
Zth=80;  %profondita
Hth=0;  %altezza della camera
T=[Xth Zth Hth].';

%Matrici di Rotazione tra Mondo e Camera

%Matrice di Rotazione 2D 
R=@(Theta)[ cos(Theta), -sin(Theta);
            sin(Theta), cos(Theta)];
%ROTAZIONE MONDO TO CAMERA        
Rw2c=[R(pi/2) [0;0]; 0 0 1];
%ROTAZIONE CAMERA TO MONDO
Rc2w=inv(Rw2c);

%Recupero gli Assi della Camera in coordinate Mondo applicando la ROTAZIONE
%attorno all' asse Y ed applicando la traslazione T
Ac=Rw2c*Aw;
Ac=Ac./4;
w=10;
h=8;
DisegnaPiano(w*Ac(:,1)./norm(Ac(:,1)),h*Ac(:,3)./norm(Ac(:,3)),T)
Vettore3D_Applicato(T,Ac(1,1),Ac(2,1),Ac(3,1),colorX)
Vettore3D_Applicato(T,Ac(1,2),Ac(2,2),Ac(3,2),colorY)
Vettore3D_Applicato(T,Ac(1,3),Ac(2,3),Ac(3,3),colorZ)

plot3(T(1),T(2),T(3),'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
Prj=[[0;0;0] T];
line(Prj(1,:),Prj(2,:),Prj(3,:),'Color','b','LineStyle','-')
%Disegno il Piano Camera
coloreC=[0 1 0];


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




%Calcolo Zw in coordinate Camera

figure(2);
title('CAMERA REFERENCE SYSTEM')
%Sistema di Riferimento Camera centrato nell'origine
plot3(0,0,0,'or','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',5);
coloreM=[1 0 0];
colorX=[1 0 0];
colorY=[0 1 0];
colorZ=[0 0 1];
%Sistema di Riferimento CAMERA centrato nell'origine
Vettore3D(Aw(1,1),Aw(2,1),Aw(3,1),colorX,2) %X in red
Vettore3D(Aw(1,2),Aw(2,2),Aw(3,2),colorY,2) %Z in green
Vettore3D(Aw(1,3),Aw(2,3),Aw(3,3),colorZ,2) %Y in blue
DisegnaPiano(Aw(:,1),Aw(:,3),[0;0;0])

%Disgeno Zw nwl sistema di riferimento della CAMERA
Ztemp=Zw-T;
Zc=Rc2w*Ztemp;
plot3(Zc(1),Zc(2),Zc(3),'or','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',5);

%Eseguo la stessa Operazione sugli Assi nel sistema di Riferimento Mondo
Ztemp=Ac;
Zc=Rc2w*Ztemp;
Vettore3D(Zc(1,1),Zc(2,1),Zc(3,1),colorX,2) %X in red
Vettore3D(Zc(1,2),Zc(2,2),Zc(3,2),colorY,2) %Z in green
Vettore3D(Zc(1,3),Zc(2,3),Zc(3,3),colorZ,2) %Y in blue

%Ruoto la Camera per vedere se si allinea al SISTEMA DI RIFERIMENTO MONDO
RuotaPiano([0;0;10],[0;10;0],Ztemp,Rw2c)
