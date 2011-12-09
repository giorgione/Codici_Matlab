% function Vettore3D(x,y,z,colore)
% Disegna un vettore in R3 di cordinate (x,y,z) applicato in nel punto 
% inziale Xo.
%
% parametri:
%
% - x: componente x del vettore che si desidera disegnare
%
% - y: componente y del vettore che si desidera disegnare
%
% - z: componente z del vettore che si desidera disegnare
%
% - colore: colore del Vettore
%
% Utilizzo:
%
% >> Vettore3D(10,10,10,[1 0 0])
% >> Vettore3D(10,0,0,[0 1 0])
% >> Vettore3D(0,10,0,[0 1 0])
% >> Vettore3D(0,0,10,[0 1 0])
%
% author: Giorgio Gemignani

function Vettore3D(x,y,z,colore,figH)
colormap(colore);
%vertice iniziale del vettore
P1=[0 0 0];

%vertice finale del vettore
P2=[x y z];

% normalizza il vettore P2
Pv=P2/norm(P2);

%Disegno il cono:

%calcola il punto dove parte la base del cono
Pstart=P2-Pv/2;

%calcola gli angoli di rotazione
[theta,phi,r]=cart2sph(P2(1),P2(2),P2(3));

%ruota il sistema di riferimento di theta gradi attorno all'asse Z
Rz=[cos(theta) -sin(theta) 0;
    sin(theta) cos(theta) 0;
    0            0        1];

%ruota il sistema di riferimento di pi/2-phi gradi attorno all'asse Y
Ry=[cos(pi/2-phi)  0  sin(pi/2-phi);
    0              1             0 ;
    -sin(pi/2-phi) 0  cos(pi/2-phi)];

%calcola il raggio del cono
raggio=norm(P2)/100;

%calcola il cono
[X,Y,Z]=cylinder([raggio 0]);
Z(2,:)=0.5;
[r c]=size(X);
%riorganizza i punti del cono in una matrice 3 x n
punti=[reshape(X,r*c,1)';reshape(Y,r*c,1)';reshape(Z,r*c,1)'];

%applica le operazioni di rotazione al cono
newpunti=Rz*Ry*punti;
 
Xnew=newpunti(1,:);Xnew=reshape(Xnew,r,c);
Ynew=newpunti(2,:);Ynew=reshape(Ynew,r,c);
Znew=newpunti(3,:);Znew=reshape(Znew,r,c);

% Disegna
figure(figH);
hold on;
P=[P1;P2];
x=P(:,1);
y=P(:,2);
z=P(:,3);
%plot3(x,y,z,'or','MarkerFaceColor','r');
line( x, y, z,'LineWidth',1.5,'Color',colore);

%trasla il cono nel verice del vettore e visualizza
surf(Xnew+Pstart(1),Ynew+Pstart(2),Znew+Pstart(3),'EdgeColor',colore);%,
%patch(Xnew(1,:)+Pstart(1),Ynew(1,:)+Pstart(2),Znew(1,:)+Pstart(3));

view(3);
box;
grid on;
