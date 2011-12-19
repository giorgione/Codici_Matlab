% function Vettore3D_Applicato(x,y,z,colore)
% Disegna un vettore in R3 di cordinate (x,y,z) applicato in nel punto 
% inziale Xo.
%
% parametri:
%
% - Xo: punto inziale (vettore colonna)
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
% >> Vettore3D_Applicato([1;1;1],10,10,10,[1 0 0])
% >> Vettore3D_Applicato([1;1;1],10,0,0,[0 1 0])
% >> Vettore3D_Applicato([1;1;1],0,10,0,[0 1 0])
% >> Vettore3D_Applicato([1;1;1],0,0,10,[0 1 0])
%
% author: Giorgio Gemignani

function Vettore3D_Applicato(Xo,x,y,z,colore,myFig)
if nargin ==6
    figure(myFig)
end
colormap(colore);
%vertice iniziale del vettore
P1=[0 0 0];

%vertice finale del vettore
P2=[x y z];

% normalizza il vettore P2
Pv=P2/norm(P2);

%Disegno il cono:

%calcola il punto dove parte la base del cono
Pstart=P2-Pv+Xo.';

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

%Matrice di Traslazione in Xo in coordinate omogenee
T=eye(4);
T(1:3,4)=Xo;

%calcola il raggio del cono
raggio=norm(P2)/100;

%calcola il cono
[X,Y,Z]=cylinder([raggio 0]);
%Z(2,:)=0.5;
[r c]=size(X);
%riorganizza i punti del cono in una matrice 3 x n
punti=[reshape(X,r*c,1)';reshape(Y,r*c,1)';reshape(Z,r*c,1)'];

%applica le operazioni di rotazione al cono
newpunti=Rz*Ry*punti;

%applico la traslazione
newpunti=[newpunti ;zeros(1,r*c)];
newpunti=T*newpunti;

Xnew=newpunti(1,:);Xnew=reshape(Xnew,r,c);
Ynew=newpunti(2,:);Ynew=reshape(Ynew,r,c);
Znew=newpunti(3,:);Znew=reshape(Znew,r,c);

% Disegna

hold on;
P=[P1+Xo.';P2+Xo.'];
x=P(:,1);
y=P(:,2);
z=P(:,3);
%plot3(x,y,z,'or','MarkerFaceColor','r');
line( x, y, z,'LineWidth',1.5,'Color',colore);

%trasla il cono nel verice del vettore e visualizza
surf(Xnew+Pstart(1),Ynew+Pstart(2),Znew+Pstart(3),'EdgeColor',colore);%,

view(3);
box;
grid on;
