% function Compas3D(Xo,Yo,Zo,x,y,z)
% Disegna un vettore in R3 di cordinate (x,y,z) applicato in (Xo,Yo,Zo)
%
% author: Giorgio Gemignani

function Compas3D_Applicato(Xo,Yo,Zo,x,y,z,colore)

P1=[0 0 0];
P1=[Xo,Yo,Zo];
%vertice finale del vettore
P2=[x y z];

% normalizza il vettore
Pv=P2/norm(P2);

%calcola il punto dove parte la base del cono
Pstart=P2-Pv;

%calcola gli angoli di rotazione
[theta,phi,r]=cart2sph(P2(1),P2(2),P2(3));

%ruota il sistema di riferimento di theta gradi attorno all'asse Z
Rz=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];

%ruota il sistema di riferimento di pi/2-phi gradi attorno all'asse Y
Ry=[cos(pi/2-phi) 0 sin(pi/2-phi);0 1 0 ;-sin(pi/2-phi) 0 cos(pi/2-phi) ];

%calcola il raggio del cono
raggio=norm(P2)/100;

%calcola il cono
[X,Y,Z]=cylinder([raggio 0]);

[r c]=size(X);
%riorganizza i punti del cono in una matrice 3 x n
punti=[reshape(X,r*c,1)';reshape(Y,r*c,1)';reshape(Z,r*c,1)'];

%applica le operazioni di rotazione al cono
newpunti=Rz*Ry*punti;

Xnew=newpunti(1,:);Xnew=reshape(Xnew,r,c);
Ynew=newpunti(2,:);Ynew=reshape(Ynew,r,c);
Znew=newpunti(3,:);Znew=reshape(Znew,r,c);

% Disegna

hold on;
P=[P1;P2];
x=P(:,1);
y=P(:,2);
z=P(:,3);
%plot3(x,y,z,'or','MarkerFaceColor','r');
line(x,y,z,'LineWidth',1.5,'Color',colore);

%trasla il cono nel verice del vettore e visualizza
surf(Xnew+Pstart(1),Ynew+Pstart(2),Znew+Pstart(3),'EdgeColor',colore);
patch(Xnew(2,:)+Pstart(1)+Xo,Ynew(2,:)+Pstart(2)+Yo,Znew(2,:)+Pstart(3)+Zo,colore)
view(3);box;
grid on;
