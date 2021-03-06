function RuotaPiano(a,b,Xo,R)
%Disegna lo span(a,b): spazio generato dai Vettori a e b
%Piano Per Xo=(x1,y2,z3)  di direttori (a,b)
% x = 	x1 	+ 	a1u 	+ 	b1t
% y = 	y2 	+ 	a2u 	+ 	b2t
% z = 	z3 	+ 	a3u 	+ 	b3 
[u,v]=meshgrid(-1:1);
[m,n]=size(u);
Points=m*n;
u=reshape(u,1,Points);
v=reshape(v,1,Points)

X=Xo(1)+a(1)*u+b(1)*v;
Y=Xo(2)+a(2)*u+b(2)*v;
Z=Xo(3)+a(3)*u+b(3)*v;
P=[X;Y;Z];
P=R*P;

X=reshape(P(1,:),m,n);
Y=reshape(P(2,:),m,n);
Z=reshape(P(3,:),m,n);
%Piano Generato dai Vettori a,b
mesh(X,Y,Z)