%function y=rettangolo(x,h)
% Valuta la finestra Gaussiana di Ampiezza h nei punti x
function F=Parzen_2d(X_tr,hn,Vn)

syms x y real;

F=0;
[d n]=size(X_tr);
Sigma=eye(2);
M=[0;0];
for I=1:n 
    X=[(x - X_tr(1,I))/hn;
       (y - X_tr(2,I))/hn];
   
    F=F+GaussianaMulti(Sigma,M,X)./Vn;
end