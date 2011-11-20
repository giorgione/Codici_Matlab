%function [J DJ]=Perceptron(a,y)
% Calcola la funzione Percettrone J(a,y) ed il suo gradiente DJ(a,y)
%
%          
%           -----          t     
% J(a,y) =    \        (-a) * y 
%           /
%           -----
%           y appartente a Y --> Insieme dei Vettori missclassified
%
% paremtri in:
%
% -a :vettore dei pesi (Spazio in cui vado a disegnare la funzione)
%
% -y :Vettore dei Pattern (Training Set)
%
% paremtri out:
%
% - J : valore della Funzione Percettorne nel punto (a,y)
%
% - DJ : valore del Gradiente della Funzione Percettorne nel punto (a,y)
%

function [J DJ]=Perceptron(a,y)
%                           t
%Calcolo la trasformazione a * y 
Jm=a.'*y;
[m,n]=size(Jm);
[M,N]=size(y);
%Memorizzo in Ymiss gli indici dei Vettori y missclassified, vettori y tali
%che:
%        t
%       a * y <= 0
Ymiss=zeros(m,n);
for I=1:m
    for J=1:n
        %errore di classificazione
        if Jm(I,J)<=0
            Ymiss(I,J)=J;
        end
    end
end

%Calcola la funzione Percettrone
%
for i=1:m
    %Estrai i vettori missclassified per ogni punto a(a1,a2)
    [R C Val]=find(Ymiss(i,:));
    if isempty(Val)
        J(i)=0;
        DJ(:,i)=zeros(M,1);%[0;0];
    else
        ym=y(:,Val);
        J(i)=sum(-a(:,i).'*ym,2);
        DJ(:,i)=sum(-ym,2);
    end
end

