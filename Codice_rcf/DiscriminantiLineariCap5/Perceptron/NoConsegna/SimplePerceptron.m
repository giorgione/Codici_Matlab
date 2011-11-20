%function [J DJ]=SimplePerceptron(a,y)
% Calcola la funzione Percettrone Semplice J(a,y) ed il suo gradiente DJ(a,y)
%
%          
%           -----       t
% J(a,y) =    \       (-a) * y
%           /
%           -----
%           y appartente a Y --> Insieme dei Vettori missclassified
%
% parametri in:
%
% - a :vettore dei pesi (Spazio in cui vado a disegnare la funzione)
%
% - y :Vettore dei Pattern (Training Set)
%
% paremtri out:
%
% - J : valore della Funzione Percettorne nel punto (a,y)
%
% - DJ : valore del Gradiente della Funzione Percettorne nel punto (a,y)
%

function [J DJ]=SimplePerceptron(a,y)
%Valuta i vettori missclassified per ogni punto dello spazio dei Pesi
% attraverso il classificatore lineare
Jm=a.'*y;
[m,n]=size(Jm);

%Memorizzo in Ymiss gli indici dei Vettori y missclassified dal
%classificatore lineare
Ymiss=zeros(m,n);
for I=1:m
    for J=1:n
        %errore di classificazione
        if Jm(I,J)<=0
            Ymiss(I,J)=J;
        end
    end
end

%Calcola la funzione Percettrone Semplice:
for i=1:m
    %Estrai i vettori missclassified per ogni punto a(a1,a2)
    [R C Val]=find(Ymiss(i,:));
    if isempty(Val)
        J(i)=0;
        DJ(:,i)=[0;0];
    else
    ym=y(:,Val)
    J(i)=-a(:,i).'*ym(:,1);
    DJ(:,i)=-ym(:,1);
    end
end
