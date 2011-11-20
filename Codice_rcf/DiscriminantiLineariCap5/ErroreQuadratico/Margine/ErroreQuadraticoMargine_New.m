%function [J DJ]=ErroreQuadraticoMargine(a,y,b)
% Calcola la funzione Percettrone J(a,y) ed il suo gradiente DJ(a,y)
% dove:
%
% -a :vettore dei pesi (Spazio in cui vado a disegnare la funzione)
%
% -y :Vettore dei Pattern (Training Set)
function [J,DJ]=ErroreQuadraticoMargine_New(a,y,b)
%Valuta i vettori missclassified per ogni punto dello spazio dei Pesi
Jm=a.'*y-b;
[m,n]=size(Jm)

%Memorizzo in Ymiss gli indici dei Vettori y missclassified
Ymiss=zeros(m,n);
for I=1:m
    for J=1:n
        %errore di classificazione
        if Jm(I,J)<=0
            Ymiss(I,J)=J;
        end
    end
end

%Calcola la funzione ErroreQuadratico
for i=1:m
    %Estrai i vettori missclassified per ogni punto a(a1,a2)
    [R C Val]=find(Ymiss(i,:));
    if isempty(Val)
        J(i)=0;
        DJ(:,i)=[0;0];
    else
    ym=y(:,Val);
    n=length(ym(1,:));
    ValJ=0;
    for j=1:n
        ymis=ym(:,j);
        ValJ=ValJ+((a(:,i).'*ymis-b(i,j))*(a(:,i).'*ymis-b(i,j)).')/norm(ymis,2)^2;
    end
    J(i)=ValJ;
    
    % 
    %Calcolo il Gradiente di F(a):
    %     t         2
    %   (a  * y - b) 
    % ----------------
    %        t 
    %       y  * y
    %
    % J(a):
    %      t          
    %    (a  * y) - b 
    %  ----------------- * y
    %       t 
    %      y  * y
    
    Coef=2*a(:,i).'*ym -b(i);
    nc=length(Coef);
    S=0;
    for I=1:nc
        S=S+Coef(I)*ym(:,I)/norm(ym(:,I));
    end
    DJ(:,i)=S;
    
    
end
end
