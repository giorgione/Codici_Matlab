%Esempio stupido
%
% Genero 2 processi X1 e X2 con n1 ed n2 stati
clc;clear;
n1=4;

n2=3;
%Calcolo la prob congiunta P(X1,X2) 
Pjoint=rand(n1,n2);
Z=sum(sum(Pjoint));
Pjoint=Pjoint/Z;

% MARGINALIZZO la Congiunta e verifico che mi restituisce la distribuzione
%della singola variabile aleatoria
%           --------
%P(x1)=     \
%               P(X1,X2) =  P(X1|X2)*P(X2)
%           /
%           ------ X2
%
% risulterà P1 == P11
P11=zeros(1,n1);
P22=zeros(1,n2);

P1=zeros(1,n1);
P2=zeros(1,n2);

for x1=1:n1
       % La somma effettiva per ogni valore assunto da X1
       for x2=1:n2
        
            P1(x1)= P1(x1)+Pjoint(x1,x2);
            P2(x2)= P2(x2)+Pjoint(x1,x2);  
            
        end
end



%Probabilità condizionale P(X1 | X2) funzione di 2 variabili:
% X1 variabile aleatoria
% X2 è generalemnte un paramentro della distribuzione
% nb:
%   Ogni Colonna di Px1_given_x2 rappresenta la Distribuzione
%   P(X1|X2=i)
Px1_given_x2=zeros(n1,n2);
for x1=1:n1
    for x2=1:n2
        Px1_given_x2(x1,x2)=Pjoint(x1,x2)/P2(x2);
    end
end

Px2_given_x1=zeros(n1,n2);
for x2=1:n2
        for x1=1:n1
            Px2_given_x1(x1,x2) = Pjoint(x1,x2)/P1(x1);
    end
end

% MARGINALIZZO la Congiunta e verifico che mi restituisce la distribuzione
%della singola variabile aleatoria
%           --------
%P(x1)=     \
%               P(X1,X2) =  P(X1|X2)*P(X2)
%           /
%           ------ X2
%
% risulterà P1 == P11   --> sum(P1)=1
for x1=1:n1
       % La somma effettiva per ogni valore assunto da X1
       for x2=1:n2
         
            P11(x1)= P11(x1)+ Px2_given_x1(x1,x2)*P1(x1);
            P22(x2)= P22(x2)+ Px1_given_x2(x1,x2)*P2(x2);
        end
end


