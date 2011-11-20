%Esempio di Modello Grafico P(x1,x2,x3) con connessione seriale
%
% x1 --> x2 --> x3
clc;clear;close all;
display('Seriale: x1 --> x2 --> x3')
Mxy=zeros(3);
Mxy(1,2)=1;
Mxy(2,3)=1;

for i=1:3
    Parents=GetParents(i,Mxy);
    display(['Pa(x' num2str(i) ')=']);
    disp(Parents);
    
    Child=GetChilds(i,Mxy);
    display(['Ch(x' num2str(i) ')=']);
    disp(Child);
end

Pjoint=Fattorizza(Mxy)

% Connessione divergente
% x1 <-- x2 --> x3
display('Divergente: x1 <-- x2 --> x3')
Mxy=zeros(3);
Mxy(2,1)=1;
Mxy(2,3)=1;

for i=1:3
    Parents=GetParents(i,Mxy);
    display(['Pa(x' num2str(i) ')=']);
    disp(Parents);
    
    Child=GetChilds(i,Mxy);
    display(['Ch(x' num2str(i) ')=']);
    disp(Child);
end
Pjoint=Fattorizza(Mxy)

% Connessione convergente
% x1 --> x2 <-- x3
display('Convergente: x1 <-- x2 --> x3')
Mxy=zeros(3);
Mxy(1,2)=1;
Mxy(3,2)=1;

for i=1:3
    Parents=GetParents(i,Mxy);
    display(['Pa(x' num2str(i) ')=']);
    disp(Parents);
    
    Child=GetChilds(i,Mxy);
    display(['Ch(x' num2str(i) ')=']);
    disp(Child);
end
Pjoint=Fattorizza(Mxy)