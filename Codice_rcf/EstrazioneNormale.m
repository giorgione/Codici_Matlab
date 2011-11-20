%function Data=EstrazioneNormale(Media,Var,n);
%
% Simula l' estrazione di n campioni da una Distribuzione Normale
% N(Media,Var)
%
% Idea: 
%
% Genero i Dati in maniera tale che:
%
%    - il 80% si trova nell' intervallo [Media-Var Media+Var]
%
%    - il 15% si trova nell' intervallo [Media-2*Var Media+2*Var]
%
%    - il 5% si trova nell' intervallo [Media-3*Var Media+3*Var]

function Data=EstrazioneNormale(Media,Var,n);

n1=round(n*0.80);
n2=round(n*0.15);
n3=round(n*.05);

a=Media-Var;
b=Media+Var;
Data1=a+(b-a)*rand(1,n1);

a=Media+Var;
b=Media+2*Var;
Data2=a+(b-a)*rand(1,ceil(n2/2));

a=Media-Var;
b=Media-2*Var;
Data2=[Data2 a+(b-a)*rand(1,floor(n2/2))];


a=Media+2*Var;
b=Media+3*Var;
Data3=a+(b-a)*rand(1,ceil(n3/2));

a=Media-2*Var;
b=Media-3*Var;
Data3=[Data3 a+(b-a)*rand(1,floor(n3/2))];

Data=[Data1 Data2 Data3];