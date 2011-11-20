% EM_Stima_1.m
%
% Stima dei parametri di un Modello Gaussiano attraverso l' Algoritmo EM
% supponendo che i dati siano Mancanti.
%
% Considero il seguente problema in R2.
%
% 1) Considero l' insieme di osservazioni estratte da una Distribuzione
%    Gaussiana in R2 avente parametri incogniti T=[u1;u2;o1;o2];:
%
%    Data=[0 1 2 x4;
%          2 0 2 4];
%
%   Il valore x4 indica che il dato è Mancante-Corrotto 
%
% 2) Risolvo il problema della Stima a Massima Verosimiglianza di T 
%    mediante l' algoritmo EM.
%
% 3) Considero la Stima iniziale dei Parametri Teta0=[0;0;1;1]
%
% 4) Applico EM andando a disegnare ad ogni passo l' Ellisse definita dai
%    parametri STIMATI.

clc;clear;close all;
% carattere di newline
s=sprintf('\n');

%siamo in R2
Teta0=[0;0;1;1];
%
syms x4 u1 u2 o1 o2 real
TetaS=[u1;u2;o1;o2];

% Matrice delle osservazioni (i-esima colonna --> osservazione xi )
Data=[0 1 2 x4;
      2 0 2 4];
Data_G=double(Data(:,1:3));
Data_B=Data(:,4);

%Disegno i 3 punti
plot(Data_G(1,1),Data_G(2,1),'or','MarkerFaceColor','r');hold on;
plot(Data_G(1,2),Data_G(2,2),'ob','MarkerFaceColor','b');
plot(Data_G(1,3),Data_G(2,3),'og','MarkerFaceColor','g');


%Definisco i settaggi dell' OTTIMIZAZIONE
options = optimset('Display','final','MaxFunEvals',500,'MaxIter',500,'TolFun',1e-4,'TolX',1e-4);


% punto iniziale dal quale inizia la ricerca del Massimo Locale della
% Funzione di Verosimiglianza dei dati completi
T_r(:,1)=rand(1,4);

% Stima Corrente dei Parmaetri
T_i=Teta0;
T(:,1)=Teta0;


%Fisso i parametri di Arresto dell' Algoritmo EM: Tolleranza e Massimo
%numero di Iterazioni.
Tol=1e-5;
Maxiter=10;
k=1;
Risultato=1;
while ( k<Maxiter && Risultato>Tol)
    
    %Disegna l' Ellisse definita dai parametri stimati
    [x,y]=DisegnaEllisse(T(1:2,k),T(3,k),T(4,k),40,0);
    plot(x,y,'Color',[1/k 1 0]);hold on
    

    %E-STEP: calcolo la funzione di Verosimiglianza dei Dati Completi
    %Calcolo la Distribuzione Marginalizzata rispetto ai DATI INCOMPLETI
    %data la STIMA CORRENTE dei PARAMETRI
    K=double(int(GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],Data_B),'x4',-inf,inf));
    M_Step=@(T) -ValoreAtteso_EM_ottimizzata(T,K,T_i,Data_G,Data_B,'n');
    
    %M-STEP: Massimizzo
    [T(:,k+1),Val(k+1),exitflag,output] = fminsearch(M_Step,[T_r(1);T_r(1);sqrt(T_r(3));sqrt(T_r(4))],options);
    
    T(3:4,k+1)=T(3:4,k+1).^2;
    
    display('Stima Corrente:')
    disp(T(:,k))
    display('Stima Raffinita:')
    disp(T(:,k+1))
    disp(s);
    
    pause;
    k=k+1;
    
    %Aggiorna la stima corrente
    T_i=T(:,k);
    
    %Aggiorna il punto iniziale dal quale parte la ricerca del Minimo
    T_r=T_i;
    
    %Calcolo lo scostamento tra 2 stime successive
    Risultato=norm(T(:,k)-T(:,k-1));
end