% Parzen_1D_4.m
% Esempio di Stima di Parzen su un set di Dati Monodimensionali 
% Applico:
%
% 1) finestra Gaussiana 
%
% Idea:
%
% for I=0:4:8
% 
%       1) Estraggo n=2^I campioni da una Distribuzione Normale N(0,1)
%
%       2) Eseguo l' algoritmo di Parzen con 4 valori di h per stimare N
%
%       3) Visualizzo i risultati.
%
% end
%
% per n-> +inf 0sservo che : 
% 
% - la Distribuzione di Probabilità Stimata Pn converge a P(x)
%
% - l' integrale di Pn -> 1
%

clc;clear; close all;
% Genero 1 campione dalla distribuzione monodimensionale Normale N(0,1)

n=10;

Media=0;
Var=1;
fun=@(x) exp( (-.5)*( ( (x-Media)/Var ).^2 ))./(sqrt(2*pi)*Var);
%punti x in cui Valutare le finestre di Parzen
t=linspace(-5,5,30);
y=fun(t);

h1=1;
for N=0:4:8
    %Fisso il numero di Campioni
    n=2^N;
    
    %Estraggo i campioni
    x=randn(1,n);
    
    hn=h1/sqrt(n);
    
    %Stima di Parzen con n campioni e 3 ampiezze distinte della Finestra
    figure;
    display(['Stima di Parzen con ' num2str(n) ' campioni'])
    for K=1:4

        h=hn/K;
        %Fisso il Volume 
        Vn=h;

        %Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
        Pn=Parzen_Window(t,x,h,Vn,'g');
        
        %disegno
        subplot(2,2,K)
        plot(t,Pn,'b');hold on;drawnow;
        plot(t,y,'r');drawnow
        
        %Verifico che l' integrale della Distribuzione Stimata sia sempre 1
        P=Parzen_Window_Simbolico(x,h,Vn);
        I=double(int(P,-10,10))
        
        legend(['Stima di Parzen (N.C=' num2str(n) ') con Gaussiana di Varianza h=' num2str(h) ' Integrale: ' num2str(I)],'D. normale')
        plot(x,fun(x),'ob','MarkerSize',5,'MarkerFaceColor','b');drawnow;


    end
    display('fine')
    pause
end

pause;
close all