% Parzen_1D_3.m
% Esempio di Stima di Parzen su un set di Dati Monodimensionali n=10
% Applico:
%
% 1) finestra Rettangolare 
%
% 2) finestra Gaussiana
clc;clear; close all;

% Considero un Set di Campioni estratti da una Distribuzione di tipo Coseno
% genero 1 campione di 10 elementi in [0 2*pi]
x=2*pi*rand(1,10);
n=10;

%punti x in cui Valutare le finestre di Parzen
t=linspace(0,2*pi,100);
fun=@(x) (cos(x)+1)/(2*pi)

%punti x in cui Valutare le finestre di Parzen
y=fun(t);

h1=0.5;
hn=h1/sqrt(n);

for K=1:3
    %Applico il metodo di Parzen per la stima di P(x) con finestra Rettangolare
    %Fisso l' ampiezza della finestra di parzen
    h=hn/K;
    %Fisso il Volume 
    Vn=h;
        
    %Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
    figure;
    hold on;
    Pn=Parzen_Window(t,x,hn,Vn,'g')/n;
    plot(t,Pn,'b')
    plot(t,y,'r')
    legend(['Stima di Parzen con finestra Rettangolare di ampiezza h=' num2str(h)],'D. da stimare')
    plot(x,fun(x),'ob','MarkerSize',5,'MarkerFaceColor','b');
     
    pause
end

pause;
close all