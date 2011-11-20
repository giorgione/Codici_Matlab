% Esempio di Stima di Parzen su un set di Dati Monodimensionali
% Applico:
%
% 1) finestra Rettangolare 
%
% 2) finestra Gaussiana
clc;clear; close all;
% Genero un campione monodimensionale di 10 elementi
x1=3*rand(1,5);
x2=6*rand(1,5);
x=[x1 x2];
n=10;

t=linspace(-1,7,1000);
for K=1:3
    %Applico il metodo di Parzen per la stima di P(x) con finestra Rettangolare
    %Fisso l' ampiezza della finestra di parzen
    
    hn=0.5/K;
    %Fisso il Volume 
    Vn=hn;
    
    Pn=Parzen_Window(t,x,hn,Vn,'r');
    figure;
    hold on
    plot(t,Pn)
    plot(x,zeros(1,10),'ob','MarkerSize',5,'MarkerFaceColor','b');
    title(['Stima di Parzen con finestra Rettangolare di ampiezza ' num2str(hn)])

    %Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
    figure;
    plot(x,zeros(1,10),'ob','MarkerSize',5,'MarkerFaceColor','b');
    hold on;
    Pn=Parzen_Window(t,x,hn,Vn,'g');
    plot(t,Pn)
    title(['Stima di Parzen con finestra Rettangolare di ampiezza ' num2str(hn)])
    pause;
end
