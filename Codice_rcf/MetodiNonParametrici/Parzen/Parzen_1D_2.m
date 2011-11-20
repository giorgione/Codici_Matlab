% Esempio di Stima di Parzen su un set di Dati Monodimensionali n=10
% Applico:
%
% 1) finestra Rettangolare 
%
% 2) finestra Gaussiana
%
% Idea:
%
%      1) Estraggo n=10 campioni da una Distribuzione Normale N(0,1)
%
%       for K=1:3
%
%           2) Eseguo l' algoritmo di Parzen con h=hn/K per stimare N
%           
%           3) Verifico che l'integrale della Distribuzione Stimata � 1
%           
%           4) Visualizzo i risultati.
%
%       end


clc;clear; close all;
% Genero 1 campione dalla distribuzione monodimensionale Normale N(0,1)
x=randn(1,10);
n=10;

Media=0;
Var=1;
fun=@(x) exp( (-.5)*( ( (x-Media)/Var ).^2 ))./(sqrt(2*pi)*Var);
%punti x in cui Valutare le finestre di Parzen
t=linspace(-5,5,100);
y=fun(t);

h1=2;
hn=h1/sqrt(n);
for K=1:3
    %Applico il metodo di Parzen per la stima di P(x) con finestra Rettangolare
    %Fisso l' ampiezza della finestra di parzen
    h=hn/K;
    %Fisso il Volume 
    Vn=h;
    
    Pn=Parzen_Window(t,x,h,Vn,'r');
    figure;
    hold on
    plot(t,Pn,'b')
    plot(t,y,'r')
    legend(['Stima di Parzen con finestra Rettangolare di ampiezza h=' num2str(h)],'D. normale')
    plot(x,fun(x),'ob','MarkerSize',5,'MarkerFaceColor','b');
    
    %Applico il metodo di Parzen per la stima di P(x) con finestra Gaussiana
    figure;
    hold on;
    Pn=Parzen_Window(t,x,h,Vn,'g');
    plot(t,Pn,'b')
    plot(t,y,'r')
    legend(['Stima di Parzen con finestra Rettangolare di ampiezza h=' num2str(h)],'D. normale')
    plot(x,fun(x),'ob','MarkerSize',5,'MarkerFaceColor','b');
    
    %Verifico che l' integrale della Distribuzione Stimata sia sempre 1
    P=Parzen_Window_Simbolico(x,h,Vn);
    I=double(int(P,-10,10))
    
    pause
end

pause;
close all