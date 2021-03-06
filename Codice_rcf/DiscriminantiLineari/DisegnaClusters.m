% function DisegnaClusters(X,Y,F)
% Disegna i punto (X,Y) con il colore specificati dal vettore F
%
% parametri in:
%
% X: ascissa del punto
%
% Y: ordinata del punto
%
% F: colore del punto i-esimo
%                               1 - blue
%                               2 - verde
%                               3 - rosso
%                               4 - magenta

function DisegnaClusters(X,Y,F)
size=length(F);
for i=1:size
    switch F(i)
        case {1}
            color='b';
        case 2
            color='g';
        case { 3, -1, 0}
            color='r';
        case 4
            color='m';
    end
            
plot(X(i),Y(i),['o' color],'MarkerSize',3,'MarkerFaceColor',color); hold on;

end

%axis([0 200 0 200 0 200])