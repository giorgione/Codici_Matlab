%Validation Gate Test:
%
% Disegno: 
%       - il Validation Gate
%       - Le Misure Generate al tempo k
%       - I Dati che passano il Test di Validita
zx = sym('zx','real');
zy = sym('zy','real');
Z=[zx;zy];

Fun=(Z-Z_predic).'*inv(S)*(Z-Z_predic)-g_sigma;
ezplot(Fun,[-10^4 10^4]); hold on

%Disegno La misura predetta
plot(Z_predic(1),Z_predic(2),'or','MarkerFaceColor','r')
pause 

%Disegno Le Misure 
plot(NOISE(1,:),NOISE(2,:),'ob','MarkerFaceColor','b') 
pause
%Dati che passano il Test di Validita
if isempty(y)==0
    plot(y(1,:),y(2,:),'og','MarkerFaceColor','g') 
    pause
end


 
