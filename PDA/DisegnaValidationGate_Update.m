%Validation Gate dopo la Fase di Update:
%
% Disegno: 
%       - il Validation Gate
%       - Le Misure Generate al tempo k
%       - I Dati che passano il Test di Validita
figure(1)
clf
ezplot(Fun,[-10^4 10^4]); hold on

%Disegno la misura predetta
plot(Z_predic(1),Z_predic(2),'or','MarkerFaceColor','r')
  
%Disegno Lo STATO STIMATO
plot( x_filter(1,t),x_filter(3,t),'ob','MarkerFaceColor','b')

%Misure che hanno passato il TEST
if isempty(y)==0
    plot(y(1,:),y(2,:),'og','MarkerFaceColor','g') 
end

figure(2)
%ezplot(Fun,[-10^4 10^4]); hold on

%Disegno la misura predetta
%plot(Z_predic(1),Z_predic(2),'or','MarkerFaceColor','r')
  
%Disegno Lo STATO STIMATO
plot( x_filter(1,t),x_filter(3,t),'or','MarkerFaceColor','r')

%Misure che hanno passato il TEST
%if isempty(y)==0
%    plot(y(1,:),y(2,:),'og','MarkerFaceColor','g') 
%end

display('Update: done')
pause 