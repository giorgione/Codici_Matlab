% fourier descriptor
close 3; close 4
%fisso la lunghezza del descrittore;
N=57;
Npunti=length(contorno(1).x(1:50:end));

c=contorno(1).x(1:50:end)+1i*contorno(1).y(1:50:end);
figure(3);hold on
plot(contorno(1).x,contorno(1).y,'g','LineWidth',2);
plot(contorno(1).x(1),contorno(1).y(1),'or','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(contorno(1).x(end),contorno(1).y(end),'ob','MarkerFaceColor','b','MarkerEdgeColor','b');
pause

C=fft(c,N);
crec=ifft(C,Npunti)
contornoR(1).x=real(crec);
contornoR(1).y=imag(crec);
figure(4);hold on
plot(contornoR(1).x,contornoR(1).y,'g','LineWidth',2);
plot(contornoR(1).x(1),contornoR(1).y(1),'or','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(contornoR(1).x(end),contornoR(1).y(end),'ob','MarkerFaceColor','b','MarkerEdgeColor','b');
