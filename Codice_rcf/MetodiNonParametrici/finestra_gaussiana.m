%function y=rettangolo(x,h)
% Valuta la finestra Gaussiana di Ampiezza h nei punti x
function y=finestra_gaussiana(x)
f=@(u) exp(- (u.^2)/2)/sqrt(2*pi);
y=f(x);