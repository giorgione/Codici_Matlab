%function y=rettangolo(x,h)
% Valuta la finestra rettangolare di Ampiezza h nei punti x
function y=rettangolo(x,h)
n=length(x);
y=zeros(1,n);
for i=1:n
    if(abs(x(i))<=h/2)
        y(i)=1;
    end
end
