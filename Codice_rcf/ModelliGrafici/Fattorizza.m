% Fattorizzazione del Probabilità congiunta del Modello Grafico specificato
% dalla Matrice di adiacenza Mxy
%             N
% P(x1,..xN)= II p( xk | Pa(xk))
%             k=1       
%
% ritorna:
%       la Probabilità congiunta fattorizzata come STRINGA

function Pjoint=Fattorizza(Mxy)
[M,N]=size(Mxy)

Pjoint='';
for i=1:N
    Parents=GetParents(i,Mxy);
    n=length(Parents);
    P='';
    for j=1:n
        if j==1 && n~=1
            P=[P sprintf('x%d,',Parents(j))];
        else
             P=[P sprintf('x%d',Parents(j))];
        end
    end
    if(i==1)
        if(n==0)
             Pjoint=[Pjoint 'P(x' num2str(i) ')' ];
        else
            Pjoint=[Pjoint 'P(x' num2str(i) ' | ' P ')' ];
        end
   else
         if(n==0)
            Pjoint=[Pjoint '*P(x' num2str(i) ')' ];
         else
            Pjoint=[Pjoint '*P(x' num2str(i) ' | ' P ')' ]; 
         end 
    end
end