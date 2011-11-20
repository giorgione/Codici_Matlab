%Calcolo il Valore Atteso della VEROSIMIGLIANZA logaritmica dei DATI
%COMPLETI rispetto ai dati INCOMPLETI 
function F=ValoreAtteso_EM_ottimizzata(T,K,A,Data_G,Data_B,method)

switch lower(method)
          case {'s','sym'}
                 %Calcolo simbolicamente l' integrale con Dati Simbolici
                 F=double(GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G))...
                + (1/K)*int(log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],Data_B))...
                *A,'x4',-inf,inf))
          case {'r','reali'}
                %Calcolo simbolicamente l' integrale con Dati Numerici
                F=double(GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G))...
                + (1/K)*int(log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],Data_B))...
                *A,'x4',-inf,inf))
         case {'n','numeric'}
                %Calcolo Numericamente l' integrale
                T_i=A;
                I=@(x) FIntegranda(T,T_i,Data_B,x,'n');
                Area=quad(I,-10,10);
                
                F=double(GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G)))...
                + (1/K)*Area;
end
 
%Verosimiglianza Logaritmica
function I=FIntegranda(T,T_i,Data_B,x,method)
switch lower(method)
          case {'s','sym'}
   
I=log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],[x;double(Data_B(2))]))...
   *GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],[x;double(Data_B(2))]);
               
          case {'n','numeric'}
              [d,N]=size(x);
              I=log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],[x;double(Data_B(2))*ones(1,N)]))...
   .*GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],[x;double(Data_B(2))*ones(1,N)]);
end
 
            