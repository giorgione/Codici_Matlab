
function F=ValoreAtteso_EM(T,T_i,K,A,Data_G,Data_B,method)

switch lower(method)
          case {'s','sym'}
                K=int(GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],Data_B),'x4',-inf,inf);
                
                F=GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G))...
                + (1/K)*int(log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],Data_B))...
                *GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],Data_B),'x4',-inf,inf);
          case {'n','numeric'}
                
                
                F=double(GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G))...
                + (1/K)*int(log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],Data_B))...
                *GaussianaMulti(diag(T_i(3:4)).^2,[T_i(1);T_i(2)],Data_B),'x4',-inf,inf))
            
                F=double(GaussianaMulti_Punti(diag(T(3:4).^2),[T(1);T(2)],double(Data_G))...
                + (1/K)*int(log(GaussianaMulti(diag(T(3:4)).^2,[T(1);T(2)],Data_B))...
                *A,'x4',-inf,inf))
          
        end
 

