function Time=EstraiSegnali(content,TypeSignal)

J=0;

N=length(content);
for I= 1:N
	  
	Time(I).name=content(I).name;
    
    
    %Numero segnali al tempo t
    nrSigns=length(content(I).signs);
    
    ns=1;
    for i=1:nrSigns
        sign  = content(I).signs(i);

        if strcmpi(sign.signTypes,TypeSignal) 
            Time(I).signs(ns)=sign;
            ns=ns+1;
             
        end   
    end 

end 

