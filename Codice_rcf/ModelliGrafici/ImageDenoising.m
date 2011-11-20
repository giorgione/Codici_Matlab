%Esempio di Image Denoising con ICM su Energia Totale di Markov Random Field
clc;clear;close all;
figure(1);
I=imread('ImageDenoise.jpg');

I=I(1:8:end,1:8:end);
imshow(I)
[m,n]=size(I);
%Trasformo i pixel in [-1 1]

for i=1:m
    for j=1:n
        if I(i,j)==0  %nero --> -1
            I1(i,j)=-1;
        else
            I1(i,j)=1; % bianco -->1
        end
    end
end
  
%Genero L'immagine Rumorosa con 10% di pixel alterati(inverto il colore)
Irs=imnoise(I,'salt & pepper',0.4);
figure(2)
imshow(Irs)
for i=1:m
    for j=1:n
        if Irs(i,j)==0  %nero --> -1
            Irumorosa(i,j)=-1;
        else
            Irumorosa(i,j)=1; % bianco -->1
        end
    end
end

%Genero il MRF
% Variabili Osservate
Y=Irumorosa;

%Variabili Latenti: Inizializzo X=Y
X=Y;
h=0;
b=1;
nn=2.1;

%Genero il MRF
% Variabili Osservate
Y=Irumorosa;

%Variabili Latenti: Inizializzo X=Y
X=Y;        %..%
h=1;     %..%
b=1;        %..%
nn=2.1;  %..%

%Energia Totale Iniziale
Eo=EnergiaTotaleRMF(X,Y,h,b,nn);
disp(['Iter 0: ' num2str(Eo)]);
Iter=0;
minE=Eo;
while Iter< 10
    
    for i=1:m
        for j=1:n

            % Calcolo il Contributo Locale di X(i,j)      
            X(i,j)=-1;
            E_=EnergiaTotaleRMF(X,Y,h,b,nn);

            % 
            X(i,j)=1;
            E=EnergiaTotaleRMF(X,Y,h,b,nn);

            % Conservo il valore di X(i,j) che produce l'energia Minimima
            if E_<E            
                 X(i,j)=-1;
                 minE=E_;
                 disp(['Pixel( ' num2str(i) ',' num2str(j) '): ' num2str(E_)])
            else            
                 X(i,j)=1;
                 minE=E;
                 disp(['Pixel( ' num2str(i) ',' num2str(j) '): ' num2str(E)])
            end

            %    disp(['Pixel ' num2str(i*n+j) 'E-: ' num2str(E_) ' E+: ' num2str(E)])

           
        end
    end
    
    Iter=Iter+1;
    disp(['Iter ' num2str( Iter) ': ' num2str(minE) ])
     % Converti X in im magine visualizzabile
     for i1=1:m
         for j1=1:n
             if X(i1,j1)==-1  %nero --> -1
                    XX(i1,j1)=0;
             else
                   XX(i1,j1)=255; % bianco -->1
            end
         end
    end

   figure(3)
   imshow(XX)
   pause;
end