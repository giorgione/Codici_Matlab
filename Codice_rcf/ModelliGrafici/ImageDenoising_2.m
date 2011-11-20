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

%Energia Totale Iniziale
Eo=EnergiaTotaleRMF(X,Y,h,b,nn);
disp(['Iter 0: ' num2str(Eo)])
Niter=20;
for I=1:Niter
    
    for i=1:m
        for j=1:n

            %Calcolo il Contributo Locale di X(i,j)
            a=X(i,j);
            E1=EnergiaLocale(i,j,X,Y,h,b,nn);

            %Sottraggo il contributo Locale Attuale All'EnergiaTotale
            Eo=Eo-E1;

            X(i,j)=-a;
            E2=EnergiaLocale(i,j,X,Y,h,b,nn);
            %Conservo il valore di X(i,j) che produce l'energia Minimima
            if(E1<E2)
                X(i,j)=a;
                Eo=Eo+E1;
                %disp(['Iter ' num2str(i*n+j) ': ' num2str(Eo)])
            else
                X(i,j)=-a; 
                Eo=Eo+E2;

            end
            XX=X;
            %Converti X in immagine visualizzabile
            for ii=1:m
                for jj=1:n
                    if X(ii,jj)==-1  %nero --> -1
                        XX(ii,jj)=0;
                    else
                        XX(ii,jj)=255; % bianco -->1
                    end
                end
            end
            figure(3)
            imshow(XX)
            disp(['Iter ' num2str(I) ' Pixel(' num2str(i) ',' num2str(j) '): ' num2str(Eo)])
            pause
        end
    end
    %disp(['Iter ' num2str(I) ': ' num2str(Eo)])
    %visualizzo
    
end
XX=X;
%Converti X in immagine visualizzabile
for i=1:m
    for j=1:n
        if X(i,j)==-1  %nero --> -1
            XX(i,j)=0;
        else
            XX(i,j)=255; % bianco -->1
        end
    end
end
figure(3)
imshow(XX)