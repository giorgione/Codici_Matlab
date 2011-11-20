% Carica le immagini
% 
% 
close all;
path='C:\Documents and Settings\Giorgio\Desktop\UniV\faces94\faces94\';
StartDir=dir(path);
numfile=size(StartDir);
Data=[];
for i=3:numfile-1
    %Apri le directory Male -Female
    path1=[path StartDir(i).name];
    DirDir=dir(path1);
    for j=3:8
        %Apri 5 Persone per Directory
        path2=[path1 '\' DirDir(j).name];
        DirImage=dir(path2);
        %Apri 5 imagini per ogni pesona
        for k=3:8
            path3=[path2 '\' DirImage(k).name];
            I=imread(path3,'JPEG');
            figure;
            
            image(I);
            title(DirImage(k).name)
        end
    end
    
       
    %Data=[Data I(:)];
end

[r c]=size(Data);
Medie=double(mean(Data,2))*ones(1,c);
B=double(Data)-Medie;
clear Data;

%Calcola le autofacce
C=B*B';
[V D]=eig(C);