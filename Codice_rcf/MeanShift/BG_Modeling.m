% Fasi di Tracking
imgdir='../Documents/DataSets/AVSS_Easy/'; 

% Lettura immagini in cui Tracciare il Background
imfiles = dir([imgdir '*.png']);
if length(imfiles) == 0
    imfiles = dir([imgdir '*.jpg']);
    if length(imfiles) == 0
        imfiles = dir([imgdir '*.bmp']);
    end
end
 
NFRAMES = length(imfiles);
 
firstframe = 2;
lastframe = length(imfiles);
 


Irgb = imread([imgdir imfiles(1).name]);
Irgb=Irgb(1:7:end,1:7:end,:);
[h w c] = size(Irgb);
Distance=@(S,M,X) exp(- (X-M).'*S*(X-M));


 

BG=zeros(h,w);
% %%%%%%%%%%
for I = 2:length(imfiles)
    disp(['Processing ' num2str(I) 'th frame of ' imgdir]);
    Irgb = imread([imgdir imfiles(I).name]);
    Irgb=Irgb(1:7:end,1:7:end,:);
    %% Analisi delle singole immagini: per ogni pixel (i,j)
    %
    %1) Cerco le regioni Candidate nelle quali il Pixel cade
    %
    %2) Cerco la Rergione R di colore + simile al pixel
    %
    %3) Classifico il Pixel
    
    for i=1:h
        for j=1:w
            PC=double([ Irgb(i,j,1); Irgb(i,j,2); Irgb(i,j,3)])./255;
            PX=[i;j];
            for k=1:AcceptedRegion
                Ds(k)=Distance(Region(k).SCov,Region(k).SMedia,PX);
            end
            
            %Cerco le regioni Candidate
            [Ri]=find(Ds>0);
            
            CandidateRegion=length(Ri);
            
            for k=Ri
                Dc(k)=Distance(Region(k).CCov,Region(k).CMedia,PC);
            end
             
            %Regione Cui e' associato il Pixel
            [Val Rj]=max(Dc);
            if(Val > 0.9)
                BG(i,j)=0;
            else
                BG(i,j)=1;
            end
        end
    end
    figure(8)
    subplot(2,1,1)
    imshow(BG)
    subplot(2,1,2)
    imshow(Irgb)
    pause;
end