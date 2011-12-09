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
Distance=@(S,X) exp(-(X.')*S*(X));
R=Irgb(:,:,1);
G=Irgb(:,:,2);
B=Irgb(:,:,3);

 
Npoints=w*h;
 
R=reshape(R,1,Npoints);
G=reshape(G,1,Npoints);
B=reshape(B,1,Npoints);


[x y]=meshgrid(1:h,1:w);

 

BG=zeros(h,w);
reshape(B,1,Npoints)
% %%%%%%%%%%
for I = 2:length(imfiles)
    disp(['Processing ' num2str(I) 'th frame of ' imgdir]);
    Irgb = imread([imgdir imfiles(I).name]);
    Irgb=Irgb(1:7:end,1:7:end,:);
    
    R=Irgb(:,:,1);
    G=Irgb(:,:,2);
    B=Irgb(:,:,3);

    R=reshape(R,1,Npoints);
    G=reshape(G,1,Npoints);
    B=reshape(B,1,Npoints);
    x=reshape(x,1,Npoints);
    y=reshape(y,1,Npoints);
    
    X=[x;y];
    
    for k=1:AcceptedRegion
        Media=repmat(Region(k).SMedia,1,Npoints);
        X=X-Media;
        D=Distance(Region(k).SCov,X);
        Ds=diag(D );
        
        [ij]=find(Ds>0.98);
            
         CandidatePixels=length(ij);
         if CandidatePixels > 0 
             X=[ R(ij); G(ij); B(ij)];

             Media=repmat(Region(k).CMedia,1,CandidatePixels);
             X=X-Media;
             D=Distance(Region(k).CCov,X);
             Dc=diag(D );

             ij1=find(Dc>0.90);

             BG(ij1(ij))=BG(ij1(ij))+1;
         end
         
    end
    
    
     
     
   
        
    figure(8)
    subplot(2,1,1)
    imshow(BG)
    subplot(2,1,2)
    imshow(Irgb)
    pause;
end