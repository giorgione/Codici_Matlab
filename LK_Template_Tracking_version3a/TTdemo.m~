% Start with clean Matlab Workspace
clear all; close all; clc

% Compile the fast Lucas Kanade c-code tracking
%mex LucasKanadeInverseAffine.c -v

% Load a Traffic Movie, frames in stack with 3th dimension time
%
% (Only differences between frames are stored, for 3x smaller .mat filesize.
% Integrate to get the approximated original movie back)
load('TTdemo_packed_movie'); Vmovie=uint8(cumsum(single(Vmovie),3)+128);

% Get the first movie frame
I = double(Vmovie(:,:,1))*(1/255);

% Show the movie frame
figure, imshow(I,[]);hold on

% Make a struct to store templates
TemplateData=struct;

% Select the coordinates of 2 templates
% rect=getrect; TempPos1=round([rect(2) rect(2)+rect(4);rect(1) rect(1)+rect(3);]);
% rect=getrect; TempPos2=round([rect(2) rect(2)+rect(4);rect(1) rect(1)+rect(3);]);
TempPos1=[086,111;100,134];
TempPos2=[054,091;224,256];
TempPos3=[093,103;285,312];
TempPos4=[154,171;250,266];

% Pad the select templates with extra boundary pixels. These boundary
% pixels are not used for the actual template tracking. But to get
% more reliable image derivatives.
b=5;
padding=[-b,b;-b,b];
TempPos1=TempPos1+padding;
TempPos2=TempPos2+padding;
TempPos3=TempPos3+padding;
TempPos4=TempPos4+padding;

%disegno i template
        hh=TempPos1(1,2)-TempPos1(1,1);
        ww=TempPos1(2,2)-TempPos1(2,1);
        rectangle('Position',[TempPos1(2,1),TempPos1(1,1),ww,hh] ,'FaceColor','r');
        
        hh=TempPos2(1,2)-TempPos2(1,1);
        ww=TempPos2(2,2)-TempPos2(2,1);
        rectangle('Position',[TempPos2(2,1),TempPos2(1,1),ww,hh], 'FaceColor','r');
        
        hh=TempPos3(1,2)-TempPos3(1,1);
        ww=TempPos3(2,2)-TempPos3(2,1);
        rectangle('Position',[TempPos3(2,1),TempPos3(1,1),ww,hh], 'FaceColor','r');
                
        hh=TempPos4(1,2)-TempPos4(1,1);
        ww=TempPos4(2,2)-TempPos4(2,1);
        rectangle('Position',[TempPos4(2,1),TempPos4(1,1),ww,hh], 'FaceColor','r');

% -Set initial parameters of the templates
% -Set padded template image.
% 
% Backwards affine Transformation Matrix is used in Lucas Kanade Tracking
% with 6 parameters
% M    = [ 1+p(1) p(3)   p(5); 
%          p(2)   1+p(4) p(6); 
%          0      0      1];
%
center=[TempPos1(1,1)+TempPos1(1,2)-1 TempPos1(2,1)+TempPos1(2,2)-1]/2;
TemplateData(1).p=[0 0 0 0 center(1) center(2)];
TemplateData(1).image=I(TempPos1(1,1):TempPos1(1,2),TempPos1(2,1):TempPos1(2,2));
% This weight function is used in the LK-Hessian and multiplied with the 
% error between in image and template. And is used to exclude unreliable pixels form
% the template tracking.
TemplateData(1).weight=im2double(imread('weight1.png'));


center=[TempPos2(1,1)+TempPos2(1,2)-1 TempPos2(2,1)+TempPos2(2,2)-1]/2;
TemplateData(2).p=[0 0 0 0 center(1) center(2)];
TemplateData(2).image=I(TempPos2(1,1):TempPos2(1,2),TempPos2(2,1):TempPos2(2,2));
TemplateData(2).weight=im2double(imread('weight2.png'));
 
center=[TempPos3(1,1)+TempPos3(1,2)-1 TempPos3(2,1)+TempPos3(2,2)-1]/2;
TemplateData(3).p=[0 0 0 0 center(1) center(2)];
TemplateData(3).image=I(TempPos3(1,1):TempPos3(1,2),TempPos3(2,1):TempPos3(2,2));
TemplateData(3).weight=im2double(imread('weight3.png'));

center=[TempPos4(1,1)+TempPos4(1,2)-1 TempPos4(2,1)+TempPos4(2,2)-1]/2;
TemplateData(4).p=[0 0 0 0 center(1) center(2)];
TemplateData(4).image=I(TempPos4(1,1):TempPos4(1,2),TempPos4(2,1):TempPos4(2,2));
TemplateData(4).weight=im2double(imread('weight4.png'));

% LK Tracking Options (other options default)
Options.TranslationIterations=30;
Options.AffineIterations=0;
Options.RoughSigma=3;
Options.FineSigma=1.5;

% Make a colormap
cmap=hot(256);

% Matrix to store squared pixel error between template and ROI in
% movieframe after template tracking.
T_error=zeros(size(Vmovie,3), length(TemplateData));

% Loop through the movie frames
for i=1:size(Vmovie,3)
    % Get the a movie frame
    I = double(Vmovie(:,:,i))*(1/255);
    
    % Do the tracking for all templates, using Lucas Kanade Inverse Affine
    for t=1:length(TemplateData)
        [TemplateData(t).p,ROIimage,T_error(i,t)]=LucasKanadeInverseAffine(I,TemplateData(t).p,TemplateData(t).image,TemplateData(t).weight,Options);
            
% % Weights update, see paper "Robust template tracking with drift correction"
%
% % Constant used for the weight function [0..1], with a lower value
% % the weight function will be less depended on the current error between
% % template and image (more average of itterations) then with a higher
% % value.
% alpha = 0.1;
%         if(i>1)
%             TemplateData(t).E(:,:,i)=(1-alpha)*TemplateData(t).E(:,:,i-1)+alpha*abs(ROIimage-TemplateData(t).image);
%         else
%             TemplateData(t).E(:,:,i)=abs(ROIimage-TemplateData(t).image);
%         end
%         if(i>5)
%             TemplateData(t).weight=double(TemplateData(t).E(:,:,i)<=median(TemplateData(t).E(:,:,4:i),3)*1.4826); 
%         end
    end
    
    % Show the movie frame
    if(i==1), figure, handle_imshow=imshow(I); hold on
        
        
    else
       for t=1:length(TemplateData)
           delete(h(t))
       end
        set(handle_imshow,'Cdata',I); 
    end

    % Show the location of the templates in the movie frame
    for t=1:length(TemplateData)
        h(t)=plot(TemplateData(t).p(6),TemplateData(t).p(5),'go','MarkerFaceColor',cmap(round(t*255/length(TemplateData))+1,:)); 
        drawnow
    end

end

% Show the squared pixel errors between template and ROI during the movie frames
figure,
subplot(2,2,1),plot(T_error(:,1)); title('Pixel^2 error template 1')
subplot(2,2,2),plot(T_error(:,2)); title('Pixel^2 error template 2')
subplot(2,2,3),plot(T_error(:,3)); title('Pixel^2 error template 3')
subplot(2,2,4),plot(T_error(:,4)); title('Pixel^2 error template 3')








