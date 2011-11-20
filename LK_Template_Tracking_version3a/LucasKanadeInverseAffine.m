function [p,I_roi,T_error]=LucasKanadeInverseAffine(I,p,I_template,Wn,Options)
% This is an Robust Affine Lucas Kanade template tracker, which performs 
% a template tracking step on a 2D image. It is called inverse because
% template and image are switched in the equations for more time 
% efficient tracking.
%
% [p,I_roi,T_error]=LucasKanadeInverseAffine(I,p,I_template,Wn,Options)
% 
% inputs,
%   I : A 2d image of type double (movie frame)
%   p : 6 parameters, which describe affine transformation
%       (Backwards) Affine Transformation Matrix is used in 
%       Lucas Kanade Tracking with 6 parameters
%       M = [ 1+p(1) p(3)   p(5); 
%             p(2)   1+p(4) p(6); 
%             0      0      1];
%       Note : The center of the image is 0,0 not a corner
%   I_template : An image of the template
%   Wn : Robust weights same size as I_template, gives the reliability
%        of each pixel. see paper "Robust template tracking with drift correction"
%   Options : A struct with Options
%       .Padding: Defines the size of the padding of the template  
%                       input image. Those padding pixels at 
%                       the boundaries are used for derivatives not 
%                       for tracking (default 5)
%       .TranslationIterations : Number of translation itterations before
%                       performing Affine (default 6)
%       .AffineIterations: Number of Affine itterations (default 6)
%       .TolP: Tollerance on parameters allowed (default 1e-5)
%       .RoughSigma: Large sigma used in the first few itterations
%                       for noise robustness (default 3)
%       .FineSigma: Sigma used in the other itterations (default 1.5)
%       .SigmaIterations: Number of itterations with rough sigma 
%                       (default 2)
%
% Outputs,
%   p : The new affine parameters 
%   I_roi : The image ROI on the found template position
%   T_error : The squared error between template and ROI
%
% Literature used, S. Baker et Al. "Lucas-Kanade 20 Years  On: A  
%  Unifying Framework"
%
% Function is written by D.Kroon University of Twente (June 2009)

% Process inputs
defaultoptions=struct('TranslationIterations',6,'AffineIterations',6,'TolP',1e-5,'FineSigma',1.5,'RoughSigma',3,'SigmaIterations',2,'Padding',5);
if(~exist('Options','var')), 
    Options=defaultoptions; 
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
         if(~isfield(Options,tags{i})),  Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))), 
        warning('LucasKanadeAffine:unknownoption','unknown options found');
    end
end

% Parameters to column vector
p=p(:);

% 3: Evaluate the Gradient of the template
[G_x1,G_y1]=image_derivatives(I_template,Options.RoughSigma);
[G_x2,G_y2]=image_derivatives(I_template,Options.FineSigma);

% Remove the padding from the derivatives
b=Options.Padding;
G_x1=G_x1((1+b):end-b,(1+b):end-b);
G_y1=G_y1((1+b):end-b,(1+b):end-b);
G_x2=G_x2((1+b):end-b,(1+b):end-b);
G_y2=G_y2((1+b):end-b,(1+b):end-b);

% Remove padding weight function
Wn=Wn((1+b):end-b,(1+b):end-b);

% Make padded x,y template imagesindices
[x,y]=ndgrid(0:size(I_template,1)-1,0:size(I_template,2)-1);
TemplateCenter=size(I_template)/2;
x_padded=x-TemplateCenter(1); y_padded=y-TemplateCenter(2);

% Remove the padding from the template
I_template=I_template((1+b):end-b,(1+b):end-b);

% Make all x,y indices
[x,y]=ndgrid(0:size(I_template,1)-1,0:size(I_template,2)-1);
% Calculate center of the template image
TemplateCenter=size(I_template)/2;
% Make center of the template image coordinates 0,0
x=x-TemplateCenter(1); y=y-TemplateCenter(2);

% 4: Evaluate Gamma(x) (Same as Jacobian in normal LK_affine)
Gamma_x=[x(:) zeros(size(x(:))) y(:) zeros(size(x(:))) ones(size(x(:))) zeros(size(x(:)))];
Gamma_y=[zeros(size(x(:))) x(:) zeros(size(x(:))) y(:) zeros(size(x(:))) ones(size(x(:)))];

% 5: Compute the modified steepest descent images gradT * Gamma(x)
gG1=zeros(numel(x),6); gG2=zeros(numel(x),6);    
for j=1:numel(x),
    Gamma=[Gamma_x(j,:);Gamma_y(j,:)];
    Gradient=[G_x1(j) G_y1(j)];
    gG1(j,:)=Gradient*Gamma;
    Gradient=[G_x2(j) G_y2(j)];
    gG2(j,:)=Gradient*Gamma;
end
% 5, for translation only
gG1t=zeros(numel(x),2); gG2t=zeros(numel(x),2);    
for j=1:numel(x),
    gG1t(j,1)=G_x1(j); gG1t(j,2)=G_y1(j);
    gG2t(j,1)=G_x2(j); gG2t(j,2)=G_y2(j);
end

% 6: Compute the modified Hessian H* using equation 58
H_mod1=zeros(6,6); H_mod2=zeros(6,6);
H_mod1t=zeros(2,2); H_mod2t=zeros(2,2);
for j=1:numel(x),
    H_mod1=H_mod1+Wn(j)*gG1(j,:)'*gG1(j,:);
    H_mod2=H_mod2+Wn(j)*gG2(j,:)'*gG2(j,:);
    H_mod1t=H_mod1t+Wn(j)*gG1t(j,:)'*gG1t(j,:);
    H_mod2t=H_mod2t+Wn(j)*gG2t(j,:)'*gG2t(j,:);
end
% Compute the inverse hessians
H_mod1_inv=inv(H_mod1); 
H_mod2_inv=inv(H_mod2);
H_mod1t_inv=inv(H_mod1t); 
H_mod2t_inv=inv(H_mod2t);

% Loop
for i=1:Options.TranslationIterations+Options.AffineIterations;
    % W(x;p)
    W_xp = [ 1+p(1) p(3) p(5); p(2) 1+p(4) p(6); 0 0 1];

    % 1: Warp I with W(x;p) to compute I(W(x;p))
    I_warped = affine_transform_2d_double(I,x,y,W_xp,true);

    % 2: Compute the error image I(W(x;p))-T(x)
    I_error= I_warped - I_template ;
    % Break if outside image
    if((p(5)>(size(I,1))-1)||(p(6)>(size(I,2)-1))||(p(5)<0)||(p(6)<0)), break; end
 
    % First itterations only do translation updates for more robustness
    % and after that Affine.
    if(i>Options.TranslationIterations)
        % Affine parameter optimalization

        % 7: Computer sum_x [gradT*Gamma(x)]^T (I(W(x;p))-T(x)]
        sum_xy=zeros(6,1);
        if(i<=Options.SigmaIterations)
            for j=1:numel(x), sum_xy=sum_xy+gG1(j,:)'*I_error(j); end
        else
            for j=1:numel(x), sum_xy=sum_xy+gG2(j,:)'*I_error(j); end
        end

        % 8: Computer delta_p using Equation 61
        if(i<=Options.SigmaIterations)
            delta_p_mod=H_mod1_inv*sum_xy;
        else
            delta_p_mod=H_mod2_inv*sum_xy;
        end
        
    else
        % Translation parameter optimalization

        % 7: Computer sum_x [gradT*Gamma(x)]^T (I(W(x;p))-T(x)]
        sum_xy=zeros(2,1);
        if(i<=Options.SigmaIterations)
            for j=1:numel(x), sum_xy=sum_xy+Wn(j)*gG1t(j,:)'*I_error(j); end
        else
            for j=1:numel(x), sum_xy=sum_xy+Wn(j)*gG2t(j,:)'*I_error(j); end
        end

        % 8: Computer delta_p using Equation 61
        if(i<=Options.SigmaIterations)
            delta_p_mod=H_mod1t_inv*sum_xy;
        else
            delta_p_mod=H_mod2t_inv*sum_xy;
        end
        delta_p_mod=[0;0;0;0;delta_p_mod(1);delta_p_mod(2)];
    end
    
    % 9: Compute Gamma(p)^-1 and Update the parameters p <- p + delta_p
    Gamma_p_inv=[ (1+p(1)) p(3) 0 0 0 0; p(2) (1+p(4)) 0 0 0 0; 0 0 (1+p(1)) p(3) 0 0; 0 0 p(2) (1+p(4)) 0 0; 0 0 0 0 (1+p(1)) p(3); 0 0 0 0 p(2) (1+p(4))];
    delta_p=Gamma_p_inv*delta_p_mod;
    p = p - delta_p;
        
    % Break if position is already good enough
    if((norm(delta_p,2)<Options.TolP)&&(i>Options.TranslationIterations)), break, end
end

% Warp to give a roi back with padding
I_roi = affine_transform_2d_double(I,x_padded,y_padded,W_xp,true);
    
T_error=sum(I_error(:).^2)/numel(I_error);
p=p';


function [Ix,Iy]=image_derivatives(I,sigma)
% Make derivatives kernels
[x,y]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
DGaussx=-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
DGaussy=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
% Filter the images to get the derivatives
Ix = imfilter(I,DGaussx,'conv');
Iy = imfilter(I,DGaussy,'conv');

function Iout=affine_transform_2d_double(Iin,x,y,M,black)
% Affine transformation function (Rotation, Translation, Resize)
% This function transforms a volume with a 3x3 transformation matrix
%
% Iout=affine_transform_2d_double(Iin,Minv,black)
%
% inputs,
%   Iin: The input image
%   Minv: The (inverse) 3x3 transformation matrix
%   black: If true pixels from outside the image are set to zero
%           if false to the nearest old pixel.
% output,
%   Iout: The transformed image
%
% example,
%   % Read image
%   I=im2double(imread('lenag2.png'))
%   % Make a transformation matrix
%   M=make_transformation_matrix([2 3],2,[1.0 1.1]);
%   % Transform the image
%   Iout=affine_transform_2d_double(I,M,0)
%   % Show the image
%   figure, imshow(Iout);
%
% Function is written by D.Kroon University of Twente (February 2009)


% Calculate the Transformed coordinates
Tlocalx =  M(1,1) * x + M(1,2) *y + M(1,3) * 1;
Tlocaly =  M(2,1) * x + M(2,2) *y + M(2,3) * 1;

% All the neighborh pixels involved in linear interpolation.
xBas0=floor(Tlocalx);
yBas0=floor(Tlocaly);
xBas1=xBas0+1;
yBas1=yBas0+1;

% Linear interpolation constants (percentages)
xCom=Tlocalx-xBas0;
yCom=Tlocaly-yBas0;
perc0=(1-xCom).*(1-yCom);
perc1=(1-xCom).*yCom;
perc2=xCom.*(1-yCom);
perc3=xCom.*yCom;

% limit indexes to boundaries
check_xBas0=(xBas0<0)|(xBas0>(size(Iin,1)-1));
check_yBas0=(yBas0<0)|(yBas0>(size(Iin,2)-1));
xBas0(check_xBas0)=0;
yBas0(check_yBas0)=0;
check_xBas1=(xBas1<0)|(xBas1>(size(Iin,1)-1));
check_yBas1=(yBas1<0)|(yBas1>(size(Iin,2)-1));
xBas1(check_xBas1)=0;
yBas1(check_yBas1)=0;

Iout=zeros([size(x) size(Iin,3)]);
for i=1:size(Iin,3);
    Iin_one=Iin(:,:,i);
    % Get the intensities
    intensity_xyz0=Iin_one(1+xBas0+yBas0*size(Iin,1));
    intensity_xyz1=Iin_one(1+xBas0+yBas1*size(Iin,1));
    intensity_xyz2=Iin_one(1+xBas1+yBas0*size(Iin,1));
    intensity_xyz3=Iin_one(1+xBas1+yBas1*size(Iin,1));
    % Make pixels before outside Ibuffer black
    if(black>0)
        intensity_xyz0(check_xBas0|check_yBas0)=0;
        intensity_xyz1(check_xBas0|check_yBas1)=0;
        intensity_xyz2(check_xBas1|check_yBas0)=0;
        intensity_xyz3(check_xBas1|check_yBas1)=0;
    end
    Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;
    Iout(:,:,i)=reshape(Iout_one, [size(x,1) size(x,2)]);
end
