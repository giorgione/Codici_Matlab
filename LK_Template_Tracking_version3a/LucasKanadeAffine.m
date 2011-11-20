function [p,I_roi,T_error]=LucasKanadeAffine(I,p,I_template,Options)
% This is an Affine Lucas Kanade template tracker, which performs 
% a template tracking step on a 2D image.
%
% [p,I_roi,T_error]=LucasKanadeAffine(I,p,I_template,Options)
% 
% inputs,
%   I : A 2d image of type double (movie frame)
%   p : 6 parameters, which describe affine transformation
%       (Backwards) Affine Transformation Matrix is used in 
%       Lucas Kanade Tracking with 6 parameters
%       M = [ 1+p(1) p(3)   p(5); 
%             p(2)   1+p(4) p(6); 
%             0      0      1];
%   I_template : An image of the template
%   Options : A struct with Options
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
defaultoptions=struct('TranslationIterations',6,'AffineIterations',6,'TolP',1e-5,'FineSigma',1.5,'RoughSigma',3,'SigmaIterations',2);
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

% Make all x,y indices
[x,y]=ndgrid(0:size(I_template,1)-1,0:size(I_template,2)-1);
% Calculate center of the template image
TemplateCenter=size(I_template)/2;
% Make center of the template image coordinates 0,0
x=x-TemplateCenter(1); y=y-TemplateCenter(2);

% Gradient of the image
[G_x1,G_y1]=image_derivatives(I,Options.RoughSigma);
[G_x2,G_y2]=image_derivatives(I,Options.FineSigma);

% Loop
for i=1:Options.TranslationIterations+Options.AffineIterations;
    % The affine matrix for template rotation and translation
    W_xp = [ 1+p(1) p(3) p(5); p(2) 1+p(4) p(6); 0 0 1];
    
    % 1: Warp I with W(x;p) to compute I(W(x;p))
    I_warped = affine_transform_2d_double(I,x,y,W_xp);

    % 2: Compute the error image T(x) - I(W(x;p))
    I_error= I_template - I_warped;

    % Break if outside image
    if((p(5)>(size(I,1))-1)||(p(6)>(size(I,2)-1))||(p(5)<0)||(p(6)<0)), break; end
    
    % 3: Warp the gradient gradI with W(x;p)
    % First itterations use a large scale derivative for more noise
    % robustness
    if(i<=Options.SigmaIterations)
        G_x_warped= affine_transform_2d_double(G_x1,x,y,W_xp);
        G_y_warped= affine_transform_2d_double(G_y1,x,y,W_xp);
    else
        G_x_warped= affine_transform_2d_double(G_x2,x,y,W_xp);
        G_y_warped= affine_transform_2d_double(G_y2,x,y,W_xp);
    end
    
    % First itterations only do translation updates for more robustness
    % and after that Affine.
    if(i>Options.TranslationIterations) 
        % Affine parameter optimalization
        
        % 4: Evaluate the Jacobian dW/ dp at (x;p)
        WP_Jacobian_x=[x(:) zeros(size(x(:))) y(:) zeros(size(x(:))) ones(size(x(:))) zeros(size(x(:)))];
        WP_Jacobian_y=[zeros(size(x(:))) x(:) zeros(size(x(:))) y(:) zeros(size(x(:))) ones(size(x(:)))];

        % 5: Compute the steepst descent image gradI * dW /dp
        I_steepest=zeros(numel(x),6);
        for j=1:numel(x),
            WP_Jacobian=[WP_Jacobian_x(j,:); WP_Jacobian_y(j,:)];
            Gradient=[G_x_warped(j) G_y_warped(j)];
            I_steepest(j,1:6)=Gradient*WP_Jacobian;
        end

        % 6: Compute the Hessian matrix using equation 11
        H=zeros(6,6);
        for j=1:numel(x), H=H+ I_steepest(j,:)'*I_steepest(j,:); end

        % 7: Computer sum_x [gradI*dW/dp]^T (T(X)-I(W(x;p))]
        sum_xy=zeros(6,1);
        for j=1:numel(x), sum_xy=sum_xy+I_steepest(j,:)'*I_error(j); end

        % 8: Computer delta_p using Equation 10
        delta_p=H\sum_xy;

        % 9: Update the parameters p <- p + delta_p
        p = p + delta_p;
    else
        % Translation parameter optimalization

        % 4: Evaluate the Jacobian dW/ dp at (x;p)
        % 5: Compute the steepst descent image gradI * dW /dp
        I_steepest(:,1)=G_x_warped(:);
        I_steepest(:,2)=G_y_warped(:);
        % 6: Compute the Hessian matrix using equation 11
        H=zeros(2,2);
        for j=1:numel(x), H=H+I_steepest(j,:)'*I_steepest(j,:); end
        % 7: Computer sum_x [gradI*dW/dp]^T (T(X)-I(W(x;p))]
        sum_xy=zeros(2,1);
        for j=1:numel(x), sum_xy=sum_xy+I_steepest(j,:)'*I_error(j); end
        % 8: Computer delta_p using Equation 10
        delta_p=H\sum_xy;
        % 9: Update the parameters p <- p + delta_p
        p(5:6) = p(5:6) + delta_p;
    end

    % Break if position is already good enough
    if((norm(delta_p,2)<Options.TolP)&&(i>Options.TranslationIterations)), break, end
end

I_roi=I_warped;
T_error=sum(I_error(:).^2)/numel(I_error);

function [Ix,Iy]=image_derivatives(I,sigma)
% Make derivatives kernels
[x,y]=ndgrid(floor(-3*sigma):ceil(3*sigma),floor(-3*sigma):ceil(3*sigma));
DGaussx=-(x./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
DGaussy=-(y./(2*pi*sigma^4)).*exp(-(x.^2+y.^2)/(2*sigma^2));
% Filter the images to get the derivatives
Ix = imfilter(I,DGaussx,'conv');
Iy = imfilter(I,DGaussy,'conv');

function Iout=affine_transform_2d_double(Iin,x,y,M)
% Affine transformation function (Rotation, Translation, Resize)
% This function transforms a volume with a 3x3 transformation matrix

%M=inv(M);

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
    Iout_one=intensity_xyz0.*perc0+intensity_xyz1.*perc1+intensity_xyz2.*perc2+intensity_xyz3.*perc3;
    Iout(:,:,i)=reshape(Iout_one, [size(x,1) size(x,2)]);
end

