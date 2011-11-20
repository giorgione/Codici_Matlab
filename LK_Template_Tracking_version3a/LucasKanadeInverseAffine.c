#include "mex.h"
#include "math.h"
#include "string.h"
#define PI 3.1415926535897931
#define E 2.7182818284590455
/*
 * This is an Robust Affine Lucas Kanade template tracker, which performs
 * a template tracking step on a 2D image. It is called inverse because
 * template and image are switched in the equations for more time
 * efficient tracking.
 *
 * [p,I_roi,T_error]=LucasKanadeInverseAffine(I,p,I_template,Options)
 *
 * inputs,
 *   I : A 2d image of type double (movie frame)
 *   p : 6 parameters, which describe affine transformation
 *       (Backwards) Affine Transformation Matrix is used in
 *       Lucas Kanade Tracking with 6 parameters
 *       M = [ 1+p(1) p(3)   p(5);
 *             p(2)   1+p(4) p(6);
 *             0      0      1];
 *       Note : The center of the image is 0,0 not a corner
 *   I_template : An image of the template
 *   Wn : Robust weights same size as I_template, gives the reliability
 *       of each pixel. see paper "Robust template tracking with drift correction"
 *   Options : A struct with Options
 *       .Padding: Defines the size of the padding of the template
 *                       input image. Those padding pixels at
 *                       the boundaries are used for derivatives not
 *                       for tracking (default 5)
 *       .TranslationIterations : Number of translation itterations before
 *                       performing Affine (default 6)
 *       .AffineIterations: Number of Affine itterations (default 6)
 *       .TolP: Tollerance on parameters allowed (default 1e-5)
 *       .RoughSigma: Large sigma used in the first few itterations
 *                       for noise robustness (default 3)
 *       .FineSigma: Sigma used in the other itterations (default 1.5)
 *       .SigmaIterations: Number of itterations with rough sigma
 *                       (default 2)
 *
 * Outputs,
 *   p : The new affine parameters
 *   I_roi : The image ROI on the found template position
 *   T_error : The squared error between template and ROI
 *
 * Literature used, S. Baker et Al. "Lucas-Kanade 20 Years  On: A
 *  Unifying Framework"
 *
 * Function is written by D.Kroon University of Twente (June 2009)
 */

struct options {
    int Padding;
    int TranslationIterations;
    int AffineIterations;
    double TolP;
    double RoughSigma;
    double FineSigma;
    int SigmaIterations;
};
void setdefaultoptions(struct options* t) {
    t->Padding=5;
    t->TranslationIterations=6;
    t->AffineIterations=6;
    t->TolP=1e-5;
    t->RoughSigma=3.0;
    t->FineSigma=1.5;
    t->SigmaIterations=2;
}

__inline double pow2(double val) { return val*val; }
__inline double pow3(double val) { return val*val*val; }
__inline double pow4(double val) { return pow2(val)*pow2(val); }

void matrix_derivatives(double *I, double sigma, double *Ix, double *Iy, int *sizeI) {
/* This function calculates a gaussian based x and y derivative of a 2D image */
    int k_min, k_max, k_size, kx, ky, ix, iy, i, j, index, index2;
    double val, valx, valy, *DGaussx, *DGaussy;
    /* Calculate kernel variables */
    k_min=(int)floor(-3*sigma);
    k_max=(int)ceil(3*sigma);
    k_size=k_max-k_min+1;
    /* Reserve memory to store kernels */
    DGaussx= (double*)malloc( k_size*k_size*sizeof(double) );
    DGaussy= (double*)malloc( k_size*k_size*sizeof(double) );
    /* Make Gaussian derivative kernels */
    for(j=0; j<k_size; j++) {
        ky = -(j+k_min);
        for(i=0; i<k_size; i++) {
            kx = -(i+k_min);
            val=-(kx/(2.0*PI*pow4(sigma)))*exp(-(pow2(kx)+pow2(ky))/(2.0*pow2(sigma)));
            DGaussx[i+j*k_size]=val; DGaussy[j+i*k_size]=val;
        }
    }
    /* Imfilter the whole image with the gaussian kernels */
    for(iy=0; iy<sizeI[1]; iy++) {
        for(ix=0; ix<sizeI[0]; ix++) {
            valx=0; valy=0;
            for(j=0; j<k_size; j++) {
                ky = iy+(j+k_min);
                for(i=0; i<k_size; i++) {
                    kx = ix+(i+k_min);
                    if((kx>=0)&&(ky>=0)&&(kx<sizeI[0])&&(ky<sizeI[1])) {
                        index = kx+ky*sizeI[0]; index2 = i+j*k_size;
                        valx+=I[index]*DGaussx[index2];
                        valy+=I[index]*DGaussy[index2];
                    }
                }
            }
            index=ix + iy*sizeI[0];
            Ix[index]=valx; Iy[index]=valy;
        }
    }
    /* Remove from kernels from memory */
    free(DGaussx); free(DGaussy);
}

void matrix_inverse(double *Min, double *Mout, int actualsize) {
    /* This function calculates the inverse of a square matrix
     *
     * matrix_inverse(double *Min, double *Mout, int actualsize)
     *
     * Min : Pointer to Input square Double Matrix
     * Mout : Pointer to Output (empty) memory space with size of Min
     * actualsize : The number of rows/columns
     *
     * Notes:
     *  - the matrix must be invertible
     *  - there's no pivoting of rows or columns, hence,
     *        accuracy might not be adequate for your needs.
     *
     * Code is rewritten from c++ template code Mike Dinolfo
     */
    /* Loop variables */
    int i, j, k;
    /* Sum variables */
    double sum, x;
    
    /*  Copy the input matrix to output matrix */
    for(i=0; i<actualsize*actualsize; i++) { Mout[i]=Min[i]; }
    
    /* Add small value to diagonal if diagonal is zero */
    for(i=0; i<actualsize; i++) {
        j=i*actualsize+i;
        if((Mout[j]<1e-12)&&(Mout[j]>-1e-12)){ Mout[j]=1e-12; }
    }
    
    /* Matrix size must be larger than one */
    if (actualsize <= 1) return;
    
    for (i=1; i < actualsize; i++) {
        Mout[i] /= Mout[0]; /* normalize row 0 */
    }
    
    for (i=1; i < actualsize; i++)  {
        for (j=i; j < actualsize; j++)  { /* do a column of L */
            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += Mout[j*actualsize+k] * Mout[k*actualsize+i];
            }
            Mout[j*actualsize+i] -= sum;
        }
        if (i == actualsize-1) continue;
        for (j=i+1; j < actualsize; j++)  {  /* do a row of U */
            sum = 0.0;
            for (k = 0; k < i; k++) {
                sum += Mout[i*actualsize+k]*Mout[k*actualsize+j];
            }
            Mout[i*actualsize+j] = (Mout[i*actualsize+j]-sum) / Mout[i*actualsize+i];
        }
    }
    for ( i = 0; i < actualsize; i++ )  /* invert L */ {
        for ( j = i; j < actualsize; j++ )  {
            x = 1.0;
            if ( i != j ) {
                x = 0.0;
                for ( k = i; k < j; k++ ) {
                    x -= Mout[j*actualsize+k]*Mout[k*actualsize+i];
                }
            }
            Mout[j*actualsize+i] = x / Mout[j*actualsize+j];
        }
    }
    for ( i = 0; i < actualsize; i++ ) /* invert U */ {
        for ( j = i; j < actualsize; j++ )  {
            if ( i == j ) continue;
            sum = 0.0;
            for ( k = i; k < j; k++ ) {
                sum += Mout[k*actualsize+j]*( (i==k) ? 1.0 : Mout[i*actualsize+k] );
            }
            Mout[i*actualsize+j] = -sum;
        }
    }
    for ( i = 0; i < actualsize; i++ ) /* final inversion */ {
        for ( j = 0; j < actualsize; j++ )  {
            sum = 0.0;
            for ( k = ((i>j)?i:j); k < actualsize; k++ ) {
                sum += ((j==k)?1.0:Mout[j*actualsize+k])*Mout[k*actualsize+i];
            }
            Mout[j*actualsize+i] = sum;
        }
    }
}

/* Get an pixel from an image */
double getintensity_mindex2(int x, int y, int sizx, int sizy, double *I) {
    return I[y*sizx+x];
}

__inline double interpolate_2d_linear_gray(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Linear interpolation variables */
    int xBas0, xBas1, yBas0, yBas1;
    double perc[4]={0, 0, 0, 0};
    double xCom, yCom, xComi, yComi;
    double color[4]={0, 0, 0, 0};
    
    /* Rounded location */
    double fTlocalx, fTlocaly;
    
    /* Determine the coordinates of the pixel(s) which will be come the current pixel */
    /* (using linear interpolation) */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    xBas1=xBas0+1; yBas1=yBas0+1;
    
    /* Linear interpolation constants (percentages) */
    xCom=Tlocalx-fTlocalx; yCom=Tlocaly-fTlocaly;
    xComi=(1-xCom); yComi=(1-yCom);
    perc[0]=xComi * yComi;
    perc[1]=xComi * yCom;
    perc[2]=xCom * yComi;
    perc[3]=xCom * yCom;
    
    if(xBas0<0) { xBas0=0; if(xBas1<0) { xBas1=0; }}
    if(yBas0<0) { yBas0=0; if(yBas1<0) { yBas1=0; }}
    if(xBas1>(Isize[0]-1)) { xBas1=Isize[0]-1; if(xBas0>(Isize[0]-1)) { xBas0=Isize[0]-1; }}
    if(yBas1>(Isize[1]-1)) { yBas1=Isize[1]-1; if(yBas0>(Isize[1]-1)) { yBas0=Isize[1]-1; }}
    
    color[0]=getintensity_mindex2(xBas0, yBas0, Isize[0], Isize[1], Iin);
    color[1]=getintensity_mindex2(xBas0, yBas1, Isize[0], Isize[1], Iin);
    color[2]=getintensity_mindex2(xBas1, yBas0, Isize[0], Isize[1], Iin);
    color[3]=getintensity_mindex2(xBas1, yBas1, Isize[0], Isize[1], Iin);
    return color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
}

__inline double interpolate_2d_cubic_gray(double Tlocalx, double Tlocaly, int *Isize, double *Iin) {
    /* Floor of coordinate */
    double fTlocalx, fTlocaly;
    /* Zero neighbor */
    int xBas0, yBas0;
    /* The location in between the pixels 0..1 */
    double tx, ty;
    /* Neighbor loccations */
    int xn[4], yn[4];

    /* The vectors */
    double vector_tx[4],vector_ty[4];
    double vector_qx[4]; /*={0.25,0.25,0.25,0.25}; */
    double vector_qy[4]; /*={0.25,0.25,0.25,0.25}; */
    double vector_b[4];
    /* Interpolated Intensity; */
	double Ipixel=0;
    /* Temporary value boundary */
    int b;
    /* Loop variable */
    int i;
    
    /* Determine of the zero neighbor */
    fTlocalx = floor(Tlocalx); fTlocaly = floor(Tlocaly);
    xBas0=(int) fTlocalx; yBas0=(int) fTlocaly;
    
    /* Determine the location in between the pixels 0..1 */
    tx=Tlocalx-fTlocalx; ty=Tlocaly-fTlocaly;

    /* Determine the t vectors */
    vector_tx[0]= 0.5; vector_tx[1]= 0.5*tx; vector_tx[2]= 0.5*pow2(tx); vector_tx[3]= 0.5*pow3(tx);
    vector_ty[0]= 0.5; vector_ty[1]= 0.5*ty; vector_ty[2]= 0.5*pow2(ty); vector_ty[3]= 0.5*pow3(ty);
    
    /* t vector multiplied with 4x4 bicubic kernel gives the to q vectors */
    vector_qx[0]= -1.0*vector_tx[1]+2.0*vector_tx[2]-1.0*vector_tx[3];
    vector_qx[1]= 2.0*vector_tx[0]-5.0*vector_tx[2]+3.0*vector_tx[3];
    vector_qx[2]= 1.0*vector_tx[1]+4.0*vector_tx[2]-3.0*vector_tx[3];
    vector_qx[3]= -1.0*vector_tx[2]+1.0*vector_tx[3];
    vector_qy[0]= -1.0*vector_ty[1]+2.0*vector_ty[2]-1.0*vector_ty[3];
    vector_qy[1]= 2.0*vector_ty[0]-5.0*vector_ty[2]+3.0*vector_ty[3];
    vector_qy[2]= 1.0*vector_ty[1]+4.0*vector_ty[2]-3.0*vector_ty[3];
    vector_qy[3]= -1.0*vector_ty[2]+1.0*vector_ty[3];
                
    /* Determine 1D neighbour coordinates */
    xn[0]=xBas0-1; xn[1]=xBas0; xn[2]=xBas0+1; xn[3]=xBas0+2;
    yn[0]=yBas0-1; yn[1]=yBas0; yn[2]=yBas0+1; yn[3]=yBas0+2;
    
    /* Clamp to image boundary if outside image */
    if(xn[0]<0) { xn[0]=0;if(xn[1]<0) { xn[1]=0;if(xn[2]<0) { xn[2]=0; if(xn[3]<0) { xn[3]=0; }}}}
    if(yn[0]<0) { yn[0]=0;if(yn[1]<0) { yn[1]=0;if(yn[2]<0) { yn[2]=0; if(yn[3]<0) { yn[3]=0; }}}}
    b=Isize[0]-1;
    if(xn[3]>b) { xn[3]=b;if(xn[2]>b) { xn[2]=b;if(xn[1]>b) { xn[1]=b; if(xn[0]>b) { xn[0]=b; }}}}
    b=Isize[1]-1; 
    if(yn[3]>b) { yn[3]=b;if(yn[2]>b) { yn[2]=b;if(yn[1]>b) { yn[1]=b; if(yn[0]>b) { yn[0]=b; }}}}
    
    /* First do interpolation in the x direction followed by interpolation in the y direction */
    for(i=0; i<4; i++)
    {
        vector_b[i] =vector_qx[0]*getintensity_mindex2(xn[0], yn[i], Isize[0], Isize[1], Iin);
        vector_b[i]+=vector_qx[1]*getintensity_mindex2(xn[1], yn[i], Isize[0], Isize[1], Iin);
        vector_b[i]+=vector_qx[2]*getintensity_mindex2(xn[2], yn[i], Isize[0], Isize[1], Iin);
        vector_b[i]+=vector_qx[3]*getintensity_mindex2(xn[3], yn[i], Isize[0], Isize[1], Iin);
        Ipixel+= vector_qy[i]*vector_b[i];
    }        
    return Ipixel;
}

void transformimage_gray(double *Iin, int *sizeI, int *sizeT, double *A, double *Iout) {
    int x, y;
    double TemplateCenter[2];
    
    /* Location of pixel which will be come the current pixel */
    double Tlocalx, Tlocaly;
    
    /* X,Y,Z coordinates of current pixel */
    double xd, yd;
    
    /* Parts of location calculation */
    double compa0, compa1, compb0, compb1;
    
    /* Variables to store 1D index */
    int index=0;
    
    /* Calculate center of template; */
    TemplateCenter[0]=((double)sizeT[0])/2.0;
    TemplateCenter[1]=((double)sizeT[1])/2.0;
    
    /* Translation part; */
    compb0= A[2];
    compb1= A[5];
    
    /* Loop through all image pixel coordinates */
    for (y=0; y<sizeT[1]; y++) {
        yd=(double)y-TemplateCenter[1];
        compa0 = A[1] *yd + compb0;
        compa1 = A[4] *yd + compb1;
        
        for (x=0; x<sizeT[0]; x++) {
            xd=(double)x-TemplateCenter[0];
            Tlocalx =  A[0] * xd + compa0;
            Tlocaly =  A[3] * xd + compa1;
            
            /* Set the current pixel value */
            Iout[index]=interpolate_2d_cubic_gray(Tlocalx, Tlocaly, sizeI, Iin); index++;
        }
    }
}

void matrixAddvectorXvectorWeighted(double *V, double wn, int lengthv, double *M) {
    int i, j; 
    for(j=0;j<lengthv;j++) { for(i=0;i<lengthv;i++) { M[i+j*lengthv]+=wn*V[j]*V[i]; } }
}

void matrix_min_matrix(double *M1, double *M2, int npixels, double *M) {
    int i; for(i=0;i<npixels;i++) { M[i]=M1[i]-M2[i]; }
}


void LucasKanadeInverseAffine(double *I, int *sizeI, double *Wn_padded, double *p_in, double *I_template_padded, int *sizeT, struct options* Options, double *p_out, double *I_roi, double *T_error) {
    double *G_x1, *G_y1, *G_x2, *G_y2, *I_template, *Wn;
    /* Loop variables */
    int i, j, k;
    /* Coordinates */
    double x, y;
    /* Temp variable */
    int b;
    /* 1D Index variables */
    int Index1, Index2, Index3;
    /* Tsize without padding */
    int sizeTc[2];
    /* Template Center */
    double TemplateCenter[2];
    /* Gamma (is like the jacobian in Lucas Kanade Affine) */
    double Gamma_x[6], Gamma_y[6];
    /* Inverse parameter Gamma */
    double *Gamma_p_inv;
    /* Steepest decent images */
    double *gG1, *gG2, *gG1t, *gG2t;
    /* Hessian */
    double *H_mod1, *H_mod2, *H_mod1t, *H_mod2t;
    double *H_mod1_inv, *H_mod2_inv, *H_mod1t_inv, *H_mod2t_inv;
    /* Number of template pixels */
    int nTPixels;
    /* Warp/ pixel transformation matrix */
    double W_xp[9];
    /* Warped image */
    double *I_warped;
    /* Error between warped image and template */
    double *I_error;
    /* . */
    double sum_xy[6];
    /* Modified and unmodified Update of affine parameters */
    double delta_p_mod[6], delta_p[6];
    /* Norm storage */
    double norm_s;
    
    
    G_x1 = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    G_y1 = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    G_x2 = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    G_y2 = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    Wn = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    
    
    /* 3: Evaluate the Gradient of the template */
    matrix_derivatives(I_template_padded, Options->RoughSigma, G_x1, G_y1, sizeT);
    matrix_derivatives(I_template_padded, Options->FineSigma, G_x2, G_y2, sizeT);
    
    /* Remove the padding from the derivatives */
    b=Options->Padding;
    sizeTc[0]=sizeT[0]-b*2; sizeTc[1]=sizeT[1]-b*2;
    I_template = (double*)malloc( sizeT[0]*sizeT[1]*sizeof(double) );
    if((sizeTc[0]<0)||(sizeTc[1]<0)){ mexErrMsgTxt("2xPadding must be smaller than template size"); }
    for(j=0;j<sizeTc[1];j++) {
        for(i=0;i<sizeTc[0];i++) {
            Index1=i+j*sizeTc[0]; Index2=(i+b)+(j+b)*sizeT[0];
            G_x1[Index1]=G_x1[Index2];
            G_y1[Index1]=G_y1[Index2];
            G_x2[Index1]=G_x2[Index2];
            G_y2[Index1]=G_y2[Index2];
            Wn[Index1]=Wn_padded[Index2];
            I_template[Index1]=I_template_padded[Index2];
        }
    }
    
    TemplateCenter[0]=((double)sizeTc[0])/2.0;
    TemplateCenter[1]=((double)sizeTc[1])/2.0;
    nTPixels=sizeTc[0]*sizeTc[1];
    
    /* memory for steepest decent images */
    gG1  = (double*)malloc( 6*nTPixels*sizeof(double) );
    gG2  = (double*)malloc( 6*nTPixels*sizeof(double) );
    gG1t  = (double*)malloc( 2*nTPixels*sizeof(double) );
    gG2t  = (double*)malloc( 2*nTPixels*sizeof(double) );
    
    /* memory for hessian and inverted hessian */
    H_mod1 = (double*)malloc( 6*6*sizeof(double) );
    H_mod2 = (double*)malloc( 6*6*sizeof(double) );
    H_mod1t = (double*)malloc( 2*2*sizeof(double) );
    H_mod2t = (double*)malloc( 2*2*sizeof(double) );
    H_mod1_inv = (double*)malloc( 6*6*sizeof(double) );
    H_mod2_inv = (double*)malloc( 6*6*sizeof(double) );
    H_mod1t_inv = (double*)malloc( 2*2*sizeof(double) );
    H_mod2t_inv = (double*)malloc( 2*2*sizeof(double) );
    for (i=0; i<36; i++) {  H_mod1[i]=0; H_mod2[i]=0; H_mod1_inv[i]=0; H_mod2_inv[i]=0;  }
    for (i=0; i<4; i++) {  H_mod1t[i]=0; H_mod2t[i]=0; H_mod1t_inv[i]=0; H_mod2t_inv[i]=0;  }
    
    /* memory for warped image and error image */
    I_warped = (double*)malloc( nTPixels*sizeof(double) );
    I_error = (double*)malloc( nTPixels*sizeof(double) );
    
    /* memory for inverse parameter gamma, and initialize to zero */
    Gamma_p_inv = (double*)malloc(36*sizeof(double) );
    for (i=0; i<36; i++) {  Gamma_p_inv[i]=0; }
    
    for(j=0;j<sizeTc[1];j++) {
        for(i=0;i<sizeTc[0];i++) {
            Index1=i+j*sizeTc[0];
            x=((double)i)-TemplateCenter[0]; y=((double)j)-TemplateCenter[1];
            
            /* 4: Evaluate Gamma(x) (Same as Jacobian in normal LK_affine) */
            Gamma_x[0]=x; Gamma_x[1]=0; Gamma_x[2]=y; Gamma_x[3]=0; Gamma_x[4]=1; Gamma_x[5]=0;
            Gamma_y[0]=0; Gamma_y[1]=x; Gamma_y[2]=0; Gamma_y[3]=y; Gamma_y[4]=0; Gamma_y[5]=1;
            
            /* 5: Compute the modified steepest descent images gradT * Gamma(x) */
            Index2=Index1*6;
            for(k=0; k<6; k++) {
                gG1[Index2+k]=Gamma_x[k]*G_x1[Index1]+Gamma_y[k]*G_y1[Index1];
                gG2[Index2+k]=Gamma_x[k]*G_x2[Index1]+Gamma_y[k]*G_y2[Index1];
            }
            
            /* 5, for translation only */
            Index3=Index1*2;
            gG1t[Index3+0]=G_x1[Index1];
            gG1t[Index3+1]=G_y1[Index1];
            gG2t[Index3+0]=G_x2[Index1];
            gG2t[Index3+1]=G_y2[Index1];
            
            /* 6: Compute the modified Hessian H* using equation 58 */
            matrixAddvectorXvectorWeighted(&gG1[Index2], Wn[Index1], 6, H_mod1);
            matrixAddvectorXvectorWeighted(&gG2[Index2], Wn[Index1], 6, H_mod2);
            matrixAddvectorXvectorWeighted(&gG1t[Index3], Wn[Index1], 2, H_mod1t);
            matrixAddvectorXvectorWeighted(&gG2t[Index3], Wn[Index1], 2, H_mod2t);
            
            /* Compute the inverse hessians */
            matrix_inverse(H_mod1, H_mod1_inv, 6);
            matrix_inverse(H_mod2, H_mod2_inv, 6);
            matrix_inverse(H_mod1t, H_mod1t_inv, 2);
            matrix_inverse(H_mod2t, H_mod2t_inv, 2);
        }
    }
    
    
    /* Copy parameters to output/working parameters */
    memcpy(p_out, p_in, 6*sizeof(double));
    
    /* Lucas Kanade Main Loop */
    for (i=0; i<(Options->TranslationIterations+Options->AffineIterations); i++) {
        /* W(x;p) */
        W_xp[0]=1+p_out[0];  W_xp[1]=p_out[2];   W_xp[2]=p_out[4];
        W_xp[3]=p_out[1];    W_xp[4]=1+p_out[3]; W_xp[5]=p_out[5];
        W_xp[6]=0;           W_xp[7]=0;          W_xp[8]=1;
        
        /* 1: Warp I with W(x;p) to compute I(W(x;p)) */
        transformimage_gray(I, sizeI, sizeTc, W_xp, I_warped);
        
        /* 2: Compute the error image I(W(x;p))-T(x) */
        matrix_min_matrix(I_warped, I_template, nTPixels, I_error);
        
        /* Break if outside image */
        if((p_out[4]>(sizeI[0]-1))||(p_out[5]>(sizeI[1]-1))||(p_out[4]<0)||(p_out[5]<0)) { break; }
        
        /* First itterations only do translation updates for more robustness */
        /* and after that Affine. */
        
        if(i>=Options->TranslationIterations) {
            /* Affine parameter optimalization */
            
            /* 7: Computer sum_x [gradT*Gamma(x)]^T (I(W(x;p))-T(x)] */
            for (j=0; j<6; j++) {sum_xy[j]=0;}
            if(i<Options->SigmaIterations) {
                for (j=0; j<nTPixels; j++) {
                    for (k=0; k<6; k++) { sum_xy[k]+=Wn[j]*gG1[k+j*6]*I_error[j]; }
                }
            }
            else {
                for (j=0; j<nTPixels; j++) {
                    for (k=0; k<6; k++) { sum_xy[k]+=Wn[j]*gG2[k+j*6]*I_error[j]; }
                }
            }
            
            /* 8: Computer delta_p using Equation 61 */
            for (j=0; j<6; j++) {delta_p_mod[j]=0;}
            if(i<Options->SigmaIterations) {
                for (j=0; j<6; j++) {
                    for (k=0; k<6; k++) { delta_p_mod[j]+=H_mod1_inv[k*6+j]*sum_xy[k]; }
                }
            }
            else {
                for (j=0; j<6; j++) {
                    for (k=0; k<6; k++) { delta_p_mod[j]+=H_mod2_inv[k*6+j]*sum_xy[k]; }
                }
            }
        }
        else {
            /* Translation parameter optimalization */
            
            /* 7: Computer sum_x [gradT*Gamma(x)]^T (I(W(x;p))-T(x)] */
            sum_xy[0]=0; sum_xy[1]=0;
            if(i<Options->SigmaIterations) {
                for (j=0; j<nTPixels; j++) {
                    sum_xy[0]+=gG1t[0+j*2]*I_error[j];
                    sum_xy[1]+=gG1t[1+j*2]*I_error[j];
                }
            }
            else {
                for (j=0; j<nTPixels; j++) {
                    sum_xy[0]+=gG2t[0+j*2]*I_error[j];
                    sum_xy[1]+=gG2t[1+j*2]*I_error[j];
                }
            }
            
            /* 8: Computer delta_p using Equation 61 */
            for (j=0; j<6; j++) {delta_p_mod[j]=0;}
            if(i<Options->SigmaIterations) {
                for (j=0; j<2; j++) {
                    for (k=0; k<2; k++) { delta_p_mod[j+4]+=H_mod1t_inv[k*2+j]*sum_xy[k]; }
                }
            }
            else {
                for (j=0; j<2; j++) {
                    for (k=0; k<2; k++) { delta_p_mod[j+4]+=H_mod2t_inv[k*2+j]*sum_xy[k]; }
                }
            }
        }
        
        /* 9: Compute Gamma(p)^-1 and Update the parameters p <- p + delta_p */
        Gamma_p_inv[0]=1+p_out[0];   Gamma_p_inv[6]=p_out[2];      Gamma_p_inv[1]=p_out[1];
        Gamma_p_inv[7]=1+p_out[3];   Gamma_p_inv[14]=1+p_out[0];   Gamma_p_inv[20]=p_out[2];
        Gamma_p_inv[15]=p_out[1];    Gamma_p_inv[21]=1+p_out[3];   Gamma_p_inv[28]=1+p_out[0];
        Gamma_p_inv[34]=p_out[2];    Gamma_p_inv[29]=p_out[1];     Gamma_p_inv[35]=1+p_out[3];

  
        for (j=0; j<6; j++) {delta_p[j]=0;}
        for (j=0; j<6; j++) {
            for (k=0; k<6; k++) { delta_p[j]+=Gamma_p_inv[k*6+j]*delta_p_mod[k]; }
        }

        for (j=0; j<6; j++) { p_out[j]-=delta_p[j];}
        
        /* Break if position is already good enough */
        norm_s=0; for (j=0; j<6; j++) {norm_s+=pow2(delta_p[j]);}
        if((norm_s<pow2(Options->TolP))&&(i>=Options->TranslationIterations)) { break; }
    }
   
    /* Warp to give a roi back with padding */
    transformimage_gray(I, sizeI, sizeT, W_xp, I_roi);
    
    norm_s=0; for(j=0; j<nTPixels; j++) {norm_s+=pow2(I_error[j]);}
    T_error[0]=norm_s/((double)nTPixels);

    /* Free memory */
    free(G_x1); free(G_y1); free(G_x2); free(G_y2);
    free(gG1); free(gG2); free(gG1t); free(gG2t);
    free(H_mod1); free(H_mod2); free(H_mod1t); free(H_mod2t);
    free(H_mod1_inv); free(H_mod2_inv); free(H_mod1t_inv); free(H_mod2t_inv);
    free(I_warped); free(I_error); free(Gamma_p_inv);
}




/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *I, *p_in, *I_template, *Wn;
    double *p_out, *I_roi, *T_error;
    double *OptionsField;
    const mwSize *dimsI, *dimsT, *dimsP;
    mwSize dimsE[2]={1, 1};
    int field_num;
    int sizeI[2];
    int sizeT[2];
    struct options Options;
    /* Check inputs */
    /* [p,I_roi,T_error]=LucasKanadeInverseAffine(I,p,I_template,Options) */
    
    /* Check number of inputs */
    if(nrhs<3) { mexErrMsgTxt("3 or 4 input variables required."); }
    
    /* Check properties of image I */
    if(mxGetNumberOfDimensions(prhs[0])!=2) { mexErrMsgTxt("Image must be 2D"); }
    if(!mxIsDouble(prhs[0])){ mexErrMsgTxt("Image must be double"); }
    dimsI = mxGetDimensions(prhs[0]);
    sizeI[0]=dimsI[0]; sizeI[1]=dimsI[1];
    
    /* Check properties of affine parameters */
    if(!mxIsDouble(prhs[1])){ mexErrMsgTxt("Parameters must be double"); }
    dimsP = mxGetDimensions(prhs[1]);
    if( mxGetNumberOfElements(prhs[1])!=6) {
        mexErrMsgTxt("Length of Parameters variable must be 6");
    }
    
    /* Check properties of template image */
    if(mxGetNumberOfDimensions(prhs[2])!=2) { mexErrMsgTxt("Template must be 2D"); }
    if(!mxIsDouble(prhs[2])){ mexErrMsgTxt("Template must be double"); }
    dimsT = mxGetDimensions(prhs[2]);
    sizeT[0]=dimsT[0]; sizeT[1]=dimsT[1];
    
    /* Check properties of Weight matrix */
    if(!mxIsDouble(prhs[3])){ mexErrMsgTxt("Weights must be double"); }
    if( mxGetNumberOfElements(prhs[2])!=mxGetNumberOfElements(prhs[3])) {
         mexErrMsgTxt("Weight matrix must be same size as image"); 
    }
    
    /* Check options */
    setdefaultoptions(&Options);
    if(nrhs==5) {
        if(!mxIsStruct(prhs[4])){ mexErrMsgTxt("Options must be structure"); }
        field_num = mxGetFieldNumber(prhs[4], "Padding");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.Padding=(int)OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "TranslationIterations");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.TranslationIterations=(int)OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "AffineIterations");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.AffineIterations=(int)OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "TolP");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.TolP=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "RoughSigma");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.RoughSigma=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "FineSigma");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.FineSigma=OptionsField[0];
        }
        field_num = mxGetFieldNumber(prhs[4], "SigmaIterations");
        if(field_num>=0) {
            OptionsField=mxGetPr(mxGetFieldByNumber(prhs[4], 0, field_num));
            Options.SigmaIterations=(int)OptionsField[0];
        }
    }
    
    /* Connect inputs */
    I = mxGetPr(prhs[0]);
    p_in = mxGetPr(prhs[1]);
    I_template = mxGetPr(prhs[2]);
    Wn = mxGetPr(prhs[3]);
    
    
    /* Connect outputs */
    plhs[0] = mxCreateNumericArray(2, dimsP, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(2, dimsT, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericArray(2, dimsE, mxDOUBLE_CLASS, mxREAL);
    p_out = mxGetPr(plhs[0]);
    I_roi = mxGetPr(plhs[1]);
    T_error = mxGetPr(plhs[2]);
    
    /* Perform the Lucaks Kanade tracking */
    LucasKanadeInverseAffine(I, sizeI, Wn, p_in, I_template, sizeT, &Options, p_out, I_roi, T_error);
}













