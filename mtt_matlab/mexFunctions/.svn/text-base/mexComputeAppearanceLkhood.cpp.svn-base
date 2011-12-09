#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "matrixwg.h"
#include "matrix.h"

#define DEBUG_ON 0

#if DEBUG_ON
#define DBGPRINT printf
#else
#define DBGPRINT
#endif

#define sqr(x)  (x)*(x)
#define TensorRef(A, r, c, h, rs, cs) *(A + r + (c + h * cs) * rs)
#define PI 3.141592

#define Bhatt

double bhattacharyya_coeff(double *a, double *b, int len)
{
    double dist = 0;
    for(int i = 0; i < len; i++) {
        dist += sqrt((*(a+i)) * (*(b+i)));
    }
    return dist;
}

double euclidean_sim(double *a, double *b, int len)
{
    double dist = 0;
    double denoma = 0;
    double denomb = 0;
    
    for(int i = 0; i < len; i++) {
        dist += ((*(a+i)) * (*(b+i)));
        denoma += ((*(a+i)) * (*(a+i)));
        denomb += ((*(b+i)) * (*(b+i)));
    }
    return dist / sqrt(denoma * denomb);
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs == 5);
    else{
		printf("Invalid number of inputs %d, parameters : (Img, Hist, colStep, wSigma, nColChannel)\n",nrhs);
		return;
	}
    
    unsigned char *pImg = (unsigned char*)mxGetData(prhs[0]);
    double *pHistIn = (double*)mxGetData(prhs[1]);
    double *pHistOut;
    
    int colStep = (int)mxGetScalar(prhs[2]);
    double wSigma = (double)mxGetScalar(prhs[3]);
    int nColChannel = (int)mxGetScalar(prhs[4]);
    int histDim; msize imSize;

    int temp, histIdx, chSize;
    double imDims, weight, totalWeight;
       
    imSize.nrows = mxGetM(prhs[0]);    imSize.ncols = mxGetN(prhs[0]) / nColChannel;
    imDims = imSize.nrows * imSize.ncols;
    
    chSize = 256 >> colStep;
    
    histDim = (int)pow((double)chSize, nColChannel);
    plhs[0] = mxCreateDoubleMatrix(1, histDim, mxREAL);
    pHistOut = (double*)mxGetData(plhs[0]);
    
    /* histogram construction */
    totalWeight = 0;
    double t1, t2;
    for(int row = 0; row < imSize.nrows; row++) {
        for(int col = 0; col < imSize.ncols; col++) {
            histIdx = 0;
            for(int ch = 0; ch < nColChannel; ch++) {
                temp = (int)TensorRef(pImg, row, col, ch, imSize.nrows, imSize.ncols);
                histIdx += pow((double)chSize, ch) * (temp >> colStep);
            }
            // add spatial weight here, if necessary
//             weight = 1.0;
            weight = exp(-(pow(row / imSize.nrows - 0.5, 2) + pow(col / imSize.ncols - 0.5, 2)) / (2 * pow(wSigma, 2)) );
            (*(pHistOut + histIdx)) += weight;
            
            totalWeight += weight;
        }
    }
    // histogram normalization
    for(int i = 0; i < histDim; i++) {
        *(pHistOut + i) /= totalWeight;
    }
    
    /* Computing probability(lkhood) of match */
    msize histInSize;
    double* lkhood;
    
    histInSize.nrows = mxGetM(prhs[1]);    histInSize.ncols = mxGetN(prhs[1]);
    
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    lkhood = (double*)mxGetData(plhs[1]);
    *lkhood = 0;
    
    if (histDim != histInSize.nrows * histInSize.ncols) {
        return;
    }
#ifdef Bhatt
    *lkhood = bhattacharyya_coeff(pHistOut, pHistIn, histDim);
#else
    *lkhood = euclidean_sim(pHistOut, pHistIn, histDim);
#endif
    
    return;
}