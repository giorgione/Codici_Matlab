#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "matrixwg.h"
#include "matrix.h"

#define DEBUG_ON 0
#define WEIGHT_FACTOR   5

#if DEBUG_ON
#define DBGPRINT printf
#else
#define DBGPRINT
#endif

#define sqr(x)  (x)*(x)

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs == 8);
    else{
		printf("Invalid number of inputs %d\n",nrhs);
		return;
	}
    
    double *gfeat = (double*)mxGetData(prhs[0]);
    int tcnt = (int)mxGetScalar(prhs[1]);
    double *pSamples = (double*)mxGetData(prhs[2]);
    double *w = (double*)mxGetData(prhs[4]);
    double *normFull = (double*)mxGetData(prhs[5]);
    int nSamples = (int)mxGetScalar(prhs[7]);
    
    const mxArray *params = prhs[6];
    const mxArray *field_array_ptr;
    
    double *prec;
    const mxArray *prec_ptr;
    
    double *prior;
    double tempdiff[10], tempmult[10], tempper[10];
    double temp;
    
    msize sgfeat, spSamples, sPrec, sTemp1, sTemp2;
    
    int ngfeat;
    field_array_ptr = mxGetField(params, 0, "ngfeat");
    if(field_array_ptr){
        ngfeat = (int)mxGetScalar(field_array_ptr);
        DBGPRINT("ngfeat = %d\n", ngfeat);
    }
    else
    {
        printf("ERROR : No field named ngfeat!!\n");
        return;
    }
    sgfeat.nrows = mxGetM(prhs[0]);    sgfeat.ncols = mxGetN(prhs[0]);    
    spSamples.nrows = mxGetM(prhs[2]);    spSamples.ncols = mxGetN(prhs[2]);
	sPrec.nrows = ngfeat-1;    sPrec.ncols = ngfeat-1;
    
    plhs[0] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
    prior = (double*)mxGetData(plhs[0]);
    
    if (tcnt >= 1) {
        sTemp1.nrows = 1; sTemp1.ncols = ngfeat-1;
        sTemp2.nrows = ngfeat-1; sTemp2.ncols = 1;
        
        for(int i = 0; i < nSamples; i++) {
            // copmute prior
            prec_ptr = mxGetCell(prhs[3], i);
            prec = (double*)mxGetData(prec_ptr);
            for(int j = 0; j < ngfeat-1; j++) {
                tempdiff[j] = *(gfeat + j) - *(pSamples + j + (i * spSamples.nrows));
            }
            MatrixMultiplication(tempdiff, prec, sTemp1, sPrec, tempmult);
            DBGPRINT("2: %f, %f\n", tempmult[0], tempmult[1]);
            MatrixMultiplication(tempmult, tempdiff, sTemp1, sTemp2, &temp);
            DBGPRINT("3: temp : %f\n", temp);
            *(prior + i) = *(w+i) * *(normFull + i) * exp(-0.5 * temp);
            DBGPRINT("4: prior : %f\n", *(prior + i));
            if(1 == *(gfeat + 2)) {
                *(prior + i) *= *(pSamples + 2 + (i * spSamples.nrows));
            }
            else {
                *(prior + i) *= (1 - *(pSamples + 2 + (i * spSamples.nrows)));
            }
            DBGPRINT("5: prior : %f\n", *(prior + i));
            
            *(prior + i) += 1e-20;
        }    
    }
    else {
        // prior on the features location.
//         assert("!!!");
        if(1 == *(gfeat + 2)) {
            temp = 0.7;
        }
        else {
            temp = 0.3;
        }
        
        for(int i = 0; i < nSamples; i++) {
            *(prior + i) = temp;
        }
    }
}