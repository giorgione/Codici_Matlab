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
#define getMatrixAddr(r, c, rsize) r + (c * rsize)
#define getTensorAddr(i, j, k, rsize, csize) i + (j * rsize) + (k * rsize * csize)
#define NUMINTERACTION 2
#define min(x, y)   (x < y)? x:y
#define max(x, y)   (x > y)? x:y
#define EPS         1e-20

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs == 4);
    else{
		printf("Invalid number of inputs %d\n",nrhs);
		return;
	}
    const mxArray *sample = prhs[0];
    const mxArray *tsample = prhs[1];
    const mxArray *prevZ = prhs[2];
    const mxArray *params = prhs[3];
    const mxArray *field_array_ptr;
    double *beta, *tbeta, *zbeta, *tmatrix;
    int nTargets, nSamples;
    
    msize sbeta, szbeta, stmatrix, stemp;
    
    field_array_ptr = mxGetField(prevZ, 0, "nTargets");
    if(field_array_ptr){
        nTargets = (int)mxGetScalar(field_array_ptr);
        DBGPRINT("nTargets = %d\n", nTargets);
    }
    else
    {
        printf("ERROR : No field named nTargets!!\n");
        return;
    }
    field_array_ptr = mxGetField(prevZ, 0, "nSamples");
    if(field_array_ptr){
        nSamples = (int)mxGetScalar(field_array_ptr);
        DBGPRINT("nSamples = %d\n", nSamples);
    }
    else
    {
        printf("ERROR : No field named nSamples!!\n");
        return;
    }
    plhs[0] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
    double *pbeta1 = (double*)mxGetData(plhs[0]);
    double *pbeta2 = (double*)mxGetData(plhs[1]);
    if (nTargets == 0) {
        for(int i = 0; i < nSamples;i++) {
            *(pbeta1 + i) = 1;
            *(pbeta2 + i) = 1;
        }        
        return;
    }
    
    field_array_ptr = mxGetField(sample, 0, "beta");
    if(field_array_ptr){
        beta = (double*)mxGetData(field_array_ptr);
        sbeta.nrows = mxGetM(field_array_ptr);
        sbeta.ncols = mxGetN(field_array_ptr);
        DBGPRINT("beta(1,1) = %d (%d, %d)\n", (*beta), sbeta.nrows, sbeta.ncols);
        if (sbeta.ncols / sbeta.nrows != NUMINTERACTION)
        {
            printf("Invalid beta matrix!!!");
            return;
        }
    }
    else
    {
        printf("ERROR : No field named beta!!\n");
        return;
    }
    field_array_ptr = mxGetField(tsample, 0, "beta");
    if(field_array_ptr){
        tbeta = (double*)mxGetData(field_array_ptr);
        stemp.nrows = mxGetM(field_array_ptr);
        stemp.ncols = mxGetN(field_array_ptr);
        DBGPRINT("beta(1,1) = %d (%d, %d)\n", (*tbeta), stemp.nrows, stemp.ncols);
        if (sbeta.ncols / sbeta.nrows != NUMINTERACTION)
        {
            printf("Invalid beta matrix!!!");
            return;
        }
    }
    else
    {
        printf("ERROR : No field named beta!!\n");
        return;
    }
    field_array_ptr = mxGetField(prevZ, 0, "beta");
    if(field_array_ptr){
        zbeta = (double*)mxGetData(field_array_ptr);
        szbeta.nrows = mxGetM(field_array_ptr);
        szbeta.ncols = mxGetN(field_array_ptr);
        DBGPRINT("zbeta(1,1) = %d (%d, %d)\n", (*zbeta), szbeta.nrows, szbeta.ncols);
        if (szbeta.ncols / szbeta.nrows != NUMINTERACTION * nSamples)
        {
            printf("Invalid zbeta matrix!!!");
            return;
        }
    }
    else
    {
        printf("ERROR : No field named zbeta!!\n");
        return;
    }
    
    field_array_ptr = mxGetField(params, 0, "betatrns");
    if(field_array_ptr){
        tmatrix = (double*)mxGetData(field_array_ptr);
        stmatrix.nrows = mxGetM(field_array_ptr);
        stmatrix.ncols = mxGetN(field_array_ptr);
        DBGPRINT("tmatrix(1,1) = %d (%d, %d)\n", (*tmatrix), stmatrix.nrows, stmatrix.ncols);
    }
    else
    {
        printf("ERROR : No field named tmatrix!!\n");
        return;
    }

    double *tzbeta;
    double temparray[NUMINTERACTION], tempres[NUMINTERACTION];
    int idInt1, idInt2;

    stemp.nrows = NUMINTERACTION;
    stemp.ncols = 1;
    
    for(int i = 0; i < nSamples;i++) {
        // start address of z.beta(:,:,:,i)
        *(pbeta1 + i) = 1;
        *(pbeta2 + i) = 1;
        
        tzbeta = zbeta + (szbeta.nrows * szbeta.nrows * NUMINTERACTION * i);
        for(int j = 0; j < nTargets; j++) {
            for(int k = j+1; k < nTargets; k++) {
                temparray[0] = *(tzbeta + getTensorAddr(j, k, 0, nTargets, nTargets));
                temparray[1] = *(tzbeta + getTensorAddr(j, k, 1, nTargets, nTargets));
//                 temparray[2] = *(tzbeta + getTensorAddr(j, k, 2, nTargets, nTargets));
                
                idInt1 = -1; idInt2 = -1;
                for(int l = 0; l < NUMINTERACTION; l++) {
                    if (*(beta + getTensorAddr(j, k, l, sbeta.nrows, sbeta.nrows)) == 1) {
#if DEBUG_ON
                        if (idInt1 > -1) {
                            DBGPRINT("ERROR MORE THAN TWO INTERACTIONS are set at the same time\n");
                        }                        
#endif
                        idInt1 = l;
                    }                    
                    if (*(tbeta + getTensorAddr(j, k, l, sbeta.nrows, sbeta.nrows)) == 1) {
#if DEBUG_ON
                        if (idInt2 > -1) {
                            DBGPRINT("ERROR MORE THAN TWO INTERACTIONS are set at the same time\n");
                        }                        
#endif
                        idInt2 = l;
                    }                    
                }
                MatrixMultiplication(tmatrix, temparray, stmatrix, stemp, tempres);
                *(pbeta1 + i) *= tempres[idInt1];
                *(pbeta2 + i) *= tempres[idInt2];
            }
        }
    } 
    
    return; 
}