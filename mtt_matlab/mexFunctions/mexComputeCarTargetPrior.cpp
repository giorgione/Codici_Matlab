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
    
    double *per = (double*)mxGetData(prhs[0]);
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
    
    msize sper, spSamples, sAper, sPrec, sTemp1, sTemp2;
    
    int nperv;
    double *Aper;
    double mh, sh, sv, normh, normv, noHprob;
    
    field_array_ptr = mxGetField(params, 0, "nperv");
    if(field_array_ptr){
        nperv = (int)mxGetScalar(field_array_ptr);
        DBGPRINT("nperv = %d\n", nperv);
    }
    else
    {
        printf("ERROR : No field named nperv!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "Aper");
    
    if(field_array_ptr){
        Aper = (double*)mxGetData(field_array_ptr);
        sAper.nrows = mxGetM(field_array_ptr);
        sAper.ncols = mxGetN(field_array_ptr);
        DBGPRINT("Aper(1,1) = %f\n", (*Aper));
    }
    else
    {
        printf("ERROR : No field named Aper!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "mh");
    if(field_array_ptr){
        mh = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("mh = %f\n", mh);
    }
    else
    {
        printf("ERROR : No field named mh!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "sh");
    if(field_array_ptr){
        sh = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("sh = %f\n", sh);
    }
    else
    {
        printf("ERROR : No field named sh!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "sv");
    if(field_array_ptr){
        sv = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("sv = %f\n", sv);
    }
    else
    {
        printf("ERROR : No field named sv!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "normh");
    if(field_array_ptr){
        normh = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("normh = %f\n", normh);
    }
    else
    {
        printf("ERROR : No field named normh!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "normv");
    if(field_array_ptr){
        normv = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("normv = %f\n", normv);
    }
    else
    {
        printf("ERROR : No field named normv!!\n");
        return;
    }
    field_array_ptr = mxGetField(params, 0, "noHprob");
    if(field_array_ptr){
        noHprob = (double)mxGetScalar(field_array_ptr);
        DBGPRINT("noHprob = %f\n", noHprob);
    }
    else
    {
        printf("ERROR : No field named noHprob!!\n");
        return;
    }
        
    sper.nrows = mxGetM(prhs[0]);    sper.ncols = mxGetN(prhs[0]);    
    spSamples.nrows = mxGetM(prhs[2]);    spSamples.ncols = mxGetN(prhs[2]);
	sPrec.nrows = nperv;    sPrec.ncols = nperv;
    
    plhs[0] = mxCreateDoubleMatrix(1, nSamples, mxREAL);
    prior = (double*)mxGetData(plhs[0]);
    
    if (*(per + nperv) == 1) {
//         temp = (normh * normv) * exp(sqr(*(per + 4) - mh) / (-2 * sqr(sh)) + (sqr(*(per + 2)) + sqr(*(per + 3)))/(-2*sqr(sv)));
        temp = (normh) * exp(sqr(*(per + 4) - mh) / (-2 * sqr(sh)));
    }
    else  {
        temp = noHprob;
    }
    for(int i = 0; i < nSamples; i++) {
        *(prior + i) = temp + 1e-20;
        DBGPRINT("%f, ", temp);
    }
    DBGPRINT("\n");
    
//     if (tcnt == 1) {
//     }
//     else 
    if(tcnt >= 1){
        sTemp1.nrows = 1; sTemp1.ncols = nperv;
        sTemp2.nrows = nperv; sTemp2.ncols = 1;
        
        for(int i = 0; i < spSamples.ncols; i++) {
            // copmute prior
            prec_ptr = mxGetCell(prhs[3], i);
            prec = (double*)mxGetData(prec_ptr);
            for(int j = 0; j < nperv; j++) {
                tempper[j] = *(per + j) - *(pSamples + j + (i * spSamples.nrows));
//                 tempdiff[j] = *(per + j) - *(pSamples + j + (i * spSamples.nrows));
            }
//             MatrixMultiplication(Aper, tempdiff, sAper, sTemp2, tempper);
//             DBGPRINT("1: %f, %f, %f, %f, %f\n", tempdiff[0], tempdiff[1], tempdiff[2], tempdiff[3], tempdiff[4]);
            MatrixMultiplication(tempper, prec, sTemp1, sPrec, tempmult);
            DBGPRINT("2: %f, %f, %f, %f, %f\n", tempmult[0], tempmult[1], tempmult[2], tempmult[3], tempmult[4]);
            MatrixMultiplication(tempmult, tempper, sTemp1, sTemp2, &temp);
            DBGPRINT("3: temp : %f\n", temp);
            *(prior + i) *= *(w+i) * *(normFull + i) * exp(-0.5 * temp);
            DBGPRINT("4: %f, %f\n", *(pSamples + 5 + (i * spSamples.nrows)), *(per + 5));
            
            if(1 == *(per + 5)) {
                *(prior + i) *= (1e-3 + *(pSamples + 5 + (i * spSamples.nrows)));
            }
            else {
                *(prior + i) *= (1 + 1e-3- *(pSamples + 5 + (i * spSamples.nrows)));
            }
            
            *(prior + i) += 1e-20;
        }    

    }
    else {
    }
}