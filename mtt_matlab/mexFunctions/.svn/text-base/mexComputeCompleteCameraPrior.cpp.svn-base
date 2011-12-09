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
// void getPIP(double x, double y, double *cam1, double *cam2, double* res, int camlen)
// {
//     double temp1[2], temp2[2];
//     
//     temp1[1] = cam1[0] * cam1[1] / (y - cam1[3]);
//     temp1[0] = (x - cam1[2]) * temp1[1] / cam1[0];
//     
//     if (camlen == 8)
//     {
//         temp2[0] = cos(cam1[4]) * temp1[0] - sin(cam1[4]) * temp1[1];
//         temp2[1] = sin(cam1[4]) * temp1[0] + cos(cam1[4]) * temp1[1];
//         
//         temp1[0] = temp2[0] + cam1[6] - cam2[6]; // transition
//         temp1[1] = temp2[1] + cam1[7] - cam2[7];
//         
//         temp2[0] = cos(cam2[4]) * temp1[0] + sin(cam2[4]) * temp1[1];
//         temp2[1] = -sin(cam2[4]) * temp1[0] + cos(cam2[4]) * temp1[1];
//     }
//     else
//     {
//         double angle = cam1[4] - cam2[4];
//         
//         temp2[0] = cos(angle) * temp1[0] + -sin(angle) * temp1[1];
//         temp2[1] = sin(angle) * temp1[0] + cos(angle) * temp1[1];
//     }
//     
//     res[0] = cam2[0] / temp2[1] * temp2[0] + cam2[2];
//     res[1] = cam2[0] / temp2[1] * cam2[1] + cam2[3];
// }
// 
// // temp = [xx(:, 2),yy(:, 2)]' - getKLTPredictions([xx(:, 1), yy(:, 1)]', cam, sample.cam);
// int getKLTdiff(double *klt_x, double *klt_y, unsigned int nklt, double minhor, double *prevCam, double *cam, double *res, int camlen)
// {
//     double temp[2];
//     int cnt = 0;
//    
// //     getPIP(*(klt_x), *(klt_y), prevCam, cam, temp);
// //     DBGPRINT("(1,1)th %f, %f\n", temp[0], temp[1]);
//     
//     for(int i = 0; i < nklt; i++) {
//         if (minhor > *(klt_y + i)){
//             continue;
//         }
//         getPIP(*(klt_x + i), *(klt_y + i), prevCam, cam, temp, camlen);
//         *(res + i * 2) = *(klt_x + i + nklt) - temp[0];
//         *(res + i * 2 + 1) = *(klt_y + i + nklt) - temp[1];
//         cnt++;
//     }
//     
//     DBGPRINT("minhor : %f, valid count %d\n", minhor, cnt);
//     return cnt;
// }
// 
// void getNextCamera(double *prevcam, double dt, double *res) {
//     for(int i = 0; i < 8; i++) {
//         *(res + i) = *(prevcam + i);
//     }
//     *(res + 6) += *(prevcam + 5) * sin(*(prevcam + 4)) * dt;
//     *(res + 7) += *(prevcam + 5) * cos(*(prevcam + 4)) * dt;
// }

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs==6){
        
    }
    else{
		printf("Invalid number of inputs %d\n",nrhs);
		return;
	}
    
    double *cam = (double*)mxGetData(prhs[0]);
    double *pSamples = (double*)mxGetData(prhs[1]);
    double *w = (double*)mxGetData(prhs[3]);
    double *norm = (double*)mxGetData(prhs[4]);
//     double *klt_x = (double*)mxGetData(prhs[5]);
//     double *klt_y = (double*)mxGetData(prhs[6]);
//     double *kltPrec = (double*)mxGetData(prhs[7]);
//     double klttrunc = (double)mxGetScalar(prhs[8]);
    double dt = (double)mxGetScalar(prhs[5]);
    
    double *prec;
    const mxArray *prec_ptr;
    
    double *prior;
    double tempdiff[10], tempmult[10], tempcam[8];
    double tempresult;
    
    msize scam, spSamples, sPrec, sTemp1, sTemp2, skltPrec;
    
    scam.nrows = mxGetM(prhs[0]);    scam.ncols = mxGetN(prhs[0]);    
    spSamples.nrows = mxGetM(prhs[1]);    spSamples.ncols = mxGetN(prhs[1]);
    sPrec.nrows = scam.nrows;    sPrec.ncols = scam.nrows;
    
    plhs[0] = mxCreateDoubleMatrix(1, spSamples.ncols, mxREAL);
    prior = (double*)mxGetData(plhs[0]);
    
    sTemp1.nrows = 1; sTemp1.ncols = scam.nrows;
    for(int i = 0; i < spSamples.ncols; i++) {
        // copmute prior
        prec_ptr = mxGetCell(prhs[2], i);
        prec = (double*)mxGetData(prec_ptr);
        // need to translate!!!!!
//         if (scam.nrows == 8) {
//             getNextCamera((pSamples + (i * spSamples.nrows)), dt, tempcam);
// //             printf("dt = %f, (%f, %f)\n", dt, tempcam[6], tempcam[7]);
//         }
        for(int j = 0; j < scam.nrows; j++) {
            tempdiff[j] = *(cam + j) - *(pSamples + j + (i * spSamples.nrows));
//             tempdiff[j] = *(cam + j) - *(tempcam + j);
//             printf("%f, ", tempdiff[j]);
        }
//         printf("\n");
        MatrixMultiplication(tempdiff, prec, sTemp1, sPrec, tempmult);
        MatrixMultiplication(tempmult, tempdiff, sTemp1, scam, &tempresult);
        
        *(prior + i) = *(w+i) * *(norm + i) * exp(-0.5 * tempresult) + 1e-20;
        
//         printf("%f, %f, %.04f\n", tempresult, *(norm + i), *(prior + i));
    }
    
    return;
//     
//     unsigned int nklt = mxGetM(prhs[5]);
//     if (nklt == 0) {
//         return;
//     }
// //     printf("%d \n", nklt);
// //     printf("%f %f \n", *(klt_x), *(klt_x + 2));
//     
//     sTemp1.nrows = 1; sTemp1.ncols = 2;
//     sTemp2.nrows = 2; sTemp2.ncols = 1;
//     
//     double *tempdiff2 = (double*)mxCalloc(nklt * 2, sizeof(double));
//     double lkhood = 0;
//     double minhor = 100000;
//     for(int i = 0; i < spSamples.ncols; i++) {
//         if (*(pSamples + 3 + (i * spSamples.nrows)) < minhor){
//             minhor = *(pSamples + 3 + (i * spSamples.nrows));
//         }
//     }
//     
//     int cnt = 0;
//     for(int i = 0; i < spSamples.ncols; i++) {
//         lkhood = 0;
//         cnt = getKLTdiff(klt_x, klt_y, nklt, minhor, (pSamples + i * (spSamples.nrows)), cam, tempdiff2, scam.nrows);    
//         for(int j = 0; j < nklt; j++) {
//             MatrixMultiplication((tempdiff2 + j * 2), kltPrec, sTemp1, skltPrec, tempmult);
//             MatrixMultiplication(tempmult, (tempdiff2 + j * 2), sTemp1, sTemp2, &tempresult);
//             if (tempresult > klttrunc) {
//                 tempresult = klttrunc;
//             }
//             lkhood -= 0.5 * tempresult * WEIGHT_FACTOR / cnt;
//         }
//         
//         *(prior + i) *= exp(lkhood);
//     }
//     
//     mxFree(tempdiff2);
}