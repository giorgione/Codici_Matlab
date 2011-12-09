#include <math.h>
#include <memory.h>
#include "mex.h"

#define DEBUG_ON 0

double *getProjection(double *state, double *cam, int camlen)
{
    double temp1[2], temp2[2];
    double *ret = (double *)mxCalloc(2, sizeof(double));
    if(camlen == 8){
        temp1[0] = state[0] - cam[6];
        temp1[1] = state[1] - cam[7];
    }
    
    temp2[0] = cos(cam[4]) *  temp1[0] + sin(cam[4]) * temp1[1];
    temp2[1] = -sin(cam[4]) * temp1[0] + cos(cam[4]) * temp1[1];
    
        
    ret[0] = cam[0] / temp2[1] * temp2[0] + cam[2];
    ret[1] = cam[0] / temp2[1] * cam[1] + cam[3];
    
    return ret;
}

double getWeightedDistance(double *x, double *y, double *w)
{
    double dist = 0, temp;
    
    for(int i = 0; i < 2; i++) {
        temp = (x[i] - y[i]);
        dist = dist + w[i] * temp * temp;
    }
    
    return dist;
}

double getObservationLkhood(double *cam, double *state, int camlen, double *obs, int nObs, double lnormGF, double *Prec, double pOut)
{
    double ret = 0;
    if (nObs > 0)    {
        if (state[2] == 1) {
            double *pred = getProjection(state, cam, camlen);
//             printf("on Image %f %f \n", pred[0], pred[1]);
            ret = -0.5 * getWeightedDistance(obs, pred, Prec) + lnormGF;
            mxFree(pred);
        }
        else {
            ret = pOut;
        }
    }
    
    return ret;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs==6);
    else{
		printf("Invalid number of inputs %d\n",nrhs);
		return;
	}
    
    double *cam = (double*)mxGetData(prhs[0]);
    double *state = (double*)mxGetData(prhs[1]);
    double *obs = (double*)mxGetData(prhs[2]);
    int nObs = mxGetNumberOfElements(prhs[2]);
    double *Prec = (double*)mxGetData(prhs[3]);
    double lnormGF = (double)mxGetScalar(prhs[4]);
    double pOut = (double)mxGetScalar(prhs[5]);
    
    double *lkhood;
    
    int camlen = mxGetNumberOfElements(prhs[0]);
//     printf("len(Xobs) %d\n", nXobs);
            
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    lkhood = (double*)mxGetData(plhs[0]);
    *lkhood = getObservationLkhood(cam, state, camlen, obs, nObs, lnormGF, Prec, pOut);
}