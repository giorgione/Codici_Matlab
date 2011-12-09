#include <math.h>
#include <memory.h>
#include "mex.h"

#define DEBUG_ON 0

double *getProjection(double *state, double *cam, int camlen)
{
// R = [cos(cam(5)), sin(cam(5)); -sin(cam(5)), cos(cam(5))];
// z(1:2) = R * z(1:2);
// 
// x(1) = cam(1) / z(2) * z(1) + cam(3);
// x(2) = cam(1) / z(2) * cam(2) + cam(4);
// x(3) = cam(1) / z(2) * z(5);
// 
// if nargout > 1
//     bb(1) = x(1) - x(3) / 4;
//     bb(2) = x(2) - x(3);
//     bb(3) = x(3) / 2;
//     bb(4) = x(3);
// end
    double temp1[2], temp2[2];
    double *ret = (double *)mxCalloc(3, sizeof(double));
    if(camlen == 8){
        temp1[0] = state[0] - cam[6];
        temp1[1] = state[1] - cam[7];
    }
    
    temp2[0] = cos(cam[4]) *  temp1[0] + sin(cam[4]) * temp1[1];
    temp2[1] = -sin(cam[4]) * temp1[0] + cos(cam[4]) * temp1[1];
    
        
    ret[0] = cam[0] / temp2[1] * temp2[0] + cam[2];
    ret[1] = cam[0] / temp2[1] * cam[1] + cam[3];
    ret[2] = cam[0] / temp2[1] * state[4];
    
    return ret;
}

double getWeightedDistance(double *x, double *y, double *w, double norm)
{
    double dist = 0, temp;
    
    for(int i = 0; i < 3; i++) {
        temp = (x[i] - y[i]);
        dist = dist + w[i] * temp * temp / (norm * norm);
    }
    
    return dist;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    // cam, state, Xobs, Yobs, params
//     int nclass = mxGetNumberOfElements(prhs[0]);
    double *cam = (double*)mxGetData(prhs[0]);
    double *state = (double*)mxGetData(prhs[1]);
    double *Xobs = (double*)mxGetData(prhs[2]);
    int nXobs = mxGetNumberOfElements(prhs[2]);
    double *Yobs = (double*)mxGetData(prhs[3]);
    int nYobs = mxGetNumberOfElements(prhs[3]);
    double *Prec1 = (double*)mxGetData(prhs[4]);
    double *Prec2 = (double*)mxGetData(prhs[5]);
    double *hdet = (double*)mxGetData(prhs[6]);
    double *nhdet = (double*)mxGetData(prhs[7]);
    double *lnormDet = (double*)mxGetData(prhs[8]);
    double *lnormKtrack = (double*)mxGetData(prhs[9]);
    
    double ImObs[3];
    double *pred;
    double *lkhood;
    double temp;
    
    int camlen = mxGetNumberOfElements(prhs[0]);
    
	if(nrhs==10);
    else{
		printf("Invalid number of inputs %d\n",nrhs);
		return;
	}
    
//     printf("len(Xobs) %d\n", nXobs);
            
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    lkhood = (double*)mxGetData(plhs[0]);
    *lkhood = 0;
    
    pred = getProjection(state, cam, camlen);
    
    if (nXobs > 0)    {
        ImObs[0] = Xobs[0] + Xobs[2] / 2;
        ImObs[1] = Xobs[1] + Xobs[3];
        ImObs[2] = Xobs[3];
        *lkhood = -0.5 * (getWeightedDistance(ImObs, pred, Prec1, ImObs[2]/128)) + *lnormDet;
        if (state[5] == 1) {
            *lkhood += log(*hdet);
        }
        else {
            *lkhood += log(*nhdet);
        }
    }
    else    {
        if (state[5] == 1) {
            *lkhood = log(1 - *hdet);
        }
        else {
            *lkhood = log(1 - *nhdet);
        }
    }
    
    temp = 0;
    for(int i = 0; i < 2; i++) {
        temp += Yobs[i];
    }
    
    if (temp != 0) {
        ImObs[0] = Yobs[0];
        ImObs[1] = Yobs[1] + Yobs[2] / 2;
        ImObs[2] = Yobs[2];
        *lkhood += -0.5 * (getWeightedDistance(ImObs, pred, Prec2, ImObs[2]/128)) + *lnormKtrack;
    }
    
//     printf("%.05f, %.05f, %.05f\n", pred[0], pred[1], pred[2]);
    
    mxFree(pred);
}