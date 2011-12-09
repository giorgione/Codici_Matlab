#include <math.h>
#include <memory.h>
#include "mex.h"
#include "matrixwg.h"

void MatrixMultiplication(double *A, double *B, msize sA, msize sB, double *C)
{
    int i = 0, j = 0, k = 0;
    double temp = 0;
    if (sA.ncols != sB.nrows) {
        printf("Invalid matrix inputs!\n");
        return;
    }    
    for(i = 0; i < sA.nrows; i++) {
        for(j = 0; j < sB.ncols; j++) {
            temp = 0;
            for(k = 0; k < sA.ncols; k++) {
                temp += *(A + (k * sA.nrows) + i) * *(B + (j * sB.nrows) + k);
            }
            *(C + (j * sA.nrows) + i) = temp;
        }
    }    
}