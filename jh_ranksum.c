/*==========================================================
 * jh_ranksum.c
 * July 2012, (c) Joerg Hipp
 *
 * This function computes statistical z-values for the 
 * Mann-Whitney U test (also called the Mann-Whitney-Wilcoxon (MWW) 
 * or Wilcoxon rank-sum test). This is a test for independent samples.
 * For dependent samples use the Wilcoxon signed-rank test (signrank.m)
 * Note, that this implementation uses an approximation that is valid only
 * for > 10 realization per group. Note, this implementation does not
 * correct for tied ranks (should be fine for continous values).
 *
 * http://en.wikipedia.org/wiki/Mann-Whitney_U
 *
 * To rebuild type "mex jh_ranksum.c"
 * 
 * Usage:
 * outMatrix = jh_ranksum(inMatrix);
 * 
 * inMatrix .. matrix with single precision input data [n_vox,n_vox,n_subjects]
 *             equal amount of subjects are assumed, while first all data from the
 *             one group and then all data from the other group.
 * outMatrix .. single precision output (statistical z values [n_vox,n_vox]
 *
 * To derive p-values in matlab: p = 2*normcdf(-abs(outMatrix));
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Functions */
void quickSort(float *arr, mwSize elements);

/*************************************************************************/

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    float *inMatrix;  /* Input*/
    float *data, *data2, w, wmean, wvar, wc, z;
    float *outMatrix; /* Output*/
    mwSize icnt, jcnt, kcnt, lcnt, idx, idx2, done, tmp;
    float nx, ny; /* number of subjects, is assubed to be equal to ny, hardcoded! */
    mwSize nDim;
    mwSize dim[2];
    
    /* Parse the input  */
    inMatrix = mxGetData(prhs[0]); /* create a pointer to the input  */
    nDim = mxGetNumberOfDimensions(prhs[0]);
    if (nDim==3) {
        dim[0] = (mxGetDimensions(prhs[0]))[0]; /* vox */
        dim[1] = (mxGetDimensions(prhs[0]))[1]; /* vox */
        dim[2] = (mxGetDimensions(prhs[0]))[2]; /* subjects */
    } else {
        mexPrintf("Error! Wrong input!\n");
        return;
    }
    mexPrintf("\nData: %i x %i voxel, %i subjects\n", dim[0],dim[1],dim[2]);
    mexEvalString("drawnow");
    
    /* Output memory allocation */
    plhs[0] = mxCreateNumericMatrix(dim[0],dim[1],mxSINGLE_CLASS,mxREAL);
    outMatrix = (float *)mxGetData(plhs[0]);
    
    data = mxCalloc(dim[2],sizeof(float)); /* memory allocation*/
    data2 = mxCalloc(dim[2],sizeof(float)); /* memory allocation*/
    nx = dim[2]/2;
    ny = nx;
    wmean = nx*(nx + ny + 1)/2;
    wvar  = nx*ny*((nx + ny + 1))/12;    
           
    for (icnt=0; icnt<dim[0]; icnt++)  {
        mexPrintf("%i of %i\n",icnt+1,dim[0]);
        mexEvalString("drawnow");
        for (jcnt=0; jcnt<dim[1]; jcnt++)  {
            /* extract data */
            for (kcnt=0; kcnt<dim[2]; kcnt++)  {
                idx = icnt+jcnt*dim[0]+kcnt*dim[1]*dim[0];
                data[kcnt]=inMatrix[idx];
                data2[kcnt]=inMatrix[idx];
            }
            /* sort data */
            quickSort(&data[0],dim[2]);
            /* extract rank */
            w=0;
            for (kcnt=0; kcnt<nx; kcnt++)  {
                for (lcnt=0; lcnt<dim[2]; lcnt++)  {
                    if (data2[kcnt]==data[lcnt]) {
                        w=w+((float)lcnt+1);
                        break;
                    }
                }
            }
            wc = w - wmean;
            idx = icnt+jcnt*dim[0];
            outMatrix[idx] = (wc-0.5*((wc>0)-(wc<0)))/sqrt(wvar);
        }
    }
    
}


/*-------------------------------------------------------------------*/
/*  quickSort
    
    This public-domain C implementation by Darel Rex Finley.
    http://alienryderflex.com/quicksort/
 
    * Returns YES if sort was successful, or NO if the nested
      pivots went too deep, in which case your array will have
      been re-ordered, but probably not sorted correctly.
  
    * This function assumes it is called with valid parameters.
  
    * Example calls:
      quickSort(&myArray[0],5); // sorts elements 0, 1, 2, 3, and 4
      quickSort(&myArray[3],5); // sorts elements 3, 4, 5, 6, and 7
*/

void quickSort(float *arr, mwSize elements) {
    
    #define  MAX_LEVELS  1000

    float piv;
    mwSize  beg[MAX_LEVELS], end[MAX_LEVELS], i=0, L, R;

    beg[0]=0; end[0]=elements;
    while (i>=0) {
        L=beg[i]; R=end[i]-1;
        if (L<R) {
            piv=arr[L]; if (i==MAX_LEVELS-1) return;
            while (L<R) {
                while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
                while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; }
            arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L; }
        else {
            i--; }}
    return; 
}

