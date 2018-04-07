#include "mex.h"
#include <math.h>
#include <string.h>

/*
 * This code is used for computing filter responses.  It computes the
 * response of a set of filters with a feature map.  
 *
 * Basic version, relatively slow but very compatible.
 */

struct thread_data {
  double *A;
  double *B;
  double *C;
  mxArray *mxC;
  const mwSize *A_dims;
  const mwSize *B_dims;
  mwSize C_dims[2];
};

// convolve A and B
void process(void *thread_arg) {
  thread_data *args = (thread_data *)thread_arg;
  double *A = args->A;
  double *B = args->B;
  double *C = args->C;
  const mwSize *A_dims = args->A_dims;
  const mwSize *B_dims = args->B_dims;
  const mwSize *C_dims = args->C_dims;
  int num_features = args->A_dims[2];

  for (int f = 0; f < num_features; f++) {
    double *dst = C;
    double *A_src = A + f*A_dims[0]*A_dims[1];      
    double *B_src = B + f*B_dims[0]*B_dims[1];
    for (int x = 0; x < C_dims[1]; x++) {
      for (int y = 0; y < C_dims[0]; y++) {
          double val = 0;
          for (int xp = 0; xp < B_dims[1]; xp++) {
              double *A_off = A_src + (x+xp)*A_dims[0] + y;
              double *B_off = B_src + xp*B_dims[0];
              for (int yp = 0; yp < B_dims[0]; yp++) {
                  val += *(A_off++) * *(B_off++);
              }
          }
          *(dst++) += val;
      }
    }
  }
}

// matlab entry point
// C = fconv(A, cell of B, start, end);
// match = fconv3d(featureM, {templateM}, 1, 1);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs != 4)
    mexErrMsgTxt("Wrong number of inputs"); 
  if (nlhs != 1)
    mexErrMsgTxt("Wrong number of outputs");

  // get A
  const mxArray *mxA = prhs[0];
  if (mxGetNumberOfDimensions(mxA) != 3 || 
      mxGetClassID(mxA) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Invalid input: A");

  // get B and start/end
  const mxArray *cellB = prhs[1];
  mwSize num_bs = mxGetNumberOfElements(cellB);  
  int start = (int)mxGetScalar(prhs[2]) - 1;
  int end = (int)mxGetScalar(prhs[3]) - 1;
  if (start < 0 || end >= num_bs || start > end)
    mexErrMsgTxt("Invalid input: start/end");
  int len = end-start+1;

  // output cell
  plhs[0] = mxCreateCellMatrix(1, len);

  // do convolutions
  thread_data td;
  const mwSize *A_dims = mxGetDimensions(mxA);
  double *A = (double *)mxGetPr(mxA);
  for (int i = 0; i < len; i++) {
    const mxArray *mxB = mxGetCell(cellB, i+start);
    td.A_dims = A_dims;
    td.A = A;
    td.B_dims = mxGetDimensions(mxB);
    td.B = (double *)mxGetPr(mxB);
    if (mxGetNumberOfDimensions(mxB) != 3 ||
        mxGetClassID(mxB) != mxDOUBLE_CLASS ||
        td.A_dims[2] != td.B_dims[2])
      mexErrMsgTxt("Invalid input: B");

    // compute size of output
    int height = td.A_dims[0] - td.B_dims[0] + 1;
    int width = td.A_dims[1] - td.B_dims[1] + 1;
    if (height < 1 || width < 1)
      mexErrMsgTxt("Invalid input: B should be smaller than A");
    td.C_dims[0] = height;
    td.C_dims[1] = width;
    td.mxC = mxCreateNumericArray(2, td.C_dims, mxDOUBLE_CLASS, mxREAL);
    td.C = (double *)mxGetPr(td.mxC);
    process((void *)&td);
    mxSetCell(plhs[0], i, td.mxC);
  }
}


