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
	mwSize C_dims[3];
	int f_min;
	int f_max;
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
	int num_features = args->A_dims[3];
	int f_min = args->f_min;
	int f_max = args->f_max;

	for (int f = f_min; f < f_max; f++){
		double *dst = C;
		double *A_src = A + f*A_dims[0]*A_dims[1]*A_dims[2];      
		double *B_src = B + f*B_dims[0]*B_dims[1]*B_dims[2];
		for (int z = 0; z < C_dims[2]; z++) {
			for (int x = 0; x < C_dims[1]; x++) {
				for (int y = 0; y < C_dims[0]; y++) {
					double val = 0;
					for (int zp = 0; zp < B_dims[2]; zp++) {
						double *A_off_z = A_src + (z+zp)*A_dims[0]*A_dims[1];
						double *B_off_z = B_src + zp*B_dims[0]*B_dims[1];
						for (int xp = 0; xp < B_dims[1]; xp++) {
							double *A_off = A_off_z + (x+xp)*A_dims[0] + y;
							double *B_off = B_off_z + xp*B_dims[0];
							for (int yp = 0; yp < B_dims[0]; yp++) {
								val += *(A_off++) * *(B_off++);
							}
						}
					}
					*(dst++) += val;
				}
			}
		}
	}
}

// matlab entry point
// C = fconv(A, cell of B, start, end);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
  if (nrhs != 4)
    mexErrMsgTxt("Wrong number of inputs"); 
  if (nlhs != 1)
    mexErrMsgTxt("Wrong number of outputs");

  // get A
  const mxArray *mxA = prhs[0];
  if (mxGetNumberOfDimensions(mxA) != 4 || 
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
    if (mxGetNumberOfDimensions(mxB) != 4 ||
        mxGetClassID(mxB) != mxDOUBLE_CLASS ||
        td.A_dims[3] != td.B_dims[3])
      mexErrMsgTxt("Invalid input: B");

    // compute size of output
    int x_dim = td.A_dims[1] - td.B_dims[1] + 1;
    int y_dim = td.A_dims[0] - td.B_dims[0] + 1;
    int z_dim = td.A_dims[2] - td.B_dims[2] + 1;
    if (x_dim < 1 || y_dim < 1 || z_dim < 1)
		mexErrMsgTxt("Invalid input: B should be smaller than A");
    td.C_dims[0] = y_dim;
    td.C_dims[1] = x_dim;
    td.C_dims[2] = z_dim;
    td.mxC = mxCreateNumericArray(3, td.C_dims, mxDOUBLE_CLASS, mxREAL);
    td.C = (double *)mxGetPr(td.mxC);
    td.f_min = 0;
    td.f_max = td.A_dims[3];
    process((void *)&td);
    mxSetCell(plhs[0], i, td.mxC);
  }
}


