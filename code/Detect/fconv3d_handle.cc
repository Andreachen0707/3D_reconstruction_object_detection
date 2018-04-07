//mex -largeArrayDims CFLAGS='-O3 -ffast-math' fconv3d_handle.cc
#include "mex.h"
#include <math.h>
#include <string.h>

// convolve A and B
void process(double* A, double* B, double* MA, bool* MB, bool* E, bool* DontCare, double* C, const mwSize* A_dims, const mwSize* B_dims, const mwSize* C_dims, double* Space, double* scalemissing) {
	double* dst = C;
	bool* isE = E;
//	int count =0;
	for (int z = 0; z < C_dims[2]; z++) {
		for (int x = 0; x < C_dims[1]; x++) {
			for (int y = 0; y < C_dims[0]; y++) {
				if (*(isE++)){
					*(dst++) = -100; // skip location without points
				}
				else{
					//for each cell
					double val = 0;
					for (int zp = 0; zp < B_dims[3]; zp++) {
						double* A_off_z = A + (z+zp)*A_dims[0]*A_dims[1]*A_dims[2];
						double* B_off_z = B + zp*B_dims[0]*B_dims[1]*B_dims[2];
						double* MA_off_z = MA + (z+zp)*A_dims[1]*A_dims[2];
						bool* MB_off_z = MB + zp*B_dims[1]*B_dims[2];
						bool* MD_off_z = DontCare + zp*B_dims[1]*B_dims[2];

						for (int xp = 0; xp < B_dims[2]; xp++) {
							double* A_off_x = A_off_z + (x+xp)*A_dims[0]*A_dims[1];
							double* B_off_x = B_off_z + xp*B_dims[0]*B_dims[1];
							double* MA_off_x = MA_off_z + (x+xp)*A_dims[1];
							bool* MB_off_x = MB_off_z + xp*B_dims[1];
							bool* MD_off_x = MD_off_z + xp*B_dims[1];

							for (int yp = 0; yp < B_dims[1]; yp++) {
								double* A_off_y = A_off_x + (y+yp)*A_dims[0];
								double* B_off_y = B_off_x + yp*B_dims[0];
								double* MA_off_y = MA_off_x + (y+yp);
								bool* MB_off_y = MB_off_x + yp;
								bool* MD_off_y = MD_off_x + yp;
								if ((*MD_off_y)){
									val += 0;
								}

								else{
/*
if (y ==10&&x == 5&&z==0){ 
//printf("x: %d, y: %d z: %d\n",xp,yp,zp);

if (*(A_off_y+B_dims[0]-1)>0.1){
	printf("non empty: %d\n",++count2);
// printf("empty dim: %f\n", *(B_off_y+B_dims[0]));
    printf("empty: %f\n", *(A_off_y+B_dims[0]-1));
//    printf("miss: %f\n", *(A_off_y+B_dims[0]-2));
	printf("x: %d, y: %d z: %d\n",xp,yp,zp);
}
//printf("Occ Source: %f, boxY: %f maskModel: %d \n",(*MA_off_y),y*Space[0]+Space[1],(*MB_off_y));
}	
*/
									if ((*MA_off_y)>0.001 && !(*MB_off_y) && (y*Space[0] +Space[1]> (*MA_off_y))) {
										val += *(B_off_y+B_dims[0]-1)*scalemissing[0];
										// this cell is don't care (occluded by souce outside box) only multiply by the information bit
									}	
									else {					
	                                    for (int f = 0; f < B_dims[0]; f++){
	                                        val += *(A_off_y++) * *(B_off_y++);
	                                	}

									}
                            
                                }
							}
						}
					}
					
					*(dst++) = val;
				}
			}
		}
	}
}

// matlab entry point
// C = fconv_occ(A, B, MA, MB, E, dontCare, SpaceY);
// A is feature matrix
// B is filter
// MA is mask for feature:  0 not occluded, >0 for occluded cell, value depth (y value) of occluder.
// MB is mask for filter:  1 for occluded cell, 0 otherwise.
// E is isEmpty:  1 for empty position, 0 otherwise.
// Space = [s, minY]
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 
	if (nrhs != 8) mexErrMsgTxt("Wrong number of inputs"); 
	if (nlhs != 1) mexErrMsgTxt("Wrong number of outputs");

	// get A
	const mxArray *mxA = prhs[0];
	if (mxGetNumberOfDimensions(mxA) != 4 || mxGetClassID(mxA) != mxDOUBLE_CLASS) mexErrMsgTxt("Invalid input: A");
	// get B
	const mxArray *mxB = prhs[1];
	if (mxGetNumberOfDimensions(mxB) != 4 || mxGetClassID(mxB) != mxDOUBLE_CLASS) mexErrMsgTxt("Invalid input: B");
	// get MA
	const mxArray *mxMA = prhs[2];
	if (mxGetNumberOfDimensions(mxMA) != 3 || mxGetClassID(mxMA) != mxDOUBLE_CLASS) mexErrMsgTxt("Invalid input: MA");
	// get MB
	const mxArray *mxMB = prhs[3];
	if (mxGetNumberOfDimensions(mxMB) != 3 || mxGetClassID(mxMB) != mxLOGICAL_CLASS) mexErrMsgTxt("Invalid input: MB");
	// get E
	const mxArray *mxE = prhs[4];
	if (mxGetNumberOfDimensions(mxE) != 3 || mxGetClassID(mxMB) != mxLOGICAL_CLASS) mexErrMsgTxt("Invalid input: E");
	
	const mxArray *mxDontCare = prhs[5];
	if (mxGetNumberOfDimensions(mxDontCare) != 3 || mxGetClassID(mxDontCare) != mxLOGICAL_CLASS) mexErrMsgTxt("Invalid input: Dontcare");

	const mxArray *mxSpace = prhs[6];
	//if (mxGetNumberOfDimensions(mxDontCare) != 3 || mxGetClassID(mxDontCare) != mxLOGICAL_CLASS) mexErrMsgTxt("Invalid input: Dontcare");

	const mxArray *mxScalemissing = prhs[7];
	
	const mwSize* A_dims = mxGetDimensions(mxA);
	double* A = (double*)mxGetPr(mxA);
	const mwSize* B_dims = mxGetDimensions(mxB);
	double* B = (double*)mxGetPr(mxB);
	const mwSize* MA_dims = mxGetDimensions(mxMA);
	double* MA = (double*)mxGetPr(mxMA);
	const mwSize* MB_dims = mxGetDimensions(mxMB);
	bool* MB = (bool*)mxGetPr(mxMB);
	const mwSize* E_dims = mxGetDimensions(mxE);
	bool* E = (bool*)mxGetPr(mxE);
	const mwSize* DontCare_dim = mxGetDimensions(mxDontCare);
	bool* DontCare = (bool*)mxGetPr(mxDontCare);
	const mwSize* Space_dim = mxGetDimensions(mxSpace);
	double* Space = (double*)mxGetPr(mxSpace);
	double* scalemissing = (double*)mxGetPr(mxScalemissing);


	if (A_dims[0] != B_dims[0]) mexErrMsgTxt("Invalid input: A and B should have the same feature dimension");
	if (A_dims[1] < B_dims[1] || A_dims[2] < B_dims[2] || A_dims[3] < B_dims[3]) mexErrMsgTxt("Invalid input: B should be smaller than A");
	for (int i=0; i<3; i++){
		if (MA_dims[i] != A_dims[i+1]) mexErrMsgTxt("Invalid input: MA should have the same 3D size as A");
		if (MB_dims[i] != B_dims[i+1]) mexErrMsgTxt("Invalid input: MB should have the same 3D size as B");
		if (DontCare_dim[i] != B_dims[i+1] ) mexErrMsgTxt("Invalid input: B should have the same 3D size as Dontcare");
	}
	
	// create output
	mwSize C_dims[3];
	C_dims[0] = A_dims[1] - B_dims[1] + 1;
	C_dims[1] = A_dims[2] - B_dims[2] + 1;
	C_dims[2] = A_dims[3] - B_dims[3] + 1;
	for (int i=0; i<3; i++){
		if (C_dims[i] != E_dims[i]) mexErrMsgTxt("Invalid input: E should have the same 3D size as convolution result");
	}
	plhs[0] = mxCreateNumericArray(3, C_dims, mxDOUBLE_CLASS, mxREAL);
	double* C = (double*)mxGetPr(plhs[0]);

	// do convolutions
	process(A, B, MA, MB, E, DontCare, C, A_dims, B_dims, C_dims,Space,scalemissing);
}


