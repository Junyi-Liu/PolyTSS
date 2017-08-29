#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mmult_accel.h"

/* Begin Extra Functions Definition */
inline int min(int a, int b) { return (a<b) ? a : b;}
inline int max(int a, int b) { return (a>b) ? a : b;}
/* End Extra Functions Definition */

// accelerator
int mmult_naive(data_t * A, data_t * B, data_t * C)
{

	// input buffers
	data_t Buf_A[M][K], Buf_B[K][N]; //, Buf_C[Tm][Tn];

	
	// C buffer initialisation
	data_t buff_out[M][N];
	init_C:for(int i=0; i<M; i++) {
		for(int j=0; j<N; j++) {
			#pragma HLS PIPELINE
			buff_out[i][j] = 0;
		 }
	}

	load_A:for(int k=0; k<K; k++) {
		for(int i=0; i<M; i++) {
			#pragma HLS PIPELINE
			Buf_A[i][k] = *(A+k*M+i);
		 }
	}

	// load B
	load_B:for(int k=0; k<K; k++) {
		for(int j=0; j<N; j++) {
			#pragma HLS PIPELINE
			Buf_B[k][j] = *(B+k*N+j);
		 }
	}	
	
	// compute loop
	for (int i = 0; i < M; i++) { // row
		for (int j = 0; j < N; j++) { // column
			for (int k = 0; k < K; k++) { // inner dot-product
				#pragma HLS PIPELINE
				buff_out[i][j] += Buf_A[i][k] * Buf_B[k][j];
			} // end inner dot-product
		} // column end
	} // row end

	// store C
	store_C:for(int i=0; i<M; i++) {
		for(int j=0; j<N; j++) {
			#pragma HLS PIPELINE
			*(C+i*N+j) = buff_out[i][j];
		 }
	}
	
	return 0;
}


