#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "mmult_accel.h"

/* Begin Extra Functions Definition */
inline int min(int a, int b) { return (a<b) ? a : b;}
inline int max(int a, int b) { return (a>b) ? a : b;}
/* End Extra Functions Definition */

// accelerator
int mmult_accel(data_t * A, data_t * B, data_t * C)
{

	// input buffers
	data_t Buf_A[Tm][Tk], Buf_B[Tk][Tn]; //, Buf_C[Tm][Tn];

	// control loop
	for (int i = 0; i < M; i=i+Tm) { // row
		for (int j = 0; j < N; j=j+Tn) { // column

			// C buffer initialisation
			data_t buff_out[Tm][Tn];
			init_C:for(int ti=0; ti<Tm; ti++) {
				for(int tj=0; tj<Tn; tj++) {
					#pragma HLS PIPELINE
					buff_out[ti][tj] = 0;
				 }
			}

			for (int k = 0; k < K; k=k+Tk) { // inner dot-product

				// load A
				load_A:for(int tk=k; tk<min(k+Tk,K); tk++) {
					for(int ti=i; ti<min(i+Tm,M); ti++) {
						#pragma HLS PIPELINE
						Buf_A[ti-i][tk-k] = *(A+tk*M+ti);
					 }
				}

				// load B
				load_B:for(int tk=k; tk<min(k+Tk,K); tk++) {
					for(int tj=j; tj<min(j+Tn,N); tj++) {
						#pragma HLS PIPELINE
						Buf_B[tk-k][tj-j] = *(B+tk*N+tj);
					 }
				}

				// tile computation
				//std::cout << "!!! Tile begin" << std::endl;
				// basic pipeline
				tile_comp:for(int ti=0; ti<min(Tm,M-i); ti++) {
					for(int tk=0; tk<min(Tk,K-k); tk++) {
						for(int tj=0; tj<min(Tn,N-j); tj++) {
							#pragma HLS PIPELINE
							#pragma HLS DEPENDENCE variable=buff_out array inter RAW false
							buff_out[ti][tj] += Buf_A[ti][tk] * Buf_B[tk][tj];
						}
					}
				}

			} // end inner dot-product

			// store C
			store_C:for(int ti=i; ti<min(i+Tm,M); ti++) {
				for(int tj=j; tj<min(j+Tn,N); tj++) {
					#pragma HLS PIPELINE
					*(C+ti*N+tj) = buff_out[ti-i][tj-j];
				 }
			}

		} // column end
	} // row end

	return 0;
}


