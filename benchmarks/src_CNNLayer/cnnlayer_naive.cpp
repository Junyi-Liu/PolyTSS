/*
 * cnnlayer_naive.cpp
 *
 *  Created on: 3 July 2017
 *      Author: jl12013
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "cnnlayer_accel.h"

/* Begin Extra Functions Definition */
inline int min(int a, int b) { return (a<b) ? a : b;}
inline int max(int a, int b) { return (a>b) ? a : b;}
/* End Extra Functions Definition */


// accelerator
int cnnlayer_naive(data_t *frame_in, data_t *weight, data_t *frame_out)
{
	// buffers
	data_t buf_fin[N][inRow][inCol];
	data_t buf_wgt[M][N][K][K];
	data_t buf_fout[M][outRow][outCol];

	
	//frame_out buffer initialisation
	init_fout:for (int tm=0; tm<M; tm++){
		for (int tr=0; tr<outRow; tr++){
			for (int tc=0; tc<outCol; tc++){
				#pragma HLS PIPELINE
				buf_fout[tm][tr][tc] = 0;
			}
		}
	}

	// read frame_in
	read_fin:for (int tn=0; tn<N; tn++){
		for (int tr=0; tr<inRow; tr++){
			for (int tc=0; tc<inCol; tc++){
				#pragma HLS PIPELINE
				buf_fin[tn][tr][tc] = *(frame_in + tn*inRow*inCol+tr*inCol+tc);
			}
		}
	}

	// read weight
	read_wgt:for (int tm=0; tm<M; tm++) {
		for (int tn=0; tn<N; tn++) {
			for(int ti=0; ti<K; ti++) {
				for(int tj=0; tj<K; tj++) {
					#pragma HLS PIPELINE
					buf_wgt[tm][tn][ti][tj] = *(weight + tm*N*K*K+tn*K*K+ti*K+tj);
				}
			}
		}
	}
	
	
	// control loop
	//std::cout << "!!! Loop begin" << std::endl;
	// basic pipeline
	loop_comp:for(int m=0; m<M; m++) {
		for(int row=0; row<outRow; row++) {
			for(int col=0; col<outCol; col++) {
				for(int n=0; n<N; n++) {
					for(int i=0; i<K; i++) {
						for(int j=0; j<K; j++) {
							#pragma HLS PIPELINE
							buf_fout[m][row][col] += buf_wgt[m][n][i][j] * buf_fin[n][ST*row+i][ST*col+j];
						}
					}
				}
			}
		}
	}

	// write frame_out
//				std::cout << "!!!Write results" << std::endl;
	write_fout:for (int tm=0; tm<M; tm++) {
		for (int tr=0; tr<outRow; tr++) {
			for (int tc=0; tc<outCol; tc++) {
				#pragma HLS PIPELINE
				*(frame_out + tm*outRow*outCol+tr*outCol+tc) = buf_fout[tm][tr][tc];
			}
		}
	}
	
	return 0;
}


