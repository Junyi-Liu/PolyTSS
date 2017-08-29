/*
 * cnnlayer_accel_mmp2.cpp
 *
 *  Created on: 5 Apr 2017
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
int cnnlayer_accel_mmp2(data_t *frame_in, data_t *weight, data_t *frame_out)
{

	// buffers
	//data_t buf_fin[Tn][ST*(Tr-1)+Ti][ST*(Tc-1)+Tj];
	data_t buf_fin[Tn][Tr][Ti][Tc][Tj]; // Ti = 3, Tj = 3
	data_t buf_wgt[Tm][Tn][Ti][Tj];

	// control loop
	for(int m=0; m<M; m+=Tm) {
		for(int row=0; row<outRow; row+=Tr) {
			for(int col=0; col<outCol; col+=Tc) {

				//frame_out buffer initialisation
				data_t buf_fout[Tm][Tr][Tc];
				init_fout:for (int tm=0; tm<Tm; tm++){
					for (int tr=0; tr<Tr; tr++){
						for (int tc=0; tc<Tc; tc++){
							#pragma HLS PIPELINE
							buf_fout[tm][tr][tc] = 0;
						}
					}
				}

				//
				for(int n=0; n<N; n+=Tn) {
					for(int i=0; i<K; i+=Ti) {
						for(int j=0; j<K; j+=Tj) {

							// read frame_in
							read_fin:for (int tn=n; tn<min(n+Tn,N); tn++){
								for (int tr=row; tr<min(row+Tr,outRow); tr++){
									for(int ti=i; ti<min(i+Ti,K); ti++) {
										for (int tc=col; tc<min(col+Tc,outCol); tc++){
											for(int tj=j; tj<min(j+Tj,K); tj++) {
											#pragma HLS PIPELINE
												buf_fin[tn-n][tr-row][ti-i][tc-col][tj-j] = *(frame_in + tn*inRow*inCol+(ST*tr+ti)*inCol+(ST*tc+tj));
											}
										}
									}
								}
							}

							// read weight
							read_wgt:for (int tm=m; tm<min(m+Tm,M); tm++) {
								for (int tn=n; tn<min(n+Tn,N); tn++) {
									for(int ti=i; ti<min(i+Ti,K); ti++) {
										for(int tj=j; tj<min(j+Tj,K); tj++) {
											#pragma HLS PIPELINE
											buf_wgt[tm-m][tn-n][ti-i][tj-j] = *(weight + tm*N*K*K+tn*K*K+ti*K+tj);
										}
									}
								}
							}

							//execute loop tile
//							std::cout << "!!! Tile begin" << std::endl;
							// unroll & pipeline
//#pragma HLS RESOURCE variable=buf_fin core=RAM_2P_BRAM
//#pragma HLS RESOURCE variable=buf_fout core=RAM_2P_BRAM
//#pragma HLS RESOURCE variable=buf_wgt core=RAM_2P_BRAM
#pragma HLS array_partition variable=buf_fin cyclic factor=3 dim=1
#pragma HLS array_partition variable=buf_fout cyclic factor=16 dim=1
#pragma HLS array_partition variable=buf_wgt cyclic factor=16 dim=1
#pragma HLS array_partition variable=buf_wgt cyclic factor=3 dim=2
							tile_comp:for(int ti=0; ti<Ti; ti++) {
								for(int tj=0; tj<Tj; tj++) {
									for (int tr=0; tr<Tr; tr++) {
										for (int tc=0; tc<Tc; tc++) {
										#pragma HLS PIPELINE
											for (int tm=0; tm<Tm; tm++) {
											#pragma HLS unroll factor=16
												data_t tmp = 0;
												for (int tn=0; tn<N; tn++) {
													 if(tr+row<outRow && tc+col<outCol && ti+i<K && tj+j<K && tm+m<M && tn+n<N)
														 tmp += buf_wgt[tm][tn][ti][tj] * buf_fin[tn][tr][ti][tc][tj];
												}
												if(tr+row<outRow && tc+col<outCol && ti+i<K && tj+j<K && tm+m<M)
													buf_fout[tm][tr][tc] += tmp;
											}
										}
									}
								}
							}


		        }}}


				// write frame_out
//				std::cout << "!!!Write results" << std::endl;
//				for (int tm=m; tm<min(m+Tm,M); tm++){
//					for (int tr=row; tr<(row+Tr,outRow); tr++){
//						int len;
//						if(col+Tc<outCol)
//							len = Tc;
//						else
//							len = outCol-col;
//						memcpy(frame_out+tm*outRow+tr*outCol, buf_fout+(tm-m)*Tr+(tr-row)*Tc, len*sizeof(data_t));
//					}
//				}
				write_fout:for (int tm=m; tm<min(m+Tm,M); tm++) {
					for (int tr=row; tr<min(row+Tr,outRow); tr++) {
						for (int tc=col; tc<min(col+Tc,outCol); tc++) {
							#pragma HLS PIPELINE
							*(frame_out + tm*outRow*outCol+tr*outCol+tc) = buf_fout[tm-m][tr-row][tc-col];
						}
					}
				}


	}}}


	return 0;
}



