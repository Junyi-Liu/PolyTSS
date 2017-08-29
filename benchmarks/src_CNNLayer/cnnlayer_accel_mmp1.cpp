/*
 * cnnlayer_accel_mmp1.cpp
 *
 *  Created on: 4 Apr 2017
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
int cnnlayer_accel_mmp1(data_t *frame_in, data_t *weight, data_t *frame_out)
{

	// buffers
	//data_t buf_fin[Tn][ST*(Tr-1)+Ti][ST*(Tc-1)+Tj];
	data_t buf_fin[Tn][ST*(Tc-1)+Tj]; // Ti = 1, Tr = 1
	//data_t buf_wgt[Tm][Tn][Ti][Tj];
	data_t buf_wgt[Tm][Tn][Tj]; // Ti = 1

	// control loop
	for(int m=0; m<M; m+=Tm) {
		for(int row=0; row<outRow; row+=Tr) {
			for(int col=0; col<outCol; col+=Tc) {

				//frame_out buffer initialisation
				//data_t buf_fout[Tm][Tr][Tc];
				data_t buf_fout[Tm][Tc]; // Tr = 1
				init_fout:for (int tm=0; tm<Tm; tm++){
					//for (int tr=0; tr<Tr; tr++){
						for (int tc=0; tc<Tc; tc++){
							#pragma HLS PIPELINE
							buf_fout[tm][tc] = 0;
						}
					//}
				}

				//
				for(int n=0; n<N; n+=Tn) {
					for(int i=0; i<K; i+=Ti) {
						for(int j=0; j<K; j+=Tj) {

							// read frame_in
							read_fin:for (int tn=n; tn<min(n+Tn,N); tn++){
								for (int tr=ST*row+i; tr<min(ST*(row+Tr-1)+i+Ti,inRow); tr++){
									for (int tc=ST*col+j; tc<min(ST*(col+Tc-1)+j+Tj,inCol); tc++){
										#pragma HLS PIPELINE
										//buf_fin[tn-n][tr-ST*row-i][tc-ST*col-j]
										buf_fin[tn-n][tc-ST*col-j] = *(frame_in + tn*inRow*inCol+tr*inCol+tc); // Ti = 1, Tr = 1
									}
								}
							}

							// read weight
							read_wgt:for (int tm=m; tm<min(m+Tm,M); tm++) {
								for (int tn=n; tn<min(n+Tn,N); tn++) {
									for(int ti=i; ti<min(i+Ti,K); ti++) {
										for(int tj=j; tj<min(j+Tj,K); tj++) {
											#pragma HLS PIPELINE
											//buf_wgt[tm-m][tn-n][ti-i][tj-j]
											buf_wgt[tm-m][tn-n][tj-j] = *(weight + tm*N*K*K+tn*K*K+ti*K+tj); // Ti = 1
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
							//for(int ti=0; ti<Ti; ti++) {
							tile_comp:for(int tj=0; tj<Tj; tj++) {
										//for (int tr=0; tr<Tr; tr++) {
										for (int tc=0; tc<Tc; tc++) {
										#pragma HLS PIPELINE
											for (int tm=0; tm<Tm; tm++) {
											#pragma HLS unroll factor=16
												data_t tmp = 0;
												for (int tn=0; tn<N; tn++) {
													 if(tc+col<outCol && tj+j<K && tm+m<M && tn+n<N)
														 tmp += buf_wgt[tm][tn][tj] * buf_fin[tn][ST*tc+tj];
												}
												if(tc+col<outCol && tj+j<K && tm+m<M)
													buf_fout[tm][tc] += tmp;
											}
										}
									//}
								}
							//}

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
							*(frame_out + tm*outRow*outCol+tr*outCol+tc) = buf_fout[tm-m][tc-col];
						}
					}
				}


	}}}


	return 0;
}





