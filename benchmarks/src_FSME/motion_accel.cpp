/*
 * motion_accel.cpp
 *
 *  Created on: 26 Aug 2017
 *      Author: jl12013
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "motion_accel.h"

/* Begin Extra Functions Definition */
inline int min(int a, int b) { return (a<b) ? a : b;}
inline int max(int a, int b) { return (a>b) ? a : b;}
/* End Extra Functions Definition */


// accelerator
int motion_accel(unsigned char * inframe, unsigned char * refframe, unsigned short * idx)
{
	// variables
	int ref,cur,by,bx,sy,sx,y,x;
	//unsigned short acc[Tsx];
	unsigned short acc[Tcur][Tby][Tbx][Tref][Tsy][Tsx];
//	unsigned short minimum[Tbx];
	unsigned short minimum[Tcur][Tby][Tbx];


	// buffers
	unsigned char  buf_in[Tcur][(Tby-1)*16+Ty][(Tbx-1)*16+Tx];
	unsigned char  buf_ref[Tref][(Tby-1)*16+Tsy+Ty][(Tbx-1)*16+Tsx+Tx];
	unsigned short vec[Tcur][Tby][Tbx];

	// control loops
	for(cur=0; cur<4; cur=cur+Tcur){//loop over four intermediate frames that must be encoded
		for(by=0; by<BY; by=by+Tby){//loop over the macro blocks in the frame
			for(bx=0; bx<BX; bx=bx+Tbx){

				// init minimum array
				init_mini:for(int tc=0; tc<min(Tcur,4-cur); tc++){
					for(int tby=0; tby<min(Tby,BY-by); tby++){
						for(int tbx=0; tbx<min(Tbx,BX-bx); tbx++){
							#pragma HLS PIPELINE
							minimum[tc][tby][tbx] = 65535; //init best match to the highest value of the data type
							vec[tc][tby][tbx] = 0;
						}
					}
				}

				for(ref=0; ref<2; ref=ref+Tref){//loop over the two reference frames
				  for(sy=0; sy<SEARCH_RANGE; sy=sy+Tsy){//search for the best match which is 32 positions in y and x
					for(sx=0; sx<SEARCH_RANGE; sx=sx+Tsx){


					  for(y=0; y<16; y=y+Ty){//test if SAD score of macroblock at this position
						for(x=0; x<16; x=x+Tx){

							int dc = min(Tcur, 4-cur);
							int dby = min(Tby, BY-by);
							int dbx = min(Tbx, BX-bx);
							int dr = min(Tref, 2-ref);
							int dsy = min(Tsy, SEARCH_RANGE-sy);
							int dsx = min(Tsx, SEARCH_RANGE-sx);
							int dy = min(Ty, 16-y);
							int dx = min(Tx, 16-x);

							// load inframe
							load_inframe:for (int tc=cur; tc<min(cur+dc,4); tc++){
								for (int row=by*16+y; row<min((by+dby-1)*16+y+dy,YIN); row++){
									for (int col=bx*16+x; col<min((bx+dbx-1)*16+x+dx,XIN); col++){
										#pragma HLS PIPELINE
										buf_in[tc-cur][row-by*16-y][col-bx*16-x] = *((unsigned char *)(inframe + tc*YIN*XIN+row*XIN+col));
									}
								}
							}

							// load refframe
							load_reframe:for (int tr=ref; tr<min(ref+dr,2); tr++){
								for (int row=by*16+sy+y; row<min((by+dby-1)*16+sy+dsy+y+dy,YREFIN); row++){
									for (int col=bx*16+sx+x; col<min((bx+dbx-1)*16+sx+dsx+x+dx,XREFIN); col++){
										#pragma HLS PIPELINE
										buf_ref[tr-ref][row-by*16-sy-y][col-bx*16-sx-x] = *((unsigned char *)(refframe + tr*YREFIN*XREFIN+row*XREFIN+col));
									}
								}
							}

							// tile begin
							tile_comp:for(int tc=0; tc<dc; tc++){
								for(int tby=0; tby<dby; tby++){
									for(int tbx=0; tbx<dbx; tbx++){

										for(int tr=0; tr<dr; tr++){
										  for(int tsy=0; tsy<dsy; tsy++){

											for(int ty=0; ty<dy; ty++){
												for(int tx=0; tx<dx; tx++){

													for(int tsx=sx; tsx<sx+dsx; tsx++){
														#pragma HLS PIPELINE
														#pragma HLS DEPENDENCE variable=acc array inter RAW false

														unsigned short tmp;
														if((y+ty==0) && (x+tx==0)) tmp = 0;
														else tmp = acc[tc][tby][tbx][tr][tsy][tsx-sx];

														acc[tc][tby][tbx][tr][tsy][tsx-sx] = tmp + abs(buf_in[tc][tby*16+ty][tbx*16+tx] - buf_ref[tr][tby*16+tsy+ty][tbx*16+tsx-sx+tx]);

														if((y+ty==15) && (x+tx==15)){
															if(acc[tc][tby][tbx][tr][tsy][tsx-sx]<minimum[tc][tby][tbx]){//check if current score is better than best match until now
																minimum[tc][tby][tbx]=acc[tc][tby][tbx][tr][tsy][tsx-sx];
																vec[tc][tby][tbx]=((ref+tr)<<15)|(((sy+tsy)<<8)|tsx);//encode reference (bit 15), encode y direction (bit 11:8), encode x direction (bit 3:0)
//																vec[tc][tby][tbx]=acc[tc][tby][tbx][tr][tsy][tsx-sx];
															}
														}

													}
												}

											}
										  }
										}
									}
								}
							}
							// tile end

						} // end x
					  } // end y

					} // end sx
				  } // end sy
				} // end ref

				// store idx
				store_idx:for(int tc=cur; tc<min(cur+Tcur,4); tc++){
					for(int tby=by; tby<min(by+Tby,BY); tby++){
						for(int tbx=bx; tbx<min(bx+Tbx,BX); tbx++){
							#pragma HLS PIPELINE
							*((unsigned short *)(idx+tc*BY*BX+tby*BX+tbx)) = vec[tc-cur][tby-by][tbx-bx]; //output the best motion vector
						}
					}
				}

			} // end bx
		} // end by
	} // end cur


	return 0;
}


