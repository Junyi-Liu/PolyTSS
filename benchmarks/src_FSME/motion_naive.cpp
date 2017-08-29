/*
 * motion_naive.cpp
 *
 *  Created on: 26 Aug 2017
 *      Author: jl12013
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "motion_accel.h"

// accelerator
int motion_naive(unsigned char * inframe, unsigned char * refframe, unsigned short * idx)
{
	// variables
	int ref,cur,by,bx,sy,sx,y,x;
	unsigned short acc;
    unsigned short minimum;

	// buffers
	unsigned char  buf_in[4][YIN][XIN];
	unsigned char  buf_ref[2][YREFIN][XREFIN];
	unsigned short vec[4][BY][BX];
	
	// load inframe
	load_inframe:for (cur=0; cur<4; cur++){
			for (int row=0; row<YIN; row++){
				for (int col=0; col<XIN; col++){
					#pragma HLS PIPELINE
					buf_in[cur][row][col] = *((unsigned char *)(inframe + cur*YIN*XIN+row*XIN+col));
				}
			}
		}
	
	
	// load refframe
	load_reframe:for (ref=0; ref<2; ref++){
			for (int row=0; row<YREFIN; row++){
				for (int col=0; col<XREFIN; col++){
					#pragma HLS PIPELINE
					buf_ref[ref][row][col] = *((unsigned char *)(refframe + ref*YREFIN*XREFIN+row*XREFIN+col));
				}
			}
		}
	
	// main loop
	for(cur=0; cur<4; cur++){//loop over four intermediate frames that must be encoded
		for(by=0; by<BY; by++){//loop over the macro blocks in the frame
		  for(bx=0; bx<BX; bx++){
			
			minimum=65535;//init best match to the highest value of the data type
			vec[cur][by][bx] = 0;
			
			for(ref=0; ref<2; ref++){//loop over the two reference frames
			  for(sy=0; sy<SEARCH_RANGE; sy++){//search for the best match which is 32 positions in y and x
				for(sx=0; sx<SEARCH_RANGE; sx++){
				  
				  acc=0;
				  for(y=0; y<16; y++){//test if SAD score of macroblock at this position
					for(x=0; x<16; x++){
						acc += abs(buf_in[cur][by*16+y][bx*16+x] - buf_ref[ref][by*16+sy+y][bx*16+sx+x]);
					}
				  }
				  if(acc<minimum){//check if current score is better than best match until now
					minimum=acc;
					vec[cur][by][bx]=(ref<<15)|((sy<<8)|sx);//encode reference (bit 15), encode y direction (bit 11:8), encode x direction (bit 3:0)
				  }
				}
			  }
			}

		  }
		}
	}
	
	// store idx
	store_idx:for(cur=0; cur<4; cur++){
		for(int tby=0; tby<BY; tby++){
			for(int tbx=0; tbx<BX; tbx++){
				#pragma HLS PIPELINE
				*((unsigned short *)(idx+cur*BY*BX+tby*BX+tbx)) = vec[cur][tby][tbx]; //output the best motion vector
			}
		}
	}

	return 0;
}


