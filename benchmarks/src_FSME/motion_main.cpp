/*
 * motion_main.cpp
 *
 *  Created on: 26 Aug 2017
 *      Author: jl12013
 */

#include <iostream>
#include <stdlib.h>
#include <stdint.h>

#include "sds_lib.h"
#include "motion_accel.h"
#include <string.h>

//#define NUM_TESTS 10

class perf_counter
{
public:
     uint64_t tot, cnt, calls;
     perf_counter() : tot(0), cnt(0), calls(0) {};
     inline void reset() { tot = cnt = calls = 0; }
     inline void start() { cnt = sds_clock_counter(); calls++; };
     inline void stop() { tot += (sds_clock_counter() - cnt); };
     inline uint64_t avg_cpu_cycles() { return (tot / calls); };
};

static void init_arrays(unsigned char *frame_in, unsigned char *reference, unsigned short *idx, unsigned char *frame_in_sw, unsigned char*reference_sw, unsigned short *idx_sw)
{
     for (int i = 0; i < 4 * YIN * XIN; i++) {
         frame_in[i] =  rand() % 256;
//         frame_in[i] =  255;
         frame_in_sw[i] = frame_in[i];
     }

     for (int i = 0; i < 2 * YREFIN * XREFIN; i++) {
    	 reference[i] =  rand() % 256;
//    	 reference[i] =  0;
    	 reference_sw[i] = reference[i];
     }

     for (int i = 0; i < 4 * BY * BX; i++) {
    	 idx[i] = 0;
         idx_sw[i] = 0;
     }
}

void motion_ref(unsigned char * inframe, unsigned char * refframe, unsigned short * idx)
{	/************************************************************************************
	* Function: void blockmatch(unsigned char inframe[], unsigned char refframe[], unsigned short idx[])
	* Input   : pointer to the sequence of four input frames that must be encoded, 
	            pointer to the two reference frames that are used for encoding 
    * Output  : array with the motion vectors each encoded in a 16 bit word
    * Procedure: perform perform block matching for eachframe on two reference frames of the sequence
    ************************************************************************************/
  int ref,cur,by,bx,sy,sx,y,x;
  unsigned short acc;
  unsigned short minimum;
  unsigned short vec;
  
  for(cur=0; cur<4; cur++){//loop over four intermediate frames that must be encoded
    for(by=0; by<BY; by++){//loop over the macro blocks in the frame
      for(bx=0; bx<BX; bx++){
        
		minimum=65535;//init best match to the highest value of the data type
		vec = 0;
	    
		for(ref=0; ref<2; ref++){//loop over the two reference frames
	      for(sy=0; sy<SEARCH_RANGE; sy++){//search for the best match which is 32 positions in y and x
            for(sx=0; sx<SEARCH_RANGE; sx++){
		      
			  acc=0;
              for(y=0; y<16; y++){//test if SAD score of macroblock at this position
			    for(x=0; x<16; x++){
			    	acc += abs(inframe[cur*YIN*XIN+(by*16+y)*XIN+bx*16+x] - refframe[ref*YREFIN*XREFIN+(by*16+sy+y)*XREFIN+bx*16+sx+x]);
			    }
              }
			  if(acc<minimum){//check if current score is better than best match until now
		        minimum=acc;
			    vec=(ref<<15)|((sy<<8)|sx);//encode reference (bit 15), encode y direction (bit 11:8), encode x direction (bit 3:0)
//		        vec=acc;
		      }
			}
          }
        }
		idx[cur*BY*BX+by*BX+bx]=vec; //output the best motion vector
      }
    }
  }
}

static int result_check(unsigned short *idx, unsigned short *idx_sw)
{
     for (int i = 0; i < 4 * BY * BX; i++) {
//	 for (int i = 0; i < 10; i++) {
          if (idx_sw[i] != idx[i]) {
               std::cout << "Mismatch: data index=" << i << ", d_sw=" << idx_sw[i]
                         << ", d_hw=" << idx[i] << std::endl;
//               return 1;
          }
     }
     return 0;
}


int motion_test(int num_tests, unsigned char *frame_in, unsigned char *reference, unsigned short *idx, unsigned char *frame_in_sw, unsigned char*reference_sw, unsigned short *idx_sw)
{
     std::cout << "Testing " << num_tests << " iterations of " << "Motion..." << std::endl;

     perf_counter hw_ctr, sw_ctr;

     for (int i = 0; i < num_tests; i++)
     {
    	  std::cout << "=== TEST: "<< i << std::endl;
          init_arrays(frame_in, reference, idx, frame_in_sw, reference_sw, idx_sw);

          std::cout << "Software running" << std::endl;
          sw_ctr.start();
          motion_ref(frame_in_sw, reference_sw, idx_sw);
          sw_ctr.stop();
          std::cout << "Software finished" << std::endl;

          std::cout << "Hardware running" << std::endl;
          hw_ctr.start();
//          motion_accel(frame_in, reference, idx);
          motion_accel_simp(frame_in, reference, idx);
          //motion_naive(frame_in, reference, idx);
          hw_ctr.stop();
          std::cout << "Hardware finished" << std::endl;

          std::cout << "Checking..." << std::endl;
          if (result_check(idx, idx_sw))
               return 1;
     }
     uint64_t sw_cycles = sw_ctr.avg_cpu_cycles();
     uint64_t hw_cycles = hw_ctr.avg_cpu_cycles();
     double speedup = (double) sw_cycles / (double) hw_cycles;

     std::cout << "Average number of CPU cycles running Motion in software: "
               << sw_cycles << std::endl;
     std::cout << "Average number of CPU cycles running Motion in hardware: "
               << hw_cycles << std::endl;
     std::cout << "Speed up: " << speedup << std::endl;

     return 0;
}


int main(int argc, char* argv[]){
     int test_passed = 0;
     unsigned char *frame_in, *reference, *frame_in_sw, *reference_sw;
	 unsigned short *idx, *idx_sw;

     // allocate memory in off-chip memory
     frame_in  = (unsigned char *)sds_alloc_non_cacheable(4 * YIN * XIN * sizeof(unsigned char));
     reference = (unsigned char *)sds_alloc_non_cacheable(2 * YREFIN * XREFIN * sizeof(unsigned char));
     idx       = (unsigned short *)sds_alloc_non_cacheable(4 * BY * BX * sizeof(unsigned short));

//     frame_in  = (unsigned char *)sds_alloc(4 * YIN * XIN * sizeof(unsigned char));
//     reference = (unsigned char *)sds_alloc(2 * YREFIN * XREFIN * sizeof(unsigned char));
//     idx       = (unsigned short *)sds_alloc(4 * BY * BX * sizeof(unsigned short));

     // allocate memory for SW
	 frame_in_sw  = (unsigned char *)sds_alloc(4 * YIN * XIN * sizeof(unsigned char));
	 reference_sw = (unsigned char *)sds_alloc(2 * YREFIN * XREFIN * sizeof(unsigned char));
     idx_sw       = (unsigned short *)sds_alloc(4 * BY * BX * sizeof(unsigned short));

     if (!frame_in || !reference || !idx || !frame_in_sw || !reference_sw || !idx_sw) {
          if (frame_in)  sds_free(frame_in);
          if (reference) sds_free(reference);
          if (idx)       sds_free(idx);
          if (frame_in_sw)  sds_free(frame_in_sw);
          if (reference_sw) sds_free(reference_sw);
          if (idx_sw)       sds_free(idx_sw);
          return 2;
     }

     // run test
     int num_tests = 1;
     for(int i=1; i<argc; i++){
    	 if(strcmp(argv[i],"-N") == 0){
    		 num_tests = atoi(argv[i+1]);
    	 }
     }
     test_passed = motion_test(num_tests, frame_in, reference, idx, frame_in_sw, reference_sw, idx_sw);

     std::cout << "TEST " << (test_passed ? "FAILED" : "PASSED") << std::endl;

     sds_free(frame_in);
     sds_free(reference);
     sds_free(idx);
     sds_free(frame_in_sw);
     sds_free(reference_sw);
     sds_free(idx_sw);

     return (test_passed ? -1 : 0);
}


