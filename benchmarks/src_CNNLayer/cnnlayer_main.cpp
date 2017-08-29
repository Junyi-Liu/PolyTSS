/*
 * cnnlayer_main.cpp
 *
 *  Created on: 3 Apr 2017
 *      Author: jl12013
 */

#include <iostream>
#include <stdlib.h>
#include <stdint.h>

#include "sds_lib.h"
#include "cnnlayer_accel.h"
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

static void init_arrays(data_t *frame_in,  data_t *weight, data_t *frame_out, data_t *frame_in_sw, data_t *weight_sw, data_t *frame_out_sw)
{
     for (int i = 0; i < N*inRow*inCol; i++) {
         frame_in[i] =  rand() % 256;
         frame_in_sw[i] = frame_in[i];
     }

     for (int i = 0; i < N*M*K*K; i++) {
    	 weight[i] =  rand() % 256;
    	 weight_sw[i] = weight[i];
     }

     for (int i = 0; i < M*outRow*outCol; i++) {
    	 frame_out[i] = 0.0;
         frame_out_sw[i] = 0.0;
     }
}

void cnnlayer_golden(data_t *frame_in,  data_t *weight, data_t *frame_out)
{

	for(int m=0; m<M; m++) {
		for(int row=0; row<outRow; row++) {
			for(int col=0; col<outCol; col++) {

				data_t rlt = 0;

				for(int n=0; n<N; n++) {
					for(int i=0; i<K; i++) {
						for(int j=0; j<K; j++) {
							rlt += weight[N*K*K*m + K*K*n + K*i + j] * frame_in[inRow*inCol*n + inCol*(ST*row+i) + (ST*col+j)];

				}}}

				frame_out[outRow*outCol*m + outCol*row + col] = rlt;
	}}}

}

static int result_check(data_t *frame_out, data_t *frame_out_sw)
{
     for (int i = 0; i < M * outRow * outCol; i++) {
          if (frame_out_sw[i] != frame_out[i]) {
               std::cout << "Mismatch: data index=" << i << ", d_sw=" << frame_out_sw[i]
                         << ", d_hw=" << frame_out[i] << std::endl;
               return 1;
          }
     }
     return 0;
}

// flush data cache
//void flush_data(char *c, int size){
//
//    const int size = 20*1024*1024;
//    char *c = (char *)malloc(size);
//
//	std::cout << "Flushing data cache ..."<< std::endl;
//    for(int i=0; i<2; i++)
//    	for(int j=0; j<size; j++)
//    		c[j] = i*j;
//    std::cout << "Finished: "<< c[0] << std::endl;
//
//    free(c);
//}

int cnnlayer_test(int num_tests, data_t *frame_in,  data_t *weight, data_t *frame_out, data_t *frame_in_sw, data_t *weight_sw, data_t *frame_out_sw)
{
     std::cout << "Testing " << num_tests << " iterations of " << " CNNLayer..." << std::endl;

     perf_counter hw_ctr, sw_ctr;

     for (int i = 0; i < num_tests; i++)
     {
    	  std::cout << "=== TEST: "<< i << std::endl;
          init_arrays(frame_in, weight, frame_out, frame_in_sw, weight_sw, frame_out_sw);

          std::cout << "Software running" << std::endl;
          sw_ctr.start();
          cnnlayer_golden(frame_in_sw, weight_sw, frame_out_sw);
          sw_ctr.stop();
          std::cout << "Software finished" << std::endl;

          std::cout << "Hardware running" << std::endl;
          hw_ctr.start();
          cnnlayer_naive(frame_in, weight, frame_out);
          hw_ctr.stop();
          std::cout << "Hardware finished" << std::endl;

          std::cout << "Checking..." << std::endl;
          if (result_check(frame_out, frame_out_sw))
               return 1;
     }
     uint64_t sw_cycles = sw_ctr.avg_cpu_cycles();
     uint64_t hw_cycles = hw_ctr.avg_cpu_cycles();
     double speedup = (double) sw_cycles / (double) hw_cycles;

     std::cout << "Average number of CPU cycles running CNNLayer in software: "
               << sw_cycles << std::endl;
     std::cout << "Average number of CPU cycles running CNNLayer in hardware: "
               << hw_cycles << std::endl;
     std::cout << "Speed up: " << speedup << std::endl;

     return 0;
}


int main(int argc, char* argv[]){
     int test_passed = 0;
     data_t *frame_in, *weight, *frame_out, *frame_in_sw, *weight_sw, *frame_out_sw;

     // allocate memory in off-chip memory
     frame_in  = (data_t *)sds_alloc_non_cacheable(N * inRow * inCol * sizeof(data_t));
     weight    = (data_t *)sds_alloc_non_cacheable(N * M * K * K * sizeof(data_t));
     frame_out = (data_t *)sds_alloc_non_cacheable(M * outRow * outCol * sizeof(data_t));

     // allocate memory for SW
	 frame_in_sw  = (data_t *)sds_alloc(N * inRow * inCol * sizeof(data_t));
	 weight_sw    = (data_t *)sds_alloc(N * M * K * K * sizeof(data_t));
     frame_out_sw = (data_t *)sds_alloc(M * outRow * outCol * sizeof(data_t));

     if (!frame_in || !weight || !frame_out || !frame_out_sw) {
          if (frame_in)  sds_free(frame_in);
          if (weight)    sds_free(weight);
          if (frame_out) sds_free(frame_out);
          if (frame_in_sw)  sds_free(frame_in_sw);
          if (weight_sw)    sds_free(weight_sw);
          if (frame_out_sw) sds_free(frame_out_sw);
          return 2;
     }

     // run test
     int num_tests = 1;
     for(int i=1; i<argc; i++){
    	 if(strcmp(argv[i],"-N") == 0){
    		 num_tests = atoi(argv[i+1]);
    	 }
     }
     test_passed = cnnlayer_test(num_tests, frame_in, weight, frame_out, frame_in_sw, weight_sw, frame_out_sw);

     std::cout << "TEST " << (test_passed ? "FAILED" : "PASSED") << std::endl;

     sds_free(frame_in);
     sds_free(weight);
     sds_free(frame_out);
     sds_free(frame_in_sw);
     sds_free(weight_sw);
     sds_free(frame_out_sw);

     return (test_passed ? -1 : 0);
}


