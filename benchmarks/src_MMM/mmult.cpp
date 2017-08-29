#include <iostream>
#include <stdlib.h>
#include <stdint.h>

#include "sds_lib.h"
#include "mmult_accel.h"
#include <string.h>

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

static void init_arrays(data_t *A, data_t *B, data_t *C, data_t *A_sw, data_t *B_sw, data_t *C_sw)
{

	 for (int j = 0; j < K; j++) {
		 for (int i = 0; i < M; i++) {
              A[j * M + i] = rand() % (10);
              A_sw[i * K + j] = A[j * M + i];
         }
    }

     for (int i = 0; i < K; i++) {
          for (int j = 0; j < N; j++) {
               B[i * N + j] = rand() % (10);
               B_sw[i * N + j] = B[i * N + j];
          }
     }

     for (int i = 0; i < M; i++) {
          for (int j = 0; j < N; j++) {
               C[i * N + j] = 0;
               C_sw[i * N + j] = 0;
          }
     }
}

void mmult_golden(data_t *A,  data_t *B, data_t *C)
{
     for (int row = 0; row < M; row++) {
          for (int col = 0; col < N; col++) {
        	  data_t result = 0;
              for (int k = 0; k < K; k++) {
            	  result += A[row*K+k] * B[k*N+col];
              }
              C[row*N+col] = result;
          }
     }
}

static int result_check(data_t *C, data_t *C_sw)
{
     for (int i = 0; i < M * N; i++) {
//	 for (int i = 0; i < 10; i++) {
          if (C_sw[i] != C[i]) {
               std::cout << "Mismatch: data index=" << i << ", d_sw=" << C_sw[i]
                         << ", d_hw=" << C[i] << std::endl;
               return 1;
          }
     }
     return 0;
}

int mmult_test(int num_tests, data_t *A,  data_t *B, data_t *C, data_t *A_sw,  data_t *B_sw, data_t *C_sw)
{
	std::cout << "Testing " << num_tests << " iterations of " << " MMM..." << std::endl;

     perf_counter hw_ctr, sw_ctr;
     
     for (int i = 0; i < num_tests; i++)
     {
    	  std::cout << "TEST: "<< i << std::endl;
          init_arrays(A, B, C, A_sw, B_sw, C_sw);

          std::cout << "Software running" << std::endl;
          sw_ctr.start();
          mmult_golden(A_sw, B_sw, C_sw);
          sw_ctr.stop();
          std::cout << "Software finished" << std::endl;

          std::cout << "Hardware running" << std::endl;
          hw_ctr.start();
		  mmult_accel(A, B, C);
          //mmult_naive(A, B, C);
          hw_ctr.stop();
          std::cout << "Hardware finished" << std::endl;

          std::cout << "Checking..." << std::endl;
          if (result_check(C, C_sw))
               return 1;
     }
     uint64_t sw_cycles = sw_ctr.avg_cpu_cycles();
     uint64_t hw_cycles = hw_ctr.avg_cpu_cycles();
     double speedup = (double) sw_cycles / (double) hw_cycles;

     std::cout << "Average number of CPU cycles running MMM in software: "
               << sw_cycles << std::endl;
     std::cout << "Average number of CPU cycles running MMM in hardware: "
               << hw_cycles << std::endl;
     std::cout << "Speed up: " << speedup << std::endl;

     return 0;
}

int main(int argc, char* argv[]){
     int test_passed = 0;
     data_t *A, *B, *C;
     data_t *A_sw, *B_sw, *C_sw;
     
     // allocate memory in off-chip memory
     A = (data_t *)sds_alloc_non_cacheable(M * K * sizeof(data_t));
     B = (data_t *)sds_alloc_non_cacheable(K * N * sizeof(data_t));
     C = (data_t *)sds_alloc_non_cacheable(M * N * sizeof(data_t));

     // allocate memory for SW
     A_sw = (data_t *)sds_alloc(M * K * sizeof(data_t));
     B_sw = (data_t *)sds_alloc(K * N * sizeof(data_t));
     C_sw = (data_t *)sds_alloc(M * N * sizeof(data_t));

     if (!A || !B || !C || !A_sw || !B_sw || !C_sw) {
          if (A) sds_free(A);
          if (B) sds_free(B);
          if (C) sds_free(C);
          if (A_sw) sds_free(A_sw);
          if (B_sw) sds_free(B_sw);
          if (C_sw) sds_free(C_sw);
          return 2;
     }
     
     // run test
     int num_tests = 1;
     for(int i=1; i<argc; i++){
    	 if(strcmp(argv[i],"-N") == 0){
    		 num_tests = atoi(argv[i+1]);
    	 }
     }

     test_passed = mmult_test(num_tests, A, B, C, A_sw, B_sw, C_sw);
     
     std::cout << "TEST " << (test_passed ? "FAILED" : "PASSED") << std::endl;

     sds_free(A);
     sds_free(B);
     sds_free(C);
     sds_free(A_sw);
     sds_free(B_sw);
     sds_free(C_sw);
     
     return (test_passed ? -1 : 0);
}


