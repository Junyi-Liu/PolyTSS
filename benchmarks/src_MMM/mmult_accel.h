#include <stdint.h>

#ifndef MMULT_ACC_H_
#define MMULT_ACC_H_

// data type
typedef uint16_t data_t;
//typedef int16_t data_t;


// loop bounds
#define M 500
#define N 300
#define K 400

// tile parameters
#define Tm 27
#define Tn 26
#define Tk 6

// Interface setting
#pragma SDS data mem_attribute(A:NON_CACHEABLE, B:NON_CACHEABLE, C:NON_CACHEABLE)
//#pragma SDS data mem_attribute(A:CACHEABLE, B:CACHEABLE, C:CACHEABLE)
#pragma SDS data sys_port(A:AFI, B:AFI, C:AFI)
#pragma SDS data zero_copy(A, B, C)
int mmult_accel (data_t * A, data_t * B, data_t * C);
//int mmult_naive(data_t * A, data_t * B, data_t * C);

#endif /* MMULT_ACC_H_ */

