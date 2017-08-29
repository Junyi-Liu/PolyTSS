/*
 * cnnlayer_accel.h
 *
 *  Created on: 3 Apr 2017
 *      Author: jl12013
 */

#include <stdint.h>

#ifndef SRC_CNNLAYER_ACCEL_H_
#define SRC_CNNLAYER_ACCEL_H_

// data type
typedef uint16_t data_t;

/********************************************
 **	Layer 1
 ********************************************/
//// input frame
//#define N 3
//#define inRow 227
//#define inCol 227
//// output frame
//#define M 48
//#define outRow 55
//#define outCol 55
//// convolution parameters
//#define K 11 // kernel window size
//#define ST 4 //stride size

/********************************************
 **	Layer 3
 ********************************************/
// input frame
#define N 256
#define inRow 16
#define inCol 16
// output frame
#define M 192
#define outRow 13
#define outCol 13
// convolution parameters
#define K 3 // kernel window size
#define ST 1 //stride size

/********************************************
 ** tile size
 ********************************************/
#define Tm 16
#define Tr 6
#define Tc 6
#define Tn 3
#define Ti 3
#define Tj 3
//#define Tm 16
//#define Tr 40
//#define Tc 46
//#define Tn 3
//#define Ti 11
//#define Tj 11


// Interface setting
//#pragma SDS data data_mover(frame_in:AXIDMA_SIMPLE, weight:AXIDMA_SIMPLE, frame_out:AXIDMA_SIMPLE)
#pragma SDS data mem_attribute(frame_in:NON_CACHEABLE, weight:NON_CACHEABLE, frame_out:NON_CACHEABLE)
//#pragma SDS data mem_attribute(frame_in:CACHEABLE, weight:CACHEABLE, frame_out:CACHEABLE)
#pragma SDS data sys_port(frame_in:AFI, weight:AFI, frame_out:AFI)
#pragma SDS data zero_copy(frame_in, weight, frame_out)

// accelerators
//int cnnlayer_accel(data_t *frame_in, data_t *weight, data_t *frame_out);
//int cnnlayer_accel_mmp1(data_t *frame_in, data_t *weight, data_t *frame_out);
//int cnnlayer_accel_mmp2(data_t *frame_in, data_t *weight, data_t *frame_out);
//int cnnlayer_accel_mmp3(data_t *frame_in, data_t *weight, data_t *frame_out);
int cnnlayer_naive(data_t *frame_in, data_t *weight, data_t *frame_out);

#endif /* SRC_CNNLAYER_ACCEL_H_ */
