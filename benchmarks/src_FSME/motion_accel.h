/*
 * motion_accel.h
 *
 *  Created on: 26 Aug 2017
 *      Author: jl12013
 */

#include <stdint.h>

#ifndef SRC_MOTION_ACCEL_H_
#define SRC_MOTION_ACCEL_H_

/********************************************
 **	Loop Bounds
 ********************************************/
#define YIN 720 //frame size of frames that must be encoded
#define XIN 1280

#define YREFIN 752 //reference frame size
#define XREFIN 1312

#define BY 45 //number of 16x16 macro blocks in the frame
#define BX 80

//Search offsets within 16 pixels of (0,0)
#define SEARCH_RANGE 32

/********************************************
 ** tile size
 ********************************************/
#define Tcur 4
#define Tby  2
#define Tbx  2
#define Tref 1
#define Tsy  32
#define Tsx  32
#define Ty   16
#define Tx   16

// Interface setting
#pragma SDS data mem_attribute(inframe:NON_CACHEABLE, refframe:NON_CACHEABLE, idx:NON_CACHEABLE)
#pragma SDS data sys_port(inframe:AFI, refframe:AFI, idx:AFI)
#pragma SDS data zero_copy(inframe, refframe, idx)

// accelerators
//int motion_accel(unsigned char * inframe, unsigned char * refframe, unsigned short * idx);
int motion_accel_simp(unsigned char * inframe, unsigned char * refframe, unsigned short * idx);
//int motion_naive(unsigned char * inframe, unsigned char * refframe, unsigned short * idx);

//inline unsigned char my_abs(unsigned char a, unsigned char b) { return (a<b) ? (b-a) : (a-b);}

#endif /* SRC_MOTION_ACCEL_H_ */
