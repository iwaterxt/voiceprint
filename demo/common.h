/*
 * common.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */



#ifndef COMMON_H_
#define COMMON_H_

#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <malloc.h>

#ifndef NULL
#define NULL   ((void *) 0)
#endif


typedef uint16_t uint16 ;
typedef uint32_t uint32 ;
typedef uint64_t uint64 ;
typedef int16_t int16 ;
typedef int32_t int32 ;
typedef int64_t int64 ;

typedef float BaseFloat;
typedef int int32 ;
typedef int MatrixIndexT ;

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865475244008443621048490
#endif

#ifndef M_LOG_2PI
#define M_LOG_2PI 1.8378770664093454835606594728112
#endif

#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559005
#endif


#ifndef COMPLEXFFT_BLOCKSIZE
#define COMPLEXFFT_BLOCKSIZE 8192
#endif

#define MEMALIGN(align, size, pp_orig) \
     (!posix_memalign(pp_orig, align, size) ? *(pp_orig) : NULL)
#define MEMALIGN_FREE(x) free(x)



typedef enum {kTrans, kNoTrans} MatrixTransposeType;

inline float Log(float x) { return logf(x); }

inline bool ISFINITE(float x) {return x != x+1.0;}

inline bool ISINF(float x) {return x == x + 1.0;}

inline bool ISNAN(float x) {return x == x;}

inline float Exp(float x) { return expf(x); }

# define SWAP8(a)  \
    {  \
		int t = ((char*)a)[0] ; ((char*)a)[0] = ((char*)a)[7] ; ((char*)a)[7] = t ;  \
			t = ((char*)a)[1] ; ((char*)a)[1] = ((char*)a)[6] ; ((char*)a)[6] = t ;  \
			t = ((char*)a)[2] ; ((char*)a)[2] = ((char*)a)[5] ; ((char*)a)[5] = t ;  \
			t = ((char*)a)[3] ; ((char*)a)[3] = ((char*)a)[4] ; ((char*)a)[4] = t ;  \
    }

# define SWAP4(a)  \
	{   \
		int t = ((char*)a)[0] ; ((char*)a)[0] = ((char*)a)[3] ; ((char*)a)[3] = t ;  \
			t = ((char*)a)[1] ; ((char*)a)[1] = ((char*)a)[2] ; ((char*)a)[2] = t ;  \
	}

# define SWAP2(a)  \
	{  \
		int t = ((char*)a)[0] ; ((char*)a)[0] = ((char*)a)[1] ; ((char*)a)[1] = t ;  \
	}


/* read wave , and store the the actual sound data */
int wav_read(const char* wave , BaseFloat** data , uint32* length) ;

void error(const char* error_message , const char* file , int32 line) ;

#endif /* COMMON_H_ */

