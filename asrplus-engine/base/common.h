/*
 * common.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */



#ifndef COMMON_H_
#define COMMON_H_
#include <cmath>
#include <cfloat>
#include <string>
#include <assert.h>
#include <stdint.h>
#include <cblas.h>

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

#ifndef RAND_MAX  /* rather dangerous non-ansi workaround */
   #define RAND_MAX ((unsigned long)(1<<30))
#endif

#define MEMALIGN(align, size, pp_orig) \
     (!posix_memalign(pp_orig, align, size) ? *(pp_orig) : NULL)
#define MEMALIGN_FREE(x) free(x)

int Rand(struct RandomState* state = NULL);

struct RandomState {
  RandomState();
  unsigned seed;
};


typedef enum {
  kTrans    = CblasTrans,
  kNoTrans = CblasNoTrans
} MatrixTransposeType;

inline float Log(float x) { return logf(x); }

inline bool ISFINITE(float x) {return x != x+1.0;}

inline bool ISINF(float x) {return x == x + 1.0;}

inline bool ISNAN(float x) {return x == x;}

inline float Exp(float x) { return expf(x); }




/*  this function is for error information  */
void error(const char* error_message , const char* file , int32 line);

std::string CharToString(const char &c);

class OptionsItf {
 public:

  virtual void Register(const std::string &name,
                bool *ptr, const std::string &doc) = 0;
  virtual void Register(const std::string &name,
                int32 *ptr, const std::string &doc) = 0;
  virtual void Register(const std::string &name,
                float *ptr, const std::string &doc) = 0;
  virtual void Register(const std::string &name,
                double *ptr, const std::string &doc) = 0;
  virtual void Register(const std::string &name,
                std::string *ptr, const std::string &doc) = 0;

  virtual ~OptionsItf() {}
};


/* read wave , and store the the actual sound data */
int wav_read(const char* wave , BaseFloat** data , uint32* length) ;


#endif /* COMMON_H_ */
