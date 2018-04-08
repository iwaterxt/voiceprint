/*
 * asr-plusplus.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef ASR_PLUS_PLUS_H_
#define ASR_PLUS_PLUS_H_
#include "common.h"

# ifdef __cplusplus
extern "C"{
# endif

struct vpconfig{
	int32 ivector_dim;
	int32 num_post;
	float min_post;
	char* final_dubm;
	char* final_ie;
	char* final_plda ;
	char* final_mean ;

};

struct mfcc_opt{

	int32 mfcc_dim;
	int32 samp_freq;
	int32 frame_shift_ms;
	int32 frame_length_ms;
};


void voice_print_lunch(const struct vpconfig* config);

void* voice_print_init();

void voice_print_start(void* engine);

void voice_print_sign(void* engine, char* key, int32 length);

void voice_print_check(void* engine, int32 length);

void voice_print_exit(void* engine);

void voice_print_put_data(void* engine, short* short_data , int32 length);

void voice_print_data_end(void* engine);


# ifdef __cplusplus
} ;
# endif

#endif /* ASR_PLUS_PLUS_H_ */
