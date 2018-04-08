/*
 * voiceprint.cc
 *
 *  Created on: Oct 10, 2017
 *      Author: tao
 */

#include "asr-plusplus.h"
#include "common.h"


# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include <time.h>
# include <string.h>
# include <sys/sysinfo.h>




void copy_str(char** dst , const char* src)
{
	size_t len = strlen(src) ;
	*dst = (char*)malloc(sizeof(char) * (len + 1)) ;
	assert(NULL != *dst) ;
	memcpy((void*)(*dst) ,(const void*)src , sizeof(char) * len) ;
	(*dst)[len] ='\0' ;
}

void free_str(char** str)
{
	if(NULL != *str)
	{
		free(*str) ;
		*str = NULL ;
	}
}

int main(){


	printf("........voice print engine start........\n");

	struct vpconfig config;

	FILE* fp;

	if((fp=fopen("./vp.config", "r")) == NULL){

		printf("error: vp.config can not be opened!\n");
		exit(1);
	}

	if(!feof(fp)){

		char contain[256];

		fscanf(fp, "%d", &config.ivector_dim);
		fscanf(fp, "%d", &config.num_post);
		fscanf(fp, "%f", &config.min_post);
		fscanf(fp, "%s", contain);
		copy_str((char**)&config.final_dubm, contain);
		fscanf(fp, "%s", contain);
		copy_str((char**)&config.final_ie, contain);
		fscanf(fp, "%s", contain);
		copy_str((char**)&config.final_plda, contain);
		fscanf(fp, "%s", contain);
		copy_str((char**)&config.final_mean, contain);

	}

	fclose(fp);

	voice_print_lunch(&config);

	void* engine = voice_print_init();

	FILE* fw = fopen("./wav_sign_old.scp", "r");
	if(fw == NULL){
		printf("error: ./wav.scp can not be opened! \n");
	}

	char key[256], wav_path[512];
	if(!feof(fw)){
		fscanf(fw, "%s", key);
		fscanf(fw, "%s", wav_path);
	}
	BaseFloat* data;
	uint32 length;
	short* short_data;
	printf("........voice print engine start sign........\n");
	while(!feof(fw)){

		voice_print_start(engine);
		/* the first step : get ready of data */
		printf("%s \n" , wav_path) ;
		wav_read(wav_path , &data , &length) ;

		short_data = (short*)malloc(sizeof(short) * length) ;

		assert(NULL != short_data) ;

		for( uint32 j = 0 ; j < length ; j++)  short_data[j] = (short)(data[j]) ;

		/* sned data into vp engine */
		voice_print_put_data(engine , short_data , length);

		/* when data is ready , tell vp engine that the process of sending data is finished */
		voice_print_data_end(engine) ;

		free(data) ;  data = NULL ;
		free(short_data) ;  short_data = NULL ;

		voice_print_sign(engine, key, length);

		/* the next utterance */
		fscanf(fw , "%s" , key) ;
		fscanf(fw , "%s" , wav_path) ;

	}

	fclose(fw) ;

	FILE* ft = fopen("./wav_check_old.scp", "r");
	if(ft == NULL){
		printf("error: ./wav_test.scp can not be opened!\n");
	}
	if(!feof(ft)){
		fscanf(ft, "%s ", key);

		fscanf(ft, "%s", wav_path);
	}

	printf("........voice print engine start check........\n");
	printf("%s ", key);
	while(!feof(ft)){


		voice_print_start(engine);
		/* the first step : get ready of data */
		//printf("%s " , wav_path) ;
		wav_read(wav_path , &data , &length) ;

		short_data = (short*)malloc(sizeof(short) * length) ;

		assert(NULL != short_data) ;

		for(uint32 j = 0 ; j < length ; j++)  short_data[j] = (short)(data[j]) ;
		/* sned data into vp engine */
		voice_print_put_data(engine , short_data , length);

		/* when data is ready , tell vp engine that the process of sending data is finished */
		voice_print_data_end(engine) ;

		free(data) ;  data = NULL ;
		free(short_data) ;  short_data = NULL ;

		voice_print_check(engine, length);

		/* the next utterance */
		fscanf(ft , "%s" , key) ;
		fscanf(ft , "%s" , wav_path) ;
		if(!feof(ft))
			printf("%s ", key);
	}

	voice_print_exit(engine);

	fclose(ft) ;
	printf("........voice print engine end check........\n");
	free_str(&config.final_dubm);
	free_str(&config.final_ie);
	free_str(&config.final_plda);
	free_str(&config.final_mean);

	return 0;
}


