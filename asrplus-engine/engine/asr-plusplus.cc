/*
 * asr-plusplus.cc
 *
 *  Created on: Sep 15, 2017
 *      Author: tao
 */

#include "engine/engine.h"
#include "engine/vpcontainer.h"
#include "engine/asr-plusplus.h"
#include <assert.h>

static VPcontainer *container_ = NULL;

void voice_print_lunch(const vpconfig* config){

	container_ = new VPcontainer ;
	if(container_ != 0){
		container_->Initialize(config);
		printf("VP container initialized successful!\n");
	}else{
		printf("VP container initialized failed!\n");
		exit(1);
	}
}

void* voice_print_init(){

	VPengine* engine = new VPengine ;

	assert(engine != 0);

	return engine ;
}

void voice_print_start(void* obj){

	VPengine* engine  = static_cast<VPengine*>(obj);

	assert(engine != 0);

	engine->initialize(container_);

}

void voice_print_sign(void* obj, char* key, int32 length){

	VPengine* engine  = static_cast<VPengine*>(obj);
	assert(engine != 0);
	std::string Key = key;
	engine->sign(Key, length);

}

void voice_print_check(void* obj, int32 length){

	VPengine* engine  = static_cast<VPengine*>(obj);
	assert(engine != 0);
	engine->check(length);

}

void voice_print_exit(void* obj){

	VPengine* engine  = static_cast<VPengine*>(obj);
	assert(engine != 0);
	engine->release();
}

void voice_print_put_data(void* obj, short* short_data , uint32 length){

	VPengine* engine  = static_cast<VPengine*>(obj);
	assert(engine != 0);
	engine->put_data(short_data, length);

}

void voice_print_data_end(void* obj){

	VPengine* engine  = static_cast<VPengine*>(obj);
	assert(engine != 0);
	engine->data_end();

}



