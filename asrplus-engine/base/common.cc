/*
 * common.cc
 *
 *  Created on: Feb 13, 2018
 *      Author: tao
 */
#include "common.h"
#include <stdio.h>
#include <mutex>
#include <stdlib.h>

int Rand(struct RandomState* state) {

  if (state) {
    return rand_r(&(state->seed));
  } else {
    return rand();
  }

}

RandomState::RandomState() {

  seed = Rand() + 27437;
}


void error(const char* error_message , const char* file , int32 line) {
	printf("%s , %s : line %u \n" , error_message , file , line) ;
	assert(NULL) ;
}

std::string CharToString(const char &c) {
  char buf[20];
  if (std::isprint(c))
    snprintf(buf, sizeof(buf), "\'%c\'", c);
  else
    snprintf(buf, sizeof(buf), "[character %d]", static_cast<int>(c));
  return (std::string) buf;
}
