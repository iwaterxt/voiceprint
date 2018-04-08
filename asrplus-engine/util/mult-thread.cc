/*
 * mult-thread.cc
 *
 *  Created on: Jan 21, 2018
 *      Author: tao
 */

# include "util/mult-thread.h"
# include <string.h>
# include <sys/time.h>
# include <stdio.h>

/* initialize */
void Mutex::initialize()
{
	if(0 != pthread_mutex_init(&mutex_ , NULL))
		error("error : can not initialize pthread mutex" , __FILE__ , __LINE__) ;
}

/* release */
void Mutex::release()
{
	if(0 != pthread_mutex_destroy(&mutex_))
		error("error : can not destroy pthread mutex" , __FILE__ , __LINE__) ;
}

void Mutex::Lock()
{
	int result = pthread_mutex_lock(&mutex_) ;
	if(0 != result)
	{
		const char* str = strerror(result) ;
		if(NULL == str)  str = "[NULL]" ;
		char err[100] ;
		sprintf(err , "error on locking pthread mutex : %s" , str) ;
		error(err , __FILE__ , __LINE__) ;
	}
}

void Mutex::Unlock()
{
	int result = pthread_mutex_unlock(&mutex_) ;
	if(0 != result)
	{
		const char* str = strerror(result) ;
		if(NULL == str)  str = "[NULL]" ;
		char err[100] ;
		sprintf(err , "error on unlocking pthread mutex : %s" , str) ;
		error(err , __FILE__ , __LINE__) ;
	}
}

bool Mutex::TryLock()
{
	int32 result = pthread_mutex_trylock(&mutex_) ;
	bool success = false ;
	switch(result)
	{
    	case 0 :
    	{
    		success = true ; break ;
    	}
    	case EBUSY :
    	{
    		success = false ; break ;
    	}
    	default :
    	{
    		const char* str = strerror(result) ;
    		if(NULL == str)  str = "[NULL]" ;
    		char err[100] ;
    		sprintf(err , "error on try-locking pthread mutex : %s" , str) ;
    		error(err , __FILE__ , __LINE__) ;
    	}
	}
	return success ;
}

/*******************************************************************************************************************/

/* initialize */
void ConditionMutex::init()
{
	Mutex::initialize() ;
	if(0 != pthread_cond_init(&cond_ , NULL))
		error("error : can not initialize pthread conditional variable" , __FILE__ , __LINE__) ;
}

/* free */
void ConditionMutex::destroy()
{
	Mutex::release() ;
	if(0 != pthread_cond_destroy(&cond_))
		error("error : can not destroy pthread conditional variable" , __FILE__ , __LINE__) ;
}

void ConditionMutex::WaitCond(unsigned int* timeout)
{
	int result ;

	if(NULL == timeout || 0 == *timeout)
	{
		result = pthread_cond_wait(&cond_ , &mutex_) ;
		if(0 != result)
		{
			const char* str = strerror(result) ;
			if(NULL == str)  str = "[NULL]" ;
			char err[100] ;
			sprintf(err , "error on waiting pthread cond : %s" , str) ;
			error(err , __FILE__ , __LINE__) ;
		}
	} else {
		struct timeval now ;
		struct timespec outtime ;
		gettimeofday(&now , NULL) ;
		outtime.tv_sec = now.tv_sec + (*timeout / 1000) ;
		outtime.tv_nsec = now.tv_usec * 1000 ;
		result = pthread_cond_timedwait(&cond_ , &mutex_ , &outtime) ;
		if(0 != result)
		{
			if(ETIMEDOUT == result)
				printf("Warning : timeout --- %s , line : %d \n" , __FILE__ , __LINE__) ;
			else {
				const char* str = strerror(result) ;
				if(NULL == str)  str = "[NULL]" ;
				char err[100] ;
				sprintf(err , "error on timed waiting pthread cond : %s" , str) ;
				error(err , __FILE__ , __LINE__) ;
			}
		}
	}
}

void ConditionMutex::SignalCond()
{
	int result = pthread_cond_signal(&cond_) ;
	if(0 != result)
	{
		const char* str = strerror(result) ;
		if(NULL == str)  str = "[NULL]" ;
		char err[100] ;
		sprintf(err , "error on signal pthread cond : %s" , str) ;
		error(err , __FILE__ , __LINE__) ;
	}
}


