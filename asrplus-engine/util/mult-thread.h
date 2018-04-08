/*
 * mult-thread.h
 *
 *  Created on: Jan 21, 2018
 *      Author: tao
 */

#ifndef MULT_THREAD_H_
#define MULT_THREAD_H_
# include "base/common.h"
# include <pthread.h>
# include <errno.h>

/* add "-lpthread" as link option */

class Mutex {
public :
	Mutex()  {}
	~Mutex()  {}

	/* initialize */
	void initialize() ;

	/* release */
	void release() ;

	void Lock() ;
	void Unlock() ;

	/*
	 * try to lock the mutex without waiting for it
	 * return : true when locking successfully , false when mutex was already locked
	 */
	bool TryLock() ;

protected :
	pthread_mutex_t mutex_ ;

} ;

/*******************************************************************************************************************/

class ConditionMutex : public Mutex
{
public :
	ConditionMutex()  {}
	~ConditionMutex()  {}

	/* initialize */
	void init() ;

	/* destroy */
	void destroy() ;

	void WaitCond(unsigned int* timeout) ;
	void SignalCond() ;

private :
	pthread_cond_t cond_ ;

} ;



#endif /* MULT_THREAD_H_ */
