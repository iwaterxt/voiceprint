/*
 * data-interface.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef DATA_INTERFACE_H_
#define DATA_INTERFACE_H_

# include "util/mult-thread.h"
# include "util/queue.h"
# include "base/common.h"

class DataInterface
{
public :
	DataInterface() ;
	~DataInterface() ;

	/* initialize */
	void initialize() ;

	/* release */
	void release() ;

	int32 input_data(int16* input , int32 acquire_length);

	int32 read(BaseFloat* output , int32 request_length) ;

  	int32 data_finished(void) ;

private :
  	Queue* queue_ ;

    ConditionMutex cond_ ;

    bool is_data_finished_ ;

    bool is_wait_producer_signal_ ;
    bool is_wait_consumer_signal_ ;

    int32 request_length_ ;
    int32 acquire_length_ ;
} ;



#endif /* DATA_INTERFACE_H_ */
