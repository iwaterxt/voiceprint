/*
 * queue.h
 *
 *  Created on: Jan 21, 2018
 *      Author: tao
 */

#ifndef QUEUE_H_
#define QUEUE_H_


# include <assert.h>
# include <malloc.h>
# include <stdio.h>
# include "base/common.h"

class Queue
{
public :
	inline Queue() : data_(NULL) , in_index_(0) , out_index_(0) , capacity_(0)
	{
		capacity_ = 2000000 ;
		data_ = (BaseFloat*)malloc(sizeof(BaseFloat) * capacity_) ;
		if(NULL == data_)
		{
			printf("error : could not allocate memory for queue ! \n") ;
			assert(NULL) ;
		}
	}

	inline ~Queue()
	{
		release() ;
	}

	/* initialize */
	inline void initialize()
	{
		in_index_ = 0 ;
		out_index_ = 0 ;
		capacity_ = 2000000 ;
		data_ = (BaseFloat*)malloc(sizeof(BaseFloat) * capacity_) ;
		if(NULL == data_)
		{
			printf("error : could not allocate memory for queue ! \n") ;
			assert(NULL) ;
		}
	}

	/* release */
	inline void release()
	{
		if(NULL != data_)
		{
			free(data_) ;
			data_ = NULL ;
		}
		in_index_ = 0 ;
		out_index_ = 0 ;
		capacity_ = 0 ;
	}

	inline bool empty() const
	{
		return in_index_ == out_index_ ;
	}

	inline bool full() const
	{
		return (in_index_ + 1) % capacity_ == out_index_ ;
	}

	inline void push(const BaseFloat& input)
	{
		if(full())
		{
			 printf("error : queue is full ! \n") ;
			 assert(NULL) ;
		}

		data_[in_index_] = input ;
		in_index_ = (in_index_ + 1) % capacity_ ;
	}

	inline BaseFloat pop(void)
	{
		if(empty())
		{
			 printf("error : queue is empty ! \n") ;
			 assert(NULL) ;
		}

		BaseFloat output = data_[out_index_] ;
		out_index_ = (out_index_ + 1) % capacity_ ;
		return output ;
	}

	inline size_t size() const
	{
		return in_index_ > out_index_ ? (in_index_ - out_index_) : (in_index_ - out_index_ + capacity_) ;
	}

	inline bool push_sequence_allowed(size_t length) const
	{
		if(empty())  return true ;

		size_t expected_destination = in_index_ + length + 1 ;
		bool allowed = true ;

		if(in_index_ > out_index_)
		{
			if(expected_destination >= capacity_ && expected_destination % capacity_ > out_index_)  allowed = false ;
		} else {
			if(expected_destination > out_index_)  allowed = false ;
		}

		return allowed ;
	}

	inline void push_sequence(BaseFloat* input , size_t length)
	{
	# if 0
		if(!push_sequence_allowed(length))
		{
			printf("error : no allowance for pushing sequence in queue ! \n") ;
			assert(NULL) ;
		}
	# endif

		for(size_t i = 0 ; i < length ; ++i)
		{
			data_[in_index_] = input[i] ;
			in_index_ = (in_index_ + 1) % capacity_ ;
		}
	}

	inline void push_sequence(int16* input , size_t length)
	{
	# if 0
		if(!push_sequence_allowed(length))
		{
			printf("error : no allowance for pushing sequence in queue ! \n") ;
			assert(NULL) ;
		}
	# endif

		for(size_t i = 0 ; i < length ; ++i)
		{
			data_[in_index_] = (BaseFloat)input[i] ;
			in_index_ = (in_index_ + 1) % capacity_ ;
		}
	}

	inline bool pop_sequence_allowed(size_t length) const
	{
		if(empty())  return false ;

		size_t expected_destination = out_index_ + length ;
		bool allowed = true ;

		if(in_index_ > out_index_)
		{
			if(expected_destination > in_index_)  allowed = false ;
		} else {
			if(length > capacity_ && (length % capacity_) > in_index_)  allowed = false ;
		}

		return allowed ;
	}

	inline void pop_sequence(BaseFloat* output , size_t length)
	{
	# if 0
		if(!pop_sequence_allowed(length))
		{
			printf("error : no allowance for popping sequence in queue ! \n") ;
			assert(NULL) ;
		}
	# endif

		for(size_t i = 0 ; i < length ; ++i)
		{
			output[i] = data_[out_index_] ;
			out_index_ = (out_index_ + 1) % capacity_ ;
		}
	}

private :
	BaseFloat* data_ ;
	size_t in_index_ ;
	size_t out_index_ ;
	size_t capacity_ ;

} ;



#endif /* QUEUE_H_ */
