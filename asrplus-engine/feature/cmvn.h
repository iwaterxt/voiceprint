/*
 * cmvn.h
 *
 *  Created on: Oct 11, 2017
 *      Author: tao
 */

#ifndef CMVN_H_
#define CMVN_H_

#include "base/common.h"

/* CMN : Cepstral Mean Normalization */
class CMN {
public :
	CMN() : row_(0) , dim_(0) , window_size_(0) , num_cached_feature_(0) , norm_(NULL) , stats_(NULL) {}
	CMN(uint32 dim , uint32 window_size) ;
	~CMN() ;

	/* initialize */
	void initialize(uint32 dim , uint32 window_size) ;

	/* release */
	void release() ;

	/* cepstral mean normalization (variance normalization not supported yet) */
	void apply_cmn(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col) ;

	void apply_cmn_offline(BaseFloat* input_feature, const uint32 num_row, const uint32 num_col);

private :
	void local_normalize(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col) ;

	uint32 row_ ;
	uint32 dim_ ;
	uint32 window_size_ ;
	uint32 num_cached_feature_ ;
	double* norm_ ;
	BaseFloat* stats_ ;

} ;

/* CMVN : Cepstral Mean and Variance Normalization */
class CMVN {
public :
	CMVN() : row_(0) , dim_(0) , window_size_(0) , num_cached_feature_(0) , norm_(NULL) , squared_mean_(NULL) , var_(NULL) , stats_(NULL) {}
	CMVN(uint32 dim , uint32 window_size) ;
	~CMVN() ;

	/* initialize */
	void initialize(uint32 dim , uint32 window_size) ;

	/* release */
	void release() ;

	/* cepstral mean and variance normalization */
	void apply_cmvn(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col) ;


private :
	void local_normalize(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col) ;

	uint32 row_ ;
	uint32 dim_ ;
	uint32 window_size_ ;
	uint32 num_cached_feature_ ;
	double* norm_ ;
	double* squared_mean_ ;
	double* var_ ;
	BaseFloat* stats_ ;

} ;



#endif /* CMVN_H_ */
