/*
 * cmvn.cc
 *
 *  Created on: Oct 11, 2017
 *      Author: tao
 */

#include "base/common.h"
#include "feature/cmvn.h"
#include <assert.h>
#include <malloc.h>
#include <string.h>


static const double var_floor = 1e-20 ;

CMN::CMN(uint32 dim , uint32 window_size) : row_(0) , dim_(0) , window_size_(0) , num_cached_feature_(0) , norm_(NULL) , stats_(NULL)
{
	initialize(dim , window_size) ;
}

CMN::~CMN()
{
	release() ;
}

/* initialize */
void CMN::initialize(uint32 dim , uint32 window_size)
{
	row_ = 0 ;
	dim_ = dim ;
	window_size_ = window_size ;
	num_cached_feature_ = 0 ;

	norm_ = NULL ;

	stats_ = (BaseFloat*)malloc(window_size_ * dim_ * sizeof(BaseFloat)) ;
	if(NULL == stats_)  error("no memory to be allocated" , __FILE__ , __LINE__) ;
	memset((void*)stats_ , 0 , window_size_ * dim_ * sizeof(BaseFloat)) ;
}

/* release */
void CMN::release()
{
	if(NULL != norm_)
	{
		free(norm_) ;
		norm_ = NULL ;
	}

	if(NULL != stats_)
	{
		free(stats_) ;
		stats_ = NULL ;
	}

	row_ = 0 ;
	dim_ = 0;
	window_size_ = 0 ;
	num_cached_feature_ = 0 ;
}

void CMN::apply_cmn(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col)
{
	assert(dim_ == num_col) ;

	/* first , fill up the history cache */
	uint32 i = 0 ;
	while(num_cached_feature_ < window_size_)
	{
		if(i < num_row)
		{
			memcpy((void*)(stats_ + num_cached_feature_ * dim_) , (const void*)(input_feature + i * dim_) , sizeof(BaseFloat) * dim_) ;
			num_cached_feature_++ ;
			i++ ;
		}else break ;
	}

	/* if the history cache not filled up yet */
	if(num_cached_feature_ !=  window_size_)
	{
		/* just compute the mean based on what have been filled so far */
		local_normalize(input_feature , num_row , num_col) ;
		return ;
	}

	/* now that having the history cache of features filled up , need to compute the full sum only once */
	if(NULL == norm_)
	{
		double* sum = (double*)malloc(dim_ * sizeof(double)) ;
		norm_ = (double*)malloc(dim_ * sizeof(double)) ;
		if(NULL == sum || NULL == norm_)  error("no memory to be allocated" , __FILE__ , __LINE__) ;

		memset((void*)sum , 0 , dim_ * sizeof(double)) ;

		for(uint32 r = 0 ; r < window_size_ ; r++)
			for(uint32 c = 0 ; c < dim_ ; c++)  sum[c] += (double)stats_[r * dim_ + c] ;

		for(uint32 c = 0 ; c < dim_ ; c++)  norm_[c] = sum[c] / num_cached_feature_ ;

		free(sum) ; sum = NULL ;
	}

	/* compute the rolling mean */
	for(; i < num_row ; ++i)
	{
		/* update rolling mean */
		for(uint32 c = 0 ; c < dim_ ; c++)  norm_[c] -= (stats_[row_ * dim_ + c] - input_feature[i * dim_ + c]) / window_size_ ;

		/* replace the oldest row with the current feature vector */
		memcpy((void*)(stats_ + row_ * dim_) , (const void*)(input_feature + i * dim_) , sizeof(BaseFloat) * dim_) ;

		/* update the oldest row_ variable */
		row_ = (row_ + 1) % window_size_ ;
	}

	/* compute the mean normalized feature */
	for(uint32 r = 0 ; r < num_row ; r++)
		for(uint32 c = 0 ; c < num_col ; c++)  input_feature[r * num_col + c] = input_feature[r * num_col + c] - (BaseFloat)norm_[c] ;
}

void CMN::apply_cmn_offline(BaseFloat* input_feature, const uint32 num_row, const uint32 num_col){

	double* mean = (double*)malloc(sizeof(double)*num_col);
	memset(mean, 0, sizeof(double)*num_col);
	for(uint32 i = 0; i < num_row; i++){
		for(uint32 j = 0; j < num_col; j++){
			mean[j] += input_feature[i * num_col + j];
		}
	}

	for(uint32 k = 0; k < num_col; k++)
		mean[k] = mean[k]/num_row ;

	for(uint32 i = 0; i < num_row; i++){
		for(uint32 j = 0; j < num_col; j++){
			input_feature[i * num_col + j] = input_feature[i * num_col + j] - mean[j];
		}
	}

}

void CMN::local_normalize(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col)
{
	assert(dim_ == num_col) ;

	double* sum = (double*)malloc(dim_ * sizeof(double)) ;
	double* norm = (double*)malloc(dim_ * sizeof(double)) ;
	if(NULL == sum || NULL == norm)  error("no memory to be allocated" , __FILE__ , __LINE__) ;

	memset((void*)sum , 0 , dim_ * sizeof(double)) ;

	for(uint32 r = 0 ; r < num_cached_feature_ ; r++)
		for(uint32 c = 0 ; c < dim_ ; c++)  sum[c] += (double)stats_[r * dim_ + c] ;

	for(uint32 c = 0 ; c < dim_ ; c++)  norm[c] = sum[c] / num_cached_feature_ ;

	/* compute the mean normalized feature */
	for(uint32 r = 0 ; r < num_row ; r++)
		for(uint32 c = 0 ; c < num_col ; c++)  input_feature[r * num_col + c] = input_feature[r * num_col + c] - (BaseFloat)norm[c] ;

	free(sum) ; sum = NULL ;
	free(norm) ; norm = NULL ;
}

CMVN::CMVN(uint32 dim , uint32 window_size) : row_(0) , dim_(0) , window_size_(0) , num_cached_feature_(0) , norm_(NULL) , squared_mean_(NULL) , var_(NULL) , stats_(NULL)
{
	initialize(dim , window_size) ;
}

CMVN::~CMVN()
{
	release() ;
}

/* initialize */
void CMVN::initialize(uint32 dim , uint32 window_size)
{
	row_ = 0 ;
	dim_ = dim ;
	window_size_ = window_size ;
	num_cached_feature_ = 0 ;

	norm_ = NULL ;
	squared_mean_ = NULL ;
	var_ = NULL ;

	stats_ = (BaseFloat*)malloc(window_size_ * dim_ * sizeof(BaseFloat)) ;
	if(NULL == stats_)  error("no memory to be allocated" , __FILE__ , __LINE__) ;
	memset((void*)stats_ , 0 , window_size_ * dim_ * sizeof(BaseFloat)) ;
}

/* release */
void CMVN::release()
{
	if(NULL != norm_)
	{
		free(norm_) ;
		norm_ = NULL ;
	}

	if(NULL != squared_mean_)
	{
		free(squared_mean_) ;
		squared_mean_ = NULL ;
	}

	if(NULL != var_)
	{
		free(var_) ;
		var_ = NULL ;
	}

	if(NULL != stats_)
	{
		free(stats_) ;
		stats_ = NULL ;
	}

	row_ = 0 ;
	dim_ = 0;
	window_size_ = 0 ;
	num_cached_feature_ = 0 ;
}

void CMVN::apply_cmvn(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col)
{
	assert(dim_ == num_col) ;

	/* first , fill up the history cache */
	uint32 i = 0 ;
	while(num_cached_feature_ < window_size_)
	{
		if(i < num_row)
		{
			memcpy((void*)(stats_ + num_cached_feature_ * dim_) , (const void*)(input_feature + i * dim_) , sizeof(BaseFloat) * dim_) ;
			num_cached_feature_++ ;
			i++ ;
		}else break ;
	}

	/* if the history cache not filled up yet */
	if(num_cached_feature_ !=  window_size_)
	{
		/* just compute the mean and variance based on what have been filled so far */
		local_normalize(input_feature , num_row , num_col) ;
		return ;
	}

	/* now that having the history cache of features filled up , need to compute the full sum or squared sum only once */
	if(NULL == norm_ || NULL == squared_mean_ || NULL == var_)
	{
		double* sum = (double*)malloc(dim_ * sizeof(double)) ;
		double* squared_sum = (double*)malloc(dim_ * sizeof(double)) ;
		norm_ = (double*)malloc(dim_ * sizeof(double)) ;
		squared_mean_ = (double*)malloc(dim_ * sizeof(double)) ;
		var_ = (double*)malloc(dim_ * sizeof(double)) ;
		if(NULL == sum || NULL == squared_sum || NULL == norm_ || squared_mean_ == NULL || var_ == NULL)  error("no memory to be allocated" , __FILE__ , __LINE__) ;

		memset((void*)sum , 0 , dim_ * sizeof(double)) ;
		memset((void*)squared_sum , 0 , dim_ * sizeof(double)) ;

		for(uint32 r = 0 ; r < window_size_ ; r++)
			for(uint32 c = 0 ; c < dim_ ; c++)
			{
				sum[c] += (double)stats_[r * dim_ + c] ;
				squared_sum[c] += pow(stats_[r * dim_ + c] , 2) ;
			}

		for(uint32 c = 0 ; c < dim_ ; c++)
		{
			norm_[c] = sum[c] / num_cached_feature_ ;
			squared_mean_[c] = squared_sum[c] / num_cached_feature_ ;
		}

		free(sum) ; sum = NULL ;
		free(squared_sum) ; squared_sum = NULL ;
	}

	/* compute the rolling (squared) mean */
	for(; i < num_row ; ++i)
	{
		/* update rolling (squared) mean */
		for(uint32 c = 0 ; c < dim_ ; c++)
		{
			norm_[c] -= (stats_[row_ * dim_ + c] - input_feature[i * dim_ + c]) / window_size_ ;
			squared_mean_[c] -= (pow(stats_[row_ * dim_ + c] , 2) - pow(input_feature[i * dim_ + c] , 2)) / window_size_ ;
		}

		/* replace the oldest row with the current feature vector */
		memcpy((void*)(stats_ + row_ * dim_) , (const void*)(input_feature + i * dim_) , sizeof(BaseFloat) * dim_) ;

		/* update the oldest row_ variable */
		row_ = (row_ + 1) % window_size_ ;
	}

	for(uint32 c = 0 ; c < dim_ ; c++)
	{
		var_[c] = squared_mean_[c] - pow(norm_[c] , 2) ;

		if(var_[c] < var_floor) var_[c] = var_floor ;

		/* switch to scale and shift */
		var_[c] = 1.0 / sqrt(var_[c]) ;

		norm_[c] *= (-1) * var_[c] ;
	}

	/* compute the mean and variance normalized feature */
	for(uint32 r = 0 ; r < num_row ; r++)
		for(uint32 c = 0 ; c < num_col ; c++)  input_feature[r * num_col + c] = (BaseFloat)(input_feature[r * num_col + c] * var_[c] + norm_[c]) ;
}

void CMVN::local_normalize(BaseFloat* input_feature , const uint32 num_row , const uint32 num_col)
{
	assert(dim_ == num_col) ;

	double* sum = (double*)malloc(dim_ * sizeof(double)) ;
	double* squared_sum = (double*)malloc(dim_ * sizeof(double)) ;
	double* norm = (double*)malloc(dim_ * sizeof(double)) ;
	double* squared_mean = (double*)malloc(dim_ * sizeof(double)) ;
	double* var = (double*)malloc(dim_ * sizeof(double)) ;
	if(NULL == sum || NULL == squared_sum || NULL == norm || squared_mean == NULL || var == NULL)  error("no memory to be allocated" , __FILE__ , __LINE__) ;

	memset((void*)sum , 0 , dim_ * sizeof(double)) ;
	memset((void*)squared_sum , 0 , dim_ * sizeof(double)) ;

	for(uint32 r = 0 ; r < num_cached_feature_ ; r++)
		for(uint32 c = 0 ; c < dim_ ; c++)
		{
			sum[c] += (double)stats_[r * dim_ + c] ;
			squared_sum[c] += pow(stats_[r * dim_ + c] , 2) ;
		}

	for(uint32 c = 0 ; c < dim_ ; c++)
	{
		norm[c] = sum[c] / num_cached_feature_ ;

		squared_mean[c] = squared_sum[c] / num_cached_feature_ ;

		var[c] = squared_mean[c] - pow(norm[c] , 2) ;

		if(var[c] < var_floor) var[c] = var_floor ;


		/* switch to scale and shift */
		var[c] = 1.0 / sqrt(var[c]) ;

		norm[c] *= (-1) * var[c] ;
	}

	/* compute the mean and variance normalized feature */
	for(uint32 r = 0 ; r < num_row ; r++)
		for(uint32 c = 0 ; c < num_col ; c++)  input_feature[r * num_col + c] = (BaseFloat)(input_feature[r * num_col + c] * var[c] + norm[c]) ;

	free(sum) ; sum = NULL ;
	free(squared_sum) ; squared_sum = NULL ;
	free(norm) ; norm = NULL ;
	free(squared_mean) ; squared_mean = NULL ;
	free(var) ; var = NULL ;
}

