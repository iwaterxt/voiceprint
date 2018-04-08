/*
 * vpcontainer.h
 *
 *  Created on: Oct 11, 2017
 *      Author: tao
 */

#ifndef VPCONTAINER_H_
#define VPCONTAINER_H_

#include "gmm/diag-gmm.h"
#include "ivector/ivector-extractor.h"
#include "ivector/plda.h"
#include "base/posterior.h"
#include "engine/asr-plusplus.h"



class VPcontainer{
public:

	VPcontainer():ivector_dim_(600), num_post_(20), min_post_(0.025),gmm_(NULL), ie_(NULL), plda_(NULL), mean_(NULL){}

	~VPcontainer();

	void Initialize(const struct vpconfig *config);

	void Release();


	int32 ivector_dim_;
	int32 num_post_;
	float min_post_;
	DiagGmm* gmm_ ;
	IvectorExtractor* ie_ ;
	Plda *plda_ ;
	Vector *mean_;




};


#endif /* VPCONTAINER_H_ */
