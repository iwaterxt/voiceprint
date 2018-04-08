/*
 * vpcontainer.cc
 *
 *  Created on: Oct 11, 2017
 *      Author: tao
 */
#include "gmm/diag-gmm.h"
#include "ivector/ivector-extractor.h"
#include "ivector/plda.h"
#include "feature/mfcc.h"
#include "engine/vpcontainer.h"
#include <fstream>

void VPcontainer::Initialize(const struct vpconfig *config){


		ivector_dim_ = config->ivector_dim ;
		num_post_ = config->num_post ;
		min_post_ = config->min_post ;

		std::string gmm_file = config->final_dubm;
		std::ifstream is(gmm_file.c_str());

		//gmm_ = (DiagGmm*)malloc(sizeof(DiagGmm));
		gmm_ = new DiagGmm;
		gmm_->Read(is, false);

		//check if read gmm sucessful
		assert(gmm_->NumGauss() != 0);

		std::string ivector_extract_file = config->final_ie ;
		std::ifstream iv(ivector_extract_file.c_str());
		ie_ = new IvectorExtractor;
		ie_->Read(iv, false);
		//check if ivector extractor read sucessful
		assert(ie_->NumGauss() != 0);
		assert(ie_->IvectorDim() == ivector_dim_);

		std::string plda_file = config->final_plda ;
		std::ifstream ip(plda_file.c_str());
		plda_ = new Plda;
		plda_->Read(ip, false);
		assert(plda_->Dim() != 0);

		std::string mean_file = config->final_mean ;
		std::ifstream im(mean_file.c_str());
		mean_ = new Vector ;
		mean_->Read(im, false);
		assert(mean_->Dim() == ie_->IvectorDim());

		std::cout<<"-----------container initialize successful!-------------"<<std::endl;
}

void VPcontainer::Release(){

	if(gmm_ != NULL){
		delete gmm_;
		gmm_ = NULL;
	}

	if(ie_ != NULL){
		delete ie_;
		ie_ = NULL;
	}

	if(plda_ != NULL){
		delete plda_;
		plda_ = NULL;
	}

	if(mean_ != NULL){
		delete mean_;
		mean_ = NULL;
	}
}

