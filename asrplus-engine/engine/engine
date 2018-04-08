/*
 * engine.cc
 *
 *  Created on: Sep 19, 2017
 *      Author: tao
 */
#include "engine/engine.h"
#include "feature/mfcc.h"
#include "feature/feature-common.h"
#include "feature/mfcc.h"
#include "util/math.h"
#include "ivector/vad.h"


	VPengine::~VPengine(){

		release();
	}

	void VPengine::release(){

		if(container_ != NULL){
			free(container_);
			container_ = NULL;
		}
		if(datainterface_ != NULL){
			free(datainterface_);
			datainterface_ = NULL;
		}

		if(mfcc_ != NULL){
			delete mfcc_;
			mfcc_ = NULL;
		}

		if(cmn_ != NULL){
			delete cmn_ ;
			cmn_ = NULL ;
		}

	}

	void VPengine::initialize(VPcontainer* container){

		container_ = container;
		assert(container_ != NULL);
		datainterface_ = (DataInterface*)malloc(sizeof(DataInterface));
		datainterface_->initialize();
		assert(datainterface_ != NULL);
		MfccOptions opt;
		mfcc_ = new Mfcc(opt);
		assert(mfcc_ != NULL);
		cmn_ = new CMN ;
		assert(cmn_ != NULL);
	}

	void VPengine::sign(std::string key, int32 length){

		std::map<std::string, int32>::iterator it;
		Vector tmp_vec ,ivector;
		Matrix feats , loglikes, feats_vad;
		int32 num_post=20;
		BaseFloat min_post=0.025;
		int32 max_count=0;
		double max_count_scale = 1.0;
		BaseFloat acoustic_weight = 1.0;
		PldaConfig pldaconfig ;
		VadEnergyOptions vad_opt;
		Vector output_voiced;

		IvectorExtractorUtteranceStats utt_states(container_->ie_->NumGauss(), container_->ie_->FeatDim(), false);
		it = count_.find(key);
		if(it != count_.end()){

			if(count_[key]==3){
				std::cout<<key <<" Already Signed!"<<std::endl;
				return;
			}else{
				tmp_vec.Resize(length);
				datainterface_->read(tmp_vec.Data(), length);
				mfcc_->Compute(tmp_vec, 1.0, &feats);
				ComputeVadEnergy(vad_opt,feats,&output_voiced);
				cmn_->apply_cmn_offline(feats.Data(), feats.NumRows(), feats.NumCols());
				int32 dim = 0;
				for(int32 i = 0; i < feats.NumRows(); i++){
					if(output_voiced(i) != 0.0){
						assert(output_voiced(i) == 1.0);
						dim++;
					}
				}

				feats_vad.Resize(dim, feats.NumCols());
				int32 index = 0;
				for(int32 i = 0; i < feats.NumRows(); i++){
					if(output_voiced(i) == 1.0){
						feats_vad.CopyRowFromVec(feats.Row(i), index);
						index++;
					}
				}

				container_->gmm_->LogLikelihoods(feats_vad, &loglikes);
				int32 T = feats_vad.NumRows();
				Posterior post(T);
				for (int32 t = 0; t < T; t++) {
				    VectorToPosteriorEntry(loglikes.Row(t), num_post, min_post, &(post[t]));
				}
				double this_t =  acoustic_weight * TotalPosterior(post);
				if (max_count > 0 && this_t > max_count) {
				    max_count_scale = max_count / this_t;
				    this_t = max_count;
				}
				ScalePosterior(acoustic_weight*max_count_scale, &post);
				utt_states.AccStats(feats_vad, post);
				ivector.Resize(container_->ie_->IvectorDim());
				ivector(0) = container_->ie_->PriorOffset();
				container_->ie_->GetIvectorDistribution(utt_states, &ivector, NULL);
				ivector(0) -= container_->ie_->PriorOffset();
				//ivector.Print(30);
				//normalize
			    BaseFloat norm = ivector.Norm(2.0);
			    BaseFloat ratio = norm / sqrt(ivector.Dim());
			    ivector.Scale(1.0 / ratio);
			    //mean
				Vector ivector_averge = database_[key];
				ivector_averge.AddVec(1.0, ivector);
				count_[key] += 1;
				database_[key] = ivector_averge ;

				if(count_[key] == 3){
					ivector_averge.Scale(1/3.0);
					//normalize
				    norm = ivector_averge.Norm(2.0);
				    ratio = norm / sqrt(ivector.Dim());
				    ivector_averge.Scale(1.0 / ratio);
				    //sub mean
				    /*ivector_averge.AddVec(-1.0, *container_->mean_);
				    //normalize
				    norm = ivector_averge.Norm(2.0);
				    ratio = norm / sqrt(ivector.Dim());
				    ivector_averge.Scale(1.0 / ratio);
	*/
					database_[key] = ivector_averge ;
				}

			}

		}else{
			tmp_vec.Resize(length);
			datainterface_->read(tmp_vec.Data(), length);
			mfcc_->Compute(tmp_vec, 1.0, &feats);
			ComputeVadEnergy(vad_opt,feats,&output_voiced);
			cmn_->apply_cmn_offline(feats.Data(), feats.NumRows(), feats.NumCols());
			int32 dim = 0;
			for(int32 i = 0; i < feats.NumRows(); i++){
				if(output_voiced(i) != 0.0){
					assert(output_voiced(i) == 1.0);
					dim++;
				}
			}

			feats_vad.Resize(dim, feats.NumCols());
			int32 index = 0;
			for(int32 i = 0; i < feats.NumRows(); i++){
				if(output_voiced(i) == 1.0){
					feats_vad.CopyRowFromVec(feats.Row(i), index);
					index++;
				}
			}

			container_->gmm_->LogLikelihoods(feats_vad, &loglikes);

			int32 T = feats_vad.NumRows();
			Posterior post(T);
		    for (int32 t = 0; t < T; t++) {

		        VectorToPosteriorEntry(loglikes.Row(t), num_post, min_post, &(post[t]));
		    }
		    double this_t =  acoustic_weight*TotalPosterior(post);
	        if (max_count > 0 && this_t > max_count) {
	            max_count_scale = max_count / this_t;
	            this_t = max_count;
	        }
	        ScalePosterior(acoustic_weight*max_count_scale, &post);
	        utt_states.AccStats(feats_vad, post);
	        ivector.Resize(container_->ie_->IvectorDim());
	        ivector(0) = container_->ie_->PriorOffset();
			container_->ie_->GetIvectorDistribution(utt_states, &ivector, NULL);
		    ivector(0) -= container_->ie_->PriorOffset();
		    //normalize
		    BaseFloat norm = ivector.Norm(2.0);
		    BaseFloat ratio = norm / sqrt(ivector.Dim());
		    ivector.Scale(1.0 / ratio);
			database_[key] = ivector ;
			count_[key] = 1 ;
		}
	}

	void VPengine::check(int32 length){

		Vector tmp_vec ,ivector;
		Matrix feats , feats_vad, loglikes;
		int32 num_post=20;
		BaseFloat min_post=0.025;
		int32 max_count=0;
		double max_count_scale = 1.0;
		BaseFloat acoustic_weight = 1.0;
		BaseFloat score=0.0, max_value=-1e20;
		std::string key ;
		PldaConfig pldaconfig ;
		VadEnergyOptions vad_opt;
		Vector output_voiced;

		IvectorExtractorUtteranceStats utt_states(container_->ie_->NumGauss(), container_->ie_->FeatDim(), false);
		tmp_vec.Resize(length);
		datainterface_->read(tmp_vec.Data(), length);
		mfcc_->Compute(tmp_vec, 1.0, &feats);
		ComputeVadEnergy(vad_opt,feats,&output_voiced);

		cmn_->apply_cmn_offline(feats.Data(), feats.NumRows(), feats.NumCols());
		int32 dim = 0;
		for(int32 i = 0; i < feats.NumRows(); i++){
			if(output_voiced(i) != 0.0){
				assert(output_voiced(i) == 1.0);
				dim++;
			}
		}
		feats_vad.Resize(dim, feats.NumCols());
		int32 index = 0;
		for(int32 i = 0; i < feats.NumRows(); i++){
			if(output_voiced(i) == 1.0){
				feats_vad.CopyRowFromVec(feats.Row(i), index);
				index++;
			}
		}
		container_->gmm_->LogLikelihoods(feats_vad, &loglikes);

		int32 T = feats_vad.NumRows();
		Posterior post(T);
	    for (int32 t = 0; t < T; t++) {
	        VectorToPosteriorEntry(loglikes.Row(t), num_post, min_post, &(post[t]));
	    }
	    double this_t =  acoustic_weight*TotalPosterior(post);
        if (max_count > 0 && this_t > max_count) {
            max_count_scale = max_count / this_t;
            this_t = max_count;
        }
        ScalePosterior(acoustic_weight*max_count_scale, &post);
        utt_states.AccStats(feats_vad, post);
        ivector.Resize(container_->ie_->IvectorDim());
        ivector(0) = container_->ie_->PriorOffset();
		container_->ie_->GetIvectorDistribution(utt_states, &ivector, NULL);
	    ivector(0) -= container_->ie_->PriorOffset();

	    //normalize
	    BaseFloat norm = ivector.Norm(2.0);
	    BaseFloat ratio = norm / sqrt(ivector.Dim());
	    ivector.Scale(1.0 / ratio);

	    //sub mean
	   /* ivector.AddVec(-1.0, *container_->mean_);

        //normalize
	    norm = ivector.Norm(2.0);
	    ratio = norm / sqrt(ivector.Dim());
	    ivector.Scale(1.0 / ratio);*/
	    //plda

		std::map<std::string, Vector>::iterator it;
		it = database_.begin();

		while(it != database_.end()){
			score = CosDistance(it->second,  ivector);
			if(max_value < score) {
				max_value = score ;
				key = it->first ;
			}
			it++;
		}

		if(max_value > threshold_){

			std::cout<<"the wave come from "<<key <<" the score is: "<<max_value<<std::endl;
		}else{
			std::cout<<"the owner of the wave is not in the database!"<<"the score is: "<<max_value<<" the key is: "<< key<<std::endl;
		}
	}

	void VPengine::put_data(short* short_data , int32 length){
		datainterface_->input_data((int16*)short_data, length);
	}

	void VPengine::data_end(){
		datainterface_->data_finished() ;
	}


