/*
 * plda.h
 *
 *  Created on: Mar 28, 2018
 *      Author: tao
 */

#ifndef PLDA_H_
#define PLDA_H_

#include "matrix/Vector.h"
#include "matrix/Matrix.h"
#include "base/io.h"

struct PldaConfig {

  bool normalize_length;
  bool simple_length_norm;
  PldaConfig(): normalize_length(true), simple_length_norm(false) { }
  void Register(OptionsItf *opts) {
    opts->Register("normalize-length", &normalize_length,
                   "If true, do length normalization as part of PLDA (see "
                   "code for details).  This does not set the length unit; "
                   "by default it instead ensures that the inner product "
                   "with the PLDA model's inverse variance (which is a "
                   "function of how many utterances the iVector was averaged "
                   "over) has the expected value, equal to the iVector "
                   "dimension.");

    opts->Register("simple-length-normalization", &simple_length_norm,
                   "If true, replace the default length normalization by an "
                   "alternative that normalizes the length of the iVectors to "
                   "be equal to the square root of the iVector dimension.");
  }
};

class Plda{

public:
	Plda() { }


	float TransformIvector(const PldaConfig &config,
	                         const Vector &ivector,
	                         int32 num_examples,
	                         Vector *transformed_ivector) const;


    float LogLikelihoodRatio(const Vector &transformed_train_ivector,
	                            int32 num_train_utts,
	                            const Vector &transformed_test_ivector)const;
    void ComputeDerivedVars();

    int32 Dim() const { return mean_.Dim(); }

    void Read(std::ifstream &is, bool binary);

protected:

	  Vector mean_;  // mean of samples in original space.
	  Matrix transform_; // of dimension Dim() by Dim(); this transform makes within-class covar unit and diagonalizes the between-class covar.
	  Vector psi_; // of dimension Dim().  The between-class (diagonal) covariance elements, in decreasing order.
	  Vector offset_;  // derived variable: -1.0 * transform_ * mean_

private:

	  float GetNormalizationFactor(const Vector &transformed_ivector,
	                                int32 num_examples) const;
};



#endif /* PLDA_H_ */
