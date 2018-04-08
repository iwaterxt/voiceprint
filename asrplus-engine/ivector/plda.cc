/*
 * plda.cc
 *
 *  Created on: Mar 28, 2018
 *      Author: tao
 */



#include "ivector/plda.h"
#include "util/math.h"
#include <fstream>

float Plda::TransformIvector(const PldaConfig &config,
        const Vector &ivector,
        int32 num_examples,
        Vector *transformed_ivector) const{

	  Vector tmp(ivector), tmp_out(ivector.Dim());
	  assert(ivector.Dim() == Dim() && transformed_ivector->Dim() == Dim());
	  float normalization_factor;
	  transformed_ivector->CopyFromVec(offset_);
	  transformed_ivector->AddMatVec(1.0, transform_, kNoTrans, ivector, 1.0);
	  if (config.simple_length_norm)
	    normalization_factor = sqrt(transformed_ivector->Dim())
	      / transformed_ivector->Norm(2.0);
	  else
	    normalization_factor = GetNormalizationFactor(*transformed_ivector,
	                                                  num_examples);
	  if (config.normalize_length)
	    transformed_ivector->Scale(normalization_factor);
	  return normalization_factor;
}


float Plda::LogLikelihoodRatio(const Vector &transformed_train_ivector,
        int32 n,
        const Vector &transformed_test_ivector)const{

      //std::cout<<"the sum value of transformed_train_ivector is: "<<transformed_train_ivector.Sum()<<std::endl;
      //std::cout<<"the sum value of transformed_test_ivector is :"<<transformed_test_ivector.Sum()<<std::endl;

	  int32 dim = Dim();
	  float loglike_given_class, loglike_without_class;
	  { // work out loglike_given_class.
	    // "mean" will be the mean of the distribution if it comes from the
	    // training example.  The mean is \frac{n \Psi}{n \Psi + I} \bar{u}^g
	    // "variance" will be the variance of that distribution, equal to
	    // I + \frac{\Psi}{n\Psi + I}.
	    Vector mean(dim);
	    Vector variance(dim);
	    for (int32 i = 0; i < dim; i++) {
	      mean(i) = n * psi_(i) / (n * psi_(i) + 1.0)
	        * transformed_train_ivector(i);
	      variance(i) = 1.0 + psi_(i) / (n * psi_(i) + 1.0);
	    }
	    float logdet = variance.SumLog();
	    Vector sqdiff(transformed_test_ivector);
	    sqdiff.AddVec(-1.0, mean);
	    sqdiff.ApplyPow(2.0);
	    variance.InvertElements();
	    loglike_given_class = -0.5 * (logdet + M_LOG_2PI * dim +
	                                  VecVec(sqdiff, variance));
	  }
	  { // work out loglike_without_class.  Here the mean is zero and the variance
	    // is I + \Psi.
	    Vector sqdiff(transformed_test_ivector); // there is no offset.
	    sqdiff.ApplyPow(2.0);
	    Vector variance(psi_);
	    variance.Add(1.0); // I + \Psi.
	    double logdet = variance.SumLog();
	    variance.InvertElements();
	    loglike_without_class = -0.5 * (logdet + M_LOG_2PI * dim +
	                                    VecVec(sqdiff, variance));
	  }
	  float loglike_ratio = loglike_given_class - loglike_without_class;
	  return loglike_ratio;
}

void Plda::ComputeDerivedVars() {
  assert(Dim() > 0);
  offset_.Resize(Dim());
  offset_.AddMatVec(-1.0, transform_, kNoTrans, mean_, 0.0);
}

void Plda::Read(std::ifstream &is, bool binary){

	  ExpectToken(is, binary, "<Plda>");
	  mean_.Read(is, binary);
	  transform_.Read(is, binary);
	  psi_.Read(is, binary);
	  ExpectToken(is, binary, "</Plda>");
	  ComputeDerivedVars();
}


float Plda::GetNormalizationFactor(const Vector &transformed_ivector,
	                                int32 num_examples) const{

	  assert(num_examples > 0);
	  Vector transformed_ivector_sq(transformed_ivector.Dim());
	  transformed_ivector_sq.CopyFromVec(transformed_ivector);
	  transformed_ivector_sq.ApplyPow(2.0);
	  Vector inv_covar(psi_);
	  inv_covar.Add(1.0 / num_examples);
	  inv_covar.InvertElements();

	  float dot_prod = VecVec(inv_covar, transformed_ivector_sq);
	  return sqrt(Dim() / dot_prod);

}
