/*
 * ivector-extractor.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef IVECTOR_EXTRACTOR_H_
#define IVECTOR_EXTRACTOR_H_

#include "base/common.h"
#include "matrix/Vector.h"
#include "matrix/Matrix.h"
#include "matrix/SpMatrix.h"
#include "base/posterior.h"
#include <vector>
#include <fstream>


class IvectorExtractor;
class IvectorExtractorUtteranceStats {
 public:
  IvectorExtractorUtteranceStats(int32 num_gauss, int32 feat_dim,
                                 bool need_2nd_order_stats):
      gamma_(num_gauss), X_(num_gauss, feat_dim) {}

  void AccStats(const Matrix &feats, const Posterior &post);

  void Scale(double scale); // Used to apply acoustic scale.

  double NumFrames() { return gamma_.Sum(); }

 protected:
  friend class IvectorExtractor;
  Vector gamma_; // zeroth-order stats (summed posteriors), dimension [I]
  Matrix X_; // first-order stats, dimension [I][D]
  std::vector<SpMatrix> S_; // 2nd-order stats, dimension [I][D][D], if
                                     // required.
};


struct IvectorExtractorOptions {
  int ivector_dim;
  int num_iters;
  bool use_weights;
  IvectorExtractorOptions(): ivector_dim(600), num_iters(2),
                             use_weights(true) { }
  void Register(OptionsItf *opts) {
    opts->Register("num-iters", &num_iters, "Number of iterations in "
                   "iVector estimation (>1 needed due to weights)");
    opts->Register("ivector-dim", &ivector_dim, "Dimension of iVector");
    opts->Register("use-weights", &use_weights, "If true, regress the "
                   "log-weights on the iVector");
  }
};



class IvectorExtractor {
 public:


  IvectorExtractor(): prior_offset_(0.0) { }

  /// Gets the distribution over ivectors (or at least, a Gaussian approximation
  /// to it).  The output "var" may be NULL if you don't need it.  "mean", and
  /// "var", if present, must be the correct dimension (this->IvectorDim()).
  /// If you only need a point estimate of the iVector, get the mean only.
  void GetIvectorDistribution(
      const IvectorExtractorUtteranceStats &utt_stats,
      Vector *mean,
      SpMatrix *var) const;

  void GetIvectorDistMean(
      const IvectorExtractorUtteranceStats &utt_stats,
      Vector *linear,
      SpMatrix *quadratic) const;

  /// Gets the linear and quadratic terms in the distribution over
  /// iVectors, that arise from the prior.  Adds to the outputs,
  /// rather than setting them.
  void GetIvectorDistPrior(
      const IvectorExtractorUtteranceStats &utt_stats,
      Vector *linear,
      SpMatrix *quadratic) const;
  /// The distribution over iVectors, in our formulation, is not centered at
  /// zero; its first dimension has a nonzero offset.  This function returns
  /// that offset.
  double PriorOffset() const { return prior_offset_; }


  int32 FeatDim() const;
  int32 IvectorDim() const;
  int32 NumGauss() const;
  bool IvectorDependentWeights() const { return w_.NumRows() != 0; }
  void Read(std::ifstream &is, bool binary);

  // Note: we allow the default assignment and copy operators
  // because they do what we want.
 protected:

  void ComputeDerivedVars();

  void ComputeDerivedVars(int32 i);

  // Imagine we'll project the iVectors with transformation T, so apply T^{-1}
  // where necessary to keep the model equivalent.  Used to keep unit variance
  // (like prior re-estimation).
  void TransformIvectors(const Matrix &T,
                         double new_prior_offset);


  /// Weight projection vectors, if used.  Dimension is [I][S]
  Matrix w_;

  /// If we are not using weight-projection vectors, stores the Gaussian mixture
  /// weights from the UBM.  This does not affect the iVector; it is only useful
  /// as a way of making sure the log-probs are comparable between systems with
  /// and without weight projection matrices.
  Vector w_vec_;

  /// Ivector-subspace projection matrices, dimension is [I][D][S].
  /// The I'th matrix projects from ivector-space to Gaussian mean.
  /// There is no mean offset to add-- we deal with it by having
  /// a prior with a nonzero mean.
  std::vector<Matrix > M_;

  /// Inverse variances of speaker-adapted model, dimension [I][D][D].
  std::vector<SpMatrix > Sigma_inv_;

  /// 1st dim of the prior over the ivector has an offset, so it is not zero.
  /// This is used to handle the global offset of the speaker-adapted means in a
  /// simple way.
  float prior_offset_;

  // Below are *derived variables* that can be computed from the
  // variables above.

  /// The constant term in the log-likelihood of each Gaussian (not
  /// counting any weight).
  Vector gconsts_;

  /// U_i = M_i^T \Sigma_i^{-1} M_i is a quantity that comes up
  /// in ivector estimation.  This is conceptually a
  /// std::vector<SpMatrix<double> >, but we store the packed-data
  /// in the rows of a matrix, which gives us an efficiency
  /// improvement (we can use matrix-multiplies).
  Matrix U_;

  /// The product of Sigma_inv_[i] with M_[i].
  std::vector<Matrix > Sigma_inv_M_;
 private:
  // var <-- quadratic_term^{-1}, but done carefully, first flooring eigenvalues
  // of quadratic_term to 1.0, which mathematically is the least they can be,
  // due to the prior term.
  static void InvertWithFlooring(const SpMatrix &quadratic_term,
                                 SpMatrix *var);
};

#endif /* IVECTOR_EXTRACTOR_H_ */
