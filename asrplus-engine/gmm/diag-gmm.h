/*
 * diag-gmm.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef DIAG_GMM_H_
#define DIAG_GMM_H_

#include "base/common.h"
#include "matrix/Matrix.h"
#include "matrix/Vector.h"
#include <fstream>

class DiagGmm {

 public:
  /// Empty constructor.
  DiagGmm() : valid_gconsts_(false) { }

  /// Returns the number of mixture components in the GMM
  int32 NumGauss() const { return weights_.Dim(); }
  /// Returns the dimensionality of the Gaussian mean vectors
  int32 Dim() const { return means_invvars_.NumCols(); }


  /// This version of the LogLikelihoods function operates on
  /// a sequence of frames simultaneously; the row index of both "data" and
  /// "loglikes" is the frame index.
  void LogLikelihoods(const Matrix &data,
                      Matrix *loglikes) const;

  int32 ComputeGconsts();

  void Read(std::ifstream &in, bool binary);



 private:
  Vector gconsts_;
  bool valid_gconsts_;   ///< Recompute gconsts_ if false
  Vector weights_;        ///< weights (not log).
  Matrix inv_vars_;       ///< Inverted (diagonal) variances
  Matrix means_invvars_;  ///< Means times inverted variance

 private:
  const DiagGmm &operator=(const DiagGmm &other);  // Disallow assignment
};


#endif /* DIAG_GMM_H_ */
