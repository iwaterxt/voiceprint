/*
 * diag-gmm.cc
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#include "gmm/diag-gmm.h"
#include "base/io.h"
#include <assert.h>

void DiagGmm::LogLikelihoods(const Matrix &data,
                             Matrix *loglikes) const {
  assert(data.NumRows() != 0);
  loglikes->Resize(data.NumRows(), gconsts_.Dim());
  loglikes->CopyRowsFromVec(gconsts_);
  if (data.NumCols() != Dim()) {
    std::cerr << "DiagGmm::ComponentLogLikelihood, dimension "
              << "mismatch " << data.NumCols() << " vs. "<< Dim();
  }
  Matrix data_sq(data);
  data_sq.ApplyPow(2.0);

  // loglikes +=  means * inv(vars) * data.
  loglikes->AddMatMat(1.0, data, kNoTrans, means_invvars_, kTrans, 1.0);
  // loglikes += -0.5 * inv(vars) * data_sq.
  loglikes->AddMatMat(-0.5, data_sq, kNoTrans, inv_vars_, kTrans, 1.0);
}


void DiagGmm::Read(std::ifstream &is, bool binary) {
//  ExpectToken(is, binary, "<DiagGMMBegin>");
  std::string token;
  ReadToken(is, binary, &token);
  // <DiagGMMBegin> is for compatibility. Will be deleted later
  if (token != "<DiagGMMBegin>" && token != "<DiagGMM>")
    std::cerr << "Expected <DiagGMM>, got " << token;
  ReadToken(is, binary, &token);
  if (token == "<GCONSTS>") {  // The gconsts are optional.
    gconsts_.Read(is, binary);
    ExpectToken(is, binary, "<WEIGHTS>");
  } else {
    if (token != "<WEIGHTS>")
      std::cerr << "DiagGmm::Read, expected <WEIGHTS> or <GCONSTS>, got "
                << token;
  }
  weights_.Read(is, binary);
  ExpectToken(is, binary, "<MEANS_INVVARS>");
  means_invvars_.Read(is, binary);
  ExpectToken(is, binary, "<INV_VARS>");
  inv_vars_.Read(is, binary);
//  ExpectToken(is, binary, "<DiagGMMEnd>");
  ReadToken(is, binary, &token);
  // <DiagGMMEnd> is for compatibility. Will be deleted later
  if (token != "<DiagGMMEnd>" && token != "</DiagGMM>")
    std::cerr << "Expected </DiagGMM>, got " << token;

  ComputeGconsts();  // safer option than trusting the read gconsts
}

int32 DiagGmm::ComputeGconsts() {
  int32 num_mix = NumGauss();
  int32 dim = Dim();
  BaseFloat offset = -0.5 * M_LOG_2PI * dim;  // constant term in gconst.
  int32 num_bad = 0;

  // Resize if Gaussians have been removed during Update()
  if (num_mix != static_cast<int32>(gconsts_.Dim()))
    gconsts_.Resize(num_mix);

  for (int32 mix = 0; mix < num_mix; mix++) {
    assert(weights_(mix) >= 0);  // Cannot have negative weights.
    BaseFloat gc = Log(weights_(mix)) + offset;  // May be -inf if weights == 0
    for (int32 d = 0; d < dim; d++) {
      gc += 0.5 * Log(inv_vars_(mix, d)) - 0.5 * means_invvars_(mix, d)
        * means_invvars_(mix, d) / inv_vars_(mix, d);
    }
    // Change sign for logdet because var is inverted. Also, note that
    // mean_invvars(mix, d)*mean_invvars(mix, d)/inv_vars(mix, d) is the
    // mean-squared times inverse variance, since mean_invvars(mix, d) contains
    // the mean times inverse variance.
    // So gc is the likelihood at zero feature value.

    if (isnan(gc)) {  // negative infinity is OK but NaN is not acceptable
      std::cerr << "At component "  << mix
                << ", not a number in gconst computation";
    }
    if (isinf(gc)) {
      num_bad++;
      // If positive infinity, make it negative infinity.
      // Want to make sure the answer becomes -inf in the end, not NaN.
      if (gc > 0) gc = -gc;
    }
    gconsts_(mix) = gc;
  }

  valid_gconsts_ = true;
  return num_bad;
}

std::istream & operator >>(std::ifstream &is, DiagGmm &gmm) {
  gmm.Read(is, false);  // false == non-binary.
  return is;
}




