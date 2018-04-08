/*
 * ivector-extractor.cc
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#include "ivector/ivector-extractor.h"
#include "base/io.h"


int32 IvectorExtractor::FeatDim() const {
  assert(!M_.empty());
  return M_[0].NumRows();
}

int32 IvectorExtractor::IvectorDim() const {
  if (M_.empty()) { return 0; }
  else { return M_[0].NumCols(); }
}

int32 IvectorExtractor::NumGauss() const {
  return static_cast<int32>(M_.size());
}


void IvectorExtractor::GetIvectorDistribution(
    const IvectorExtractorUtteranceStats &utt_stats,
    Vector *mean,
    SpMatrix *var) const {

    Vector linear(IvectorDim());
    SpMatrix quadratic(IvectorDim());
    GetIvectorDistMean(utt_stats, &linear, &quadratic);
    GetIvectorDistPrior(utt_stats, &linear, &quadratic);
    if (var != NULL) {
      var->CopyFromSp(quadratic);
      var->Invert(); // now it's a variance.

      // mean of distribution = quadratic^{-1} * linear...
      mean->AddSpVec(1.0, *var, linear, 0.0);
    } else {
      quadratic.Invert();
      mean->AddSpVec(1.0, quadratic, linear, 0.0);
    }
}

void IvectorExtractor::GetIvectorDistMean(
    const IvectorExtractorUtteranceStats &utt_stats,
    Vector *linear,
    SpMatrix *quadratic) const {
  int32 N = NumGauss();
  for (int32 i = 0; i < N; i++) {
    double gamma = utt_stats.gamma_(i);
    if (gamma != 0.0) {
      Vector x(utt_stats.X_, i); // == \gamma(i) \m_i
      // next line: a += \gamma_i \M_i^T \Sigma_i^{-1} \m_i
      linear->AddMatVec(1.0, Sigma_inv_M_[i], kTrans, x, 1.0);
    }
  }
  Vector q_vec( IvectorDim()*(IvectorDim()+1)/2);
  q_vec.CopyFromData(quadratic->Data(), IvectorDim()*(IvectorDim()+1)/2);
  q_vec.AddMatVec(1.0, U_, kTrans, utt_stats.gamma_, 1.0);
  quadratic->CopyFromVec(q_vec);
}

void IvectorExtractor::GetIvectorDistPrior(
    const IvectorExtractorUtteranceStats &utt_stats,
    Vector *linear,
    SpMatrix *quadratic) const {

  (*linear)(0) += prior_offset_; // the zero'th dimension has an offset mean.
  /// The inverse-variance for the prior is the unit matrix.
  quadratic->AddToDiag(1.0);
}

void IvectorExtractor::ComputeDerivedVars() {
  std::cout << "Computing derived variables for iVector extractor"<<std::endl;
  gconsts_.Resize(NumGauss());
  for (int32 i = 0; i < NumGauss(); i++) {
    double var_logdet = -Sigma_inv_[i].LogPosDefDet();
    gconsts_(i) = -0.5 * (var_logdet + FeatDim() * M_LOG_2PI);
    // the gconsts don't contain any weight-related terms.
  }
  U_.Resize(NumGauss(), IvectorDim() * (IvectorDim() + 1) / 2);
  Sigma_inv_M_.resize(NumGauss());

  for (int32 i = 0; i < NumGauss(); i++){
	  ComputeDerivedVars(i);
  }

  std::cout << "Done."<< std::endl;
}

void IvectorExtractor::ComputeDerivedVars(int32 i) {
  SpMatrix temp_U(IvectorDim());
  // temp_U = M_i^T Sigma_i^{-1} M_i
  temp_U.AddMat2Sp(1.0, M_[i], kTrans, Sigma_inv_[i], 0.0);
  Vector temp_U_vec(IvectorDim() * (IvectorDim() + 1) / 2);
  temp_U_vec.CopyFromData(temp_U.Data(), IvectorDim() * (IvectorDim() + 1) / 2);

  U_.CopyRowFromVec(temp_U_vec, i);

  Sigma_inv_M_[i].Resize(FeatDim(), IvectorDim());
  Sigma_inv_M_[i].AddSpMat(1.0, Sigma_inv_[i], M_[i], kNoTrans, 0.0);
}


void IvectorExtractor::InvertWithFlooring(const SpMatrix &inverse_var,
                                          SpMatrix *var) {
  SpMatrix dbl_var(inverse_var);
  int32 dim = inverse_var.NumRows();
  Vector s(dim);
  Matrix P(dim, dim);
  // Solve the symmetric eigenvalue problem, inverse_var = P diag(s) P^T.
  inverse_var.Eig(&s, &P);
  s.ApplyFloor(1.0);
  s.InvertElements();
  var->AddMat2Vec(1.0, P, kNoTrans, s, 0.0);
}



void IvectorExtractor::Read(std::ifstream &is, bool binary) {
  ExpectToken(is, binary, "<IvectorExtractor>");
  ExpectToken(is, binary, "<w>");
  w_.Read(is, binary);
  ExpectToken(is, binary, "<w_vec>");
  w_vec_.Read(is, binary);
  ExpectToken(is, binary, "<M>");
  int32 size;
  ReadBasicType(is, binary, &size);
  assert(size > 0);
  M_.resize(size);
  for (int32 i = 0; i < size; i++)
    M_[i].Read(is, binary);
  ExpectToken(is, binary, "<SigmaInv>");
  Sigma_inv_.resize(size);
  for (int32 i = 0; i < size; i++)
    Sigma_inv_[i].Read(is, binary);
  ExpectToken(is, binary, "<IvectorOffset>");
  ReadBasicType(is, binary, &prior_offset_);
  ExpectToken(is, binary, "</IvectorExtractor>");
  ComputeDerivedVars();
}


void IvectorExtractorUtteranceStats::AccStats(
    const Matrix &feats,
    const Posterior &post) {
  typedef std::vector<std::pair<int32, BaseFloat> > VecType;
  int32 num_frames = feats.NumRows(),
      num_gauss = X_.NumRows(),
      feat_dim = feats.NumCols();
  assert(X_.NumCols() == feat_dim);
  assert(feats.NumRows() == static_cast<int32>(post.size()));
  bool update_variance = (!S_.empty());
  SpMatrix outer_prod(feat_dim);
  for (int32 t = 0; t < num_frames; t++) {
    Vector frame(feats, t);
    const VecType &this_post(post[t]);
    if (update_variance) {
      outer_prod.SetZero();
      outer_prod.AddVec2(1.0, frame);
    }
    for (VecType::const_iterator iter = this_post.begin();
         iter != this_post.end(); ++iter) {
      int32 i = iter->first; // Gaussian index.
      assert(i >= 0 && i < num_gauss &&
                   "Out-of-range Gaussian (mismatched posteriors?)");
      BaseFloat weight = iter->second;
      gamma_(i) += weight;
      X_.AddVec2Row(i, weight, frame);
      if (update_variance)
        S_[i].AddSp(weight, outer_prod);
    }
  }
}


