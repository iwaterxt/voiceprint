/*
 * math.cc
 *
 *  Created on: Feb 7, 2018
 *      Author: tao
 */

#include <assert.h>
#include <iostream>

#include "util/math.h"
#include "base/cblas-warppers.h"


int32 RoundUpToNearestPowerOfTwo(int32 n) {
  assert(n > 0);
  n--;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  return n+1;
}


void ComplexImExp(BaseFloat x, BaseFloat *a_re, BaseFloat *a_im) {
  *a_re = std::cos(x);
  *a_im = std::sin(x);
}

void ComplexMul(const BaseFloat &a_re, const BaseFloat &a_im,
                                            BaseFloat *b_re, BaseFloat *b_im) {
  BaseFloat tmp_re = (*b_re * a_re) - (*b_im * a_im);
  *b_im = *b_re * a_im + *b_im * a_re;
  *b_re = tmp_re;
}

void ComplexAddProduct(const BaseFloat &a_re, const BaseFloat &a_im,
                                                   const BaseFloat &b_re, const BaseFloat &b_im,
                                                   BaseFloat *c_re, BaseFloat *c_im) {
  *c_re += b_re*a_re - b_im*a_im;
  *c_im += b_re*a_im + b_im*a_re;
}

void Factorize(int32 m, std::vector<int32> *factors) {
  // Splits a number into its prime factors, in sorted order from
  // least to greatest,  with duplication.  A very inefficient
  // algorithm, which is mainly intended for use in the
  // mixed-radix FFT computation (where we assume most factors
  // are small).
  assert(factors != NULL);
  assert(m >= 1);  // Doesn't work for zero or negative numbers.
  factors->clear();
  int32 small_factors[10] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };

  // First try small factors.
  for (int32 i = 0; i < 10; i++) {
    if (m == 1) return;  // We're done.
    while (m % small_factors[i] == 0) {
      m /= small_factors[i];
      factors->push_back(small_factors[i]);
    }
  }
  // Next try all odd numbers starting from 31.
  for (int32 j = 31;; j += 2) {
    if (m == 1) return;
    while (m % j == 0) {
      m /= j;
      factors->push_back(j);
    }
  }
}

void ComplexFftRecursive (BaseFloat *data, int nffts, int N,
                          const int *factor_begin,
                          const int *factor_end, bool forward,
                          Vector *tmp_vec) {
  if (factor_begin == factor_end) {
    assert(N == 1);
    return;
  }

  {  // an optimization: compute in smaller blocks.
    // this block of code could be removed and it would still work.
    MatrixIndexT size_perblock = N * 2 * sizeof(BaseFloat);
    if (nffts > 1 && size_perblock*nffts > COMPLEXFFT_BLOCKSIZE) {  // can break it up...
      // Break up into multiple blocks.  This is an optimization.  We make
      // no progress on the FFT when we do this.
      int block_skip = COMPLEXFFT_BLOCKSIZE / size_perblock;  // n blocks per call
      if (block_skip == 0) block_skip = 1;
      if (block_skip < nffts) {
        int blocks_left = nffts;
        while (blocks_left > 0) {
          int skip_now = std::min(blocks_left, block_skip);
          ComplexFftRecursive(data, skip_now, N, factor_begin, factor_end, forward, tmp_vec);
          blocks_left -= skip_now;
          data += skip_now * N*2;
        }
        return;
      } // else do the actual algorithm.
    } // else do the actual algorithm.
  }

  int P = *factor_begin;
  assert(P > 1);
  int Q = N / P;


  if (P > 1 && Q > 1) {  // Do the rearrangement.   C.f. eq. (8) below.  Transform
    // (a) to (b).
    BaseFloat *data_thisblock = data;
    if (tmp_vec->Dim() < (MatrixIndexT)N) tmp_vec->Resize(N);
    BaseFloat *data_tmp = tmp_vec->Data();
    for (int thisfft = 0; thisfft < nffts; thisfft++, data_thisblock+=N*2) {
      for (int offset = 0; offset < 2; offset++) {  // 0 == real, 1 == im.
        for (int p = 0; p < P; p++) {
          for (int q = 0; q < Q; q++) {
            int aidx = q*P + p, bidx = p*Q + q;
            data_tmp[bidx] = data_thisblock[2*aidx+offset];
          }
        }
        for (int n = 0;n < P*Q;n++) data_thisblock[2*n+offset] = data_tmp[n];
      }
    }
  }

  {  // Recurse.
    ComplexFftRecursive(data, nffts*P, Q, factor_begin+1, factor_end, forward, tmp_vec);
  }

  int exp_sign = (forward ? -1 : 1);
  BaseFloat rootN_re, rootN_im;  // Nth root of unity.
  ComplexImExp(static_cast<BaseFloat>(exp_sign * M_2PI / N), &rootN_re, &rootN_im);

  BaseFloat rootP_re, rootP_im;  // Pth root of unity.
  ComplexImExp(static_cast<BaseFloat>(exp_sign * M_2PI / P), &rootP_re, &rootP_im);

  {  // Do the multiplication
    // could avoid a bunch of complex multiplies by moving the loop over data_thisblock
    // inside.
    if (tmp_vec->Dim() < (MatrixIndexT)(P*2)) tmp_vec->Resize(P*2);
    BaseFloat *temp_a = tmp_vec->Data();

    BaseFloat *data_thisblock = data, *data_end = data+(N*2*nffts);
    for (; data_thisblock != data_end; data_thisblock += N*2) {  // for each separate fft.
      BaseFloat qd_re = 1.0, qd_im = 0.0;  // 1^(q'/N)
      for (int qd = 0; qd < Q; qd++) {
        BaseFloat pdQ_qd_re = qd_re, pdQ_qd_im = qd_im;  // 1^((p'Q+q') / N) == 1^((p'/P) + (q'/N))
                                              // Initialize to q'/N, corresponding to p' == 0.
        for (int pd = 0; pd < P; pd++) {  // pd == p'
          {  // This is the p = 0 case of the loop below [an optimization].
            temp_a[pd*2] = data_thisblock[qd*2];
            temp_a[pd*2 + 1] = data_thisblock[qd*2 + 1];
          }
          {  // This is the p = 1 case of the loop below [an optimization]
            // **** MOST OF THE TIME (>60% I think) gets spent here. ***
            ComplexAddProduct(pdQ_qd_re, pdQ_qd_im,
                              data_thisblock[(qd+Q)*2], data_thisblock[(qd+Q)*2 + 1],
                              &(temp_a[pd*2]), &(temp_a[pd*2 + 1]));
          }
          if (P > 2) {
            BaseFloat p_pdQ_qd_re = pdQ_qd_re, p_pdQ_qd_im = pdQ_qd_im;  // 1^(p(p'Q+q')/N)
            for (int p = 2; p < P; p++) {
              ComplexMul(pdQ_qd_re, pdQ_qd_im, &p_pdQ_qd_re, &p_pdQ_qd_im);  // p_pdQ_qd *= pdQ_qd.
              int data_idx = p*Q + qd;
              ComplexAddProduct(p_pdQ_qd_re, p_pdQ_qd_im,
                                data_thisblock[data_idx*2], data_thisblock[data_idx*2 + 1],
                                &(temp_a[pd*2]), &(temp_a[pd*2 + 1]));
            }
          }
          if (pd != P-1)
            ComplexMul(rootP_re, rootP_im, &pdQ_qd_re, &pdQ_qd_im);  // pdQ_qd *= (rootP == 1^{1/P})
          // (using 1/P == Q/N)
        }
        for (int pd = 0; pd < P; pd++) {
          data_thisblock[(pd*Q + qd)*2] = temp_a[pd*2];
          data_thisblock[(pd*Q + qd)*2 + 1] = temp_a[pd*2 + 1];
        }
        ComplexMul(rootN_re, rootN_im, &qd_re, &qd_im);  // qd *= rootN.
      }
    }
  }
}



void ComplexFft(Vector *v, bool forward, Vector *tmp_in) {
  assert(v != NULL);

  if (v->Dim()<=1) return;
  assert(v->Dim() % 2 == 0);  // complex input.
  int N = v->Dim() / 2;
  std::vector<int> factors;
  Factorize(N, &factors);
  int *factor_beg = NULL;
  if (factors.size() > 0)
    factor_beg = &(factors[0]);
  Vector tmp;  // allocated in ComplexFftRecursive.
  ComplexFftRecursive(v->Data(), 1, N, factor_beg, factor_beg+factors.size(), forward, (tmp_in?tmp_in:&tmp));
}



void RealFft (Vector *v, bool forward) {
  assert(v != NULL);
  MatrixIndexT N = v->Dim(), N2 = N/2;
  assert(N%2 == 0);
  if (N == 0) return;

  if (forward) ComplexFft(v, true);

  BaseFloat *data = v->Data();
  BaseFloat rootN_re, rootN_im;  // exp(-2pi/N), forward; exp(2pi/N), backward
  int forward_sign = forward ? -1 : 1;
  ComplexImExp(static_cast<BaseFloat>(M_2PI/N *forward_sign), &rootN_re, &rootN_im);
  BaseFloat kN_re = -forward_sign, kN_im = 0.0;  // exp(-2pik/N), forward; exp(-2pik/N), backward
  // kN starts out as 1.0 for forward algorithm but -1.0 for backward.
  for (MatrixIndexT k = 1; 2*k <= N2; k++) {
    ComplexMul(rootN_re, rootN_im, &kN_re, &kN_im);

    BaseFloat Ck_re, Ck_im, Dk_re, Dk_im;
    // C_k = 1/2 (B_k + B_{N/2 - k}^*) :
    Ck_re = 0.5 * (data[2*k] + data[N - 2*k]);
    Ck_im = 0.5 * (data[2*k + 1] - data[N - 2*k + 1]);
    // re(D_k)= 1/2 (im(B_k) + im(B_{N/2-k})):
    Dk_re = 0.5 * (data[2*k + 1] + data[N - 2*k + 1]);
    // im(D_k) = -1/2 (re(B_k) - re(B_{N/2-k}))
    Dk_im =-0.5 * (data[2*k] - data[N - 2*k]);
    // A_k = C_k + 1^(k/N) D_k:
    data[2*k] = Ck_re;  // A_k <-- C_k
    data[2*k+1] = Ck_im;
    // now A_k += D_k 1^(k/N)
    ComplexAddProduct(Dk_re, Dk_im, kN_re, kN_im, &(data[2*k]), &(data[2*k+1]));

    MatrixIndexT kdash = N2 - k;
    if (kdash != k) {
      // Next we handle the index k' = N/2 - k.  This is necessary
      // to do now, to avoid invalidating data that we will later need.
      // The quantities C_{k'} and D_{k'} are just the conjugates of C_k
      // and D_k, so the equations are simple modifications of the above,
      // replacing Ck_im and Dk_im with their negatives.
      data[2*kdash] = Ck_re;  // A_k' <-- C_k'
      data[2*kdash+1] = -Ck_im;
      // now A_k' += D_k' 1^(k'/N)
      // We use 1^(k'/N) = 1^((N/2 - k) / N) = 1^(1/2) 1^(-k/N) = -1 * (1^(k/N))^*
      // so it's the same as 1^(k/N) but with the real part negated.
      ComplexAddProduct(Dk_re, -Dk_im, -kN_re, kN_im, &(data[2*kdash]), &(data[2*kdash+1]));
    }
  }

  {  // Now handle k = 0.
    // In simple terms: after the complex fft, data[0] becomes the sum of real
    // parts input[0], input[2]... and data[1] becomes the sum of imaginary
    // pats input[1], input[3]...
    // "zeroth" [A_0] is just the sum of input[0]+input[1]+input[2]..
    // and "n2th" [A_{N/2}] is input[0]-input[1]+input[2]... .
    BaseFloat zeroth = data[0] + data[1],
        n2th = data[0] - data[1];
    data[0] = zeroth;
    data[1] = n2th;
    if (!forward) {
      data[0] /= 2;
      data[1] /= 2;
    }
  }

  if (!forward) {
    ComplexFft(v, false);
    v->Scale(2.0);  // This is so we get a factor of N increase, rather than N/2 which we would
    // otherwise get from [ComplexFft, forward] + [ComplexFft, backward] in dimension N/2.
    // It's for consistency with our normal FFT convensions.
  }
}


void ComputeDctMatrix(Matrix *M) {
  MatrixIndexT K = M->NumRows();
  MatrixIndexT N = M->NumCols();

  assert(K > 0);
  assert(N > 0);
  BaseFloat normalizer = std::sqrt(1.0 / static_cast<BaseFloat>(N));  // normalizer for
  // X_0.
  for (MatrixIndexT j = 0; j < N; j++) (*M)(0, j) = normalizer;
  normalizer = std::sqrt(2.0 / static_cast<BaseFloat>(N));  // normalizer for other
   // elements.
  for (MatrixIndexT k = 1; k < K; k++)
    for (MatrixIndexT n = 0; n < N; n++)
      (*M)(k, n) = normalizer
          * std::cos( static_cast<double>(M_PI)/N * (n + 0.5) * k );
}


void ComputeLifterCoeffs(BaseFloat Q, Vector *coeffs) {
  // Compute liftering coefficients (scaling on cepstral coeffs)
  // coeffs are numbered slightly differently from HTK: the zeroth
  // index is C0, which is not affected.
  for (int32 i = 0; i < coeffs->Dim(); i++)
    (*coeffs)(i) = 1.0 + 0.5 * Q * sin (M_PI * i / Q);
}


BaseFloat CosDistance(Vector& v1, Vector& v2){

	assert(v1.Dim() == v2.Dim());
	BaseFloat M_v1 = v1.Norm(2.0);
	BaseFloat M_v2 = v2.Norm(2.0);
	BaseFloat Dot = VecVec(v1, v2);

	return Dot/(M_v1*M_v2);
}


void HouseBackward(MatrixIndexT dim, const BaseFloat *x, BaseFloat *v, BaseFloat *beta) {
  assert(dim > 0);
  // To avoid overflow, we first compute the max of x_ (or
  // one if that's zero, and we'll replace "x" by x/max(x_i)
  // below.  The householder vector is anyway invariant to
  // the magnitude of x.  We could actually avoid this extra loop
  // over x if we wanted to be a bit smarter, but anyway this
  // doesn't dominate the O(N) performance of the algorithm.
  BaseFloat s; // s is a scale on x.
  {
    BaseFloat max_x = std::numeric_limits<BaseFloat>::min();
    for (MatrixIndexT i = 0; i < dim; i++)
      max_x = std::max(max_x, (x[i] < 0 ? -x[i] : x[i]));
    s = 1.0 / max_x;
  }
  BaseFloat sigma = 0.0;
  v[dim-1] = 1.0;
  for (MatrixIndexT i = 0; i + 1  < dim; i++) {
    sigma += (x[i] * s) * (x[i] * s);
    v[i] = x[i] * s;
  }
  assert(ISFINITE(sigma) && "Tridiagonalizing matrix that is too large or has NaNs.");
  if (sigma == 0.0) *beta = 0.0;
  else {
    BaseFloat x1 = x[dim-1] * s, mu = std::sqrt(x1 * x1 + sigma);
    if (x1 <= 0) {
      v[dim-1] = x1 - mu;
    } else {
      v[dim-1] = -sigma / (x1 + mu);
      assert(ISFINITE(v[dim-1]));
    }
    BaseFloat v1 = v[dim-1];
    BaseFloat v1sq = v1 * v1;
    *beta = 2 * v1sq / (sigma + v1sq);
    BaseFloat inv_v1 = 1.0 / v1;
    if (ISINF(inv_v1)) {
      // can happen if v1 is denormal.
      assert(v1 == v1 && v1 != 0.0);
      for (MatrixIndexT i = 0; i < dim; i++) v[i] /= v1;
    } else {
      cblas_Xscal(dim, inv_v1, v, 1);
    }
    if (ISNAN(inv_v1)) {
      std::cerr << "NaN encountered in HouseBackward"<<std::endl;
    }
  }
}

void Givens(BaseFloat a, BaseFloat b, BaseFloat *c, BaseFloat *s) {
  if (b == 0) {
    *c = 1;
    *s = 0;
  } else {
    if (std::abs(b) > std::abs(a)) {
      BaseFloat tau = -a / b;
      *s = 1 / std::sqrt(1 + tau*tau);
      *c = *s * tau;
    } else {
      BaseFloat tau = -b / a;
      *c = 1 / std::sqrt(1 + tau*tau);
      *s = *c * tau;
    }
  }
}

void QrStep(MatrixIndexT n,
            BaseFloat *diag,
            BaseFloat *off_diag,
            Matrix *Q) {
  assert(n >= 2);
  // below, "scale" could be any number; we introduce it to keep the
  // floating point quantities within a good range.
  BaseFloat   d = (diag[n-2] - diag[n-1]) / 2.0,
      t = off_diag[n-2],
      inv_scale = std::max(std::max(std::abs(d), std::abs(t)),
                           std::numeric_limits<BaseFloat>::min()),
      scale = 1.0 / inv_scale,
      d_scaled = d * scale,
      off_diag_n2_scaled = off_diag[n-2] * scale,
      t2_n_n1_scaled = off_diag_n2_scaled * off_diag_n2_scaled,
      sgn_d = (d > 0.0 ? 1.0 : -1.0),
      mu = diag[n-1] - inv_scale * t2_n_n1_scaled /
      (d_scaled + sgn_d * std::sqrt(d_scaled * d_scaled + t2_n_n1_scaled)),
      x = diag[0] - mu,
      z = off_diag[0];
  assert(ISFINITE(x));
  BaseFloat *Qdata = (Q == NULL ? NULL : Q->Data());
  MatrixIndexT Qstride = (Q == NULL ? 0 : Q->NumCols()),
      Qcols = (Q == NULL ? 0 : Q->NumCols());
  for (MatrixIndexT k = 0; k < n-1; k++) {
    BaseFloat c, s;
    Givens(x, z, &c, &s);
    // Rotate dimensions k and k+1 with the Givens matrix G, as
    // T <== G^T T G.
    // In 2d, a Givens matrix is [ c s; -s c ].  Forget about
    // the dimension-indexing issues and assume we have a 2x2
    // symmetric matrix [ p q ; q r ]
    // We ask our friends at Wolfram Alpha about
    // { { c, -s}, {s, c} } * { {p, q}, {q, r} } * { { c, s}, {-s, c} }
    // Interpreting the result as [ p', q' ; q', r ]
    //    p' = c (c p - s q) - s (c q - s r)
    //    q' = s (c p - s q) + c (c q - s r)
    //    r' = s (s p + c q) + c (s q + c r)
    BaseFloat p = diag[k], q = off_diag[k], r = diag[k+1];
    // p is element k,k; r is element k+1,k+1; q is element k,k+1 or k+1,k.
    // We'll let the compiler optimize this.
    diag[k] = c * (c*p - s*q) - s * (c*q - s*r);
    off_diag[k] = s * (c*p - s*q) + c * (c*q - s*r);
    diag[k+1] = s * (s*p + c*q) + c * (s*q + c*r);

    // We also have some other elements to think of that
    // got rotated in a simpler way: if k>0,
    // then element (k, k-1) and (k+1, k-1) get rotated.  Here,
    // element k+1, k-1 will be present as z; it's the out-of-band
    // element that we remembered from last time.  This is
    // on the left as it's the row indexes that differ, so think of
    // this as being premultiplied by G^T.  In fact we're multiplying
    // T by in some sense the opposite/transpose of the Givens rotation.
    if (k > 0) { // Note, in rotations, going backward, (x,y) -> ((cx - sy), (sx + cy))
      BaseFloat &elem_k_km1 = off_diag[k-1],
          elem_kp1_km1 = z; // , tmp = elem_k_km1;
      elem_k_km1 = c*elem_k_km1 - s*elem_kp1_km1;
      // The next line will set elem_kp1_km1 to zero and we'll never access this
      // value, so we comment it out.
      // elem_kp1_km1 = s*tmp + c*elem_kp1_km1;
    }
    if (Q != NULL)
      cblas_Xrot(Qcols, Qdata + k*Qstride, 1,
                 Qdata + (k+1)*Qstride, 1, c, -s);
    if (k < n-2) {
      // Next is the elements (k+2, k) and (k+2, k-1), to be rotated, again
      // backwards.
      BaseFloat &elem_kp2_k = z,
          &elem_kp2_kp1 = off_diag[k+1];
      // Note: elem_kp2_k == z would start off as zero because it's
       // two off the diagonal, and not been touched yet.  Therefore
      // we eliminate it in expressions below, commenting it out.
      // If we didn't do this we should set it to zero first.
      elem_kp2_k =  - s * elem_kp2_kp1; // + c*elem_kp2_k
      elem_kp2_kp1 =  c * elem_kp2_kp1; // + s*elem_kp2_k (original value).
      // The next part is from the algorithm they describe: x = t_{k+1,k}
      x = off_diag[k];
    }
  }
}


void QrInternal(MatrixIndexT n,
                BaseFloat *diag,
                BaseFloat *off_diag,
                Matrix *Q) {
  assert(Q == NULL || Q->NumCols() == n); // We may
  // later relax the condition that Q->NumCols() == n.

  MatrixIndexT counter = 0, max_iters = 500 + 4*n, // Should never take this many iters.
      large_iters = 100 + 2*n;
  BaseFloat epsilon = (pow(2.0, sizeof(BaseFloat) == 4 ? -23.0 : -52.0));

  for (; counter < max_iters; counter++) { // this takes the place of "until
                                           // q=n"... we'll break out of the
                                           // loop when we converge.
    if (counter == large_iters ||
        (counter > large_iters && (counter - large_iters) % 50 == 0)) {
      std::cout << " Warning! Took " << counter
                 << " iterations in QR (dim is " << n << "), doubling epsilon.";
      Vector d(n);
      d.CopyFromData(diag, n);
      Vector o(n-1);
      o.CopyFromData(off_diag, n-1);
      epsilon *= 2.0;
    }
    for (MatrixIndexT i = 0; i+1 < n; i++) {
      if (std::abs(off_diag[i]) <= epsilon *
          (std::abs(diag[i]) + std::abs(diag[i+1])))
        off_diag[i] = 0.0;
    }
    // The next code works out p, q, and npq which is n - p - q.
    // For the definitions of q and p, see Golub and Van Loan; we
    // partition the n dims into pieces of size (p, n-p-q, q) where
    // the part of size q is diagonal and the part of size n-p-p is
    // "unreduced", i.e. has no zero off-diagonal elements.
    MatrixIndexT q = 0;
    // Note: below, "n-q < 2" should more clearly be "n-2-q < 0", but that
    // causes problems if MatrixIndexT is unsigned.
    while (q < n && (n-q < 2 || off_diag[n-2-q] == 0.0))
      q++;
    if (q == n) break; // we're done.  It's diagonal.
    assert(n - q >= 2);
    MatrixIndexT npq = 2; // Value of n - p - q, where n - p - q must be
    // unreduced.  This is the size of "middle" band of elements.  If q != n,
    // we must have hit a nonzero off-diag element, so the size of this
    // band must be at least two.
    while (npq + q < n && (n-q-npq-1 < 0 || off_diag[n-q-npq-1] != 0.0))
      npq++;
    MatrixIndexT p = n - q - npq;
    { // Checks.
      for (MatrixIndexT i = 0; i+1 < npq; i++)
        assert(off_diag[p + i] != 0.0);
      for (MatrixIndexT i = 0; i+1 < q; i++)
        assert(off_diag[p + npq - 1 + i] == 0.0);
      if (p > 1) // Something must have stopped npq from growing further..
        assert(off_diag[p-1] == 0.0); // so last off-diag elem in
      // group of size p must be zero.
    }

    if (Q != NULL) {
      // Do one QR step on the middle part of Q only.
      // Qpart will be a subset of the rows of Q.
      Matrix Qpart(npq, Q->NumCols());
      Qpart.CopyFromMat(*Q, p, npq, 0, Q->NumCols());
      QrStep(npq, diag + p, off_diag + p, &Qpart);
    } else {
      QrStep(npq, diag + p, off_diag + p,
             static_cast<Matrix*>(NULL));
    }
  }
  if (counter == max_iters) {
    std::clog << "Failure to converge in QR algorithm. "
               << "Exiting with partial output.";
  }
}


BaseFloat VecVec(const Vector &ra,
            const Vector &rb) {
  MatrixIndexT adim = ra.Dim();
  assert(adim == rb.Dim());
  const BaseFloat *a_data = ra.Data();
  const BaseFloat *b_data = rb.Data();
  BaseFloat sum = 0.0;
  for (MatrixIndexT i = 0; i < adim; i++)
    sum += a_data[i]*b_data[i];
  return sum;
}
