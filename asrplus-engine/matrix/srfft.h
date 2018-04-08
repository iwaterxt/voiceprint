/*
 * srfft.h
 *
 *  Created on: Oct 10, 2017
 *      Author: tao
 */

#ifndef SRFFT_H_
#define SRFFT_H_

#include "base/common.h"
#include <vector>


class SplitRadixComplexFft {
 public:
  typedef MatrixIndexT Integer;

  // N is the number of complex points (must be a power of two, or this
  // will crash).  Note that the constructor does some work so it's best to
  // initialize the object once and do the computation many times.
  SplitRadixComplexFft(Integer N);

  // Copy constructor
  SplitRadixComplexFft(const SplitRadixComplexFft &other);

  // Does the FFT computation, given pointers to the real and
  // imaginary parts.  If "forward", do the forward FFT; else
  // do the inverse FFT (without the 1/N factor).
  // xr and xi are pointers to zero-based arrays of size N,
  // containing the real and imaginary parts
  // respectively.
  void Compute(BaseFloat *xr, BaseFloat *xi, bool forward) const;

  // This version of Compute takes a single array of size N*2,
  // containing [ r0 im0 r1 im1 ... ].  Otherwise its behavior is  the
  // same as the version above.
  void Compute(BaseFloat *x, bool forward);


  // This version of Compute is const; it operates on an array of size N*2
  // containing [ r0 im0 r1 im1 ... ], but it uses the argument "temp_buffer" as
  // temporary storage instead of a class-member variable.  It will allocate it if
  // needed.
  void Compute(BaseFloat *x, bool forward, std::vector<BaseFloat> *temp_buffer) const;

  ~SplitRadixComplexFft();

 protected:
  // temp_buffer_ is allocated only if someone calls Compute with only one Real*
  // argument and we need a temporary buffer while creating interleaved data.
  std::vector<BaseFloat> temp_buffer_;
 private:
  void ComputeTables();
  void ComputeRecursive(BaseFloat *xr, BaseFloat *xi, Integer logn) const;
  void BitReversePermute(BaseFloat *x, Integer logn) const;

  Integer N_;
  Integer logn_;  // log(N)

  Integer *brseed_;
  // brseed is Evans' seed table, ref:  (Ref: D. M. W.
  // Evans, "An improved digit-reversal permutation algorithm ...",
  // IEEE Trans. ASSP, Aug. 1987, pp. 1120-1125).
  BaseFloat **tab_;       // Tables of butterfly coefficients.

  // Disallow assignment.
  SplitRadixComplexFft &operator =(const SplitRadixComplexFft &other);
};


class SplitRadixRealFft: private SplitRadixComplexFft {
 public:
  SplitRadixRealFft(MatrixIndexT N):  // will fail unless N>=4 and N is a power of 2.
      SplitRadixComplexFft (N/2), N_(N) { }

  // Copy constructor
  SplitRadixRealFft(const SplitRadixRealFft &other):
      SplitRadixComplexFft(other), N_(other.N_) { }
  /// If forward == true, this function transforms from a sequence of N real points to its complex fourier
  /// transform; otherwise it goes in the reverse direction.  If you call it
  /// in the forward and then reverse direction and multiply by 1.0/N, you
  /// will get back the original data.
  /// The interpretation of the complex-FFT data is as follows: the array
  /// is a sequence of complex numbers C_n of length N/2 with (real, im) format,
  /// i.e. [real0, real_{N/2}, real1, im1, real2, im2, real3, im3, ...].
  void Compute(BaseFloat *x, bool forward);

  /// This is as the other Compute() function, but it is a const version that
  /// uses a user-supplied buffer.
  void Compute(BaseFloat *x, bool forward, std::vector<BaseFloat> *temp_buffer) const;

 private:
  // Disallow assignment.
  SplitRadixRealFft &operator =(const SplitRadixRealFft &other);
  uint32 N_;
};



#endif /* SRFFT_H_ */
