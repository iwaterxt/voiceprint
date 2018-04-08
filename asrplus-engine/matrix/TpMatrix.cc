/*
 * TpMatrix.cc
 *
 *  Created on: Feb 18, 2018
 *      Author: tao
 */

#include "matrix/TpMatrix.h"
#include "base/common.h"
#include <stdexcept>
#include <string.h>



void TpMatrix::Resize(MatrixIndexT nRows){

	  if (data_ != NULL) {
		  delete []data_;
		  rows_ = 0;
	  }
	  Init(nRows);
	  this->SetZero();

}

void TpMatrix::SetZero(){
	memset(data_, 0, SizeInBytes());
}

void TpMatrix::Cholesky(const SpMatrix &orig) {
  assert(orig.NumRows() == this->NumRows());
  MatrixIndexT n = this->NumRows();
  this->SetZero();
  BaseFloat *data = this->data_, *jdata = data;  // start of j'th row of matrix.
  const BaseFloat *orig_jdata = orig.Data(); // start of j'th row of matrix.
  for (MatrixIndexT j = 0; j < n; j++, jdata += j, orig_jdata += j) {
    BaseFloat *kdata = data; // start of k'th row of matrix.
    BaseFloat d(0.0);
    for (MatrixIndexT k = 0; k < j; k++, kdata += k) {
      BaseFloat s = cblas_Xdot(k, kdata, 1, jdata, 1);
      // (*this)(j, k) = s = (orig(j, k) - s)/(*this)(k, k);
      jdata[k] = s = (orig_jdata[k] - s)/kdata[k];
      d = d + s*s;
    }
    // d = orig(j, j) - d;
    d = orig_jdata[j] - d;

    if (d >= 0.0) {
      // (*this)(j, j) = std::sqrt(d);
      jdata[j] = std::sqrt(d);
    } else {
      std::cout << "Cholesky decomposition failed. Maybe matrix "
          "is not positive definite. Throwing error";
      throw std::runtime_error("Cholesky decomposition failed.");
    }
  }
}


void TpMatrix::Init(MatrixIndexT r) {
  if (r == 0) {
    rows_ = 0;
    data_ = NULL;
    return;
  }
  size_t size = ((static_cast<size_t>(r) * static_cast<size_t>(r + 1)) / 2);

  if (static_cast<size_t>(static_cast<MatrixIndexT>(size)) != size) {
    std::cout << " Warning Allocating packed matrix whose full dimension does not fit "
               << "in MatrixIndexT: not all code is tested for this case.";
  }

  data_ = new BaseFloat[size];
  if ( data_ != NULL) {
    this->rows_ = r;
  } else {
    throw std::bad_alloc();
  }
}
