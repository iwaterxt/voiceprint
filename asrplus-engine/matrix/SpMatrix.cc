/*
 * SpMatrix.cc
 *
 *  Created on: Feb 2, 2018
 *      Author: tao
 */

#include "base/common.h"
#include "matrix/SpMatrix.h"
#include "matrix/TpMatrix.h"
#include "util/math.h"
#include "base/io.h"
#include <assert.h>
#include <string.h>
#include <limits>
#include <vector>

bool SpMatrix::IsTridiagonal(BaseFloat cutoff) const {
  MatrixIndexT R = this->NumRows();
  BaseFloat max_abs_2diag = 0.0, max_abs_offdiag = 0.0;
  for (MatrixIndexT i = 0; i < R; i++)
    for (MatrixIndexT j = 0; j <= i; j++) {
      if (j+1 < i)
        max_abs_offdiag = std::max(max_abs_offdiag,
                                   std::abs((*this)(i, j)));
      else
        max_abs_2diag = std::max(max_abs_2diag,
                                 std::abs((*this)(i, j)));
    }
  return (max_abs_offdiag <= cutoff * max_abs_2diag);
}


void SpMatrix::Resize(MatrixIndexT nRows){

	  if (data_ != NULL) {
		  delete []data_;
		  rows_ = 0;
	  }
	  Init(nRows);
	  this->SetZero();

}

void SpMatrix::Invert(BaseFloat *logdet, BaseFloat *det_sign, bool need_inverse) {
	 int32   result;
	  int32   rows = static_cast<int>(this->rows_);
	  int32*  p_ipiv = new int32[rows];
	  BaseFloat *p_work;  // workspace for the lapack function
	  void *temp;
	  if ((p_work = static_cast<BaseFloat*>(
	          MEMALIGN(16, sizeof(BaseFloat) * rows, &temp))) == NULL) {
	    delete[] p_ipiv;
	    throw std::bad_alloc();
	  }

	  memset(p_work, 0, sizeof(BaseFloat) * rows); // gets rid of a probably
	  // spurious Valgrind warning about jumps depending upon uninitialized values.



	  // NOTE: Even though "U" is for upper, lapack assumes column-wise storage
	  // of the data. We have a row-wise storage, therefore, we need to "invert"
	  clapack_Xsptrf(&rows, this->data_, p_ipiv, &result);


	  assert(result >= 0 && "Call to CLAPACK ssptrf_ called with wrong arguments");

	  if (result > 0) {  // Singular...
	    if (det_sign) *det_sign = 0;
	    if (logdet) *logdet = -std::numeric_limits<BaseFloat>::infinity();
	    if (need_inverse) std::cerr << "CLAPACK stptrf_ : factorization failed";
	  } else {  // Not singular.. compute log-determinant if needed.
	    if (logdet != NULL || det_sign != NULL) {
	      BaseFloat prod = 1.0, log_prod = 0.0;
	      int sign = 1;
	      for (int i = 0; i < (int)this->rows_; i++) {
	        if (p_ipiv[i] > 0) {  // not a 2x2 block...
	          // if (p_ipiv[i] != i+1) sign *= -1;  // row swap.
	          BaseFloat diag = (*this)(i, i);
	          prod *= diag;
	        } else {  // negative: 2x2 block. [we are in first of the two].
	          i++;  // skip over the first of the pair.
	          // each 2x2 block...
	          BaseFloat diag1 = (*this)(i, i), diag2 = (*this)(i-1, i-1),
	              offdiag = (*this)(i, i-1);
	          BaseFloat thisdet = diag1*diag2 - offdiag*offdiag;
	          // thisdet == determinant of 2x2 block.
	          // The following line is more complex than it looks: there are 2 offsets of
	          // 1 that cancel.
	          prod *= thisdet;
	        }
	        if (i == (int)(this->rows_-1) || fabs(prod) < 1.0e-10 || fabs(prod) > 1.0e+10) {
	          if (prod < 0) { prod = -prod; sign *= -1; }
	          log_prod += Log(std::abs(prod));
	          prod = 1.0;
	        }
	      }
	      if (logdet != NULL) *logdet = log_prod;
	      if (det_sign != NULL) *det_sign = sign;
	    }
	  }
	  if (!need_inverse) {
	    delete [] p_ipiv;
	    MEMALIGN_FREE(p_work);
	    return;  // Don't need what is computed next.
	  }
	  // NOTE: Even though "U" is for upper, lapack assumes column-wise storage
	  // of the data. We have a row-wise storage, therefore, we need to "invert"
	  clapack_Xsptri(&rows, this->data_, p_ipiv, p_work, &result);

	  assert(result >=0 &&
	               "Call to CLAPACK ssptri_ called with wrong arguments");

	  if (result != 0) {
	    std::cerr << "CLAPACK ssptrf_ : Matrix is singular";
	  }

	  delete [] p_ipiv;
	  MEMALIGN_FREE(p_work);
}

void SpMatrix::SetZero(){
	memset(data_, 0, SizeInBytes());
}

BaseFloat SpMatrix::LogPosDefDet() const {
  TpMatrix chol(this->NumRows());
  double det = 0.0;
  double diag;
  chol.Cholesky(*this);  // Will throw exception if not +ve definite!

  for (MatrixIndexT i = 0; i < this->NumRows(); i++) {
    diag = static_cast<double>(chol(i, i));
    det += Log(diag);
  }
  return static_cast<BaseFloat>(2*det);
}


void SpMatrix::Scale(BaseFloat alpha) {
  size_t nr = rows_,
      sz = (nr * (nr + 1)) / 2;
  cblas_Xscal(sz, alpha, data_, 1);
}


void SpMatrix::AddToDiag(BaseFloat r) {
  BaseFloat *ptr = data_;
  for (MatrixIndexT i = 2; i <= rows_+1; i++) {
    *ptr += r;
    ptr += i;
  }
}


void SpMatrix::CopyFromSp(const SpMatrix &orig) {
  assert(NumRows() == orig.NumRows());
  if (sizeof(BaseFloat) == sizeof(BaseFloat)) {
    memcpy(data_, orig.Data(), SizeInBytes());
  } else {
    BaseFloat *dst = data_;
    const BaseFloat *src = orig.Data();
    size_t nr = NumRows(),
        size = (nr * (nr + 1)) / 2;
    for (size_t i = 0; i < size; i++, dst++, src++)
      *dst = *src;
  }
}

void SpMatrix::CopyFromVec(const Vector& v){
	int32 nr = static_cast<size_t>(rows_);
	int32 size = (nr*(nr+1))/2;
	const BaseFloat* data_v = v.Data();
	assert(size == v.Dim());
	for(int32 i = 0; i < size; i++)
		data_[i] = data_v[i];
}

void SpMatrix::AddMat2Vec(const BaseFloat alpha, const Matrix &M,
                MatrixTransposeType transM, const Vector &v,
                const BaseFloat beta){

	  this->Scale(beta);
	  assert((transM == kNoTrans && this->NumRows() == M.NumRows() &&
	                M.NumCols() == v.Dim()) ||
	               (transM == kTrans && this->NumRows() == M.NumCols() &&
	                M.NumRows() == v.Dim()));

	  if (transM == kNoTrans) {
	    const BaseFloat *Mdata = M.Data(), *vdata = v.Data();
	    BaseFloat *data = this->data_;
	    MatrixIndexT dim = this->NumRows(), mcols = M.NumCols(),
	        mstride = M.NumCols();
	    for (MatrixIndexT col = 0; col < mcols; col++, vdata++, Mdata += 1)
	      cblas_Xspr(dim, *vdata*alpha, Mdata, mstride, data);
	  } else {
	    const BaseFloat *Mdata = M.Data(), *vdata = v.Data();
	    BaseFloat *data = this->data_;
	    MatrixIndexT dim = this->NumRows(), mrows = M.NumRows(),
	        mstride = M.NumCols();
	    for (MatrixIndexT row = 0; row < mrows; row++, vdata++, Mdata += mstride)
	      cblas_Xspr(dim, *vdata*alpha, Mdata, 1, data);
	  }
}

void SpMatrix::AddMat2Sp(
    const BaseFloat alpha, const Matrix &M,
    MatrixTransposeType transM, const SpMatrix &A, const BaseFloat beta) {
  if (transM == kNoTrans) {
    assert(M.NumCols() == A.NumRows() && M.NumRows() == this->rows_);
  } else {
    assert(M.NumRows() == A.NumRows() && M.NumCols() == this->rows_);
  }
  Vector tmp_vec(A.NumRows());
  BaseFloat *tmp_vec_data = tmp_vec.Data();
  SpMatrix tmp_A;
  const BaseFloat *p_A_data = A.Data();
  BaseFloat *p_row_data = this->Data();
  MatrixIndexT M_other_dim = (transM == kNoTrans ? M.NumCols() : M.NumRows()),
      M_same_dim = (transM == kNoTrans ? M.NumRows() : M.NumCols()),
      M_stride = M.NumCols(), dim = this->NumRows();
  assert(M_same_dim == dim);

  const BaseFloat *M_data = M.Data();

  if (this->Data() <= A.Data() + A.SizeInBytes() &&
      this->Data() + this->SizeInBytes() >= A.Data()) {
    // Matrices A and *this overlap. Make copy of A
    tmp_A.Resize(A.NumRows());
    tmp_A.CopyFromSp(A);
    p_A_data = tmp_A.Data();
  }

  if (transM == kNoTrans) {
    for (MatrixIndexT r = 0; r < dim; r++, p_row_data += r) {
      cblas_Xspmv(A.NumRows(), 1.0, p_A_data, M.RowData(r), 1, 0.0, tmp_vec_data, 1);
      cblas_Xgemv(transM, r+1, M_other_dim, alpha, M_data, M_stride,
                  tmp_vec_data, 1, beta, p_row_data, 1);
    }
  } else {
    for (MatrixIndexT r = 0; r < dim; r++, p_row_data += r) {
      cblas_Xspmv(A.NumRows(), 1.0, p_A_data, M.Data() + r, M.NumCols(), 0.0, tmp_vec_data, 1);
      cblas_Xgemv(transM, M_other_dim, r+1, alpha, M_data, M_stride, tmp_vec_data, 1, beta, p_row_data, 1);
    }
  }
}

void SpMatrix::AddVec2(const BaseFloat alpha, const Vector &v) {

  assert(v.Dim() == this->NumRows());
  BaseFloat *data = this->data_;
  const BaseFloat *v_data = v.Data();
  MatrixIndexT nr = this->rows_;
  for (MatrixIndexT i = 0; i < nr; i++)
    for (MatrixIndexT j = 0; j <= i; j++, data++)
      *data += alpha * v_data[i] * v_data[j];
}

void SpMatrix::AddSp(const BaseFloat alpha, const Matrix &rMa) {
  assert(rows_ == rMa.NumRows());
  size_t nr = rows_,
      sz = (nr * (nr + 1)) / 2;
  cblas_Xaxpy(sz, alpha, rMa.Data(), 1, data_, 1);
}

void SpMatrix::Tridiagonalize(Matrix *Q) {
  MatrixIndexT n = this->NumRows();
  assert(Q == NULL || (Q->NumRows() == n &&
                             Q->NumCols() == n));
  if (Q != NULL) Q->SetUnit();
  BaseFloat *data = this->Data();
  BaseFloat *qdata = (Q == NULL ? NULL : Q->Data());
  MatrixIndexT qstride = (Q == NULL ? 0 : Q->NumCols());
  Vector tmp_v(n-1), tmp_p(n);
  BaseFloat beta, *v = tmp_v.Data(), *p = tmp_p.Data(), *w = p, *x = p;
  for (MatrixIndexT k = n-1; k >= 2; k--) {
    MatrixIndexT ksize = ((k+1)*k)/2;
    // ksize is the packed size of the lower-triangular matrix of size k,
    // which is the size of "all rows previous to this one."
    BaseFloat *Arow = data + ksize; // In Golub+Van Loan it was A(k+1:n, k), we
    // have Arow = A(k, 0:k-1).
    HouseBackward(k, Arow, v, &beta); // sets v and beta.
    cblas_Xspmv(k, beta, data, v, 1, 0.0, p, 1); // p = beta * A(0:k-1,0:k-1) v
    BaseFloat minus_half_beta_pv = -0.5 * beta * cblas_Xdot(k, p, 1, v, 1);
    cblas_Xaxpy(k, minus_half_beta_pv, v, 1, w, 1); // w = p - (beta p^T v/2) v;
    // this relies on the fact that w and p are the same pointer.
    // We're doing A(k, k-1) = ||Arow||.  It happens that this element
    // is indexed at ksize + k - 1 in the packed lower-triangular format.
    data[ksize + k - 1] = std::sqrt(cblas_Xdot(k, Arow, 1, Arow, 1));
    for (MatrixIndexT i = 0; i + 1 < k; i++)
      data[ksize + i] = 0; // This is not in Golub and Van Loan but is
    // necessary if we're not using parts of A to store the Householder
    // vectors.
    // We're doing A(0:k-1,0:k-1) -= (v w' + w v')
    cblas_Xspr2(k, -1.0, v, 1, w, 1, data);
    if (Q != NULL) { // C.f. Golub, Q is H_1 .. H_n-2... in this
      // case we apply them in the opposite order so it's H_n-1 .. H_1,
      // but also Q is transposed so we really have Q = H_1 .. H_n-1.
      // It's a double negative.
      // Anyway, we left-multiply Q by each one.  The H_n would each be
      // diag(I + beta v v', I) but we don't ever touch the last dims.
      // We do (in Matlab notation):
      // Q(0:k-1,:) = (I - beta v v') * Q, i.e.:
      // Q(:,0:i-1) += -beta v (v' Q(:,0:k-1)v .. let x = -beta Q(0:k-1,:)^T v.
      cblas_Xgemv(kTrans, k, n, -beta, qdata, qstride, v, 1, 0.0, x, 1);
      // now x = -beta Q(:,0:k-1) v.
      // The next line does: Q(:,0:k-1) += v x'.
      cblas_Xger(k, n, 1.0, v, 1, x, 1, qdata, qstride);
    }
  }
}

void SpMatrix::Qr(Matrix *Q) {
  assert(this->IsTridiagonal());
  // We envisage that Q would be square but we don't check for this,
  // as there are situations where you might not want this.
  assert(Q == NULL || Q->NumRows() == this->NumRows());
  // Note: the first couple of lines of the algorithm they give would be done
  // outside of this function, by calling Tridiagonalize().

  MatrixIndexT n = this->NumRows();
  Vector diag(n), off_diag(n-1);
  for (MatrixIndexT i = 0; i < n; i++) {
    diag(i) = (*this)(i, i);
    if (i > 0) off_diag(i-1) = (*this)(i, i-1);
  }
  QrInternal(n, diag.Data(), off_diag.Data(), Q);
  // Now set *this to the value represented by diag and off_diag.
  this->SetZero();
  for (MatrixIndexT i = 0; i < n; i++) {
    (*this)(i, i) = diag(i);
    if (i > 0) (*this)(i, i-1) = off_diag(i-1);
  }
}


void SpMatrix::Eig(Vector *s, Matrix *P) const {
  MatrixIndexT dim = this->NumRows();
  assert(s->Dim() == dim);
  assert(P == NULL || (P->NumRows() == dim && P->NumCols() == dim));

  SpMatrix A(*this); // Copy *this, since the tridiagonalization

  A.Tridiagonalize(P); // Tridiagonalizes.
  A.Qr(P); // Diagonalizes.
  if(P) P->Transpose();
  s->CopyDiagFromSp(A);
}

void SpMatrix::Init(MatrixIndexT r) {
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


void SpMatrix::Read(std::istream& is, bool binary) {


  std::ostringstream specific_error;
  MatrixIndexT pos_at_start = is.tellg();

  const char *my_token =  (sizeof(BaseFloat) == 4 ? "FP" : "DP");
  const char *new_format_token = "[";
  bool is_new_format = false;//added by hxu

  int32 size;
  MatrixIndexT num_elems;

  std::string token;
  ReadToken(is, binary, &token);
  if (token != my_token) {
    if(token != new_format_token) {
      specific_error << ": Expected token " << my_token << ", got " << token;
      goto bad;
    }
    //new format it is
    is_new_format = true;
  }
  if(!is_new_format) {
    ReadBasicType(is, binary, &size);  // throws on error.
    if ((MatrixIndexT)size != this->NumRows()) {
      assert(size>=0);
      this->Resize(size);
    }
    num_elems = ((size+1)*(MatrixIndexT)size)/2;
    if (!binary) {
      for (MatrixIndexT i = 0; i < num_elems; i++) {
        ReadBasicType(is, false, data_+i);  // will throw on error.
      }
    } else {
      if (num_elems)
        is.read(reinterpret_cast<char*>(data_), sizeof(BaseFloat)*num_elems);
    }
    if (is.fail()) goto bad;
    return;
  }
  else {
    std::vector<BaseFloat> data;
    while(1) {
      int32 num_lines = 0;
      int i = is.peek();
      if (i == -1) { specific_error << "Got EOF while reading matrix data"; goto bad; }
      else if (static_cast<char>(i) == ']') {  // Finished reading matrix.
        is.get();  // eat the "]".
        i = is.peek();
        if (static_cast<char>(i) == '\r') {
          is.get();
          is.get();  // get \r\n (must eat what we wrote)
        }// I don't actually understand what it's doing here
        else if (static_cast<char>(i) == '\n') { is.get(); } // get \n (must eat what we wrote)

        if (is.fail()) {
          std::cout << "After end of matrix data, read error."<<std::endl;
          // we got the data we needed, so just warn for this error.
        }
        //now process the data:
        num_lines = int32(sqrt(data.size()*2));

        assert(data.size() == num_lines*(num_lines+1)/2);

        this->Resize(num_lines);

        //std::cout<<data.size()<<' '<<num_lines<<'\n';

        for(int32 i = 0; i < data.size(); i++) {
          data_[i] = data[i];
        }
        return;
        //std::cout<<"here!!!!!hxu!!!!!"<<std::endl;
      }
      else if ( (i >= '0' && i <= '9') || i == '-' ) {  // A number...
        BaseFloat r;
        is >> r;
        if (is.fail()) {
          specific_error << "Stream failure/EOF while reading matrix data.";
          goto bad;
        }
        data.push_back(r);
      }
      else if (isspace(i)) {
        is.get();  // eat the space and do nothing.
      } else {  // NaN or inf or error.
        std::string str;
        is >> str;
        if (!strcasecmp(str.c_str(), "inf") ||
            !strcasecmp(str.c_str(), "infinity")) {
          data.push_back(std::numeric_limits<BaseFloat>::infinity());
          std::cout << "Reading infinite value into matrix."<<std::endl;
        } else if (!strcasecmp(str.c_str(), "nan")) {
          data.push_back(std::numeric_limits<BaseFloat>::quiet_NaN());
          std::cout << "Reading NaN value into matrix."<<std::endl;
        } else {
          specific_error << "Expecting numeric matrix data, got " << str;
          goto bad;
        }
      }
    }
  }
bad:
  std::cerr << "Failed to read packed matrix from stream. " << specific_error.str()
            << " File position at start is "
            << pos_at_start << ", currently " << is.tellg()<<std::endl;
}

void SpMatrix::Print(int32 num){
	assert(num <= rows_);
	for(int32 i = 0; i < num; i++)
		std::cout<<data_[i]<<" ";
	std::cout<<std::endl;
}
