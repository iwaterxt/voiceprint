/*
 * SpMatrix.h
 *
 *  Created on: Feb 2, 2018
 *      Author: tao
 */

#ifndef SPMATRIX_H_
#define SPMATRIX_H_

#include <assert.h>
#include <stddef.h>
#include <iostream>
#include "base/common.h"
#include "matrix/Matrix.h"
#include "matrix/Vector.h"
#include "base/cblas-warppers.h"

class Matrix;
class Vector;

class SpMatrix{

public:
	SpMatrix() : data_(NULL), rows_(0) {}

	SpMatrix(MatrixIndexT r):
	      data_(NULL) {  Resize(r);  }

	SpMatrix(const SpMatrix &orig) : data_(NULL) {
	    Resize(orig.rows_);
	    CopyFromSp(orig);
	}

	~SpMatrix(){
		if(data_ != NULL){
			delete []data_;
			data_ = NULL;
		}
		rows_ = 0;
	}

	inline const BaseFloat* Data() const {return data_;}

	inline BaseFloat* Data() { return data_; }

	inline MatrixIndexT NumRows() const { return rows_; }

	inline MatrixIndexT NumCols() const { return rows_; }

	bool IsTridiagonal(BaseFloat cutoff = 1.0e-05) const;

	void Resize(MatrixIndexT dim);

	void Invert(BaseFloat *logdet = NULL, BaseFloat*det_sign= NULL,
	              bool inverse_needed = true);
	void SetZero();

	BaseFloat LogPosDefDet() const;

	void Scale(BaseFloat alpha);

	void AddToDiag(BaseFloat value);

	void CopyFromSp(const SpMatrix &orig);

	void CopyFromVec(const Vector& v);

	void AddMat2Vec(const BaseFloat alpha, const Matrix &M,
	                MatrixTransposeType transM, const Vector &v,
	                const BaseFloat beta = 0.0);
	void AddMat2Sp(const BaseFloat alpha, const Matrix &M,
            MatrixTransposeType transM, const SpMatrix &A,
            const BaseFloat beta = 0.0);
	void AddVec2(const BaseFloat alpha, const Vector &v);

	void AddSp(const BaseFloat alpha, const SpMatrix &Ma) {
		  assert(rows_ == Ma.NumRows());
		  size_t nr = rows_,
		      sz = (nr * (nr + 1)) / 2;
		  cblas_Xaxpy(sz, alpha, Ma.Data(), 1, data_, 1);
	}

	void AddSp(const BaseFloat alpha, const Matrix &rMa) ;

	void Tridiagonalize(Matrix *Q);

	void Eig(Vector *s, Matrix *P = NULL) const;

	void Qr(Matrix *Q);

	void Read(std::istream &is, bool binary=false);

	size_t SizeInBytes() const {
	    size_t nr = static_cast<size_t>(rows_);
	    return ((nr * (nr+1)) / 2) * sizeof(BaseFloat);
	}

	  inline BaseFloat operator() (MatrixIndexT r, MatrixIndexT c) const {
	    // if column is less than row, then swap these as matrix is stored
	    // as upper-triangular...  only allowed for const matrix object.
	    if (static_cast<MatrixIndexT>(c) >
	        static_cast<MatrixIndexT>(r))
	      std::swap(c, r);
	    // c<=r now so don't have to check c.
	    assert(static_cast<MatrixIndexT>(r) <
	                 static_cast<MatrixIndexT>(this->rows_));
	    return *(this->data_ + (r*(r+1)) / 2 + c);
	    // Duplicating code from PackedMatrix.h
	  }

	  inline BaseFloat &operator() (MatrixIndexT r, MatrixIndexT c) {
	    if (static_cast<MatrixIndexT>(c) >
	        static_cast<MatrixIndexT>(r))
	      std::swap(c, r);
	    // c<=r now so don't have to check c.
	    assert(static_cast<MatrixIndexT>(r) <
	                 static_cast<MatrixIndexT>(this->rows_));
	    return *(this->data_ + (r * (r + 1)) / 2 + c);
	    // Duplicating code from PackedMatrix.h
	  }

	  void Print(int num);

protected:

	BaseFloat* data_;
	int32 rows_;
private:
	void Init(MatrixIndexT dim);

};





#endif /* SPMATRIX_H_ */
