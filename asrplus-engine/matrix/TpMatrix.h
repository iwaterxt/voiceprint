/*
 * TpMatrix.h
 *
 *  Created on: Feb 18, 2018
 *      Author: tao
 */

#ifndef TPMATRIX_H_
#define TPMATRIX_H_

#include <assert.h>
#include "base/common.h"
#include "matrix/SpMatrix.h"

class SpMatrix;

class TpMatrix{

public:
	TpMatrix():data_(NULL), rows_(0){}

	TpMatrix(MatrixIndexT r):data_(NULL) {  Resize(r);  }

	~TpMatrix(){
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

	void Resize(MatrixIndexT r);

	void SetZero();


	size_t SizeInBytes() const {
	    size_t nr = static_cast<size_t>(rows_);
	    return ((nr * (nr+1)) / 2) * sizeof(BaseFloat);
	}


	void Cholesky(const SpMatrix& orig);

	  BaseFloat operator() (MatrixIndexT r, MatrixIndexT c) const {
	    if (static_cast<MatrixIndexT>(c) >
	        static_cast<MatrixIndexT>(r)) {
	      assert(static_cast<MatrixIndexT>(c) <
	                   static_cast<MatrixIndexT>(this->rows_));
	      return 0;
	    }
	    assert(static_cast<MatrixIndexT>(r) <
	                 static_cast<MatrixIndexT>(this->rows_));
	    // c<=r now so don't have to check c.
	    return *(this->data_ + (r*(r+1)) / 2 + c);
	    // Duplicating code from PackedMatrix.h
	  }

	  BaseFloat &operator() (MatrixIndexT r, MatrixIndexT c) {
	    assert(static_cast<MatrixIndexT>(r) <
	                 static_cast<MatrixIndexT>(this->rows_));
	    assert(static_cast<MatrixIndexT>(c) <=
	                 static_cast<MatrixIndexT>(r) &&
	                 "you cannot access the upper triangle of TpMatrix using "
	                 "a non-const matrix object.");
	    return *(this->data_ + (r*(r+1)) / 2 + c);
	    // Duplicating code from PackedMatrix.h
	  }

protected:

	BaseFloat* data_;
	int32 rows_;
private:
	void Init(MatrixIndexT dim);

};



#endif /* TPMATRIX_H_ */
