/*
 *  matrix.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <assert.h>
#include <string.h>
#include "base/common.h"
#include "matrix/Vector.h"
#include "matrix/SpMatrix.h"


class Vector;
class SpMatrix;

class Matrix{

public:

	Matrix():data_(NULL), rows_(0), cols_(0){}


    Matrix(BaseFloat *data, MatrixIndexT cols, MatrixIndexT rows, MatrixIndexT stride) :
	    data_(data), rows_(rows), cols_(cols) {}

	Matrix(int32 rows, int32 cols): data_(0) , rows_(rows),cols_(cols){

		data_ = new BaseFloat[rows_*cols_];
		memset(data_, 0, sizeof(BaseFloat)*rows_*cols_);
	}


	Matrix(const Matrix& M, MatrixTransposeType trans = kNoTrans);


	~Matrix(){
		if(data_ != 0){
			delete []data_;
			data_ = NULL;
		}
		rows_ = 0;
		cols_ = 0;
	}

	inline int32 NumRows()const {return rows_;}

	inline int32 NumCols()const {return cols_;}

	void Resize(int32 numrows, int32 numcols);

	void SetZero();

	BaseFloat* RowData(int32 rows);

	void Init(int32 numrows, int32 numcols);

	void ApplyPow(float power);

	void Scale(BaseFloat scale);

	void CopyColFromVec(const Vector &rv, const MatrixIndexT col);

	void CopyRowFromVec(const Vector &rv, const MatrixIndexT row);

	void CopyRowsFromVec(const Vector &rv);

	void CopyFromMat(const Matrix& M, MatrixTransposeType trans = kNoTrans);

	void CopyFromMat(const Matrix& M, MatrixIndexT row_offet, MatrixIndexT numrows, MatrixIndexT col_offset, MatrixIndexT numcols);

	void CopyFromSp(const SpMatrix &M);


	void AddMatMat(const BaseFloat alpha,
            const Matrix& A,
            MatrixTransposeType transA,
            const Matrix& B,
            MatrixTransposeType transB,
            const BaseFloat beta);

	void AddSpMat(const BaseFloat alpha,
	                const SpMatrix& A,
	                const Matrix& B, MatrixTransposeType transB,
	                const BaseFloat beta) ;
	void AddVec2Row(MatrixIndexT row, const BaseFloat alpha, Vector &v);

	void SetRandn();

	void SetUnit();

	void Transpose();

	inline BaseFloat* Data()const {return data_;}

	Vector Row(int32 i) const ;

	inline BaseFloat* RowData(int32 i) const {return data_ + i * cols_;}

	inline const BaseFloat operator() (MatrixIndexT r, MatrixIndexT c) const {
	    assert(static_cast<int32>(r) <
	                          static_cast<int32>(rows_) &&
	                          static_cast<int32>(c) <
	                          static_cast<int32>(cols_));
	    return *(data_ + r * cols_ + c);
	}

	inline  BaseFloat& operator() (MatrixIndexT r, MatrixIndexT c)  {
	    assert(static_cast<int32>(r) <
	                          static_cast<int32>(rows_) &&
	                          static_cast<int32>(c) <
	                          static_cast<int32>(cols_));
	    return *(data_ + r * cols_ + c);
	}


	void Read(std::istream &is, bool binary=false);

	Matrix &operator = ( Matrix &other) {
	    if (Matrix::NumRows() != other.NumRows() ||
	        Matrix::NumCols() != other.NumCols())
	      Resize(other.NumRows(), other.NumCols());
	    Matrix::CopyFromMat(other);
	    return *this;
	}

	Matrix &operator = (const Matrix &other) {
	    if (Matrix::NumRows() != other.NumRows() ||
	        Matrix::NumCols() != other.NumCols())
	      Resize(other.NumRows(), other.NumCols());
	    Matrix::CopyFromMat(other);
	    return *this;
	  }

	  inline BaseFloat*  Data_workaround() const {
	    return data_;
	  }

protected:

	BaseFloat* data_;
	int32 rows_;
	int32 cols_;

};


#endif /* MATRIX_H_ */
