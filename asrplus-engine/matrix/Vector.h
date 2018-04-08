/*
 *  vector.h
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include "base/common.h"
#include "matrix/SpMatrix.h"
#include "matrix/Matrix.h"
#include <assert.h>
#include <stddef.h>


class Matrix;
class SpMatrix;

class Vector{

public:
	Vector():data_(NULL), dim_(0){}

	Vector(int dim):data_(NULL), dim_(dim)
	{
		data_ = new BaseFloat[dim];
		memset(data_, 0, sizeof(BaseFloat)*dim_);
	}

	Vector(const Matrix& M, int32 r);

	Vector(const Vector &other);

	~Vector(){
		if(data_ != NULL){
			delete []data_;
			data_ = NULL;
		}
		dim_ = 0;
	}

	inline int32 Dim() const {return dim_;}


	inline const BaseFloat* Data() const {return data_;}

	inline BaseFloat* Data() { return data_; }

	Vector Range(const MatrixIndexT o,const MatrixIndexT l) const ;

	void Init(int32 dim);

	BaseFloat Sum() const;

	BaseFloat Max() const;

	BaseFloat Min() const;

	BaseFloat Norm(BaseFloat p)const;

	void Add(BaseFloat c);

	void MulElements(const Vector &v) ;

	void Resize(int32 dim);

	void AddSpVec(const BaseFloat alpha, const SpMatrix &M, const Vector &v, const BaseFloat beta);

	void ApplyFloor(BaseFloat floor);

	void ApplyLog();

	void ApplyPow(BaseFloat pow);

	void Scale(BaseFloat value);

	float SumLog() const;
	//*this = *this + alpha * rv
	void AddVec(const BaseFloat alpha, const Vector &v);

	void SetZero();

	void InvertElements();

	void MulElements(Vector& vec);

	void Read(std::istream &is, bool binary);

	void AddMatVec( const BaseFloat alpha, const Matrix &M, MatrixTransposeType trans, const Vector &v,  const BaseFloat beta);

	void CopyFromData(BaseFloat* data, int32 dim);

	void CopyFromVec(const Vector& v);

	void CopyFromVecRange(const Vector& v);

	void CopyFromVecRange(const Vector& v, int32 length);

	void CopyColFromMat(const Matrix& M, int32 col);


	void CopyDiagFromSp(SpMatrix& M);

	BaseFloat ApplySoftMax();

	void Print(int32 num)const;

	inline BaseFloat operator() (int32 i) const {
	    assert(static_cast<int32>(i) <
	                 static_cast<int32>(dim_));
	    return *(data_ + i);
	}

	inline BaseFloat& operator() (int32 i) {
	    assert(static_cast<int32>(i) <
	                 static_cast<int32>(dim_));
	    return *(data_ + i);
	}

	Vector &operator = (const Vector &other) {
	    Resize(other.Dim());
	    this->CopyFromVec(other);
	    return *this;
	}

	Vector &operator = ( Vector &other) {
	    Resize(other.Dim());
	    this->CopyFromVec(other);
	    return *this;
	}




protected:
	BaseFloat* data_ ;
	int32 dim_;
};



#endif /* VECTOR_H_ */
