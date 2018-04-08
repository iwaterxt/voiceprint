/*
 *  vector.cc
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */
#include "base/common.h"
#include "matrix/SpMatrix.h"
#include "matrix/Vector.h"
#include "matrix/Matrix.h"
#include "base/io.h"
#include <assert.h>
#include <cstring>
#include <limits>
#include <vector>


Vector::Vector(const Matrix& M, int32 r){

	BaseFloat* data = M.Data();
	dim_ = M.NumCols();
	data_ = new BaseFloat[dim_];
	for(int32 i = 0; i < dim_ ; i++)
		data_[i] = data[r*dim_ + i];
}

Vector::Vector(const Vector &v) { //  (cannot be explicit)
  data_ = new BaseFloat[v.Dim()];
  dim_ = v.Dim();
  this->CopyFromVec(v);
}

Vector Vector::Range(const MatrixIndexT o,const MatrixIndexT l) const {
	assert(o < dim_&& o >=0);
	assert(l > 0 && l <= dim_-o);
    Vector tmp(l) ;
    tmp.CopyFromData(data_+o, l);
    return tmp;
}

void Vector::Init(int32 dim){

	if(this->data_ != NULL){
		delete []this->data_;
		this->data_ = NULL;
	}

	this->data_ = new BaseFloat[dim];
	this->dim_ = dim ;
	memset(data_, 0, sizeof(BaseFloat)*dim_);

}


BaseFloat Vector::Sum() const{

	  double sum = 0.0;
	  for (MatrixIndexT i = 0; i < dim_; i++) { sum += data_[i]; }
	  return sum;
}

BaseFloat Vector::Max() const{
	BaseFloat max = -1e20;
	for (MatrixIndexT i = 0; i < dim_; i++){
		if(data_[i] > max) max = data_[i];

	}
	return max;
}

BaseFloat Vector::Min() const{
	BaseFloat min = 1e20;
	for (MatrixIndexT i = 0; i < dim_; i++){
		if(data_[i] < min) min = data_[i];

	}
	return min;
}


BaseFloat Vector::Norm(BaseFloat p) const {
  assert(p >= 0.0);
  BaseFloat sum = 0.0;
  if (p == 0.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      if (data_[i] != 0.0) sum += 1.0;
    return sum;
  } else if (p == 1.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum += std::abs(data_[i]);
    return sum;
  } else if (p == 2.0) {
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum += data_[i] * data_[i];
    return std::sqrt(sum);
  } else if (p == std::numeric_limits<BaseFloat>::infinity()){
    for (MatrixIndexT i = 0; i < dim_; i++)
      sum = std::max(sum, std::abs(data_[i]));
    return sum;
  } else {
    BaseFloat tmp;
    bool ok = true;
    for (MatrixIndexT i = 0; i < dim_; i++) {
      tmp = pow(std::abs(data_[i]), p);
      if (tmp == HUGE_VAL) // HUGE_VAL is what pow returns on error.
        ok = false;
      sum += tmp;
    }
    tmp = pow(sum, static_cast<BaseFloat>(1.0/p));
    assert(tmp != HUGE_VAL); // should not happen here.
    if (ok) {
      return tmp;
    } else {
      BaseFloat maximum = Vector::Max(), minimum = Vector::Min(),
          max_abs = std::max(maximum, -minimum);
      assert(max_abs > 0); // Or should not have reached here.
      Vector tmp(*this);
      tmp.Scale(1.0 / max_abs);
      return tmp.Norm(p) * max_abs;
    }
  }
}

void Vector::Add(BaseFloat c) {
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] += c;
  }
}

void Vector::MulElements(const Vector &v) {
  assert(dim_ == v.dim_);
  for (MatrixIndexT i = 0; i < dim_; i++) {
    data_[i] *= v.data_[i];
  }
}

void Vector::Resize(int32 dim){

	  if (this->data_ != NULL) {
	    if (this->dim_ == dim) {
	      return;
	    } else {
	      delete []data_ ;
	      this->data_ = NULL;
	      this->dim_ = 0;
	    }
	  }
	  Init(dim);
}

void Vector::AddSpVec(const BaseFloat alpha,
                                 const SpMatrix &M,
                                 const Vector &v,
                                 const BaseFloat beta) {
  assert(M.NumRows() == v.dim_ && dim_ == v.dim_);
  assert(&v != this);
  cblas_Xspmv(alpha, M.NumRows(), M.Data(), v.Data(), 1, beta, data_, 1);
}


void Vector::ApplyFloor(BaseFloat floor){
	for (MatrixIndexT i = 0; i < dim_; i++){if(data_[i] < floor) data_[i] = floor;}
}

void Vector::ApplyLog(){
	for (MatrixIndexT i = 0; i < dim_; i++){data_[i] = Log(data_[i]);}
}

void Vector::ApplyPow(BaseFloat power){


		  if (power == 1.0) return;
		  if (power == 2.0) {
		    for (MatrixIndexT i = 0; i < dim_; i++)
		      data_[i] = data_[i] * data_[i];
		  } else if (power == 0.5) {
		    for (MatrixIndexT i = 0; i < dim_; i++) {
		      if (!(data_[i] >= 0.0))
		        std::cerr << "Cannot take square root of negative value "
		                  << data_[i];
		      data_[i] = std::sqrt(data_[i]);
		    }
		  } else {
		    for (MatrixIndexT i = 0; i < dim_; i++) {
		      data_[i] = pow(data_[i], power);
		      if (data_[i] == HUGE_VAL) {  // HUGE_VAL is what errno returns on error.
		        std::cerr << "Could not raise element "  << i << " to power "
		                  << power << ": returned value = " << data_[i];
		      }
		    }
		  }

}

void Vector::Scale(BaseFloat scale){
	if(scale == 1.0) return;
	for (MatrixIndexT i = 0; i < dim_; i++){

		data_[i] = data_[i] * scale ;
	}
}

float Vector::SumLog() const{
	  double sum_log = 0.0;
	  double prod = 1.0;
	  for (MatrixIndexT i = 0; i < dim_; i++) {
	    prod *= data_[i];
	    // Possible future work (arnab): change these magic values to pre-defined
	    // constants
	    if (prod < 1.0e-10 || prod > 1.0e+10) {
	      sum_log += Log(prod);
	      prod = 1.0;
	    }
	  }
	  if (prod != 1.0) sum_log += Log(prod);
	  return sum_log;
}

void Vector::AddVec(const float alpha,
                               const Vector &v) {
  assert(dim_ == v.dim_);
  assert(&v != this);
  cblas_Xaxpy(dim_, alpha, v.Data(), 1, data_, 1);
}


void Vector::SetZero(){
	for (MatrixIndexT i = 0; i < dim_; i++){data_[i] = 0;}
}

void Vector::InvertElements(){

	  for (MatrixIndexT i = 0; i < dim_; i++) {
		assert(data_[i] != 0.0);
	    data_[i] = static_cast<BaseFloat>(1 / data_[i]);
	  }

}

void Vector::MulElements(Vector& vec){
	assert(dim_ == vec.Dim());
	BaseFloat* vec_data = vec.Data();
	for (MatrixIndexT i = 0; i < dim_; i++) {

		data_[i] = data_[i] * vec_data[i];
	}
}

void Vector::Read(std::istream &is, bool binary){


	  std::ostringstream specific_error;
	  MatrixIndexT pos_at_start = is.tellg();

	  if (binary) {

	    const char *my_token =  (sizeof(BaseFloat) == 4 ? "FV" : "DV");
	    std::string token;
	    ReadToken(is, binary, &token);
	    if (token != my_token) {
	      specific_error << ": Expected token " << my_token << ", got " << token;
	      goto bad;
	    }
	    int32 size;
	    ReadBasicType(is, binary, &size);  // throws on error.
	    if ((MatrixIndexT)size != this->Dim())  this->Resize(size);
	    if (size > 0)
	      is.read(reinterpret_cast<char*>(this->data_), sizeof(BaseFloat)*size);
	    if (is.fail()) {
	      specific_error << "Error reading vector data (binary mode); truncated "
	          "stream? (size = " << size << ")";
	      goto bad;
	    }
	    return;
	  } else {  // Text mode reading; format is " [ 1.1 2.0 3.4 ]\n"
	    std::string s;
	    is >> s;
	    // if ((s.compare("DV") == 0) || (s.compare("FV") == 0)) {  // Back compatibility.
	    //  is >> s;  // get dimension
	    //  is >> s;  // get "["
	    // }
	    if (is.fail()) { specific_error << "EOF while trying to read vector."; goto bad; }
	    if (s.compare("[]") == 0) { Resize(0); return; } // tolerate this variant.
	    if (s.compare("[")) { specific_error << "Expected \"[\" but got " << s; goto bad; }
	    std::vector<BaseFloat> data;
	    while (1) {
	      int i = is.peek();
	      if (i == '-' || (i >= '0' && i <= '9')) {  // common cases first.
	        BaseFloat r;
	        is >> r;
	        if (is.fail()) { specific_error << "Failed to read number."; goto bad; }
	        if (! std::isspace(is.peek()) && is.peek() != ']') {
	          specific_error << "Expected whitespace after number."; goto bad;
	        }
	        data.push_back(r);
	        // But don't eat whitespace... we want to check that it's not newlines
	        // which would be valid only for a matrix.
	      } else if (i == ' ' || i == '\t') {
	        is.get();
	      } else if (i == ']') {
	        is.get();  // eat the ']'
	        this->Resize(data.size());
	        for (size_t j = 0; j < data.size(); j++)
	          this->data_[j] = data[j];
	        i = is.peek();
	        if (static_cast<char>(i) == '\r') {
	          is.get();
	          is.get();  // get \r\n (must eat what we wrote)
	        } else if (static_cast<char>(i) == '\n') { is.get(); } // get \n (must eat what we wrote)
	        if (is.fail()) {
	          std::cout << "After end of vector data, read error."<<std::endl;
	          // we got the data we needed, so just warn for this error.
	        }
	        return;  // success.
	      } else if (i == -1) {
	        specific_error << "EOF while reading vector data.";
	        goto bad;
	      } else if (i == '\n' || i == '\r') {
	        specific_error << "Newline found while reading vector (maybe it's a matrix?)";
	        goto bad;
	      } else {
	        is >> s;  // read string.
	        if (!strcasecmp(s.c_str(), "inf") ||
	            !strcasecmp(s.c_str(), "infinity")) {
	          data.push_back(std::numeric_limits<BaseFloat>::infinity());
	          std::cout << "Reading infinite value into vector."<<std::endl;
	        } else if (!strcasecmp(s.c_str(), "nan")) {
	          data.push_back(std::numeric_limits<BaseFloat>::quiet_NaN());
	          std::cout << "Reading NaN value into vector."<<std::endl;
	        } else {
	          specific_error << "Expecting numeric vector data, got " << s;
	          goto  bad;
	        }
	      }
	    }
	  }
	  // we never reach this line (the while loop returns directly).
	bad:
	  std::cout << "Failed to read vector from stream.  " << specific_error.str()
	            << " File position at start is "
	            << pos_at_start<<", currently "<<is.tellg()<<std::endl;

}

void Vector::AddMatVec( const BaseFloat alpha, const Matrix &M, MatrixTransposeType trans, const Vector &v,  const BaseFloat beta){

	  assert((trans == kNoTrans && M.NumCols() == v.dim_ && M.NumRows() == dim_)
	               || (trans == kTrans && M.NumRows() == v.dim_ && M.NumCols() == dim_));
	  assert(&v != this);
	  cblas_Xgemv(trans, M.NumRows(), M.NumCols(), alpha, M.Data(), M.NumCols(),
	              v.Data(), 1, beta, data_, 1);
}


void Vector::CopyFromData(BaseFloat* data, int32 dim){

	  assert(Dim() == dim);
	  if (data_ != data) {
	    memcpy(this->data_, data, dim_ * sizeof(BaseFloat));
	  }

}

void Vector::CopyFromVec(const Vector& v){
	  assert(Dim() == v.Dim());
	  if (data_ != v.data_) {
	    memcpy(this->data_, v.data_, dim_ * sizeof(BaseFloat));
	  }
}

void Vector::CopyFromVecRange(const Vector& v){

	  if (data_ != v.data_){
	  	memcpy(this->data_, v.data_, dim_ * sizeof(BaseFloat));
	  }
}

void Vector::CopyFromVecRange(const Vector& v, int32 length){

	  if (data_ != v.data_){
	  	memcpy(this->data_, v.data_, length * sizeof(BaseFloat));
	  }
}

void Vector::CopyColFromMat(const Matrix& M, int32 col){
	assert(dim_ == M.NumRows());
	for(MatrixIndexT i = 0; i < dim_; i++)
		data_[i] = M(i, col);
}

void Vector::CopyDiagFromSp(SpMatrix& M){
	  assert(dim_ == M.NumCols());
	  for (MatrixIndexT i = 0; i < dim_; i++)
	    data_[i] = M(i, i);

}

BaseFloat Vector::ApplySoftMax(){

	  BaseFloat max = this->Max(), sum = 0.0;
	  for (MatrixIndexT i = 0; i < dim_; i++) {
	    sum += (data_[i] = Exp(data_[i] - max));
	  }
	  this->Scale(1.0 / sum);
	  return max + Log(sum);

}

void Vector::Print(int32 num) const{

	for(int32 i = 0; i < num; i++)
		std::cout<<data_[i]<<" ";
	std::cout<<std::endl;
}
