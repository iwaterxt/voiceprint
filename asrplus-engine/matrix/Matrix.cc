/*
 *  matrix.cc
 *
 *  Created on: Sep 11, 2017
 *      Author: tao
 */

#include "matrix/Matrix.h"
#include "util/math.h"
#include "base/io.h"
#include <cstring>
#include <assert.h>
#include <vector>




Matrix::Matrix(const Matrix& M, MatrixTransposeType trans){

	data_ = new BaseFloat[M.NumRows()*M.NumCols()];
	rows_ = M.NumRows();
	cols_ = M.NumCols();

	if(trans == kNoTrans){
		this->CopyFromMat(M);
		rows_ = M.NumRows();
		cols_ = M.NumCols();
	}else{
		BaseFloat* M_data = M.Data();
		rows_ = M.NumCols();
		cols_ = M.NumRows();

		for(int32 i = 0; i < rows_; i++ )
			for(int32 j = 0; j < cols_; j++){
				data_[i * cols_ + j] = M_data[j * rows_ + i];
			}
	}
}


void Matrix::Resize(int32 numrows, int32 numcols){

	  if (this->data_ != NULL) {
	    if (numrows == rows_
	        && numcols == cols_) {
	        this->SetZero();
	        return;
	    }
	    else{
	    	delete []this->data_;
	    	this->data_ = NULL;
	    	rows_ = 0 ;
	    	cols_ = 0 ;
	    }

	  }
	  Init(numrows, numcols);
	  this->SetZero();

}

void Matrix::SetZero(){

	std::memset(data_, 0, sizeof(BaseFloat)*rows_*cols_);
}
Vector Matrix::Row(int32 i) const
{
	assert(i < rows_);
	Vector v(cols_);
	v.CopyFromData(data_ + i*cols_, cols_);
	return v;
}

BaseFloat* Matrix::RowData(int32 rows){
	assert(rows < rows_);
	return data_ + rows * cols_;
}

void Matrix::Init(int32 numrows, int32 numcols){

	if(numrows * numcols == 0){
		this->data_ =NULL;
		this->rows_ = 0 ;
		this->cols_ = 0;
		return;
	}
	assert(numrows * numcols != 0);
	this->data_ = new BaseFloat[numrows * numcols];
	if(this->data_ == NULL){
		throw std::bad_alloc();
	}
	this->rows_ = numrows;
	this->cols_ = numcols;
}

void Matrix::ApplyPow(float power){

	  for (MatrixIndexT i = 0; i < rows_; i++) {
	    for(MatrixIndexT j = 0; j < cols_ ; j++){
	    	data_[i*cols_ + j] = pow(data_[i*cols_ + j], power);
	    }
	  }

}

void Matrix::Scale(BaseFloat scale){

	if(scale == 1.0) return;
	for (MatrixIndexT i = 0; i < rows_; i++) {
	    for(MatrixIndexT j = 0; j < cols_ ; j++){
	    	data_[i*cols_ + j] = data_[i*cols_ + j] * scale;
	    }
	}

}

void Matrix::CopyColFromVec(const Vector &rv, const MatrixIndexT col){

	  assert(rv.Dim() == rows_ &&
	               static_cast<MatrixIndexT>(col) <
	               static_cast<MatrixIndexT>(cols_));

	  const BaseFloat *rv_data = rv.Data();
	  BaseFloat *col_data = data_ + col;

	  for (MatrixIndexT r = 0; r < rows_; r++)
	    col_data[r * cols_] = rv_data[r];

}

void Matrix::CopyRowFromVec(const Vector &rv, const MatrixIndexT row){
	  assert(rv.Dim() == cols_ &&
	               static_cast<MatrixIndexT>(row) <
	               static_cast<MatrixIndexT>(rows_));

	  const BaseFloat *rv_data = rv.Data();
	  BaseFloat *row_data = RowData(row);

	  std::memcpy(row_data, rv_data, cols_ * sizeof(BaseFloat));
}

void Matrix::CopyRowsFromVec(const Vector &rv){
	  if (rv.Dim() == rows_*cols_) {
	      // one big copy operation.
	      const BaseFloat *rv_data = rv.Data();
	      std::memcpy(data_, rv_data, sizeof(BaseFloat)*rows_*cols_);

	  } else if (rv.Dim() == cols_) {
	    const BaseFloat *rv_data = rv.Data();
	    for (MatrixIndexT r = 0; r < rows_; r++)
	      std::memcpy(RowData(r), rv_data, sizeof(BaseFloat)*cols_);
	  } else {
	    std::cerr << "Wrong sized arguments"<<std::endl;
	  }
}

void Matrix::CopyFromMat(const Matrix& M, MatrixTransposeType trans){

	  if (static_cast<const void*>(M.Data()) ==
	      static_cast<const void*>(this->Data())) {
	    // CopyFromMat called on same data.  Nothing to do (except sanity checks).
	    assert(trans == kNoTrans && M.NumRows() == NumRows() &&
	                 M.NumCols() == NumCols() );
	    return;
	  }
	  if (trans == kNoTrans) {
	    assert(rows_ == M.NumRows() && cols_ == M.NumCols());
	    for (MatrixIndexT i = 0; i < rows_; i++)
	      this->CopyRowFromVec(M.Row(i), i);//TODO
	  } else {
	    assert(cols_ == M.NumRows() && rows_ == M.NumCols());
	    BaseFloat *this_data = data_;
	    const BaseFloat *other_data = M.Data();
	    for (MatrixIndexT i = 0; i < rows_; i++)
	      for (MatrixIndexT j = 0; j < cols_; j++)
	        this_data[i * cols_ + j] = other_data[j * cols_ + i];
	  }
}

void Matrix::CopyFromMat(const Matrix& M, MatrixIndexT row_offet, MatrixIndexT numrows, MatrixIndexT col_offset, MatrixIndexT numcols){

	int32 M_rows = M.NumRows();
	int32 M_cols = M.NumCols();
	assert(rows_ == numrows && cols_ == numcols);
	assert(row_offet >= 0 && row_offet < M_rows);
	assert(numrows > 0 && numrows <= M_rows - row_offet);

	assert(col_offset >= 0 && col_offset < M_cols);
	assert(numcols > 0 && numcols <= M_cols - col_offset);

	BaseFloat* data_M = M.Data();
	for(int r = 0; r < numrows;  r++)
		for(int c = 0; c < numcols; c++){

			data_[r*numcols + c] = data_M[(row_offet+r)*M_cols + (col_offset + c)];
		}
}

void Matrix::CopyFromSp(const SpMatrix &M){

	  assert(rows_ == M.NumRows());
	  assert( cols_ == rows_);
	  MatrixIndexT num_rows = rows_, stride = cols_;
	  const float *Mdata = M.Data();
	  float *row_data = data_, *col_data = data_;
	  for (MatrixIndexT i = 0; i < num_rows; i++) {
	    cblas_scopy(i+1, Mdata, 1, row_data, 1); // copy to the row.
	    cblas_scopy(i, Mdata, 1, col_data, stride); // copy to the column.
	    Mdata += i+1;
	    row_data += stride;
	    col_data += 1;
	  }
}


void Matrix::AddMatMat(const BaseFloat alpha,
        const Matrix& A,
        MatrixTransposeType transA,
        const Matrix& B,
        MatrixTransposeType transB,
        const BaseFloat beta){
	  assert((transA == kNoTrans && transB == kNoTrans && A.cols_ == B.rows_ && A.rows_ == rows_ && B.cols_ == cols_)
	               || (transA == kTrans && transB == kNoTrans && A.rows_ == B.rows_ && A.cols_ == rows_ && B.cols_ == cols_)
	               || (transA == kNoTrans && transB == kTrans && A.cols_ == B.cols_ && A.rows_ == rows_ && B.rows_ == cols_)
	               || (transA == kTrans && transB == kTrans && A.rows_ == B.cols_ && A.cols_ == rows_ && B.rows_ == cols_));
	  assert(&A !=  this && &B != this);
	  if (rows_ == 0) return;
	  cblas_Xgemm(alpha, transA, A.data_, A.rows_, A.cols_, A.cols_,
	              transB, B.data_, B.cols_, beta, data_, rows_, cols_, cols_);
}

void Matrix::AddSpMat(const BaseFloat alpha,
                const SpMatrix& A,
                const Matrix& B, MatrixTransposeType transB,
                const BaseFloat beta) {
    Matrix M(A.NumRows(), A.NumRows());
    M.CopyFromSp(A);
    AddMatMat(alpha, M, kNoTrans, B, transB, beta);
}

void Matrix::AddVec2Row(MatrixIndexT row, const BaseFloat alpha, Vector &v){

	assert(cols_ == v.Dim());
	BaseFloat* data = v.Data();
	for(int32 i = 0; i < cols_; i++)
		data_[row*cols_ + i] += alpha * data[i];
}

void Matrix::SetRandn(){

	for (int32 r=0; r<rows_; r++) {
     	for (int32 c=0; c<cols_; c++) {
     		srand(time(0));
     		int32 r= rand();
       	 	(*this)(r,c) = r * RandGauss(); // 0-mean Gauss with given std_dev
      	}
    }
}

void Matrix::SetUnit() {
  SetZero();
  for (MatrixIndexT row = 0; row < std::min(rows_, cols_); row++)
    (*this)(row, row) = 1.0;
}

void Matrix::Transpose() {
  if (this->rows_ != this->cols_) {
    Matrix tmp(*this, kTrans);
    Resize(this->cols_, this->rows_);
    this->CopyFromMat(tmp);
  } else {
	  MatrixIndexT M = rows_;
	  for (MatrixIndexT i = 0;i < M;i++)
	    for (MatrixIndexT j = 0;j < i;j++) {
	      BaseFloat &a = (*this)(i, j), &b = (*this)(j, i);
	      std::swap(a, b);
	    }
  }
}

void Matrix::Read(std::istream &is, bool binary){

	  // now assume add == false.
	  MatrixIndexT pos_at_start = is.tellg();
	  std::ostringstream specific_error;

	  if (binary) {  // Read in binary mode.

	    const char *my_token =  (sizeof(BaseFloat) == 4 ? "FM" : "DM");

	    std::string token;
	    ReadToken(is, binary, &token);
	    if (token != my_token) {
	      specific_error << ": Expected token " << my_token << ", got " << token;
	      goto bad;
	    }
	    int32 rows, cols;
	    ReadBasicType(is, binary, &rows);  // throws on error.
	    ReadBasicType(is, binary, &cols);  // throws on error.
	    if ((MatrixIndexT)rows != this->rows_ || (MatrixIndexT)cols != this->cols_) {
	      this->Resize(rows, cols);
	    }
	    if (rows*cols!=0) {
	      is.read(reinterpret_cast<char*>(this->Data()),
	              sizeof(BaseFloat)*rows*cols);
	      if (is.fail()) goto bad;
	    } else {
	      for (MatrixIndexT i = 0; i < (MatrixIndexT)rows; i++) {
	        is.read(reinterpret_cast<char*>(this->RowData(i)), sizeof(BaseFloat)*cols);
	        if (is.fail()) goto bad;
	      }
	    }
	    if (is.eof()) return;
	    if (is.fail()) goto bad;
	    return;
	  } else {  // Text mode.
	    std::string str;
	    is >> str; // get a token
	    if (is.fail()) { specific_error << ": Expected \"[\", got EOF"; goto bad; }
	    // if ((str.compare("DM") == 0) || (str.compare("FM") == 0)) {  // Back compatibility.
	    // is >> str;  // get #rows
	    //  is >> str;  // get #cols
	    //  is >> str;  // get "["
	    // }
	    if (str == "[]") { Resize(0, 0); return; } // Be tolerant of variants.
	    else if (str != "[") {
	      specific_error << ": Expected \"[\", got \"" << str << '"';
	      goto bad;
	    }
	    // At this point, we have read "[".
	    std::vector<std::vector<BaseFloat>* > data;
	    std::vector<BaseFloat> *cur_row = new std::vector<BaseFloat>;
	    while (1) {
	      int i = is.peek();
	      if (i == -1) { specific_error << "Got EOF while reading matrix data"; goto cleanup; }
	      else if (static_cast<char>(i) == ']') {  // Finished reading matrix.
	        is.get();  // eat the "]".
	        i = is.peek();
	        if (static_cast<char>(i) == '\r') {
	          is.get();
	          is.get();  // get \r\n (must eat what we wrote)
	        } else if (static_cast<char>(i) == '\n') { is.get(); } // get \n (must eat what we wrote)
	        if (is.fail()) {
	          std::cout << "After end of matrix data, read error."<<std::endl;
	          // we got the data we needed, so just warn for this error.
	        }
	        // Now process the data.
	        if (!cur_row->empty()) data.push_back(cur_row);
	        else delete(cur_row);
	        cur_row = NULL;
	        if (data.empty()) { this->Resize(0, 0); return; }
	        else {
	          int32 num_rows = data.size(), num_cols = data[0]->size();
	          this->Resize(num_rows, num_cols);
	          for (int32 i = 0; i < num_rows; i++) {
	            if (static_cast<int32>(data[i]->size()) != num_cols) {
	              specific_error << "Matrix has inconsistent #cols: " << num_cols
	                             << " vs." << data[i]->size() << " (processing row"
	                             << i << ")";
	              goto cleanup;
	            }
	            for (int32 j = 0; j < num_cols; j++)
	              (*this)(i, j) = (*(data[i]))[j];
	            delete data[i];
	            data[i] = NULL;
	          }
	        }
	        return;
	      } else if (static_cast<char>(i) == '\n' || static_cast<char>(i) == ';') {
	        // End of matrix row.
	        is.get();
	        if (cur_row->size() != 0) {
	          data.push_back(cur_row);
	          cur_row = new std::vector<BaseFloat>;
	          cur_row->reserve(data.back()->size());
	        }
	      } else if ( (i >= '0' && i <= '9') || i == '-' ) {  // A number...
	        BaseFloat r;
	        is >> r;
	        if (is.fail()) {
	          specific_error << "Stream failure/EOF while reading matrix data.";
	          goto cleanup;
	        }
	        cur_row->push_back(r);
	      } else if (isspace(i)) {
	        is.get();  // eat the space and do nothing.
	      } else {  // NaN or inf or error.
	        std::string str;
	        is >> str;
	        if (!strcasecmp(str.c_str(), "inf") ||
	            !strcasecmp(str.c_str(), "infinity")) {
	          cur_row->push_back(std::numeric_limits<BaseFloat>::infinity());
	          std::cout << "Reading infinite value into matrix."<<std::endl;
	        } else if (!strcasecmp(str.c_str(), "nan")) {
	          cur_row->push_back(std::numeric_limits<BaseFloat>::quiet_NaN());
	          std::cout << "Reading NaN value into matrix.";
	        } else {
	          specific_error << "Expecting numeric matrix data, got " << str;
	          goto cleanup;
	        }
	      }
	    }
	    // Note, we never leave the while () loop before this
	    // line (we return from it.)
	 cleanup: // We only reach here in case of error in the while loop above.
	    if(cur_row != NULL)
	      delete cur_row;
	    for (size_t i = 0; i < data.size(); i++)
	      if(data[i] != NULL)
	        delete data[i];
	    // and then go on to "bad" below, where we print error.
	  }
	bad:
	  std::cerr << "Failed to read matrix from stream.  " << specific_error.str()
	            << " File position at start is "
	            << pos_at_start << ", currently " << is.tellg()<<std::endl;

}




