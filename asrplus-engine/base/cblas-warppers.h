/*
 * cblas-warppers.h
 *
 *  Created on: Feb 12, 2018
 *      Author: tao
 */

#ifndef CBLAS_WARPPERS_H_
#define CBLAS_WARPPERS_H_

#include "base/common.h"
#include <cblas.h>
#include <lapacke.h>


inline void cblas_Xcopy(const int N, const float *X, const int incX, float *Y,
                        const int incY) {
  cblas_scopy(N, X, incX, Y, incY);
}



inline float cblas_Xasum(const int N, const float *X, const int incX) {
  return cblas_sasum(N, X, incX);
}


inline void cblas_Xrot(const int N, float *X, const int incX, float *Y,
                       const int incY, const float c, const float s) {
  cblas_srot(N, X, incX, Y, incY, c, s);
}

inline float cblas_Xdot(const int N, const float *const X,
                        const int incX, const float *const Y,
                        const int incY) {
  return cblas_sdot(N, X, incX, Y, incY);
}

inline void cblas_Xaxpy(const int N, const float alpha, const float *X,
                        const int incX, float *Y, const int incY) {
  cblas_saxpy(N, alpha, X, incX, Y, incY);
}

inline void cblas_Xscal(const int N, const float alpha, float *data,
                        const int inc) {
  cblas_sscal(N, alpha, data, inc);
}

inline void cblas_Xspmv(const float alpha, const int num_rows, const float *Mdata,
                        const float *v, const int v_inc,
                        const float beta, float *y, const int y_inc) {
  cblas_sspmv(CblasRowMajor, CblasLower, num_rows, alpha, Mdata, v, v_inc, beta, y, y_inc);
}

inline void cblas_Xtpmv(MatrixTransposeType trans, const float *Mdata,
                        const int num_rows, float *y, const int y_inc) {
  cblas_stpmv(CblasRowMajor, CblasLower, static_cast<CBLAS_TRANSPOSE>(trans),
              CblasNonUnit, num_rows, Mdata, y, y_inc);
}


inline void cblas_Xtpsv(MatrixTransposeType trans, const float *Mdata,
                        const int num_rows, float *y, const int y_inc) {
  cblas_stpsv(CblasRowMajor, CblasLower, static_cast<CBLAS_TRANSPOSE>(trans),
              CblasNonUnit, num_rows, Mdata, y, y_inc);
}

// x = alpha * M * y + beta * x
inline void cblas_Xspmv(MatrixIndexT dim, float alpha, const float *Mdata,
                        const float *ydata, MatrixIndexT ystride,
                        float beta, float *xdata, MatrixIndexT xstride) {
  cblas_sspmv(CblasRowMajor, CblasLower, dim, alpha, Mdata,
              ydata, ystride, beta, xdata, xstride);
}

// Implements  A += alpha * (x y'  + y x'); A is symmetric matrix.
inline void cblas_Xspr2(MatrixIndexT dim, float alpha, const float *Xdata,
                        MatrixIndexT incX, const float *Ydata, MatrixIndexT incY,
                          float *Adata) {
  cblas_sspr2(CblasRowMajor, CblasLower, dim, alpha, Xdata,
              incX, Ydata, incY, Adata);
}


// Implements  A += alpha * (x x'); A is symmetric matrix.
inline void cblas_Xspr(MatrixIndexT dim, float alpha, const float *Xdata,
                       MatrixIndexT incX, float *Adata) {
  cblas_sspr(CblasRowMajor, CblasLower, dim, alpha, Xdata, incX, Adata);
}


// sgemv,dgemv: y = alpha M x + beta y.
inline void cblas_Xgemv(MatrixTransposeType trans, MatrixIndexT num_rows,
                        MatrixIndexT num_cols, float alpha, const float *Mdata,
                        MatrixIndexT stride, const float *xdata,
                        MatrixIndexT incX, float beta, float *ydata, MatrixIndexT incY) {
  cblas_sgemv(CblasRowMajor, static_cast<CBLAS_TRANSPOSE>(trans), num_rows,
              num_cols, alpha, Mdata, stride, xdata, incX, beta, ydata, incY);
}



// sgbmv, dgmmv: y = alpha M x +  + beta * y.
inline void cblas_Xgbmv(MatrixTransposeType trans, MatrixIndexT num_rows,
                        MatrixIndexT num_cols, MatrixIndexT num_below,
                        MatrixIndexT num_above, float alpha, const float *Mdata,
                        MatrixIndexT stride, const float *xdata,
                        MatrixIndexT incX, float beta, float *ydata, MatrixIndexT incY) {
  cblas_sgbmv(CblasRowMajor, static_cast<CBLAS_TRANSPOSE>(trans), num_rows,
              num_cols, num_below, num_above, alpha, Mdata, stride, xdata,
              incX, beta, ydata, incY);
}

inline void cblas_Xgemm(const float alpha,
                        MatrixTransposeType transA,
                        const float *Adata,
                        MatrixIndexT a_num_rows, MatrixIndexT a_num_cols, MatrixIndexT a_stride,
                        MatrixTransposeType transB,
                        const float *Bdata, MatrixIndexT b_stride,
                        const float beta,
                        float *Mdata,
                        MatrixIndexT num_rows, MatrixIndexT num_cols,MatrixIndexT stride) {
  cblas_sgemm(CblasRowMajor, static_cast<CBLAS_TRANSPOSE>(transA),
              static_cast<CBLAS_TRANSPOSE>(transB),
              num_rows, num_cols, transA == kNoTrans ? a_num_cols : a_num_rows,
              alpha, Adata, a_stride, Bdata, b_stride,
              beta, Mdata, stride);
}

inline void cblas_Xsymm(const float alpha,
                        MatrixIndexT sz,
                        const float *Adata,MatrixIndexT a_stride,
                        const float *Bdata,MatrixIndexT b_stride,
                        const float beta,
                        float *Mdata, MatrixIndexT stride) {
  cblas_ssymm(CblasRowMajor, CblasLeft, CblasLower, sz, sz, alpha, Adata,
              a_stride, Bdata, b_stride, beta, Mdata, stride);
}

// ger: M += alpha x y^T.
inline void cblas_Xger(MatrixIndexT num_rows, MatrixIndexT num_cols, float alpha,
                       const float *xdata, MatrixIndexT incX, const float *ydata,
                       MatrixIndexT incY, float *Mdata, MatrixIndexT stride) {
  cblas_sger(CblasRowMajor, num_rows, num_cols, alpha, xdata, 1, ydata, 1,
             Mdata, stride);
}


// syrk: symmetric rank-k update.
// if trans==kNoTrans, then C = alpha A A^T + beta C
// else C = alpha A^T A + beta C.
// note: dim_c is dim(C), other_dim_a is the "other" dimension of A, i.e.
// num-cols(A) if kNoTrans, or num-rows(A) if kTrans.
// We only need the row-major and lower-triangular option of this, and this
// is hard-coded.
inline void cblas_Xsyrk (
    const MatrixTransposeType trans, const MatrixIndexT dim_c,
    const MatrixIndexT other_dim_a, const float alpha, const float *A,
    const MatrixIndexT a_stride, const float beta, float *C,
    const MatrixIndexT c_stride) {
  cblas_ssyrk(CblasRowMajor, CblasLower, static_cast<CBLAS_TRANSPOSE>(trans),
              dim_c, other_dim_a, alpha, A, a_stride, beta, C, c_stride);
}


/// matrix-vector multiply using a banded matrix; we always call this
/// with b = 1 meaning we're multiplying by a diagonal matrix.  This is used for
/// elementwise multiplication.  We miss some of the arguments out of this
/// wrapper.


inline void cblas_Xsbmv1(
    const MatrixIndexT dim,
    const float *A,
    const float alpha,
    const float *x,
    const float beta,
    float *y) {
  cblas_ssbmv(CblasRowMajor, CblasLower, dim, 0, alpha, A,
              1, x, 1, beta, y, 1);
}


/// This is not really a wrapper for CBLAS as CBLAS does not have this; in future we could
/// extend this somehow.


inline void mul_elements(
    const MatrixIndexT dim,
    const float *a,
    float *b) { // does b *= a, elementwise.
  float c1, c2, c3, c4;
  MatrixIndexT i;
  for (i = 0; i + 4 <= dim; i += 4) {
    c1 = a[i] * b[i];
    c2 = a[i+1] * b[i+1];
    c3 = a[i+2] * b[i+2];
    c4 = a[i+3] * b[i+3];
    b[i] = c1;
    b[i+1] = c2;
    b[i+2] = c3;
    b[i+3] = c4;
  }
  for (; i < dim; i++)
    b[i] *= a[i];
}

void inline clapack_Xsptrf(int32 *num_rows, float *Mdata,
                           int32 *ipiv, int32 *result) {
  ssptrf_(const_cast<char *>("U"), num_rows, Mdata, ipiv, result);
}

void inline clapack_Xsptri(int32 *num_rows, float *Mdata,
                           int32 *ipiv, float *work, int32 *result) {
  ssptri_(const_cast<char *>("U"), num_rows, Mdata, ipiv, work, result);
}

#endif /* CBLAS_WARPPERS_H_ */
