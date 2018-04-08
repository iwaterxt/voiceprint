/*
 * math.h
 *
 *  Created on: Feb 7, 2018
 *      Author: tao
 */

#ifndef MATH_H_
#define MATH_H_

#include <cmath>
#include <limits>
#include <stdlib.h>
#include <vector>
#include "base/common.h"
#include "matrix/Vector.h"
#include "matrix/Matrix.h"



inline float RandUniform(struct RandomState* state = NULL) {
  return static_cast<float>((Rand(state) + 1.0) / (RAND_MAX+2.0));
}

inline float RandGauss(struct RandomState* state = NULL) {
  return static_cast<float>(sqrtf (-2 * Log(RandUniform(state))));
}

int32 RoundUpToNearestPowerOfTwo(int32 n);


void ComplexImExp(BaseFloat x, BaseFloat *a_re, BaseFloat *a_im);

void ComplexMul(const BaseFloat &a_re, const BaseFloat &a_im, BaseFloat *b_re, BaseFloat *b_im);

void ComplexAddProduct(const BaseFloat &a_re, const BaseFloat &a_im,
                                                   const BaseFloat &b_re, const BaseFloat &b_im,
                                                   BaseFloat *c_re, BaseFloat *c_im);

void Factorize(int32 m, std::vector<int32> *factors);

void ComplexFft (Vector *v, bool forward, Vector *tmp_work = NULL);


void RealFft (Vector *v, bool forward);

void ComputeDctMatrix(Matrix *M);

void ComputeLifterCoeffs(BaseFloat Q, Vector *coeffs);

BaseFloat CosDistance(Vector& v1, Vector& v2);

void HouseBackward(MatrixIndexT dim, const BaseFloat *x, BaseFloat *v, BaseFloat *beta);

void Givens(BaseFloat a, BaseFloat b, BaseFloat *c, BaseFloat *s);

void QrStep(MatrixIndexT n,
            BaseFloat *diag,
            BaseFloat *off_diag,
            Matrix *Q);

void QrInternal(MatrixIndexT n,
                BaseFloat *diag,
                BaseFloat *off_diag,
                Matrix *Q);

BaseFloat VecVec(const Vector &ra,
            const Vector &rb);

#endif /* MATH_H_ */
