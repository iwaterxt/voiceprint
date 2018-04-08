/*
 * test.cc
 *
 *  Created on: Feb 22, 2018
 *      Author: tao
 */




/*
 * test.cc
 *
 *  Created on: Feb 22, 2018
 *      Author: tao
 */

#include <stdio.h>
#include "matrix/Matrix.h"
#include "matrix/Vector.h"
#include "base/common.h"



int main_test(){



	Matrix a(10, 10);
	Vector b(10);
	Vector c(10);

	b.Add(2.0);
	a.SetUnit();
	a.Scale(2.0);
	c.AddVec(1.0, b);
	c.AddMatVec(1.0, a, kNoTrans, b, 1.0);

	printf("%f\n", c(0));


	return 0;
}
