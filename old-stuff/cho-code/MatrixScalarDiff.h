#ifndef MATRIXSCALARDIFF_H
#define MATRIXSCALARDIFF_H

#include <cstdio>
#include <cstdlib>
#include "MatrixUtil.h"

using namespace std;

class MatrixScalarDiff
{
	public:
		double value;
		mat_t diff;

		MatrixScalarDiff() {}
		MatrixScalarDiff(double value_, mat_t diff_) {
			value = value_;
			diff = diff_;

		}
		~MatrixScalarDiff() {}

		void print(FILE *fp) const {
			printf("value: %lf\n", value);
			printf("gradient: \n");
			diff.print(fp);
		}

		MatrixScalarDiff& operator= ( const MatrixScalarDiff &x) {
			value = x.value;
			diff = x.diff;
			return(*this);
		}

};

MatrixScalarDiff operator- (const MatrixScalarDiff &lhs, const MatrixScalarDiff &rhs)
{
	if ( (lhs.diff.m != rhs.diff.m) || (lhs.diff.n != rhs.diff.n) )
	{
		printf("ERROR!! input dimension mismatrch for MatrixScalarDiff + MatrixScalarDiff\n");
		exit(1);
	}
	return MatrixScalarDiff(lhs.value-rhs.value, lhs.diff-rhs.diff);
}


MatrixScalarDiff operator+ (const MatrixScalarDiff &lhs, const MatrixScalarDiff &rhs)
{
	if ( (lhs.diff.m != rhs.diff.m) || (lhs.diff.n != rhs.diff.n) )
	{
		printf("ERROR!! input dimension mismatrch for MatrixScalarDiff + MatrixScalarDiff\n");
		exit(1);
	}
	return MatrixScalarDiff(lhs.value+rhs.value, lhs.diff+rhs.diff);
}

MatrixScalarDiff operator* (const double &lhs, const MatrixScalarDiff &rhs)
{
	return MatrixScalarDiff(lhs*rhs.value, lhs*rhs.diff);
}


MatrixScalarDiff operator* (const MatrixScalarDiff &lhs, const MatrixScalarDiff &rhs)
{
	if ( (lhs.diff.m != rhs.diff.m) || (lhs.diff.n != rhs.diff.n) )
	{
		printf("ERROR!! input dimension mismatrch for MatrixScalarDiff * MatrixScalarDiff\n");
		exit(1);
	}
	return MatrixScalarDiff(lhs.value*rhs.value, rhs.value*lhs.diff+lhs.value*rhs.diff);
}

MatrixScalarDiff sqrt(const MatrixScalarDiff &X)
{
	return MatrixScalarDiff(sqrt(X.value), 0.5*(1.0/sqrt(X.value))*X.diff);
}

#endif
