#ifndef OPTIMIZATIONUTIL_H
#define OPTIMIZATIONUTIL_H

#include <cstdio>
#include <cstdlib>
#include "MatrixUtil.h"
#include "MatrixScalarDiff.h"

using namespace std;

// Solve H*vec(X) =vec(B)
int CGSolve(RDiff &H, mat_t &B, mat_t &X, mat_t &Xt, double tol = 1e-4)
{
	int m = B.m; 
	int n = B.n;
	mat_t R(m,n,MAT_ZERO), P(m,n,MAT_ZERO), HX_result(m,n,MAT_ZERO);
	int maxiter = m;

	printf("START CG Solver\n");
	// Initial from X = 0
	X = zeros(m,n);
	R = B;
	P = R;
	double R_norm = pow(fnorm(R),2);
	double initial = R_norm;
	printf("initial rnorm: %lf\n", R_norm);
	if ( R_norm < 1e-10 )
		return 0;
	int iter;
	for ( iter = 0; iter<maxiter ; iter++)
	{
		X.setname("X");
		R.setname("R");
		P.setname("P");
		HX_result = GradientVectorProduct(H, P);
		HX_result = (HX_result+transpose(HX_result))/2;
		printf("Hessian vector product: %s\n", HX_result.getname().c_str());
//		printf("norm of Hessian vector product: %lf\n", fnorm(HX_result));
		double alpha = pow(fnorm(R),2)/innerproduct(HX_result, P);
		X = X + alpha*P;
		R = R - alpha*HX_result;
		double R_norm_now = pow(fnorm(R), 2);
		printf("CG iter %d residual %lf\n", iter, R_norm_now);
		if ( R_norm_now < tol*initial)
			break;
		double beta = R_norm_now/R_norm;
		R_norm = R_norm_now;
		P = R + beta*P;
	}
	return (iter+1);
}

#endif
