#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "MatrixUtil.h"
#include "RDiff.h"
#include "MatrixScalarDiff.h"
#include "OptimizationUtil.h"

using namespace std;

/*****
 * This is the program to solve min_X logdet(X) + trace(S*X) + ||X||_F^2 st X>0
 *****///


int main(int argc, char* argv[]){
	int p = 500;
	double lambda = 1000;
	srand(0);
	mat_t xx(p,p,MAT_RANDOM);
	mat_t S = xx*transpose(xx)+eye(p);
	S.setname("S");
//	S.print(stdout);
	S.save_ascii("S_test");
//	mat_t X = eye(p);
	mat_t X(p,p,MAT_RANDOM);
	X = X*transpose(X)+eye(p);
	X.setname("X");
	X.save_ascii("X_test");

//	X.print(stdout);
	// Gradient descent
	int maxiter = 100;
	for ( int iter=0 ; iter<maxiter ; iter++ )
	{
//		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) -logdet(RDiff(X)) ;
		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) + lambda*trace(transpose(RDiff(X))*RDiff(X))-logdet(RDiff(X));
//		result.print(stdout);
		double nowvalue = result.value;
		double gnorm = fnorm(result.diff);
		mat_t truegradient = inv(X)-transpose(inv(X)); 
		printf("iter %d obj val: %.13lf, gradient norm = %.13lf Xnorm = %.13lf , truegradientnorm %.13lf\n", iter, nowvalue, pow(gnorm,2), fnorm(X), fnorm(truegradient));
		if ( gnorm < 1e-5)
			break;
		printf("%s\n", result.diff.val_s.val.c_str());

		vector<mat_t> const_list(1);
		const_list[0] = S;
		printf("gradient: %s\n", result.diff.val_s.val.c_str());
		RDiff *hessian = parsestring(result.diff.val_s.val, X,const_list );
		mat_t D;
//		G.setname("G");
		CGSolve(*hessian, result.diff, D, X);
//		D.print(stdout);

		double alpha = 1;
		for ( int lineiter = 0 ; lineiter<20 ; lineiter++, alpha/=2 )
		{
			mat_t newX = X - alpha*D;

			double newvalue = (-1)*logdet(newX) + trace(S*newX) + lambda*trace(transpose(newX)*newX);
			printf("newvlaue: %lf\n", newvalue);
			if ( newvalue < result.value )
			{
				X = (newX+transpose(newX))/2;
				X.setname("X");
				break;
			}
		}
	}
	return 0;
}

