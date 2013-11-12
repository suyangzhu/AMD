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


using namespace std;

/*****
 * This is the program to solve min_X logdet(X) + trace(S*X) + ||X||_F^2 st X>0
 * *///


int main(int argc, char* argv[]){
	int p = 4;
	double lambda = 1;
	mat_t xx(p,p,MAT_RANDOM);
	mat_t S = xx*transpose(xx)+eye(p);
	S.setname("S");
	S.print(stdout);
//	mat_t X = eye(p);
	mat_t X(p,p,MAT_RANDOM);
	X = X*transpose(X)+eye(p);
	X.setname("X");

	X.print(stdout);
	// Gradient descent
	int maxiter = 1;
	for ( int iter=0 ; iter<maxiter ; iter++ )
	{
//		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) -logdet(RDiff(X)) ;
		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) + lambda*trace(transpose(RDiff(X))*RDiff(X))-logdet(RDiff(X));
//		result.diff.print(stdout);
		double nowvalue = result.value;
		result.diff.print(stdout);

		vector<mat_t> const_list(1);
		const_list[0] = S;
//		string tmp = "((1.00*(X+X))-inv(X))";
//		RDiff *hessian = parsestring(tmp, X,const_list );
		RDiff *hessian = parsestring(result.diff.val_s.val, X,const_list );
		mat_t D(p,p,MAT_RANDOM);
		D.setname("D");
		D.print(stdout);
		mat_t E = GradientVectorProduct(*hessian, D);
		E.print(stdout);

/*
		double gnorm = fnorm(result.diff);
		printf("iter %d obj val: %lf, gradient norm = %lf\n", iter, nowvalue, gnorm);
		if ( gnorm < 1e-2)
			break;
		printf("%s\n", result.diff.val_s.val.c_str());
		double alpha = 1;
		for ( int lineiter = 0 ; lineiter<20 ; lineiter++, alpha/=2 )
		{
			mat_t newX = X - alpha*result.diff;

			double newvalue = (-1)*logdet(newX) + trace(S*newX) + lambda*trace(transpose(newX)*newX);
			printf("newvlaue: %lf\n", newvalue);
			if ( newvalue < result.value )
			{
				X = (newX+transpose(newX))/2;
				X.setname("X");
				break;
			}
		}

	*/
	}
	return 0;
}

