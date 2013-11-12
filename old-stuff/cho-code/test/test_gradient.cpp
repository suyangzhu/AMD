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
	int p = 100;
	double lambda = 1000;
	mat_t xx(p,p,MAT_RANDOM);
	mat_t S = xx*transpose(xx)+eye(p);
	S.setname("S");
//	S.print(stdout);
	mat_t X = eye(p);
	X.setname("X");

	mat_t U(p,p,MAT_EMPTY);
	mat_t eigvalues(p,1,MAT_EMPTY);
	// Gradient descent
	int maxiter = 200;
	for ( int iter=0 ; iter<maxiter ; iter++ )
	{
//		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) -logdet(RDiff(X)) ;
		MatrixScalarDiff result = trace(RDiff(S,OP_CONST)*RDiff(X)) + lambda*trace(transpose(RDiff(X))*RDiff(X))-logdet(RDiff(X)) ;
//		printf("ADRESULT: \n");
//		X.print(stdout);
//		result.diff.print(stdout);
		double nowvalue = result.value;
		double gnorm = fnorm(result.diff);
		printf("iter %d obj val: %lf, gradient norm = %lf\n", iter, nowvalue, gnorm);
		if ( gnorm < 1e-2)
			break;
		printf("%s\n", result.diff.val_s.val.c_str());
		double alpha = 1;
		for ( int lineiter = 0 ; lineiter<20 ; lineiter++, alpha/=2 )
		{
			mat_t newX = X - alpha*result.diff;

			// projection
/*			eig(newX, U, eigvalues);
			for ( int i=0 ; i<p ; i++ )
				if (eigvalues.val[i] < 0 )
					eigvalues.val[i] = 0;
*/
//			double newvalue = (-1)*logdet(newX) + trace(S*newX) ;
			double newvalue = (-1)*logdet(newX) + trace(S*newX) + lambda*trace(transpose(newX)*newX);
			printf("newvlaue: %lf\n", newvalue);
			if ( newvalue < result.value )
			{
//				printf("dsafdsaf\n");
//
//				X.print(stdout);
//				newX.print(stdout);
				X = (newX+transpose(newX))/2;
				X.setname("X");
				break;
			}
		}
	}
//X.print(stdout);
	return 0;
}

