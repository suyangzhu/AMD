#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <cassert>
#include "MatrixUtil.h"
#include "Diff.h"
#include "RDiff.h"
#include "MatrixScalarDiff.h"

using namespace std;

int test_fnorm(mat_t &X, mat_t &S)
{
	// F'(X) = d/dx (sqrt(tr(X^T X)))
  //std::pair<mat_t,double> foo = sqrt(trace(transpose(RDiff(X))*RDiff(X)));
	MatrixScalarDiff result = sqrt(trace(transpose(RDiff(X))*RDiff(X)));
  // result.evaluate();
	mat_t trueresult = X/fnorm(X);
	printf("result: \n");
	result.print(stdout);
	mat_t diff = result.diff - trueresult;
	if ( fnorm(diff) < 1e-10 )
		return 1;
	else
		return 0;
}

int test_logdet(mat_t &X, mat_t &S)
{
	// f(X) = logdet(X) + tr(S*X)
	MatrixScalarDiff result = logdet(RDiff(X))+trace(RDiff(S,OP_CONST)*RDiff(X));
	mat_t trueresult = transpose(inv(X))+transpose(S);
	printf("result: \n");
	result.print(stdout);

	mat_t diff = result.diff - trueresult;
	if ( fnorm(diff) < 1e-10)
		return 1;
	else
		return 0;
}

int test_schur2(mat_t &X, mat_t &A)
{
	// f(X) = logdet(X.*X')
	MatrixScalarDiff result = logdet(schur(RDiff(X),transpose(RDiff(X))));
	mat_t trueresult = transpose(schur(X,transpose(inv(schur(X,transpose(X))))+inv(schur(X,transpose(X)))));
//	trueresult.print(stdout);
	printf("result: \n");
	result.print(stdout);
	mat_t diff = result.diff - trueresult;
	if ( fnorm(diff) < 1e-10)
		return 1;
	else
		return 0;
}

/*
int test_schur(mat_t &X, mat_t &A)
{
	mat_t result = logdet(schur(RDiff(A,OP_CONST),RDiff(X)));
	mat_t trueresult = schur(A,transpose(inv(schur(A,X))));
//	trueresult.print(stdout);
//	result.print(stdout);
	mat_t diff = result - trueresult;
	if ( fnorm(diff) < 1e-10)
		return 1;
	else
		return 0;
}

int test_logdet(mat_t &X) 
{
	mat_t result = logdet(RDiff(X));
	mat_t trueresult = transpose(inv(X));
	mat_t diff = result - trueresult;
	if ( fnorm(diff) < 1e-10 )
		return 1;
	else
		return 0;
}

int test_frobenius(mat_t &X) 
{
	mat_t result = trace(RDiff(X)*transpose(RDiff(X)));
	mat_t trueresult = X+X;
	mat_t diff = result - trueresult;
	if ( fnorm(diff) < 1e-10 )
		return 1;
	else
		return 0;
}


int test_traceAX(mat_t &X, mat_t &A) 
{
	mat_t result = trace(RDiff(A,OP_CONST)*RDiff(X));
	mat_t trueresult = transpose(A);
	mat_t diff = result - trueresult;
	if ( fnorm(diff) < 1e-10 )
		return 1;
	else
		return 0;
}
*/

int main(int argc, char* argv[]){

/*	mat_t A;
	A.load_ascii(argv[1]);
	RDiff a(A);
	RDiff b(A);
	RDiff d = transpose(b);
	RDiff c = d*a;
	c.print(stdout);
	printf("asdf\n");
	c.value.print(stdout);
	mat_t result = trace(c);
	result.print(stdout);
*/
	mat_t A(3,3,MAT_RANDOM), B(3,3,MAT_RANDOM);
	A.setname("A");
	B.setname("B");
	mat_t C = zeros(3,3);

//	printf("input:\n");
//	A.print(stdout);
//	printf("\n");
//	A.load_ascii(argv[1]);
//	B.load_ascii(argv[1]);
// 	assert(test_logdet(A));
//	assert(test_traceAX(A,B));
	assert(test_schur2(A,B));
	assert(test_logdet(A,B));
	assert(test_fnorm(A,B));
	printf("ALL test pass\n");
/*	RDiff a(A);
	RDiff b(A);
	RDiff d = transpose(b);
	mat_t I3 = mat_t(3);
	RDiff e(A);
	RDiff f = e+I3;
	RDiff g = inv(f);
	RDiff h = g*d;
	RDiff y = h*a;
	y.value.print(stdout);
	mat_t result = trace(y);*/
//	mat_t result = logdet(y);
//	mat_t I = eye(3);
//	mat_t result = trace(inv(RDiff(A)+eye(3))*transpose(RDiff(A))*RDiff(A));
//	printf("Gradient output: \n");
//	result.print(stdout);

//	mat_t result2 = trace(RDiff(A)*transpose(RDiff(A)));
//	result2.print(stdout);

//	mat_t result = logdet(RDiff(A));
//	result.print(stdout);

	return 0;
}

