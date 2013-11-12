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
#include "Diff.h"

using namespace std;

Diff test_logdet(Diff &x) {
	return logdet(x);
}

Diff test_111(Diff &x, mat_t &A, mat_t &B) {
	return A*transpose(inv(x))*B;
}

Diff test_invX(Diff &x){
	return inv(x);
}

Diff test_XTAX(Diff &x, mat_t &A) {
	return transpose(x)*A*x;
}

Diff test_XTX(Diff &x) {
	return transpose(x)*x;
}

Diff test_XX(Diff &x) {
	return x*x;
}

Diff test_AXTB(Diff &x, mat_t &A, mat_t &B) {
	return trace(A*(transpose(x))*B);
}

Diff test_transpose(Diff& x){
	return transpose(x);
}

Diff test_fun_2(Diff& x, mat_t& A, mat_t& B) {
	return trace(A*x*B);
}

Diff test_fun_1(Diff& x, mat_t& S) {
	return trace(x*S);
}


Diff test_fun_0(Diff& x) {
	return logdet(x);
}

Diff test_fun(Diff& x, mat_t& S) {
	return (logdet(x))+(trace(S*x));
}

int main(int argc, char* argv[]){
/*
	mat_t A;
	A.load_ascii(argv[1]);
	mat_t B;
	B.load_ascii(argv[2]);
	mat_t x;
	x.load_ascii(argv[3]);
	*/
	//B=A+A;
	

	// test matrix operations
	/*
	mat_t B = A;
	mat_t C = A+B;

	mat_t D = A*C;
	printf("A: \n");
	A.print(stdout);
	printf("C: \n");
	C.print(stdout);
	printf("D: \n");
	D.print(stdout);
	double deter = B.det();
	printf("det: %lf\n", deter);
	mat_t E = inv(D);
	E.print(stdout);
*/
	// test diff
//	Diff B; 

//	Diff X(A);
//	Diff result = test_fun_0(X); 
//	result.value.print(stdout);
//	result.diff.print(stdout);

	// Test logdet+trace
/*	Diff X(A);
	printf("A: \n");
	A.print(stdout);
	printf("B:\n");
	B.print(stdout);
	Diff result = test_fun(X, B); 
	result.value.print(stdout);
	result.diff.print(stdout);
*/

	// Test trace(A*X*B)

	/*
	Diff X(x);
	x.print(stdout);
	A.print(stdout);
	B.print(stdout);
	Diff result = test_fun_2(X,A,B);
	result.value.print(stdout);
	result.diff.print(stdout);
*/

	// Test f(X)=AX^T B
/*	A.print(stdout);
	B.print(stdout);
	Diff X(x);
	x.print(stdout);
	Diff result = test_AXTB(X, A, B);
	result.value.print(stdout);
	result.diff.print(stdout);
*/

	// Test f(X) = XX
/*	Diff X(x);
	x.print(stdout);
	Diff result = test_XX(X);
	result.value.print(stdout);
	result.diff.print(stdout);
*/
	// Test f(X) = X^T X
/*	Diff X(x);
	x.print(stdout);
	Diff result = test_XTX(X);
	result.value.print(stdout);
	result.diff.print(stdout);
*/
	// Test f(X) = X^T A X
/*	Diff X(x);
	A.print(stdout);
	x.print(stdout);
	Diff result = test_XTAX(X,A);
	result.value.print(stdout);
	result.diff.print(stdout);
*/
	// Test f(X) = X^{-1}
/*	Diff X(x);
	x.print(stdout);
	Diff result = test_invX(X);
	result.value.print(stdout);
	result.diff.print(stdout);*/
/*
	Diff X(x);
	x.print(stdout);
	A.print(stdout);
	B.print(stdout);
	Diff result = test_111(X, A, B);
	result.value.print(stdout);
	result.diff.print(stdout);
*/
	mat_t x;
	x.load_ascii(argv[1]);
	Diff X(x);
	x.print(stdout);
	Diff result = test_logdet(X);
	result.value.print(stdout);
	result.diff.print(stdout);
	return 0;
}

