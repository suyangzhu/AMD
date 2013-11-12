#ifndef Diff_H
#define Diff_H

#include "MatrixUtil.h"
using namespace std;

class Diff  
{

	public:
		mat_t value, diff;
		int in_m, in_n;
		int out_m, out_n;
	
		Diff(){}
		~Diff(){}

		// f(X) = X
		Diff(mat_t &x_):value(x_), diff((x_.m)*(x_.n)){
			in_m = x_.m;
			in_n = x_.n;
			out_m = x_.m;
			out_n = x_.n;
		}

		// Copy from an existing Diff object
		Diff(const Diff &x_) {
			in_m = x_.in_m;
			in_n = x_.in_n;
			out_m = x_.out_m;
			out_n = x_.out_n;
			value = x_.value;
			diff = x_.diff;
		}

		mat_t getVal() const  { return value; }

		mat_t getDiff() const  { return diff; }

		Diff& operator= ( const Diff &x) {
			in_m = x.in_m;
			in_n = x.in_n;
			out_m = x.out_m;
			out_n = x.out_n;
			value = x.getVal();
			diff = x.getDiff();
			return(*this);
		}

		Diff& operator+= ( const Diff &x) {
			if ( (in_m!=x.in_m) || (in_n!=x.in_n) || (out_m != x.out_m) || (out_n != x.out_n) )
			{
				printf("ERROR on += on Diff: dimension mismatch");
				exit(1);
			}
			diff = diff + x.getDiff();
			value = value + x.getVal();
			return (*this);
		}

//		Diff& operator*= (const Diff &x) {
//		}
};

Diff logdet(Diff X) {
	if ( X.out_m != X.out_n)
	{
		printf("ERROR!! logdet with out_m=%d, out_n=%d\n", X.out_m, X.out_n);
		exit(1);
	}
	Diff y;
	y.in_m = X.in_m;
	y.in_n = X.in_n;
	y.out_m = 1;
	y.out_n = 1;
	y.value = logdet(X.value);
	mat_t tmpvec1 = inv(X.value);
	mat_t tmpvec = reshape(tmpvec1,1,(X.out_m)*(X.out_n));
	mat_t aaa = tmpvec*(X.diff) ;
	y.diff = transpose(reshape(aaa, X.in_m, X.in_n));
	return y;
}

Diff trace(Diff X) {
	if ( X.out_m != X.out_n)
	{
		printf("ERROR using trace: matrix dimension (%d,%d)\n", X.out_m, X.out_n);
		exit(1);
	}

	Diff y;
	y.in_m = X.in_m;
	y.in_n = X.in_n;
	y.out_m = 1;
	y.out_n = 1;
	y.value = trace(X.value);
	mat_t vecI = form_vecI(X.out_m);
	mat_t tmpaa = transpose(X.diff);
	mat_t tmp2 = tmpaa*vecI;
	y.diff = transpose(reshape(tmp2, X.in_m, X.in_n));
	return y;
}

Diff operator* (const mat_t &lhs, const Diff &rhs) 
{
	Diff result;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	result.out_m = lhs.m;
	result.out_n = rhs.out_n;
	result.value = lhs*(rhs.value);
	result.diff = kron(lhs,mat_t(lhs.n))*rhs.diff;
	return result;
}

Diff operator* (const Diff &lhs, const mat_t &rhs) 
{
	Diff result;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	result.out_m = lhs.out_m;
	result.out_n = rhs.n;
	result.value = (lhs.value)*rhs;
	result.diff = kron(mat_t(lhs.out_m),transpose(rhs))*lhs.diff;
	return result;
}


Diff operator+(const Diff &lhs, const Diff &rhs) 
{
	Diff result(lhs);
	result+=rhs;
	return result;
}

Diff transpose(Diff X)
{
	Diff result;
	result.in_m = X.in_m;
	result.in_n = X.in_n;
	result.out_m = X.out_n;
	result.out_n = X.out_m;
	result.value = transpose(X.value);
	result.diff = box(eye(X.out_n),eye(X.out_m))*X.diff;
	return result;
}

Diff operator* (const Diff &lhs, const Diff &rhs) 
{
	if ( (lhs.in_m != rhs.in_m) || (lhs.in_n != rhs.in_n))
	{
		printf("ERROR!! * in Diff, input dimension mismatch\n");
		exit(1);
	}
	if ( lhs.out_n != rhs.out_m )
	{
		printf("ERROR!! * in Diff, output dimension error. A(%d %d) * B(%d %d)\n", lhs.out_m, lhs.out_n, rhs.out_m, rhs.out_n);
		exit(1);
	}
	Diff result;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	result.out_m = lhs.out_m;
	result.out_n = rhs.out_n;
	result.value = lhs.value*rhs.value;
	result.diff = kron(mat_t(lhs.out_m),transpose(rhs.value))*lhs.diff + kron(lhs.value,mat_t(rhs.out_n))*rhs.diff;
	return result;
}

Diff inv(Diff &X)
{
	if ( X.out_m != X.out_n) 
	{
		printf("ERROR inv(Diff), output dimension (%d %d)\n", X.out_m, X.out_n);
		exit(1);
	}
	Diff result;
	result.in_m = X.in_m;
	result.in_n = X.in_n;
	result.out_m = X.out_m;
	result.out_n = X.out_n;
	result.value = inv(X.value);
	result.diff = (-1)*kron(inv(X.value), transpose(inv(X.value)))*X.diff;
	return result;
}

Diff operator* (const double lhs, const Diff &rhs) 
{
	Diff result(rhs);
	result.value = lhs*rhs.value;
	result.diff = lhs*rhs.diff;
	return result;
}

Diff operator* (const Diff &lhs, const double rhs) 
{
	Diff result(lhs);
	result.value = rhs*lhs.value;
	result.diff = rhs*lhs.diff;
	return result;
}


#endif
