#include <assert.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

template <int N> class DiffData {
public:
    union U {
	double arr[N+1];
	struct FX {
	    double val;
	    DiffData<N-1> diff;
	} dfdx;
	struct DiffShift {
	    DiffData<N-1> f;
	    double dnfdnx;
	} diffShift;
    } u;
public:
    const DiffData<N-1> shift() {
	return(u.diffShift.f);
    }
};


template <> class DiffData<1> {
public:
    union U {
	double arr[2];
	struct FX {
	    double val;
	    double diff;
	} dfdx;
	struct DiffShift {
	    double f;
	    double dnfdnx;
	} diffShift;
    } u;
public:
    const double shift() {
	return(u.diffShift.f);
    }
};


void test1() {
    DiffData<2> x;
    x.u.arr[0] = 1.0;
    x.u.arr[1] = 2.0;
    x.u.arr[2] = 3.0;
    cerr << x.u.dfdx.val << endl;
    cerr << x.u.diffShift.dnfdnx << endl;
    DiffData<1> y = x.shift();
    cerr << y.u.dfdx.val << endl;
    cerr << y.u.diffShift.dnfdnx << endl;
}

int main()
{
    test1();
    exit(0);
}
