#include <assert.h>
#include <iostream>
#include <math.h>
#include "Diff2.H"

using namespace std;

void newton1();
void newton2();

int main()
{
  newton1();
  newton2();
  exit(0);
}

template <typename T>
T f(T& x) {
    T y = 1.0 / (1.0 + x) - x;
    //    T y = exp( - x * x / 2.0 )/sqrt(2.0*M_PI);
    return(y);
}

void newton1() {
    Diff<1> x(1.0,1.0);
    cout << "Newton's algorithm to solve x = 1/(1+x), i.e. f(x) = 0 with f(x) = 1/(1+x) - x." << endl;
    for (int i=0;i<10;i++) {
	Diff<1> y = f(x);
	cout << "x = " << x.getVal() << " f(x) = " << y.getVal() << endl;
	double z = x.getVal() - y.getVal() / y.getDiff();
	x.setVal( z );
    }
}


template <typename T>
T g(T& x) {
    T y = x / (1.0 + x * x);
    //    T y = exp( - x * x / 2.0 )/sqrt(2.0*M_PI);
    return(y);
}


void newton2() {
    Diff<2> x(0.5,1.0);
    cout << "Newton's algorithm to find critical points of x/(1+x*x)." << endl;
    for (int i=0;i<10;i++) {
	Diff<2> y = g(x);
	cout << "x = " << x.get(0) << " f(x) = " << y.get(0) << endl;
	double z = x.get(0) - y.get(1) / y.get(2);
	x.setVal( z );
    }
}
