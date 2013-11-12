#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "Diff2.H"
#include "DiffMath.H"

using namespace std;

template <typename T> T func(T& x);

void test1();
void test2();
void test3();

void normSample();
template <typename T>
T normCDF(T& x);

int main()
{
    srand48(time(0));
    normSample();
    // test1();
    //  test3();
  exit(0);
}

template <typename T>
T func(T& x) {
    // compute y = x^8
    // T y = 1.0 / x;
    T y = (sqrt(M_PI) / 2.0) * erf( x );
    //    T y = cbrt( 1.0 + x\ );
    return(y);
}

void test1() {
    Diff<1> x(1.0,1.0);
    Diff<1> y = func(x);
    cout << " x = " << x.getVal() 
	 << " f(x) = " << y.getVal() 
	 << " f'(x) = " << y.getDiff() 
	 << endl;
}

void test2() {
    Diff<2> x(1.0,Diff<1>(1.0,0.0));
    Diff<2> y = func(x);
    cout << " x = " << x.getVal() 
	 << " f(x) = " << y.getVal() 
	 << " f'(x) = " << y.getDiff().getVal() 
	 << " f''(x) = " << y.getDiff().getDiff() 
	 << endl;
}

void test3() {
    Diff<10> x(0.0,1.0);
    Diff<10> y = func(x);
    cout << " x = " << x.getVal() 
	 << " f(x) = " << y.get(0) 
	 << " f'(x) = " << y.get(1)
	 << " f''(x) = " << y.get(2)
	 << " f'''(x) = " << y.get(3)
	 << " f''''(x) = " << y.get(4)
	 << " f'''''(x) = " << y.get(5)
	 << " f''''''(x) = " << y.get(6)
	 << " f'''''''(x) = " << y.get(7)
	 << endl;
}

template <typename T>
T normCDF(T& x) {
    T y = 0.5 * ( 1.0 + erf ( x / sqrt(2.0) ) );
    return(y);
}


void normSample() {
    double y = drand48();
    double z = y;
    // solve for phi(x) = y
    for (int i=0;i<10;i++) {
	Diff<1> xFunc(z,1.0);
	Diff<1> phiFunc = normCDF(xFunc);
	z = z - (phiFunc.getVal()-y) / phiFunc.getDiff();
	cout << " z = " << z << "\ty = " << y << "\tN(z) = " 
	     << phiFunc.getVal() << endl;
    }
}


/*
void test2() {
    Diff<Diff<double> > x(1.0);
    Diff<Diff<double> > y = func(x);
    cout << " x = " << x.getValue() 
	 << " f(x) = " << y.getValue() 
	 << " f'(x) = " << y.getDerivative().getValue()
	 << " f''(x) = " << y.getDerivative().getDerivative()
	 << endl;
    x = 2.0;
    y = func(x);
    cout << " x = " << x.getValue()
	 << " f(x) = " << y.getValue() 
	 << " f'(x) = " << y.getDerivative().getValue()
	 << " f''(x) = " << y.getDerivative().getDerivative()
	 << endl;
}

void test1() {
    Diff<double> x(1.0);
    Diff<double> y = func(x);
    cout << " x = " << x.getValue() 
	 << " f(x) = " << y.getValue() 
	 << " f'(x) = " << y.getDerivative()
	 << endl;
    x = 2.0;
    y = func(x);
    cout << " x = " << x.getValue()
	 << " f(x) = " << y.getValue() 
	 << " f'(x) = " << y.getDerivative()
	 << endl;
}
*/
