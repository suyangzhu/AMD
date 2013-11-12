#ifndef Diff2_H
#define Diff2_H

#include <iostream>
using namespace std;

#include "boost/operators.hpp"

template <int n> class Diff;

template <> class Diff<1> : 
  boost::arithmetic< Diff<1>,
    boost::arithmetic< Diff<1>, double > >
{
private:
    double val, diff;
public:
    typedef double DiffT;
    Diff() : val(), diff() {}
    Diff(double val_, double diff_) : val(val_), diff(diff_) {}
    Diff(const Diff<1>& x) : val(x.getVal()), diff(x.getDiff()) {}
    Diff(const double& x) : val(x), diff() {}
    void shift(const double val_, const Diff<1>& x) { 
	val = val_; 
	diff = x.getVal();

    }
    double shift() const {
	return val;
    }


    void print() const {
	cout << val << " " << diff << endl;
    }

    double getVal() const { return val; }
    double get(int k) const
    {
	if (k==0) {
	    return(val);
	} else {
	    if (k==1) {
		return(diff);
	    }
	}
	// If we come here the user has asked for a
	// derivative that is not available.
	assert(0);
	return(0.0);
    }
    void setVal(double val_) { val = val_; }
    void set(const int k, const double x)
    {
	if (k==0) {
	    val = x;
	} else {
	    if (k==1) {
		diff = x;
	    } else {
		assert(0);
	    }
	}
    }
    void setDiff(double diff_) { diff = diff_; }
    double getDiff() const { return diff; }
    Diff<1>& operator= (const Diff<1>& x) {
	val  = x.getVal();
	diff = x.getDiff();
	return(*this);
    }

    Diff<1>& operator= (const double& x) {
	val  = x;
	diff = 0.0;
	return(*this);
    }

    // (f(x)*g(x))' = f'(x) * g(x) + g'(x) * f(x)
    Diff<1>& operator *= (const Diff<1>& y) {
	diff  = val * y.getDiff() + y.getVal() * diff;
	val  *= y.getVal();
	return(*this);
    }

    // (f(x)/g(x))' = f'(x) / g(x) - f(x) * g'(x) / g(x)*g(x) 
    //              = (f'(x) * g(x) - f(x) * g'(x) ) / g(x)^2
    Diff<1>& operator /= (const Diff<1>& y) {

	diff  =( y.getVal() * diff - val * y.getDiff() ) 
	    / ( y.getVal() * y.getVal() ) ;
	val  /= y.getVal();
	return(*this);
    }

    Diff<1>& operator += (const Diff<1>& y) {
	diff += y.getDiff();
	val  += y.getVal();
	return(*this);
    }

    Diff<1>& operator -= (const Diff<1>& y) {
	diff -= y.getDiff();
	val  -= y.getVal();
	return(*this);
    }

    Diff<1>& operator += (const double& y) {
	val  += y;
	return(*this);
    }
    Diff<1>& operator -= (const double& y) {
	val  -= y;
	return(*this);
    }
    Diff<1>& operator *= (const double& y) {
	val  *= y;
	diff *= y;
	return(*this);
    }
    Diff<1>& operator /= (const double& y) {
	val  /= y;
	diff /= y;
	return(*this);
    }

};

template <int n> class Diff : 
  boost::arithmetic<Diff<n>,
      boost::arithmetic< Diff<n>, double > >
{
private:
    double val;
    Diff<n-1> diff;
public:
    typedef Diff<n-1> DiffT;
    Diff() : val(), diff() {}
    Diff(double val_, Diff<n-1> diff_) : val(val_), diff(diff_) {}
    Diff(const Diff<n>& x) : val(x.getVal()), diff(x.getDiff()) {}
    Diff(const double& x) : val(x), diff() {}
    void shift(const double val_, const Diff<n>& x) { 
	val = val_; 
	diff.shift(x.getVal(),x.getDiff());
    }
    Diff<n-1> shift() const {
	Diff<n-1> tmp;
	tmp.shift(val,diff);
	return tmp;
    }

    void print() const {
	cout << val << " ";
	diff.print();
    }

    double getVal() const { return val; }
    double get(const int k) const {
	if (k==0) {
	    return(val);
	}
	return(diff.get(k-1));
    }
    void setVal(double val_) { val = val_; }
    void set(const int k, const double x) {
	if (k==0) {
	    val = x;
	} else {
	    diff.set(k-1,x);
	}
    }

    void setDiff(Diff<n-1> diff_) { diff = diff_; }
    Diff<n-1> getDiff() const { return diff; }
    Diff<n>& operator= (const Diff<n>& x) {
	val  = x.getVal();
	diff = x.getDiff();
	return(*this);
    }

    Diff<n>& operator= (const double& x) {
	val  = x;
	diff = Diff<n-1>();
	return(*this);
    }

    // (f(x)*g(x))' = f'(x) * g(x) + g'(x) * f(x)
    Diff<n>& operator *= (const Diff<n>& y) {
	double nchoosek[n+1];
	double fk[n+1];
	double yk[n+1];
	for (int k=0;k<n+1;k++) {
	    nchoosek[k] = 0.0;
	    fk[k] = get(k);
	    yk[k] = y.get(k);
	}
	nchoosek[0] = 1.0;

	val *= y.getVal();
	// compute binomial coefficients and non-recursive
	// implementation of higher order derivatives
	// (f*g)^(k) = \sum_i (i choose k) f^(i) g^(k-i)
	for (int k=0;k<n;k++) {
	    for (int i=k+1;i>0;i--) {
		nchoosek[i] += nchoosek[i-1];
	    }
	    double kderiv = 0.0;
	    for (int i=0;i<=k+1;i++) {
		kderiv += nchoosek[i] * fk[i] * yk[k+1-i];
	    }
	    set(k+1,kderiv);
	}
	return(*this);
    }
    /*
    // (f(x)*g(x))' = f'(x) * g(x) + g'(x) * f(x)
    Diff<n>& operator *= (const Diff<n>& y) {
	Diff<n-1> xS,yS;
	xS.shift(val,diff);
	yS.shift(y.getVal(),y.getDiff());
	diff  = xS * y.getDiff() + yS * diff;
	val  *= y.getVal();
	return(*this);
    }
    */

    // (f(x)/g(x))' = f'(x) / g(x) - f(x) * g'(x) / g(x)*g(x) 
    //              = (f'(x) * g(x) - f(x) * g'(x) ) / g(x)^2
    Diff<n>& operator /= (const Diff<n>& y) {
	Diff<n-1> xS,yS;
	xS.shift(val,diff);
	yS.shift(y.getVal(),y.getDiff());
	diff  =( yS * diff - xS * y.getDiff() ) 
	    / ( yS * yS ) ;
	val  /= y.getVal();
	return(*this);
    }

    Diff<n>& operator += (const Diff<n>& y) {
	diff += y.getDiff();
	val  += y.getVal();
	return(*this);
    }
    
    Diff<n>& operator -= (const Diff<n>& y) {
	diff -= y.getDiff();
	val  -= y.getVal();
	return(*this);
    }
    Diff<n>& operator += (const double y) {
	val  += y;
	return(*this);
    }
    Diff<n>& operator -= (const double y) {
	val  -= y;
	return(*this);
    }
    Diff<n>& operator *= (const double y) {
	val  *= y;
	diff *= y;
	return(*this);
    }
    Diff<n>& operator /= (const double y) {
	val  /= y;
	diff /= y;
	return(*this);
    }

}; 



#endif
