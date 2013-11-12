#include "boost/operators.hpp"
#include <algorithms>

template <int n> class Diff;

template <int n> class Diff : 
  boost::arithmetic< Diff<n>,
    boost::arithmetic< Diff<n>, double > >
{
private:
    double val[n+1];
public:
    // copy constructors
    Diff() {
	fill(&val[0],&val[n+1],0.0);
    }
    Diff(const double val_[n+1]) {
	copy(&val_[0],&val_[n+1],&val[0]);
    }
    Diff(double* val_) {
	copy(&val_[0],&val_[n+1],&val[0]);
    }

    Diff(double val_) : val[0](val_) {
	fill(&val[1],&val[n+1],0.0);
    }
    Diff(double val_, double diff_) : val[0](val_), val[1](diff_) {
	fill(&val[2],&val[n+1],0.0);
    }

    const double* getData() const {
	return(&val[0]);
    }

    Diff(const Diff<n>& x) {
	const double* y = x.getData();
	copy(y,y+n+1,val);
    }

    void shift(const Diff<n+1>& x) {
	const double* y = x.getData();
	copy(y,y+n+1,val);
    }


    double get(int k) const;
    double getVal() const { return val; }
    void setVal(double val_) { val = val_; }
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
    Diff() : val(), diff() {}
    Diff(double val_, Diff<n-1> diff_) : val(val_), diff(diff_) {}
    Diff(const Diff<n>& x) : val(x.getVal()), diff(x.getDiff()) {}
    Diff(const double& x) : val(x), diff() {}
    void shift(const double val_, const Diff<n>& x) { 
	val = val_; 
	diff.shift(x.getVal(),x.getDiff());
    }
    double get(int k) const;
    double getVal() const { return val; }
    void setVal(double val_) { val = val_; }
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
	Diff<n-1> xS,yS;
	xS.shift(val,diff);
	yS.shift(y.getVal(),y.getDiff());
	diff  = xS * y.getDiff() + yS * diff;
	val  *= y.getVal();
	return(*this);
    }

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

template <int n>  double Diff<n>::get(int k) const
{
    if (k==0) {
	return(val);
    }
    return(diff.get(k-1));
}

double Diff<1>::get(int k) const
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

