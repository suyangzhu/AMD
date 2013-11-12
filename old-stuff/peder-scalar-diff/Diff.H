#include "boost/operators.hpp"

template <class T> class Diff;
template <typename T> class DS;

template <typename T> class Diff : 
  boost::multipliable<Diff<T>, 
    boost::multipliable<Diff<T>, double,
      boost::addable<Diff<T> > > > {
private:
    double value_;
    T derivative_;
public:
    inline double getValue() const { return value_; }
    inline T getDerivative() const { return derivative_; }
    inline void setValue(double value) { 
	value_ = value;
	derivative_ = 1.0;
    }
    inline void setValue(double value,T derivative) { 
	value_ = value;
	derivative_ = derivative;
    }

    // default constructors initializes f(x) = 0, and f'(x)=0
    Diff() : value_(), derivative_(0) {};
    // Diff(value) sets x = value and f(x) = x, f'(x) = 1
    explicit Diff(const double value) : value_(value){ derivative_ = 1.0; };
    Diff(const double value, const T derivative) : value_(value), derivative_(derivative) { };

    // choose f(x) and df/dx
    Diff(const Diff<T>& y) : value_(y.getValue()), 
			     derivative_(y.getDerivative()) {};
    
    inline Diff<T>& operator= (const double& x) {
	value_ = x;
	derivative_ = T();
	return(*this);
    }

    inline Diff<T>& operator= (const Diff<T>& x) {
	value_ = x.getValue();
	derivative_ = x.getDerivative();
	return(*this);
    }

    // multiplication of constant:
    // d/dx { v * f(x) )  = v * df/dx
    inline Diff<T>& operator *= (const double& v) {
	value_      *= v;
	derivative_ *= v;
	return(*this);
    }

    // multiplication of functions:
    // d/dx { y(x) * f(x) )  = y(x) * df/dx + dy/dx * f(x)
    inline Diff<T>& operator *= (const Diff<T>& y);

    // multiplication of functions:
    // d/dx { y(x) + f(x) )  = y'(x) + f'(x)
    template <typename T2>
    inline Diff<T>& operator += (const Diff<T2>& y) {
	value_      += y.getValue();
	derivative_ += y.getDerivative();
	return(*this);
    }
    
}; 

template <> class DS<double> {
public:
    double get(Diff<double> f) {
	return f.getValue();
    }
};

template <typename T> class DS<Diff<T> > {
public:
    Diff<T> get(Diff<Diff<T> > f) {
	DS<T> x;
	return Diff<T>(f.getValue(), x.get(f.getDerivative()));
    }
};


template <typename T> inline 
Diff<T>& Diff<T>::operator *= (const Diff<T>& y) {
    DS<double> dshifter();
    double x = dshifter.get(Diff<double>());
    //    T yDs = dshifter.get(y);
    derivative_ = 
	value_ * y.getDerivative() 
	+ y.value_ * derivative_;
    value_      *= y.getValue();
    return(*this);
}



