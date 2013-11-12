#ifndef DiffMath_H
#define DiffMath_H

#include <cmath>

// (e^{f(x)})' = f'(x) e^f(x)
template <typename T> T exp(const T& f) {
    T y( exp( f.getVal() ),
	 f.getDiff() * 
	 exp( f.shift() ) 
	 );
    return(y);
}

// (log(f(x))' = f'(x) / f(x)
template <typename T> T log(const T& f) {
    T y( log( f.getVal() ),
	 f.getDiff() / f.shift() 
	 );
    return(y);
}

// (f(x)^a)' = a f'(x) (f(x))^{a-1}   if a != 0
//           = 0                    if a == 0 
template <typename T> T pow(const T& f, double a) {
    if (a==0.0) return(T());
    T z( pow( f.getVal(), a ), 
	 a * f.getDiff() * pow( f.shift(), a-1 )
	 );
    return(z);
}

// (sqrt(f(x))' = 1/(2*f'(x)* sqrt(f(x)))
template <typename T> T sqrt(const T& f) {
    T z( sqrt( f.getVal() ), 
	 0.5 / ( f.getDiff() * sqrt( f.shift() ) )
	 );
    return(z);
}


// (f(x)^{1/3})' = 1/(3*f'(x)* (f(x)^2/3))
template <typename T> T cbrt(const T& f) {
    T z( cbrt( f.getVal() ), 
	 (1.0/3.0) / ( f.getDiff() * cbrt( f.shift() ) * cbrt( f.shift() ) )
	 );
    return(z);
}


// (sqrt(f(x))' = 1/(2*f'(x)* sqrt(f(x)))
template <typename T> T erf(const T& f) {
    T z( erf( f.getVal() ), 
	 (2.0 / sqrt(M_PI)) * exp ( -1.0 *  f.shift() * f.shift() )
	 );
    return(z);
}

// (sqrt(f(x))' = 1/(2*f'(x)* sqrt(f(x)))
template <typename T> T erfc(const T& f) {
    return(1.0 - erf(f));
}



// (sin(f(x))' = f'(x) cos(f(x))
template <typename T> T sin(const T& f) {
    T y( sin( f.getVal() ),
	 f.getDiff() * 
	 cos( f.shift() ) 
	 );
    return(y);
}


// (cos(f(x))' = -f'(x) sin(f(x))
template <typename T> T cos(const T& f) {
    T y( cos( f.getVal() ),
	 -1.0 * f.getDiff() * 
	 sin( f.shift() ) 
	 );
    return(y);
}

// (tan(f(x))' = f'(x) / cos^2(f(x))
template <typename T> T tan(const T& f) {
    typename T::DiffT z( cos( f.shift() ) );
    T y( tan( f.getVal() ),
	 f.getDiff() / ( z * z )
	 );
    return(y);
}



#endif
