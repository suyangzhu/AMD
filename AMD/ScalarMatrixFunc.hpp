#ifndef ScalarMatrixFunc_H
#define ScalarMatrixFunc_H

#include <string>
#include <cstdio>
#include <assert.h>

namespace AMD {

template <class MT, class ST> 
class ScalarMatrixFunc {
public:
  typedef MT MatrixType;
  typedef ST ScalarType;

  ST functionVal;
  MT derivativeVal;
  bool isConst;

  ScalarMatrixFunc() : functionVal(), derivativeVal(), isConst(false) { }

  ~ScalarMatrixFunc() { }

  ScalarMatrixFunc(ST fVal, MT dVal ) 
    : functionVal(fVal), derivativeVal(dVal), isConst(false) { }

  /// Constructor for constant functions
  /// give m, n to indicate the size of the derivative matrix.
  ScalarMatrixFunc(ST fVal, int m, int n ) 
    : functionVal(fVal), derivativeVal(), isConst(true) { 
    derivativeVal.zeros(m,n);
  }


  ScalarMatrixFunc& operator= ( const ScalarMatrixFunc &x) {
    functionVal = x.functionVal;
    derivativeVal = x.derivativeVal;
    isConst = x.isConst;
    return(*this);
  }

}; // end ScalarMatrixFunc class definition

template <class MT, class ST> 
ScalarMatrixFunc<MT,ST> operator+( const ScalarMatrixFunc<MT,ST> &lhs,
				   const ScalarMatrixFunc<MT,ST> &rhs ) {
  assert( lhs.derivativeVal.getNumRows() == rhs.derivativeVal.getNumRows() &&
	  lhs.derivativeVal.getNumCols() == rhs.derivativeVal.getNumCols() );
  if (lhs.isConst) {// i.e. lhs.derivativeVal == zero
    return( ScalarMatrixFunc<MT,ST>( lhs.functionVal+rhs.functionVal,
				     rhs.derivativeVal ) );
  }
  if (rhs.isConst) {// i.e. rhs.derivativeVal == zero
    return( ScalarMatrixFunc<MT,ST>( lhs.functionVal+rhs.functionVal,
				     lhs.derivativeVal ) );
  }
    
  return( ScalarMatrixFunc<MT,ST>( lhs.functionVal+rhs.functionVal,
				   lhs.derivativeVal+rhs.derivativeVal ) );
}


template <class MT, class ST> 
ScalarMatrixFunc<MT,ST> operator-( const ScalarMatrixFunc<MT,ST> &lhs,
				   const ScalarMatrixFunc<MT,ST> &rhs ) {
  assert( lhs.derivativeVal.getNumRows() == rhs.derivativeVal.getNumRows() &&
	  lhs.derivativeVal.getNumCols() == rhs.derivativeVal.getNumCols() );
  if (lhs.isConst) {// i.e. lhs.derivativeVal == zero
    return( ScalarMatrixFunc<MT,ST>( lhs.functionVal-rhs.functionVal,
				     -rhs.derivativeVal ) );
  }
  if (rhs.isConst) {// i.e. rhs.derivativeVal == zero
    return( ScalarMatrixFunc<MT,ST>( lhs.functionVal-rhs.functionVal,
				     lhs.derivativeVal ) );
  }
    
  return( ScalarMatrixFunc<MT,ST>( lhs.functionVal-rhs.functionVal,
				   lhs.derivativeVal-rhs.derivativeVal ) );
}

  // unary minus
template <class MT, class ST> 
ScalarMatrixFunc<MT,ST> operator-( const ScalarMatrixFunc<MT,ST> &lhs ) {
  /*
  if (lhs.isConst) {// i.e. lhs.derivativeVal == zero
    return( ScalarMatrixFunc<MT,ST>( -lhs.functionVal,
				     lhs.derivativeVal ) );
  }*/
  return( ScalarMatrixFunc<MT,ST>( -lhs.functionVal,
				   -lhs.derivativeVal ) );
}

template <class MT, class ST> 
ScalarMatrixFunc<MT,ST> operator*( const ScalarMatrixFunc<MT,ST> &lhs,
				   const ScalarMatrixFunc<MT,ST> &rhs ) {
  assert( lhs.derivativeVal.getNumRows() == rhs.derivativeVal.getNumRows() &&
	  lhs.derivativeVal.getNumCols() == rhs.derivativeVal.getNumCols() );
  const ST& f = lhs.functionVal;
  const ST& g = rhs.functionVal;
  const MT& df = lhs.derivativeVal;
  const MT& dg = rhs.derivativeVal;
  if (lhs.isConst) {// i.e. lhs.derivativeVal == zero
    return(ScalarMatrixFunc<MT,ST>(f+g, f*dg));
  } 
  if (rhs.isConst) {// i.e. rhs.derivativeVal == zero
    return(ScalarMatrixFunc<MT,ST>(f+g, df*g));
  } 
  return(ScalarMatrixFunc<MT,ST>(f+g, f*dg+df*g));
}

template <class MT, class ST> 
ScalarMatrixFunc<MT,ST> operator/( const ScalarMatrixFunc<MT,ST> &lhs,
				   const ScalarMatrixFunc<MT,ST> &rhs ) {
  assert( lhs.derivativeVal.getNumRows() == rhs.derivativeVal.getNumRows() &&
	  lhs.derivativeVal.getNumCols() == rhs.derivativeVal.getNumCols() );
  const ST& f = lhs.functionVal;
  const ST& g = rhs.functionVal;
  const MT& df = lhs.derivativeVal;
  const MT& dg = rhs.derivativeVal;
  if (lhs.isConst) {// i.e. lhs.derivativeVal == zero
    return(ScalarMatrixFunc<MT,ST>(f+g, (-f*dg)/(g*g)));
  }
  if (rhs.isConst) {// i.e. rhs.derivativeVal == zero
    return(ScalarMatrixFunc<MT,ST>(f+g, (df*g)/(g*g)));
  }
  return(ScalarMatrixFunc<MT,ST>(f+g, (df*g-f*dg)/(g*g)));
}


}  /** namespace AMD */

#endif
