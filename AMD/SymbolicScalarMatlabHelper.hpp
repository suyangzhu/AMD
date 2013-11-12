#ifndef AMD_SYMBOLIC_SCALAR_MATLAB_HPP
#define AMD_SYMBOLIC_SCALAR_MATLAB_HPP

#include "SymbolicScalarMatlab.hpp"

namespace AMD {

  /** 
   * Compute sqrt(a).  If a is trace(X) then the return value will hold
   * sqrt(trace(X)).
   * @param[in] a The symbolic argument passed to sqrt.
   * @return A new SymbolicScalarMatlab.
   */
  SymbolicScalarMatlab sqrt(const SymbolicScalarMatlab& a) {
    return SymbolicScalarMatlab("sqrt("+detail::removeParenthesis(a.symbol)+")");
  }

  /** 
   * Adds two scalars.
   * @param[in] a The first scalar argument.
   * @param[in] b The second scalar argument.
   * @return The sum a+b.
   */
  SymbolicScalarMatlab operator+(const SymbolicScalarMatlab& a, 
				 const SymbolicScalarMatlab& b) {
    return SymbolicScalarMatlab("("+a.symbol+"+"+b.symbol+")");
  }

  /** 
   * Subtracts two scalars.
   * @param[in] a The first scalar argument.
   * @param[in] b The second scalar argument.
   * @return The difference a-b.
   */
  SymbolicScalarMatlab operator-(const SymbolicScalarMatlab& a, 
				 const SymbolicScalarMatlab& b) {
    return SymbolicScalarMatlab("("+a.symbol+"-"+b.symbol+")");
  }

  /** 
   * Computes the unary minus operation.
   * @param[in] a The first scalar argument.
   * @return The negative value -a.
   */
  SymbolicScalarMatlab operator-(const SymbolicScalarMatlab& a) {
    return SymbolicScalarMatlab("(-"+a.symbol+")");
  }

  /** 
   * Multiplies two scalars.
   * @param[in] a The first scalar argument.
   * @param[in] b The second scalar argument.
   * @return The product a*b.
   */
  SymbolicScalarMatlab operator*(const SymbolicScalarMatlab& a, 
				 const SymbolicScalarMatlab& b) {
    return SymbolicScalarMatlab("("+a.symbol+"*"+b.symbol+")");
  }

  /** 
   * Ration between two scalars.
   * @param[in] a The first scalar argument.
   * @param[in] b The second scalar argument.
   * @return The quotient a/b.
   */
  SymbolicScalarMatlab operator/(const SymbolicScalarMatlab& a, 
				 const SymbolicScalarMatlab& b) {
    return SymbolicScalarMatlab("("+a.symbol+"/"+b.symbol+")");
  }

} /** namespace AMD */

#endif /** AMD_SYMBOLIC_SCALAR_MATLAB_HPP */
