#ifndef AMD_UTILITY_HPP
#define AMD_UTILITY_HPP

#include <string>

namespace AMD { namespace detail {

template <typename T> T* pointerCopy(T x) {
  T* tmp = new T;
  (*tmp) = x;
  return(tmp);
}

std::string removeParenthesis(std::string str) {
  std::string retVal = str;

  if ('(' == str[0] && 
      ')' == str[str.length()-1]) {
    retVal=str.substr(1,str.length()-2);
  }
  
  return retVal;
}

} } /** namespace AMD::detail */

#endif /** AMD_UTILITY_HPP */
