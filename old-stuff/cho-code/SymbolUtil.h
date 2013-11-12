#ifndef SYMBOLUTIL
#define SYMBOLUTIL
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <vector>
#include <cmath>
#include <omp.h>
#include <assert.h>


using namespace std;
class symbol_t;

class symbol_t{
	public:
		string val;

		symbol_t() {}
		symbol_t(string x_) {val = x_; }
		
		~symbol_t(){}
	
		void print(FILE *fp) const {
			fprintf(fp, "%s\n", val.c_str());
		}

		symbol_t& operator=(const symbol_t &rhs){ 
			val = rhs.val;
			return (*this);
		}

		symbol_t& operator-=(const symbol_t &other)  
		{
			val = "("+val+"-"+other.val+")";
			return (*this);
		}

		symbol_t& operator+=(const symbol_t &other)  
		{
			val = "("+val+"+"+other.val+")";
			return (*this);
		}

		symbol_t& operator*=(const symbol_t &other) 
		{
			val = "("+val+"*"+other.val+")";
			return (*this);
		}
};

symbol_t operator-(const symbol_t &lhs, const symbol_t &rhs) 
{
	symbol_t result("("+lhs.val+"-"+rhs.val+")");
	return result;
}

symbol_t operator+(const symbol_t &lhs, const symbol_t &rhs) 
{
	symbol_t result("("+lhs.val+"+"+rhs.val+")");
	return result;
}

symbol_t operator*(const symbol_t &lhs, const symbol_t &rhs) 
{
	symbol_t result("("+lhs.val+"*"+rhs.val+")");
	return result;
}

symbol_t operator/(const symbol_t &lhs, const symbol_t &rhs) 
{
	symbol_t result("("+lhs.val+"/"+rhs.val+")");
	return result;
}

symbol_t transpose(symbol_t A)
{
	symbol_t result(A.val+"'");
	return result;
}


symbol_t inv(const symbol_t &A)
{
	symbol_t result("inv("+A.val+")");
	return result;
}
  
symbol_t logdet(symbol_t X)
{
	symbol_t result("logdet("+X.val+")");
	return result;
}

symbol_t trace(symbol_t X)
{
	symbol_t result("trace("+X.val+")");
	return result;
}

symbol_t kron(symbol_t A, symbol_t B)
{
	symbol_t result("kron("+A.val+","+B.val+")");
	return result;
}

symbol_t box(symbol_t A, symbol_t B)
{
	symbol_t result("box("+A.val+","+B.val+")");
	return result;
}

symbol_t eye()
{
	symbol_t result("I");
	return result;
}

symbol_t zeros()
{
	symbol_t result("0");
	return result;
}

symbol_t fnorm(symbol_t &X)
{
	symbol_t result("||"+X.val+"||_F");
	return result;
}

symbol_t schur(symbol_t A, symbol_t B)
{
	symbol_t result("(A.*B)");
	return result;
}

symbol_t det(symbol_t x)
{
	return symbol_t("det("+x.val+")");
}

#endif
