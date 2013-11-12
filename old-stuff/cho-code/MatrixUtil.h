#ifndef MATRIXUTIL
#define MATRIXUTIL

/**
 * @file This file contains a sample implementation of a matrix class for 
 *       single node code. This should be replaced with somthing that is more
 *       heavy duty and possibly parallel (multi-threaded and/or distributed)
 *       
 * @note This is a dense matrix!
 */       
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
#include "SymbolUtil.h"
#include <limits>

#define MALLOC(type, size) (type*)malloc(sizeof(type)*(size))
static const int MAT_EMPTY = 0; /** 0x0 matrix */
static const int MAT_ZERO = 1; /** mxn mtrix with all zeros */
static const int MAT_IDENTITY = 2;
static const int MAT_RANDOM = 3;

extern "C" int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *A, const int *LDA, double *B, const int *LDB, double *beta, double *C, int *LDC); 
extern "C" int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern "C" void dgetri_(int* N, double *A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
extern "C" void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern "C" void dsyevd_(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );

using namespace std;
class mat_t;
class mat_zero;

class mat_t{
	public:
		int m,n;
		double *val;
		symbol_t val_s;
		bool is_zero;
		bool is_identity;
		bool mem_alloc_by_me;

    /**
     * Deep copy construction from one matrix to another 
     */ 
		void copyfrom(const mat_t& X) {
			m = X.m;
			n = X.n;
			if (NULL!=val) free(val);

			if (NULL!=X.val) {
				val = (double*)malloc(sizeof(double)*m*n); // replace with "new"
				for (int i=0 ; i<(m*n); i++) val[i] = X.val[i]; // replace with bcopy
			}

			val_s = X.val_s;
			is_zero = X.is_zero;
			is_identity = X.is_identity;
			mem_alloc_by_me = true;
		}

    /**
     * Constructing a 1x1 matrix 
     */
		mat_t (double x) : m(1), 
                       n(1), 
                       val((double *)malloc(sizeof(double))),
                       mem_alloc_by_me(true) {
      val[0] = x; 
    }

    /**
     * Create an empty matrix --- cannot be used until initilization.
     */ 
		mat_t():mem_alloc_by_me(false){val = NULL; }

    /**
     * Copy constructor.
     */ 
		mat_t(const mat_t& A) { copyfrom(A); }
		
    /**
     * Consutrction of an mxn matrix. Memory is allocated here.
     */ 
		mat_t(int m_, int n_, int type) {
			m = m_;
			n = n_;
			is_zero = false;
			is_identity = false;
			if ( type == MAT_ZERO) {
				is_zero = true;
				val_s.val = "0";
				val = NULL;
			} else if ( type == MAT_IDENTITY) {
				if ( m != n) {
					printf("ERROR!!! Create a %d by %d matrix to be identity.\n", m, n);
					exit(1);
				}
				val = (double *)malloc(sizeof(double)*m*n);
				for ( int i=0 ; i<n*m ; i++ )
					val[i] = 0;
				for ( int i=0 ; i<n ; i++ )
					val[i*n+i] = 1;
				is_identity = true;
				val_s.val = "I";
			} else {
				val = (double *)malloc(sizeof(double)*m*n);
				if ( type == MAT_EMPTY)
				{
					for ( int i=0 ; i<m*n ; i++ )
						val[i] = 0;
					val_s.val = "E";
				}
				if ( type == MAT_RANDOM )
				{
					for ( int i=0 ; i<m*n ; i++ )
						val[i] = ((double)rand() / RAND_MAX) - 0.5;
					val_s.val = "R";
				}
			}
		}

		string getname() { return val_s.val; }

		void setname(string s) { val_s.val = s; }

		void print(FILE *fp, bool print_symbol=true) const {
			if (true==print_symbol)
        printf("dim: %d %d symbol: %s\n", m, n, val_s.val.c_str());

			if ( is_zero == false && is_identity == false) {
				for(int i = 0; i<m; i++ ) {
					for ( int j=0 ; j<n ; j++ )
						fprintf(fp, "%.14lf ", val[j*m+i]);
					fprintf(fp, "\n");
				}
			}
		}

		void save_ascii(char *filename) {
			FILE *fp = fopen(filename, "w"); 
      print (fp, false /* dont print symbol */);
			fclose(fp);
		}

		void load_ascii(char *filename, string s) {
			FILE *fp = fopen(filename, "r");
			int t = fscanf(fp, "%d", &m);
			t = fscanf(fp, "%d", &n);
			val = (double *)malloc(sizeof(double)*m*n);
			for ( int i=0 ; i<m ; i++ )
				for ( int j=0 ; j<n ; j++ )
					t = fscanf(fp, "%lf", &(val[j*m+i]));
			val_s.val = s;
			is_zero = false;
			is_identity = false;
		}

		double det() {
			if ( m!=n ) {
				printf("ERROR: not a square matrix\n");
				return -1;
			}

			if ( is_zero == true) return 0;
			if ( is_identity ==  true) return 1;

  		char uplo = 'u';
			int info;
			double *tmp = (double *)malloc(sizeof(double)*m*n);
			memcpy(tmp, val, sizeof(double)*m*n);
			dpotrf_(&uplo, &m, tmp, &n, &info);
			double result = 1;
			for (int i=0 ; i<n ; i++ ) result *= tmp[i*n+i];
			return result*result;
		}

		void free(void *ptr) {if(!ptr) ::free(ptr);}
		
    /** Destructor */
		~mat_t() { if(NULL != val) free(val); }

		mat_t& operator=(const mat_t &rhs){ 
			copyfrom(rhs);
			return (*this);
		}

		mat_t& operator-=(const mat_t &other) {
			if ( other.m != n || other.n != n) {
				printf("ERROR!! matrix minus size mismatch. (%d x %d)*(%d x %d)\n", m, n, other.m, other.n);
				exit(1);
			}

			if ( other.is_zero == true) { }
			else if ( is_zero == true ) {
				copyfrom(other);
				for ( int i=0 ; i<m*n ; i++ )
					val[i] = -val[i];
				val_s.val = "-"+other.val_s.val;
			}
			else {
				for ( int i=0 ; i<m*n ; i++ )
					val[i] -= other.val[i];
				val_s -= other.val_s;
				is_identity = false;
			}
			return (*this);
		}

		mat_t& operator+=(const mat_t &other)  
		{
			if ( other.m != n || other.n != n)
			{
				printf("ERROR!! matrix plus size mismatch. (%d x %d)*(%d x %d)\n", m, n, other.m, other.n);
				exit(1);
			}

			if ( other.is_zero == true)
			{

			}
			else if ( is_zero == true )
			{
				copyfrom(other);
			}
			else
			{
				for ( int i=0 ; i<m*n ; i++ )
					val[i] += other.val[i];
				val_s += other.val_s;
				is_identity = false;
			}

			return (*this);
		}

		mat_t& operator*=(const mat_t &other) 
		{
			if ( other.m != n)
			{
				printf("ERROR!! matrix multiplication size mismatch. (%d x %d)*(%d x %d)\n", m, n, other.m, other.n);
				exit(1);
			}
			if ( is_zero == true || other.is_zero == true)
			{
				is_zero = true;
				is_identity = false;
				if ( val)
					free(val);
				val_s.val = "0";
			}
			else if ( other.is_identity == true )
			{
			}
			else if ( is_identity == true)
			{
				copyfrom(other);
			}
			else
			{
				double *original_val = val;
				val = (double *)malloc(sizeof(double)*m*(other.n));
				char trans = 'n';
				double alpha = 1.0;
				int LDA = 2;
				int LDB = 2;
				double beta = 0.0;
				int LDC = 2;
				int m_ = m;
				int n_ = other.n;
				int k_ = n;
		
				dgemm_(&trans,&trans,&m_,&n_,&k_,&alpha, original_val, &m_, other.val, &k_, &beta, val, &m_);
				n = other.n;
				
				val_s *= other.val_s;
				is_identity = false;
				
				free(original_val);
			}
			return (*this);
		}

		mat_t& operator/=(double other) {
			if ( is_zero == true)
				return (*this);
			const register double size=m*n;
			for (int i=0 ; i<size ; i++)
				val[i]/=other;
			val_s.val = "("+val_s.val+"/const)";
			is_identity = false;
			return (*this);
		}
};

mat_t operator-(const mat_t &lhs, const mat_t &rhs) 
{
	mat_t result(lhs);
	result -= rhs;
	return result;
}

mat_t operator+(const mat_t &lhs, const mat_t &rhs) 
{
	mat_t result(lhs);
	result += rhs;
	return result;
}

mat_t operator*(const mat_t &lhs, const mat_t &rhs) 
{
	mat_t result(lhs);
	result *= rhs;
	return result;
}

mat_t operator/(const mat_t &lhs, const double &rhs) 
{
	mat_t result(lhs);
	result /= rhs;
	return result;
}

mat_t operator*(const double lhs, const mat_t &rhs) {
	if (rhs.is_zero == true)
		return mat_t(rhs.m, rhs.n, MAT_ZERO);
	mat_t result(rhs);
	int total = result.m*result.n;
	for ( int i=0 ; i<total ; i++ )
		result.val[i] *= lhs;
	char tmp[100];
	sprintf(tmp, "(%lf*", lhs);
	string str(tmp);
	result.val_s.val = str+result.val_s.val+")";
	result.is_identity = false;
	return result;
}

mat_t operator*(const mat_t &lhs, const double rhs) 
{
	if ( lhs.is_zero == true)
		return mat_t(lhs.m, lhs.n, MAT_ZERO);
	mat_t result(lhs);
	int total = result.m*result.n;
	for ( int i=0 ; i<total ; i++ )
		result.val[i] *= rhs;
	result.val_s.val = "(const*"+result.val_s.val+")";
	result.is_identity = false;
	return result;
}


mat_t transpose(mat_t A)
{
	if ( A.is_zero == true)
	{
		return mat_t(A.n, A.m, MAT_ZERO);
	}
	else if ( A.is_identity == true)
	{
		return mat_t(A.n, A.m, MAT_IDENTITY);
	}
	
	mat_t result(A.n, A.m, MAT_EMPTY);
	for (int i=0 ; i< result.m ; i++ )
		for ( int j=0 ; j<result.n ; j++ )
			result.val[j*result.m+i] = A.val[i*A.m+j]; 
	result.val_s.val = A.val_s.val + "'";
	return result;
}


mat_t inv(const mat_t &A)
{
	if ( A.is_zero == true)
	{
		printf("ERROR!!! inverse of a 0 matrix\n");
		exit(1);
	}
	if ( A.is_identity == true)
	{
		return mat_t(A.m, A.n, MAT_IDENTITY);
	}
	int *ipiv = (int *)malloc(sizeof(int)*(A.m+1));
	mat_t result = A;
	int lwork = A.m*A.m;
	double *work = (double *)malloc(sizeof(double)*lwork);
	int info;

	dgetrf_(&result.m,&result.m,result.val,&result.m,ipiv,&info);

	dgetri_(&result.m,result.val, &result.m, ipiv, work, &lwork, &info );
	free(work);
	result.val_s.val = "inv("+A.val_s.val+")";
	return result;
}

mat_t reshape(mat_t &X, int  m_, int n_)
{
	if (X.m*X.n != m_*n_)
	{
		printf("ERROR!! reshape size error. From (%d,%d) to (%d,%d)", X.m, X.n, m_, n_);
		exit(1);
	}
	if ( X.is_zero == true)
	{
		return mat_t(m_, n_, MAT_ZERO);
	}
	mat_t result(m_,n_, MAT_EMPTY);
	memcpy(result.val, X.val, sizeof(double)*m_*n_);
	result.val_s = X.val_s;
	return result;
}

double logdet(mat_t X)
{
	int m = X.m, n=X.n;
	if ( X.m!=X.n )
	{
		printf("ERROR: not a square matrix\n");
		return -1;
	}
//	if ( X.is_zero == true)
//		return 0;
	if ( X.is_identity == true)
		return 0;
	char uplo = 'u';
	int info;
	double *tmp = (double *)malloc(sizeof(double)*m*n);
	memcpy(tmp, X.val, sizeof(double)*m*n);
	dpotrf_(&uplo, &m, tmp, &n, &info);
	double result = 0;
	if ( info> 0)
	{
		printf("Matrix Not Positive Definite!\n");
		return -1.0/0.0; // -1.0*std::numeric_limits<double>::infinity();
	}
	for (int i=0 ; i<n ; i++ )
		result += 2*log(tmp[i*n+i]);
	free(tmp);
	return result;
//	return mat_t(result);
}

double trace(mat_t X)
{
	if ( X.is_zero )
		return 0;
	if ( X.is_identity)
		return X.n;
	double result = 0;
	int d = X.m;
	for ( int i=0, idx = 0 ; i<d ; i++, idx+=(d+1))
		result += X.val[idx];
	return result;
//	return mat_t(result);
}

mat_t form_vecI(int d)
{
	mat_t result(d*d,1,MAT_EMPTY);
	for ( int i=0, ind=0 ; i<d ; i++, ind+=(d+1))
		result.val[ind] = 1;
	result.val_s.val = "vec(I)";
	result.is_zero = false;
	return result;
}

mat_t A_times_spvec(mat_t X, mat_t v)
{
	if ( (v.m != X.n) || (v.n != 1) )
	{
		printf("ERROR! A_times_spvec, A (%d,%d), x(%d,%d)", X.m, X.n, v.m, v.n);
		exit(1);
	}

	if ( X.is_zero == true || v.is_zero == true)
	{
		return mat_t(X.m, 1, MAT_ZERO);
	}
	if ( X.is_identity == true)
	{
		return mat_t(v);
	}
	mat_t result(X.m,1, MAT_EMPTY);
	for ( int i=0 ; i<v.m ; i++ )
		if ( v.val[i] != 0)
		{
			double tmpv = v.val[i];
			for ( int j=0, idx=i*(X.m) ; j<X.m ; j++, idx++)
				result.val[j] += tmpv*(X.val[idx]);
		}
	result.val_s = X.val_s*v.val_s;
	return result;
}

mat_t kron(mat_t A, mat_t B)
{
	if ( A.is_zero || B.is_zero )
	{
		return mat_t(A.m*B.m, A.n*B.n, MAT_ZERO);
	}
	if ( A.is_identity && B.is_identity )
	{
		return mat_t(A.m*B.m, A.n*B.n, MAT_IDENTITY);
	}
	mat_t result(A.m*B.m, A.n*B.n, MAT_EMPTY);
	for ( int i=0 ; i<A.m ; i++ )
		for (int j=0 ; j<A.n ; j++ )
			for ( int p=0 ; p<B.m ; p++ )
				for ( int q=0 ; q<B.n ; q++ )
				{
					int xidx = i*B.m+p;
					int yidx = j*B.n+q;
					result.val[yidx*A.m*B.m+xidx] = A.val[j*A.m+i]*B.val[q*B.m+p];
				}
	result.val_s.val = "kron("+A.val_s.val+","+B.val_s.val+")";
	return result;
}

mat_t box(mat_t A, mat_t B)
{
	if ( A.is_zero || B.is_zero )
	{
		return mat_t(A.m*B.m, A.n*B.n, MAT_ZERO);
	}
/*	if ( A.is_identity && B.is_identity )
	{
		return mat_t(A.m*B.m, A.n*B.n, MAT_IDENTITY);
	}*/
	mat_t result(A.m*B.m, A.n*B.n, MAT_EMPTY);
	for ( int i=0 ; i<A.m ; i++ )
		for (int j=0 ; j<A.n ; j++ )
			for ( int p=0 ; p<B.m ; p++ )
				for ( int q=0 ; q<B.n ; q++ )
				{
					int xidx = i*B.m+p;
					int yidx = j*B.n+q;
					result.val[yidx*A.m*B.m+xidx] = A.val[q*A.m+i]*B.val[j*B.m+p];
				}

	result.val_s.val = "box("+A.val_s.val+","+B.val_s.val+")";
	return result;
}

mat_t eye(int d)
{
	mat_t result(d, d, MAT_IDENTITY);
	return result;
}

double fnorm(mat_t &X)
{
	if ( X.is_zero )
		return 0;
	if ( X.is_identity )
		return 1;
	double result = 0;
	for ( int i=0 ; i<(X.m*X.n) ; i++ )
		result += (X.val[i])*(X.val[i]);
	return sqrt(result);
}

mat_t schur(mat_t A, mat_t B)
{
	if ( A.m != B.m || A.n != B.n )
	{
		printf("ERROR!! Matrix schur product A(%d,%d) dot B(%d,%d)\n", A.m, A.n, B.m, B.n);
		exit(1);
	}
	if ( A.is_zero || B.is_zero )
		return mat_t(A.m, A.n, MAT_ZERO);
	if ( A.is_identity && B.is_identity )
		return mat_t(A.m, A.n, MAT_IDENTITY);
	mat_t result(A);
	for ( int i=0 ; i<A.m*A.n ; i++ )
		result.val[i] *= B.val[i];
	result.val_s.val = "("+A.val_s.val + ".*"+B.val_s.val+")";
	return result;
}

mat_t zeros(int m_, int n_)
{
	mat_t result(m_, n_, MAT_ZERO);

//	result.val_s.val = "0";
	return result;
}

mat_t diag(mat_t &w)
{
	if (w.m != 1 && w.n!=1)
	{
		printf("ERROR! input of diag() is a %d by %d matrix.\n", w.m, w.n);
		exit(1);
	}
	int p = w.m;
	if ( w.n > p) p=w.n;
	mat_t result(p,p,MAT_EMPTY);
	for ( int i=0 ; i<p ; i++ )
		result.val[i*p+i] = w.val[i];
	return result;
}

void eig(mat_t &X, mat_t &U, mat_t &w) {
	if ( X.m != X.n )
	{
		printf("ERROR!! input eig with (%d,%d) matrix.\n", X.m, X.n);
		exit(1);
	}
	U = X;
	int n = X.m, info, lwork=1+6*n+2*n*n, liwork = 3+5*n;
	int* iwork = (int *)malloc(sizeof(int)*liwork);
	double* work = (double *)malloc(sizeof(double)*lwork);
	char *jobz="V";
	char *uplo = "U";

	mat_t tmp(n,1,MAT_EMPTY);
	dsyevd_(jobz, uplo, &n, U.val,&n,tmp.val, work, &lwork, iwork, &liwork, &info);
	if ( info>0 )
	{
		printf("FAIL SOLVING eig.\n");
		exit(1);
	}
	w = tmp;
}

double innerproduct(mat_t &X, mat_t &Y) {
	if ( X.m != Y.m || X.n != Y.n ) {
		printf("ERROR: innerproduct dimension mismatch: X(%d,%d), Y(%d,%d)\n", X.n, X.n, Y.m, Y.n);
		exit(1);
	}
	if ( X.is_zero || Y.is_zero) return 0;
	if ( X.is_identity ) return trace(Y);
	if ( Y.is_identity ) return trace(X);
	double result = 0;
	for ( int i=0 ; i<(X.m*X.n) ; i++ )
		result += (X.val[i])*(Y.val[i]);
	return result;
}

#endif
