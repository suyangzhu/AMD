#ifndef RDiff_H
#define RDiff_H

/**
 * A placeholder for both binary and unary matrix-matrix operations.
 * Representing something as f(X) implies that the function returns a scalar
 * Representing something as F(X) implies that the function returns a matrix
 */ 

#include <cstdio>
#include <cstdlib>
#include "MatrixUtil.h"
#include "MatrixScalarDiff.h"
#define OP_LEAF 0 /** The value is the operand itself */
#define OP_CONST 1  /** Derivative is zero */
#define OP_MULT 2 /** F(X)*G(X) */
#define OP_TRANSPOSE 3 /** F(X)' */
#define OP_ADD 4 /** F(X)+G(X) */
#define OP_INV 5 /** F(X)^{-1} */
#define OP_SCHUR 6 /** F(X).*G(X) Schur or Hadamard or Elementwise product */
#define OP_AgX 7 /** A*g(X) diff wrt X */
#define OP_FXgX 8 /** F(X)*g(X) */
#define OP_MULT_CONSTSCALAR 9 /** a*F(X) */
#define OP_CONSTSCALAR 10 /** a */
#define ONLY_SCALAR 11 /**  */
#define OP_MINUS 12 /** F(X)-G(X) */

using namespace std;

class RDiff {
	public:
		int op;  /** One of the operators defined above */
		mat_t value; /** F(X) */
		mat_t diff; /** A matrix needed to compute F'(X) */
		mat_t constMatrix; // For A*g(X)
		double constscalar; // For a*g(X)
		double gX_value;
		const RDiff *lchild;
		const RDiff *rchild;
		int in_m;
		int in_n;
		int get_symbolic;
		
		RDiff() {lchild = NULL; rchild=NULL; get_symbolic = 0;}

		RDiff(mat_t &X, int flag) {
			init(X, flag);
			get_symbolic = 0;
		}

		RDiff(double scalar) {
			constscalar = scalar;
			op = ONLY_SCALAR;
		}

		RDiff(mat_t &X) {
			init(X, OP_LEAF);
			get_symbolic = 0;
		}
		
		RDiff(const mat_t &value_, const mat_t &diff_, const mat_t &constMatrix_, int flag)
		{
			value = value_;
			lchild = NULL;
			rchild = NULL;
			in_m = diff_.m;
			in_n = diff_.n;
			op = flag;
			diff = diff_;
			constMatrix = constMatrix_;
			get_symbolic = 0;
		}

		void init(mat_t &X, int flag)
		{
			value = X;
			lchild = NULL;
			rchild = NULL;
			in_m = X.m;
			in_n = X.n;
			op = flag;
		}

		~RDiff() {
			if (!lchild)
				delete lchild;
			if (!rchild)
				delete rchild;
		}

		RDiff(const RDiff *A, const RDiff *B, int op_type)
		{
			if ( (A->in_m != B->in_m) || (A->in_n != B->in_n) )
			{
				printf("ERROR! RDIFF type mismatch: (%d,%d) and (%d,%d)\n", A->in_m, A->in_n, B->in_m, B->in_n);
				exit(1);
			}
			if ( op_type == OP_MULT)
			{
				value = (A->value)*(B->value);
				lchild = A;
				rchild = B;
			}
			op = op_type;
			in_m = A->in_m;
			in_n = A->in_n;
		}


		void print(FILE *fp) const {
			if ( op == OP_LEAF )
				value.print(fp);
			else
			{
				if ( op == OP_MULT )
					fprintf(fp,"operator: *\n");
				else if ( op == OP_TRANSPOSE )
					fprintf(fp,"operator: transpose\n");
				printf("value: ");
				value.print(fp);
				fprintf(fp,"{\n");
				if ( lchild != NULL)
					lchild->print(fp);
				if (rchild != NULL)
					rchild->print(fp);
				fprintf(fp,"}\n");
			}
		}

		RDiff& operator= ( const RDiff &x) {
			value = x.value;
			op = x.op;
			lchild = x.lchild;
			rchild = x.rchild;
			in_m = x.in_m; 
			in_n = x.in_n;
			constscalar = x.constscalar;
			return(*this);
		}
};


RDiff operator* (const RDiff &lhs, const RDiff &rhs)
{
	RDiff result(&lhs, &rhs, OP_MULT);
	return result;
}

RDiff operator* (const double lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs*rhs.value;
	result.lchild = &(rhs);
	result.rchild = NULL;
	result.op = OP_MULT_CONSTSCALAR;
	result.constscalar = lhs;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}

RDiff operator* (const RDiff &rhs, const double lhs)
{
	RDiff result;
	result.value = lhs*rhs.value;
	result.lchild = &(rhs);
	result.rchild = NULL;
	result.op = OP_MULT_CONSTSCALAR;
	result.constscalar = lhs;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}

RDiff operator/ (const RDiff &rhs, const double lhs )
{
	RDiff result;
	result.value = (1/lhs)*rhs.value;
	result.lchild = &(rhs);
	result.rchild = NULL;
	result.op = OP_MULT_CONSTSCALAR;
	result.constscalar = 1/lhs;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}

RDiff operator* (const mat_t &lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs*(rhs.value);
	result.lchild = NULL;
	result.rchild = &(rhs);
	result.op = OP_MULT;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}

RDiff operator* (const RDiff &lhs, const mat_t &rhs)
{
	RDiff result;
	result.value = (lhs.value)*rhs;
	result.lchild = &(lhs);
	result.rchild = NULL;
	result.op = OP_MULT;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	return result;
}


RDiff operator+ (const RDiff &lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs.value + rhs.value;
	result.lchild = &(lhs);
	result.rchild = &(rhs);
	result.op = OP_ADD;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	return result;
}


RDiff operator+ (const RDiff &lhs, const mat_t &rhs)
{
	RDiff result;
	result.value = lhs.value + rhs;
	result.lchild = &(lhs);
	result.rchild = NULL;
	result.op = OP_ADD;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	return result;
}

RDiff operator+ (const mat_t &lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs + rhs.value;
	result.lchild = NULL;
	result.rchild = &(rhs);
	result.op = OP_ADD;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}

RDiff operator- (const RDiff &lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs.value - rhs.value;
	result.lchild = &(lhs);
	result.rchild = &(rhs);
	result.op = OP_MINUS;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	return result;
}


RDiff operator- (const RDiff &lhs, const mat_t &rhs)
{
	RDiff result;
	result.value = lhs.value - rhs;
	result.lchild = &(lhs);
	result.rchild = NULL;
	result.op = OP_MINUS;
	result.in_m = lhs.in_m;
	result.in_n = lhs.in_n;
	return result;
}

RDiff operator- (const mat_t &lhs, const RDiff &rhs)
{
	RDiff result;
	result.value = lhs - rhs.value;
	result.lchild = NULL;
	result.rchild = &(rhs);
	result.op = OP_MINUS;
	result.in_m = rhs.in_m;
	result.in_n = rhs.in_n;
	return result;
}


RDiff transpose(const RDiff &X)
{
	RDiff result;
	result.op = OP_TRANSPOSE; 
	result.value = transpose(X.value);
	result.lchild = &(X);
	result.rchild = NULL;
	result.in_m = X.in_m;
	result.in_n = X.in_n;
	return result;
}

RDiff inv(const RDiff &X)
{
	RDiff result;
	result.op = OP_INV; 
	result.value = inv(X.value);
	result.lchild = &(X);
	result.rchild = NULL;
	result.in_m = X.in_m;
	result.in_n = X.in_n;
	return result;
}


/**
 * All RDiff's are held symbolically without being evaluated as the RDiff are
 * *huge* tensors. The only operation performed on these are mat-vec's in the
 * form of CompGradient. These mat-vecs are triggered when a scalr function is
 * encountered.
 * vec(result) = F'(X)*vec(initial)
 * result == matrix.
 */ 
void CompGradient(const RDiff *root, 
                  mat_t *initial, 
                  mat_t *result, 
                  int transpose_flag) {
	vector<const RDiff *> RDiff_stack(100);
	vector<mat_t *> R_stack(100);

	int stack_size = 1;
	RDiff_stack[0] = root;
	R_stack[0] = initial;
	while (stack_size > 0 ) {
		const RDiff *current_rdiff = RDiff_stack[stack_size-1];
		mat_t *current_R = R_stack[stack_size-1];

		stack_size--;
		if ( current_rdiff->op == OP_CONST)
		{
		}
		else if ( current_rdiff->op == OP_SCHUR )
		{
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = schur(transpose(current_rdiff->rchild->value),(*current_R));
			RDiff_stack[stack_size] = current_rdiff->lchild;
			stack_size++;
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = schur((*current_R),transpose(current_rdiff->lchild->value));
			RDiff_stack[stack_size] = current_rdiff->rchild;

			stack_size++;

		}
		else if ( current_rdiff->op == OP_MULT )
		{
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = (current_rdiff->rchild->value)*(*current_R);
			RDiff_stack[stack_size] = current_rdiff->lchild;
			stack_size++;
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = (*current_R)*(current_rdiff->lchild->value);
			RDiff_stack[stack_size] = current_rdiff->rchild;

			stack_size++;
		}
		else if ( current_rdiff->op == OP_MULT_CONSTSCALAR )
		{
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = current_rdiff->constscalar*(*current_R);
			RDiff_stack[stack_size] = current_rdiff->lchild;
			stack_size++;
		}
		else if (current_rdiff->op == OP_TRANSPOSE )
		{
			// F(X) = G(X)^T
			// \nabla F(X) = I_n box I_m
			// Compute (I_n box I_m)^T vec(A) = (I_m box I_n) vec(A) 
			// = vec(A^T)
//			printf("Deal with transpose\n");
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = transpose(*current_R);
			RDiff_stack[stack_size] = current_rdiff->lchild;
//			printf("A:\n");
//			R_stack[stack_size]->print(stdout);
			stack_size++;
		}
		else if (current_rdiff->op == OP_INV )
		{
			// F(X) = inv(G(X))
			// Compute -G(X)^(-1) A G(X)^(-1)
//			printf("Deal with inverse\n");
			R_stack[stack_size] = new mat_t;
			(*R_stack[stack_size]) = (-1)*(current_rdiff->value)*(*current_R)*(current_rdiff->value); // is there faster way?
			RDiff_stack[stack_size] = current_rdiff->lchild;
//			printf("A:\n");
//			R_stack[stack_size]->print(stdout);
			stack_size++;
		}
		else if (current_rdiff->op == OP_ADD )
		{
//			printf("Deal with Add\n");
			if (current_rdiff->lchild != NULL )
			{
				R_stack[stack_size] = new mat_t;
				(*R_stack[stack_size]) = (*current_R);
				RDiff_stack[stack_size] = current_rdiff->lchild;
				stack_size++;
			}
			if ( current_rdiff->rchild != NULL)
			{
				R_stack[stack_size] = new mat_t;
				(*R_stack[stack_size]) = (*current_R);
				RDiff_stack[stack_size] = current_rdiff->rchild;
				stack_size++;
			}
		}
		else if (current_rdiff->op == OP_MINUS )
		{
//			printf("Deal with Add\n");
			if (current_rdiff->lchild != NULL )
			{
				R_stack[stack_size] = new mat_t;
				(*R_stack[stack_size]) = (*current_R);
				RDiff_stack[stack_size] = current_rdiff->lchild;
				stack_size++;
			}
			if ( current_rdiff->rchild != NULL)
			{
				R_stack[stack_size] = new mat_t;
				(*R_stack[stack_size]) = (-1)*(*current_R);
				RDiff_stack[stack_size] = current_rdiff->rchild;
				stack_size++;
			}
		}

		else if (current_rdiff->op == OP_LEAF )
		{
//			mat_t AA = transpose(*current_R);
			if ( transpose_flag )
				(*result) = (*result) + transpose(*current_R);
			else
				(*result) = (*result) + (*current_R);
		}
		else if (current_rdiff->op == OP_AgX) 
		{
//			printf("Deal with Ag(X)\n");
			if ( transpose_flag)
				(*result) = (*result) + schur(current_rdiff->constMatrix, (*current_R))*current_rdiff->diff;
			else
				(*result) = (*result) + transpose(schur(current_rdiff->constMatrix, (*current_R))*current_rdiff->diff);
		}
		else if (current_rdiff->op == OP_FXgX )
		{
//			printf("Deal with F(X)g(X)\n");
			if ( transpose_flag )
				(*result) = (*result) + schur(current_rdiff->constMatrix, (*current_R))*current_rdiff->diff;
			else
				(*result) = (*result) + transpose(schur(current_rdiff->constMatrix, (*current_R))*current_rdiff->diff);
			if ( current_rdiff->lchild != NULL )
			{
				R_stack[stack_size] = new mat_t;
				(*R_stack[stack_size]) = (*current_R)*(current_rdiff->gX_value); // is there faster way?
				RDiff_stack[stack_size] = current_rdiff->lchild;
//				printf("A:\n");
				R_stack[stack_size]->print(stdout);
				stack_size++;
			}
		}
	}

}


MatrixScalarDiff trace(const RDiff &X)
{
	if ( X.value.m != X.value.n )
	{
		printf("ERROR! trace of a rectangular matrix (%d,%d)\n", X.value.n, X.value.n);
		exit(1);
	}
	MatrixScalarDiff result(trace(X.value), zeros(X.in_m, X.in_n));
//	mat_t result = zeros(X.in_m, X.in_n);
	mat_t tmp = eye(X.value.m);
	CompGradient(&X, &tmp, &(result.diff), 1); // modify to eye(X.value.m)
	printf("FINISH TRACE\n");
	return result;
}

MatrixScalarDiff logdet(const RDiff &X)
{
	if ( X.value.m != X.value.n )
	{
		printf("ERROR! logdet of a rectangular matrix (%d,%d)\n", X.value.n, X.value.n);
		exit(1);
	}

	MatrixScalarDiff result(logdet(X.value), zeros(X.in_m, X.in_n));
//	mat_t result = zeros(X.in_m, X.in_n);
	mat_t initial = inv(X.value);

	CompGradient(&X, &initial, &(result.diff), 1);
	return result;
}

RDiff schur(const RDiff &G, const RDiff &H)
{
	RDiff result;
	result.op = OP_SCHUR; 
	result.value = schur(G.value, H.value);
	result.lchild = &(G);
	result.rchild = &(H);
	result.in_m = G.in_m;
	result.in_n = G.in_n;
	return result;
}

RDiff operator* (const mat_t &lhs, const MatrixScalarDiff &rhs)
{
	RDiff result(lhs*(rhs.value), rhs.diff, lhs, OP_AgX);
	return result;
}

RDiff operator* (const RDiff &lhs, const MatrixScalarDiff &rhs)
{
	RDiff result(lhs.value*(rhs.value), rhs.diff, lhs.value, OP_FXgX);
	result.lchild = (&lhs) ;
	result.gX_value = rhs.value;
	return result;
}

int do_one_operation(vector<int> &operator_stack, int &ostack_len, vector<RDiff *> &rdiff_stack, int &rstack_len)
{
		if ( ostack_len == 0 )
			return 0;
		if ( operator_stack[ostack_len-1] == '*' )
		{
			RDiff *lhs = rdiff_stack[rstack_len-2];
			RDiff *rhs = rdiff_stack[rstack_len-1];
			rstack_len-=2;
			if ( lhs->op == ONLY_SCALAR && rhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff(lhs->constscalar*rhs->constscalar);
				rstack_len++;
			}
			else if ( lhs->op == ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				printf("NMAKE SURE HERE!!!!!!!!!!!!!!!!!!!!\n");

				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = lhs->constscalar*(*rhs);
				rstack_len++;
			}
			else if ( lhs->op != ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = rhs->constscalar*(*lhs);
				rstack_len++;
			}
			else
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = (*lhs)*(*rhs);
				rstack_len++;
			}
		}
		else if ( operator_stack[ostack_len-1] == '/' )
		{
			RDiff *lhs = rdiff_stack[rstack_len-2];
			RDiff *rhs = rdiff_stack[rstack_len-1];
			rstack_len-=2;
			if ( rhs->op == ONLY_SCALAR  && lhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff(lhs->constscalar*rhs->constscalar);
				rstack_len++;
			}
			else if ( lhs->op == ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				printf("ERROR in parser: constscalar/mat_t!!!\n");
				exit(1);
			}
			else if ( lhs->op != ONLY_SCALAR && rhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len]= new RDiff;
				(*rdiff_stack[rstack_len]) = (*lhs)*(1/rhs->constscalar);
				rstack_len++;
			}
			else if ( lhs->op != ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				printf("ERROR in parser: mat_t/mat_t\n");
				exit(1);
			}
		}
		else if ( operator_stack[ostack_len-1] == '+' )
		{
			printf("+!\n");
			RDiff *lhs = rdiff_stack[rstack_len-2];
			RDiff *rhs = rdiff_stack[rstack_len-1];
			rstack_len-=2;
			if ( lhs->op == ONLY_SCALAR && rhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff(lhs->constscalar+rhs->constscalar);
				rstack_len++;
			}
			else if ( lhs->op != ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = (*lhs)+(*rhs);
				rstack_len++;
			}
			else
			{
				printf("ERROR in parser: mat_t + scalar.\n");
				exit(1);
			}
		}
		else if ( operator_stack[ostack_len-1] == '-')
		{
			RDiff *lhs = rdiff_stack[rstack_len-2];
			RDiff *rhs = rdiff_stack[rstack_len-1];
			rstack_len-=2;
			if ( lhs->op == ONLY_SCALAR && rhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff(lhs->constscalar-rhs->constscalar);
				rstack_len++;
			}
			else if ( lhs->op != ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = (*lhs) - (*rhs);
				rstack_len++;
			}
			else
			{
				printf("ERROR in parser: mat_t - scalar.\n");
				exit(1);
			}
		}
		else if ( operator_stack[ostack_len-1] == 'i')
		{
			printf("gogo i\n");
			RDiff *lhs = rdiff_stack[rstack_len-1];
			rstack_len-=1;
			if ( lhs->op == ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff(1/lhs->constscalar);
				rstack_len++;
			}
			else
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = inv(*lhs);
				rstack_len++;
			}
		}
		else if ( operator_stack[ostack_len-1] == '.')
		{
			RDiff *lhs = rdiff_stack[rstack_len-2];
			RDiff *rhs = rdiff_stack[rstack_len-1];
			rstack_len-=2;
			if ( lhs->op != ONLY_SCALAR && rhs->op != ONLY_SCALAR)
			{
				rdiff_stack[rstack_len] = new RDiff;
				(*rdiff_stack[rstack_len]) = schur(*lhs, *rhs);
				rstack_len++;
			}
		}
		else
			return 0;
		ostack_len--;
	return 1;
}


int do_operations(vector<int> &operator_stack, int &ostack_len, vector<RDiff *> &rdiff_stack, int &rstack_len)
{
	while (1)
	{
		int noper = do_one_operation(operator_stack, ostack_len, rdiff_stack, rstack_len);
		if ( noper == 0)
			break;
	}
	return 1;
}


RDiff *parsestring(string str, mat_t &var, vector<mat_t> &const_list)
{
	int len = str.length();
	vector<int> operator_stack(100);
	vector<RDiff *> rdiff_stack(100);
	int ostack_len = 0, rstack_len=0;
	for ( int i=0 ; i<len ; i++ )
	{
//		printf("now deal with %c\n", str[i]);
		if ( str[i]=='(' || str[i] == '+' || str[i]== '-' || str[i]=='*' || str[i]=='/')
			operator_stack[ostack_len++] = str[i];
		else if ( str[i] == '\'') // deal with transpose
		{
			RDiff *now = rdiff_stack[rstack_len-1];
			int count = 1;
			for (i++; i<len ; i++, count++ )
				if ( str[i] != '\'')
					break;

			if ( count%2==1)
			{
			if ( now->op != ONLY_SCALAR)
			{
				rdiff_stack[rstack_len-1] = new RDiff;
				(*rdiff_stack[rstack_len-1]) = transpose(*now);
			}
			}
			i--;
		}
		else if ( str[i] == ')')
		{
			do_operations(operator_stack, ostack_len, rdiff_stack, rstack_len);
			if ( operator_stack[ostack_len-1] != '(')
			{
				printf("brackets misatch %c!!\n", operator_stack[ostack_len-1]);
				exit(1);
			}
			ostack_len--;
			if ( operator_stack[ostack_len-1] == 'i')
			{
				do_one_operation(operator_stack, ostack_len, rdiff_stack, rstack_len);
			}
		}
		else if ( str[i]-'0' >= 0  && str[i]-'0'<10 ) // Deal with numbers
		{
			int begin = i;
			for (; (str[i]-'0'>=0 && str[i]-'0'<10) || str[i] == '.'; i++);
			string tmpnum(str, begin, i-begin);
			double num = atof(tmpnum.c_str());
			i--;
			rdiff_stack[rstack_len] = new RDiff(num);
			rstack_len++;
			do_operations(operator_stack, ostack_len, rdiff_stack, rstack_len);
		}
		else if ( str[i] == '.' && str[i+1]=='*')
		{
			operator_stack[ostack_len++] = '.';
			i++;
		}
		else if ( str[i]-'a'>=0 && str[i]-'a'<26) // DEAL with operation
		{
			string s = "";
			for (; str[i]-'a'>=0 && str[i]-'a'<26 ; i++ )
				s += str[i];
			i--;
			if ( s == "inv" )
				operator_stack[ostack_len++] = 'i';
			else
			{
				printf("ERROR in parser: can't find %s\n", s.c_str());
				exit(1);
			}
		}
		else if ( str[i]-'A'>=0 && str[i]-'A' < 26 ) // DEAL with matrix
		{
			string tmp = var.getname();
			if (str[i] == tmp[0])
			{
				rdiff_stack[rstack_len] = new RDiff(var);
				rstack_len++;
			}
			else
			{
				int id;
				for ( id=0 ; id<const_list.size() ; id++ )
				{
					string tmp = const_list[id].getname();
					if ( tmp[0] == str[i] )
						break;
				}
				if ( id>= const_list.size() )
				{
					printf("ERROR in parser: CANNOT find symbol %c in the list.\n",str[i] );
					exit(1);
				}
				rdiff_stack[rstack_len] = new RDiff(const_list[id], OP_CONST );
				rstack_len++;
			}
			do_operations(operator_stack, ostack_len, rdiff_stack, rstack_len);
		}
	}

	do_operations(operator_stack, ostack_len, rdiff_stack, rstack_len);
	if ( rstack_len >1)
	{
		printf("ERROR in parser: rdiff stack has size > 1.\n");
		exit(1);
	}
	else if  ( rstack_len ==0 )
	{
		printf("ERROR in parser: rdiff stack has size 0 .\n");
		exit(1);
	}
	return rdiff_stack[rstack_len-1];
}

mat_t GradientVectorProduct(const RDiff &X, mat_t &D)
{
	if ( X.in_m != D.m || X.in_n != D.n )
	{
		printf("ERROR! Hessian vector product dimension wrong X.in(%d %d), D(%d %d)\n", X.in_m, X.in_n, D.m, D.n);
		exit(1);
	}

	mat_t result = zeros(X.in_m, X.in_n);

	CompGradient(&(X), &(D), &(result), 0);
	return result;
}

#endif
