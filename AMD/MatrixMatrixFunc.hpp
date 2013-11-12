#ifndef MatrixMatrixFunc_H
#define MatrixMatrixFunc_H

#include <string>
#include <cstdio>
#include <assert.h>
#include "boost/shared_ptr.hpp"
#include "utility.hpp"
#include "ScalarMatrixFunc.hpp"


namespace AMD {

// forward declaration
template <class MT, class ST> class MatrixMatrixFunc;

/// Callback function for differentiation involving constant matrices
template <class MT,class ST>
void constOp(MT* result, MT* current, 
	     MT* left, MT* right,
	     const MatrixMatrixFunc<MT,ST>* node, 
	     int& transposeFlag,
	     bool& identityCurrentFlag, 
	     bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL == node->leftChild &&
	  NULL == node->rightChild &&
	  node->isConst );
}

/// Callback function for differentiation involving a variable matrix
/// on the leaf node.
template <class MT, class ST>
void varOp(MT* result, MT* current, 
	   MT* left, MT* right,
	   const MatrixMatrixFunc<MT,ST>* node, 
	   int& transposeFlag,
	   bool& identityCurrentFlag, 
	   bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL == node->leftChild &&
	  NULL == node->rightChild &&
	  "var" == node->opName &&
	  NULL!=result && // result and current must be valid
	  NULL!=current );
  if (zeroResultFlag) {
    zeroResultFlag = false;
    if (transposeFlag) {
      (*result) = transpose(*current);
    } else {
      (*result) = (*current);
    }
  } else {
    if (transposeFlag) {
      (*result) = (*result)+transpose(*current);
    } else {
      (*result) = (*result)+(*current);
    }
  }
}


/// Callback function for differentiation involving operator+
template <class MT, class ST>
void plusOp(MT* result, MT* current, 
	    MT* left, MT* right,
	    const MatrixMatrixFunc<MT,ST>* node, 
	    int& transposeFlag,
	    bool& identityCurrentFlag, 
	    bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL != node->leftChild &&
	  NULL != node->rightChild &&
	  "+" == node->opName &&
	  NULL!=current && // current, left and right must be present 
	  NULL!=left && 
	  NULL!=right );
  (*left)  = (*current);
  (*right) = (*current);
  if (transposeFlag) {
    transposeFlag=3; // both left and right should inherit transpose
  }
}

/// Callback function for differentiation involving operator-
template <class MT, class ST>
void minusOp(MT* result, MT* current, 
	     MT* left, MT* right,
	     const MatrixMatrixFunc<MT,ST>* node, 
	     int& transposeFlag,
	     bool& identityCurrentFlag, 
	     bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL != node->leftChild &&
	  NULL != node->rightChild &&
	  "-" == node->opName &&
	  NULL!=current && // current, left and right must be present 
	  NULL!=left && 
	  NULL!=right );
  (*left)  = (*current);
  (*right) = -(*current);
  if (transposeFlag) {
    transposeFlag=3; // both left and right should inherit transpose
  }
}

/// Callback function for differentiation involving operator*
template <class MT, class ST>
void timesOp(MT* result, MT* current, 
	     MT* left, MT* right,
	     const MatrixMatrixFunc<MT,ST>* node, 
	     int& transposeFlag,
	     bool& identityCurrentFlag, 
	     bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL != node->leftChild &&
	  NULL != node->rightChild &&
	  "*" == node->opName &&
	  NULL!=current && // current, left and right must be present 
	  NULL!=left && 
	  NULL!=right );
  if (identityCurrentFlag) { // avoid superfluous multiplication
    transposeFlag = 0;
    if ("transpose" == node->rightChild->opName) {
      // if right is R^T then get R from it's left child
      (*left) = node->rightChild->leftChild->val;
    } else {
      (*left) = node->rightChild->val;
      transposeFlag |= 1; // set left transpose on
    }
    if ("transpose"==node->leftChild->opName) {
      // if right is R^T then get R from it's left child
      (*right) = node->leftChild->leftChild->val;
    } else {
      (*right) = node->leftChild->val;
      transposeFlag |= 2; // set right transpose on
    }
    identityCurrentFlag = false;
  } else {
    if (transposeFlag) {
      // use A^T*B^T = (B*A)^T to reduce the numbder of trans
      //(*left) = transpose(*current) * transpose(node->rightChild->val);
      (*left) = node->rightChild->val * (*current);
      //(*right) = transpose(node->leftChild->val) * transpose(*current);
      (*right) = (*current) * node->leftChild->val;
      transposeFlag = 3;
    } else {
      (*left) = (*current) * transpose(node->rightChild->val);
      (*right) = transpose(node->leftChild->val) * (*current);
      transposeFlag = 0;
    }
  }
}


/// Callback function for differentiation involving operator*
template <class MT, class ST>
void transposeOp(MT* result, MT* current, 
		 MT* left, MT* right,
		 const MatrixMatrixFunc<MT,ST>* node, 
		 int& transposeFlag,
		 bool& identityCurrentFlag, 
		 bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL != node->leftChild &&
	  NULL == node->rightChild &&
	   "transpose" == node->opName &&
	  NULL!=current && // current, left and right must be present 
	  NULL!=left );

  (*left) = (*current);
  if (!identityCurrentFlag) {
    if (transposeFlag) {
      transposeFlag = 0;
    } else {
      transposeFlag = 3;  // transpose all living children
    }
  }
}

/// Callback function for differentiation involving operator*
template <class MT, class ST>
void invOp(MT* result, MT* current, 
	   MT* left, MT* right,
	   const MatrixMatrixFunc<MT,ST>* node, 
	   int& transposeFlag,
	   bool& identityCurrentFlag, 
	   bool& zeroResultFlag) {
  assert( NULL != node && // check node type
	  NULL != node->leftChild &&
	  NULL == node->rightChild &&
	  "inv" == node->opName &&
	  NULL!=current && // current, left and right must be present 
	  NULL!=left );
  // TO DO: If we had the pointer to the parent matrix we wouldn't have
  // to do the inversion
  const MT& tmp = node->val;
  if (!identityCurrentFlag) {
    if (transposeFlag) {
      (*left) = - tmp * (*current) * tmp;
    } else {
      (*left) = - tmp * transpose(*current) * tmp;
    }
  } else {
    (*left) = - tmp * tmp;
  }
  transposeFlag = 3;
  identityCurrentFlag = false;
}



template <class MT, class ST> 
class MatrixMatrixFunc {
public:
  typedef MT MatrixType;
  typedef ST ScalarType;

  MT val;
  /// Used to compute derivative after evaluation tree is recorded.
  void (*callBackFunc)( MatrixType*, MatrixType*, 
			MatrixType*, MatrixType*, 
			const MatrixMatrixFunc<MT,ST>*,
			int&, bool&, bool& );

  std::string opName;
  bool isConst;

  MatrixMatrixFunc* leftChild;
  MatrixMatrixFunc* rightChild;
  //  const MatrixMatrixFunc* leftChild;
  //  const MatrixMatrixFunc* rightChild;
  int varNumRows; // number of rows in matrix variable
  int varNumCols; // number of cols in matrix variable

  MatrixMatrixFunc() : val(), callBackFunc(NULL), isConst(true), 
		       leftChild(NULL), rightChild(NULL),
		       varNumRows(0), varNumCols(0)
  { opName = "none"; }

  MatrixMatrixFunc(MT val) 
    : val(val), isConst(true), leftChild(NULL),  rightChild(NULL),
      varNumRows(0), varNumCols(0) { 
    opName = "const"; 
    callBackFunc = constOp<MatrixType,ScalarType>;
  }

  MatrixMatrixFunc(MatrixType &mat, bool isVariable) 
    : val(mat), callBackFunc(NULL), isConst(true), leftChild(NULL),  rightChild(NULL),
      varNumRows(0), varNumCols(0) {
    callBackFunc = constOp<MatrixType,ScalarType>;
    opName = "const";
    isConst = true;
    if (isVariable) {
      varNumRows = mat.getNumRows();
      varNumCols = mat.getNumCols();
      callBackFunc = varOp<MatrixType,ScalarType>;
      opName = "var";
      isConst = false;
    }
  }
  
  ~MatrixMatrixFunc() {
    if (NULL!=leftChild) {
      //delete leftChild;
    }
    if (NULL!=rightChild) {
      //delete rightChild;
    }
  }

  MatrixMatrixFunc& operator= ( const MatrixMatrixFunc &x) {
    val = x.val;
    leftChild = x.leftChild;
    rightChild = x.rightChild;
    opName = x.opName;
    isConst = x.isConst;
    callBackFunc = x.callBackFunc;
    varNumRows = x.varNumRows;
    varNumCols = x.varNumCols;
    return(*this);
  }

  void print(FILE *fp) const {
    val.print();
    if (NULL==leftChild && NULL==rightChild) {
      fprintf(fp,":%s",opName.c_str());
    } else {
      fprintf(fp,":%s(",opName.c_str());
      if (NULL!=leftChild) {
	leftChild->print(fp);
      } 
      fprintf(fp,",");
      if (NULL!=rightChild) {
	rightChild->print(fp);
      } 
      fprintf(fp,")");
    }
  }

  void println(FILE *fp) const {
    print(fp);
    fprintf(fp,"\n");
  }

  void print() const {
    print(stdout);
  }

  void println() const {
    println(stdout);
  }

  int numRows() {
    return(val.numRows());
  }

  int numCols() {
    return(val.numCols());
  }

  virtual MatrixType value() {
    return val;
  }

  // initial and result must point to existing MatrixTypes
  void gradientVec(MatrixType* initial, 
		   MatrixType* result, 
		   int transposeFlag, 
		   bool identityInitialFlag,
		   bool& zeroResultFlag) const {
    assert( NULL != initial && NULL != result &&
	    isConst || ( result->getNumRows() == varNumRows &&
			 result->getNumCols() == varNumCols ) );
    // If the function is constant then the derivative is zero -- so
    // we don't need to do any computation
    if (!isConst) { 
      MatrixType left;
      if (leftChild) left = leftChild->val;
      MatrixType right;
      if (rightChild) right = rightChild->val;
      MatrixType* leftPtr  = &left;
      MatrixType* rightPtr = &right;

      callBackFunc(result, initial, leftPtr, rightPtr, this,
		   transposeFlag,identityInitialFlag, zeroResultFlag);
      if (NULL!=leftChild) {
	int leftFlag = transposeFlag & 1;
	leftChild->gradientVec( leftPtr, 
				result, 
				leftFlag, 
				identityInitialFlag,
				zeroResultFlag );
      }
      if (NULL!=rightChild) {
	int rightFlag = transposeFlag & 2;
	rightChild->gradientVec( rightPtr, 
				 result, 
				 rightFlag, 
				 identityInitialFlag,
				 zeroResultFlag );
      }
    } // end if (!isConst) 

  }

};

template <class MT, class ST>
MatrixMatrixFunc<MT,ST> operator+ (const MatrixMatrixFunc<MT,ST> &lhs, 
				   const MatrixMatrixFunc<MT,ST> &rhs)
{
  assert( lhs.isConst || rhs.isConst || 
	  (lhs.varNumRows==rhs.varNumRows && 
	   lhs.varNumCols==rhs.varNumCols ) );
  MatrixMatrixFunc<MT,ST> result;
  result.varNumRows = lhs.varNumRows | rhs.varNumRows;  // should be 0 or a size
  result.varNumCols = lhs.varNumCols | rhs.varNumCols;
  result.val = lhs.val + rhs.val;
  result.leftChild = detail::pointerCopy(lhs);
  result.rightChild = detail::pointerCopy(rhs);
  result.opName = "+";
  result.isConst = lhs.isConst && rhs.isConst;
  result.callBackFunc = plusOp<MT,ST>;
  return(result);
}

template <class MT, class ST>
MatrixMatrixFunc<MT,ST> operator- (const MatrixMatrixFunc<MT,ST> &lhs, 
				   const MatrixMatrixFunc<MT,ST> &rhs)
{
  assert( lhs.isConst || rhs.isConst || 
	  (lhs.varNumRows==rhs.varNumRows && 
	   lhs.varNumCols==rhs.varNumCols ) );
  MatrixMatrixFunc<MT,ST> result;
  result.varNumRows = lhs.varNumRows | rhs.varNumRows;  // should be 0 or a size
  result.varNumCols = lhs.varNumCols | rhs.varNumCols;
  result.val = lhs.val - rhs.val;
  result.leftChild = detail::pointerCopy(lhs);
  result.rightChild = detail::pointerCopy(rhs);
  result.opName = "-";
  result.isConst = lhs.isConst && rhs.isConst;
  result.callBackFunc = minusOp<MT,ST>;
  return(result);
}

template <class MT, class ST>
MatrixMatrixFunc<MT,ST> operator* (const MatrixMatrixFunc<MT,ST> &lhs, 
				   const MatrixMatrixFunc<MT,ST> &rhs)
{
  assert( lhs.isConst || rhs.isConst || 
	  (lhs.varNumRows==rhs.varNumRows && 
	   lhs.varNumCols==rhs.varNumCols ) );
  MatrixMatrixFunc<MT,ST> result;
  result.varNumRows = lhs.varNumRows | rhs.varNumRows;  // should be 0 or a size
  result.varNumCols = lhs.varNumCols | rhs.varNumCols;
  result.val = lhs.val * rhs.val;
  result.leftChild = detail::pointerCopy(lhs);
  result.rightChild = detail::pointerCopy(rhs);
  //result.leftChild  = &(lhs);
  //result.rightChild = &(rhs);
  result.opName = "*";
  result.isConst = lhs.isConst && rhs.isConst;
  result.callBackFunc = timesOp<MT,ST>;
  //   fprintf(stdout,"exit operator *\n");
  return(result);
}

template <class MT, class ST>
MatrixMatrixFunc<MT,ST> transpose(const MatrixMatrixFunc<MT,ST> &lhs)
{
  //  fprintf(stdout,"enter operator *\n");
  MatrixMatrixFunc<MT,ST> result;
  result.varNumRows = lhs.varNumRows;
  result.varNumCols = lhs.varNumCols;
  result.opName = "transpose";
  if (lhs.opName != result.opName) {
    result.val = transpose(lhs.val);
    result.leftChild = detail::pointerCopy(lhs);
    result.rightChild = NULL;
    result.isConst = lhs.isConst;
    result.callBackFunc = transposeOp<MT,ST>;
  } else { // 2 transposes cancel each other out
    assert(NULL!=lhs.leftChild);
    result = (*lhs.leftChild);
  }
  //   fprintf(stdout,"exit operator *\n");
  return(result);
}


template <class MT, class ST>
MatrixMatrixFunc<MT,ST> inv(const MatrixMatrixFunc<MT,ST> &lhs)
{
  //  fprintf(stdout,"enter operator *\n");
  MatrixMatrixFunc<MT,ST> result;
  result.varNumRows = lhs.varNumRows;
  result.varNumCols = lhs.varNumCols;
  result.opName = "inv";
  if (lhs.opName != result.opName) {
    result.val = inv(lhs.val);
    result.isConst = lhs.isConst;
    result.callBackFunc = invOp<MT,ST>;
    result.rightChild = NULL;
    result.leftChild = detail::pointerCopy(lhs);
  } else { // 2 inverses cancel each other out
    assert(NULL!=lhs.leftChild);
    result = (*lhs.leftChild);
  }
  //   fprintf(stdout,"exit operator *\n");
  return(result);
}

template <class MT, class ST>
ScalarMatrixFunc<MT,ST> trace(const MatrixMatrixFunc<MT,ST> &lhs) {
  // matrix must be square in order to compute trace
  assert( lhs.val.getNumRows() == lhs.val.getNumCols() );
  int n = lhs.val.getNumRows();
  MT tmp;
  MT eye = tmp.eye(n);
  MT res = tmp.zeros(lhs.varNumRows,lhs.varNumCols);
  MT* initPtr = &(eye);
  MT* resPtr = &(res);
  bool zeroFlag = true;
  lhs.gradientVec(initPtr,resPtr,false,true,zeroFlag);
  if (zeroFlag) { 
    ScalarMatrixFunc<MT,ST> result(trace(lhs.val),lhs.varNumRows,lhs.varNumCols);
    // pass on knowledge that function is constant
    return(result);
  } else {
    ScalarMatrixFunc<MT,ST> result(trace(lhs.val),res);
    return(result);
  }
}

template <class MT, class ST>
ScalarMatrixFunc<MT,ST> logdet(const MatrixMatrixFunc<MT,ST> &lhs) {
  // matrix must be square in order to compute trace
  assert( lhs.val.getNumRows() == lhs.val.getNumCols() );
  int n = lhs.val.getNumRows();
  MT tmp;
  MT res = tmp.zeros(lhs.varNumRows,lhs.varNumCols);
  MT* resPtr = &(res);
  bool transposeFlag = true;
  MT init;
  if ("transpose"==lhs.opName) { // det(X^T) == det(X)
    return(logdet((*lhs.leftChild)));
    //init = inv(lhs.leftChild->val);
    //transposeFlag = false;
  }
  if ("inv" == lhs.opName) { // logdet(X^{-1)) == - logdet(X)
    return(-logdet((*lhs.leftChild)));
  }
  //else {
  //  init = inv(lhs.val);
  //  transposeFlag = true;
  //}
  init = inv(lhs.val);
  MT* initPtr = &(init);
  bool zeroFlag = true;
  lhs.gradientVec(initPtr,resPtr,transposeFlag,false,zeroFlag);
  if (zeroFlag) { 
    ScalarMatrixFunc<MT,ST> result(logdet(lhs.val),lhs.varNumRows,lhs.varNumCols);
    // pass on knowledge that function is constant
    return(result);
  } else {
    ScalarMatrixFunc<MT,ST> result(logdet(lhs.val),res);
    return(result);
  }
}

}  /** namespace AMD */

#endif
