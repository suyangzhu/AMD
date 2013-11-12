#ifndef MatrixMatrixFunc2_H
#define MatrixMatrixFunc2_H

#include <string>
#include <cstdio>
#include <assert.h>
#include "boost/shared_ptr.hpp"
#include "utility.hpp"
#include "ScalarMatrixFunc.hpp"

// definitions of enum, operators and callBackFunctions
#include "MatrixMatrixFuncHelper.hpp" 

namespace AMD {

template <class MT, class ST> 
class MatrixMatrixFunc2 {

public:
  typedef MT MatrixType;
  typedef ST ScalarType;
  typedef void (*CallBackFuncType)( boost::shared_ptr<MatrixType>, 
				      boost::shared_ptr<MatrixType>,
				      boost::shared_ptr<MatrixType>, 
				      boost::shared_ptr<MatrixType>,
				      const MatrixMatrixFunc2<MT,ST>*,
				      int&, bool&, bool& );

  // Once recorded the matrix should never be changed - so a pointer
  // is safe
  boost::shared_ptr<MT> matrixPtr; 
  /// Used to compute derivative after evaluation tree is recorded.
  CallBackFuncType callBackFunc;
  OpType opNum;
  bool isConst;
  int varNumRows; // number of rows in matrix variable
  int varNumCols; // number of cols in matrix variable

  MatrixMatrixFunc2* leftChild;
  MatrixMatrixFunc2* rightChild;


  MatrixMatrixFunc2() : matrixPtr(), callBackFunc(NULL), opNum(NONE), 
			isConst(true), varNumRows(0), varNumCols(0), 
			leftChild(NULL), rightChild(NULL)
  { }

  // Makes an expensive copy of matrix -- avoid this constructor
  // if your matrices are large.
  MatrixMatrixFunc2(MT matrix, bool isVariable=false) 
    : matrixPtr(), leftChild(NULL), rightChild(NULL) { 
    boost::shared_ptr<MT> copy (new MT(matrix));
    matrixPtr = copy;
    setVariableType(isVariable);
  }



  MatrixMatrixFunc2(boost::shared_ptr<MT> _matrixPtr, bool isVariable=false) 
    : matrixPtr(), leftChild(NULL), rightChild(NULL) { 
    matrixPtr = _matrixPtr;
    setVariableType(isVariable);
  }

  // set up variables to determine whether matrix is constant or a variable
  void setVariableType( bool isVariable ) {
    // constant and variable matrices must be leaf nodes
    assert( NULL == leftChild && NULL == rightChild );

    callBackFunc = constOp<MatrixType,ScalarType>;
    opNum = CONST;
    isConst = true;
    varNumRows = 0;
    varNumCols = 0;

    if (isVariable) {
      varNumRows = matrixPtr->getNumRows();
      varNumCols = matrixPtr->getNumCols();
      callBackFunc = varOp<MatrixType,ScalarType>;
      opNum = VAR;
      isConst = false;
    }
  }

  // Need to recursively delete information
  ~MatrixMatrixFunc2() {
    if (NULL!=leftChild) {
      delete leftChild;
    }
    if (NULL!=rightChild) {
      delete rightChild;
    }
  }

  // need to recursively copy information
  MatrixMatrixFunc2& operator= ( const MatrixMatrixFunc2 &x) {
    deepCopy(*this, x);
    return(*this);
  }

  void shallowCopy( MatrixMatrixFunc2 &copy, const MatrixMatrixFunc2 &x) { 
    copy.matrixPtr = x.matrixPtr;
    copy.opNum = x.opNum;
    copy.isConst = x.isConst;
    copy.callBackFunc = x.callBackFunc;
    copy.varNumRows = x.varNumRows;
    copy.varNumCols = x.varNumCols;
    copy.leftChild = NULL;
    copy.rightChild = NULL;
  }

  void deepCopy( MatrixMatrixFunc2 &copy, const MatrixMatrixFunc2 &x) { 
    if (copy.leftChild) {
      delete copy.leftChild;
      copy.leftChild = NULL;
    }
    if (copy.rightChild) {
      delete copy.rightChild;
      copy.rightChild = NULL;
    }
    shallowCopy(copy,x);
    if (NULL != x.leftChild) {
      copy.leftChild = new MatrixMatrixFunc2<MT,ST>;
      deepCopy(*(copy.leftChild), *(x.leftChild));
    }
    if (NULL != x.rightChild) {
      copy.rightChild = new MatrixMatrixFunc2<MT,ST>;
      deepCopy(*(copy.rightChild), *(x.rightChild));
    }
  }

  void binOpSet( boost::shared_ptr<MT> resultPtr,
		 OpType _opNum,
		 CallBackFuncType cbf,
		 const MatrixMatrixFunc2<MT,ST> &lhs, 
		 const MatrixMatrixFunc2<MT,ST> &rhs ) {
    unaryOpSet(resultPtr,_opNum,cbf,lhs);
    rightChild = new MatrixMatrixFunc2<MT,ST>;
    deepCopy(*rightChild,rhs);

    isConst = lhs.isConst && rhs.isConst;
    varNumRows = lhs.varNumRows | rhs.varNumRows;  // should be 0 or a size
    varNumCols = lhs.varNumCols | rhs.varNumCols;
  }

  void unaryOpSet( boost::shared_ptr<MT> resultPtr,
		   OpType _opNum,
		   CallBackFuncType cbf,
		   const MatrixMatrixFunc2<MT,ST> &lhs ) {
    varNumRows = lhs.varNumRows;  // should be 0 or a size
    varNumCols = lhs.varNumCols;
    matrixPtr = resultPtr;
    opNum = _opNum;
    isConst = lhs.isConst;
    callBackFunc = cbf;

    leftChild = new MatrixMatrixFunc2<MT,ST>;
    deepCopy(*leftChild,lhs);
  }


  void print(FILE *fp) const {
    matrixPtr->print();
    if (NULL==leftChild && NULL==rightChild) {
      fprintf(fp,":%s",opName[opNum].c_str());
    } else {
      fprintf(fp,":%s(",opName[opNum].c_str());
      if (NULL!=leftChild) {
	leftChild->print(fp);
      } 
      if (NULL!=rightChild) {
	fprintf(fp,",");
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
    return(matrixPtr->numRows());
  }

  int numCols() {
    return(matrixPtr->numCols());
  }

  virtual MatrixType value() {
    return(*matrixPtr);
  }

  // initial and result must point to existing MatrixTypes
  void gradientVec( boost::shared_ptr<MatrixType> initial, 
		    boost::shared_ptr<MatrixType> result, 
		    int transposeFlag, 
		    bool identityInitialFlag,
		    bool& zeroResultFlag ) const {
    assert( initial.use_count()>=1 && // must be a valid pointer
	    result.use_count()>=1  &&
	    isConst || ( result->getNumRows() == varNumRows &&
			 result->getNumCols() == varNumCols ) );
    // If the function is constant then the derivative is zero -- so
    // we don't need to do any computation
    if (!isConst) { 
      boost::shared_ptr<MatrixType> leftPtr(new MatrixType);
      if (NULL != leftChild) *leftPtr  = *(leftChild->matrixPtr);
      boost::shared_ptr<MatrixType> rightPtr(new MatrixType);
      if (NULL != rightChild) *rightPtr  = *(rightChild->matrixPtr);

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



}  /** namespace AMD */


#endif
