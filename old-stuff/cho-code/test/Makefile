SRCDIR=../

CXX=g++

#LIBS=-L/home/pkambadu/Software/lib -llapack -lcblas  -lm
ifeq ($(shell uname -s),Darwin)
  LIBS=-framework vecLib 
endif
ifeq ($(shell uname -s),Linux)
  LIBS=-llapack -lblas -lgfortran
endif


FLAGS_OPT=-O3 -fPIC -I${SRCDIR}

ALL: go_test_parser go_test_hessian go_test_newton go_test_gradient

LIB_HEADERS=${SRCDIR}/MatrixUtil.h \
						${SRCDIR}/MatrixScalarDiff.h \
						${SRCDIR}/SymbolUtil.h \
						${SRCDIR}/RDiff.h \
						${SRCDIR}/OptimizationUtil.h 
						#${SRCDIR}/FwdDiff.h 

go_test: test.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

go_test_reverse: test_reverse.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

go_test_gradient: test_gradient.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

go_test_parser: test_parser.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

go_test_hessian: test_hessian.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

go_test_newton: test_newton.cpp ${LIB_HEADERS}
	${CXX} ${FLAGS_OPT} -o $@ $< ${LIBS} 

.PHONY: clean

clean:
	rm -rf go_*
