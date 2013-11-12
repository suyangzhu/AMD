# change CFLAGS to -Wall -g to compile for 
# debug.
# CFLAGS = -Wall -O
CFLAGS = -Wall -g
CC = g++

OBJ_DIR = ../obj
EXE_DIR = ../exec

${OBJ_DIR}/%.o : %.C
	${CC} ${CFLAGS} -c $^ -o $@

${EXE_DIR}/test.exe : ${OBJ_DIR}/Diff2.o ${OBJ_DIR}/DiffMath.o ${OBJ_DIR}/test.o 
	g++ ${CFLAGS} ${OPT} -o $@ $?

${EXE_DIR}/newton.exe : ${OBJ_DIR}/Diff2.o ${OBJ_DIR}/DiffMath.o ${OBJ_DIR}/newton.o 
	g++ ${CFLAGS} ${OPT} -o $@ $?

${EXE_DIR}/testBinomial.exe : ${OBJ_DIR}/testBinomial.o 
	g++ ${CFLAGS} ${OPT} -o $@ $?

test : ${EXE_DIR}/test.exe
	echo "executable in $?."

newton : ${EXE_DIR}/newton.exe
	echo "executable in $?."

binomial : ${EXE_DIR}/testBinomial.exe
	echo "executable in $?."

clean :
	rm -f ${OBJ_DIR}/*.o
	rm -f ${EXE_DIR}/*.exe
	rm -f *~

