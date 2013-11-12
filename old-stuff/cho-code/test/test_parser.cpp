#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "MatrixUtil.h"
#include "RDiff.h"
#include "MatrixScalarDiff.h"


using namespace std;

/*****
 * *///


int main(int argc, char* argv[]){
	
	FILE *fin = fopen(argv[1],"r");
	char tmp[1000];
	fscanf(fin, "%s", tmp);
	mat_t X;
	X.load_ascii(argv[2], "X");
	vector<mat_t> const_list(1);
	const_list[0].load_ascii(argv[3], "S");
	X.print(stdout);
	const_list[0].print(stdout);
	printf("Begin parse %s\n", tmp);
	RDiff *gg = parsestring(tmp, X, const_list);
	MatrixScalarDiff result = logdet(*gg);
	result.diff.print(stdout);
//	gg->value.print(stdout);

	return 0;
}

