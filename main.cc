#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include "def.h"
#include "util.h"
#include "afn.h"

// -----------------------------------------------------------------------------
void usage() 						// usage of the package
{
	printf("\n"
		"--------------------------------------------------------------------\n"
		" Usage of the Package for Internal c-k-AFN Search:                  \n"
		"--------------------------------------------------------------------\n"
		"    -alg   (integer)   options of algorithms (0 - 6)\n"
		"    -n     (integer)   number of data  objects\n"
		"    -qn    (integer)   number of query objects\n"
		"    -d     (integer)   dimensionality\n"
		"    -L	    (integer)   number of projection\n"
		"    -M     (integer)   number of candidates\n"
		"    -c     (real)      approximation ratio (c > 1)\n"
		"    -ds    (string)    address of data  set\n"
		"    -qs    (string)    address of query set\n"
		"    -ts    (string)    address of truth set\n"
		"    -op    (string)    output path\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" The Options of Algorithms (-alg) are:                              \n"
		"--------------------------------------------------------------------\n"
		"    0 - Ground-Truth\n"
		"        Params: -alg 0 -n -qn -d -ds -qs -ts\n"
		"\n"
		"    1 - Linear Scan\n"
		"        Params: -alg 1 -n -qn -d -ds -qs -ts -op\n"
		"\n"
		"    2 - QDAFN\n"
		"        Params: -alg 2 -n -qn -d -L -M -c -ds -qs -ts -op\n"
		"\n"
		"    3 - Drusilla Select\n"
		"        Params: -alg 3 -n -qn -d -L -M -ds -qs -ts -op\n"
		"\n"
		"    4 - RQALSH\n"
		"        Params: -alg 4 -n -qn -d -c -ds -qs -ts -op\n"
		"\n"
		"    5 - RQALSH*\n"
		"        Params: -alg 5 -n -qn -d -L -M -c -ds -qs -ts -op\n"
		"\n"
		"    6 - ML_RQALSH\n"
		"        Params: -alg 6 -n -qn -d -c -ds -qs -ts -op\n"
		"\n"
		"--------------------------------------------------------------------\n"
		" Author: Qiang HUANG  (huangq2011@gmail.com)                        \n"
		"--------------------------------------------------------------------\n"
		"\n\n\n");
}

// -----------------------------------------------------------------------------
int main(int nargs, char** args)
{
	srand(6);						// set the random seed
	// usage();

	char   data_set[200];			// address of data  set
	char   query_set[200];			// address of query set
	char   truth_set[200];			// address of truth set
	char   out_path[200];			// output path

	int    alg    = -1;				// option of algorithm
	int    n      = -1;				// cardinality
	int    qn     = -1;				// query number
	int    d      = -1;				// dimensionality
	int    L      = -1;				// number of projection
	int    M      = -1;				// number of candidates
	float  ratio  = -1.0f;			// approximation ratio
	float  *data  = NULL;			// data set
	float  *query = NULL;			// query set
	Result *R     = NULL;			// k-NN ground truth
	int    cnt    = 1;

	while (cnt < nargs) {
		if (strcmp(args[cnt], "-alg") == 0) {
			alg = atoi(args[++cnt]);
			printf("alg       = %d\n", alg);
			assert(alg >= 0);
		}
		else if (strcmp(args[cnt], "-n") == 0) {
			n = atoi(args[++cnt]);
			printf("n         = %d\n", n);
			assert(n > 0);
		}
		else if (strcmp(args[cnt], "-qn") == 0) {
			qn = atoi(args[++cnt]);
			printf("qn        = %d\n", qn);
			assert(qn > 0);
		}
		else if (strcmp(args[cnt], "-d") == 0) {
			d = atoi(args[++cnt]);
			printf("d         = %d\n", d);
			assert(d > 0);
		}
		else if (strcmp(args[cnt], "-L") == 0) {
			L = atoi(args[++cnt]);
			printf("L         = %d\n", L);
			assert(L >= 0);
		}
		else if (strcmp(args[cnt], "-M") == 0) {
			M = atoi(args[++cnt]);
			printf("M         = %d\n", M);
			assert(M >= 0);
		}
		else if (strcmp(args[cnt], "-c") == 0) {
			ratio = (float) atof(args[++cnt]);
			printf("c         = %.1f\n", ratio);
			assert(ratio > 1.0f);
		}
		else if (strcmp(args[cnt], "-ds") == 0) {
			strncpy(data_set, args[++cnt], sizeof(data_set));
			printf("data_set  = %s\n", data_set);
		}
		else if (strcmp(args[cnt], "-qs") == 0) {
			strncpy(query_set, args[++cnt], sizeof(query_set));
			printf("query_set = %s\n", query_set);
		}
		else if (strcmp(args[cnt], "-ts") == 0) {
			strncpy(truth_set, args[++cnt], sizeof(truth_set));
			printf("truth_set = %s\n", truth_set);
		}
		else if (strcmp(args[cnt], "-op") == 0) {
			strncpy(out_path, args[++cnt], sizeof(out_path));
			printf("out_path  = %s\n", out_path);

			int len = (int) strlen(out_path);
			if (out_path[len - 1] != '/') {
				out_path[len] = '/';
				out_path[len + 1] = '\0';
			}
			create_dir(out_path);
		}
		else {
			usage();
			exit(1);
		}
		++cnt;
	}
	printf("\n");

	// -------------------------------------------------------------------------
	//  read data set, query set, and truth set (optional)
	// -------------------------------------------------------------------------
	data = new float[n * d];
	if (read_bin_data(n, d, true, data_set, data)) exit(1);

	query = new float[qn * d];
	if (read_bin_data(qn, d, false, query_set, query)) exit(1);

	if (alg > 0) {
		R = new Result[qn * MAXK];
		if (read_ground_truth(qn, truth_set, R)) exit(1);
	}

	// -------------------------------------------------------------------------
	//  methods
	// -------------------------------------------------------------------------
	switch (alg) {
	case 0:
		ground_truth(n, qn, d, (const float*) data, (const float*) query, 
			truth_set);
		break;
	case 1:
		linear_scan(n, qn, d, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	case 2:
		qdafn(n, qn, d, L, M, ratio, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	case 3:
		drusilla_select(n, qn, d, L, M, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	case 4:
		rqalsh(n, qn, d, ratio, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	case 5:
		rqalsh_star(n, qn, d, L, M, ratio, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	case 6:
		ml_rqalsh(n, qn, d, ratio, (const float*) data, (const float*) query, 
			(const Result*) R, out_path);
		break;
	default:
		printf("Parameters Error!\n");
		usage();
		break;
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] data;
	delete[] query;
	if (alg > 0) delete[] R;

	return 0;
}
