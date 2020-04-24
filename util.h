#ifndef __UTIL_H
#define __UTIL_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include <unistd.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

#include "def.h"
#include "pri_queue.h"

struct Result;
class  MaxK_List;

extern timeval g_start_time;		// global parameter: start time
extern timeval g_end_time;			// global parameter: end time

extern float    g_runtime;			// global parameter: running time
extern float    g_ratio;			// global parameter: overall ratio
extern float    g_recall;			// global parameter: recall
extern uint64_t g_memory;			// global parameter: memory usage

// -----------------------------------------------------------------------------
//  uitlity functions
// -------------------------------------------------------------------------
void create_dir(					// create directory
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_txt_data(					// read data (text) from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *fname,					// address of data
	float **data);						// data (return)

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	const char *fname,					// address of data
	float **data);						// data (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_inner_product(			// calc inner product (data type is float)
	int   dim,							// dimension
	const float *p1,					// 1st point
	const float *p2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int   k,							// top-k value
	const Result *R,					// ground truth results 
	const Result *result);				// results returned by algorithms

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data set
	const float **query,				// query set
	const char *truth_set);				// address of truth set

// -----------------------------------------------------------------------------
void k_fn_search(					// k-FN search
	int   n, 							// cardinality
	int   qn,							// query number
	int   d, 							// dimensionality
	int   k,							// top-k value
	const float **data,					// data objects
	const float **query,				// query objects
	Result **result);					// k-MIP results (return)

#endif
