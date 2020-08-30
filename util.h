#pragma once

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

extern timeval g_start_time;		// global param: start time
extern timeval g_end_time;			// global param: end time

extern float g_indextime;			// global param: indexing time (seconds)
extern float g_memory;				// global param: estimated memory usage (MB)

extern float g_runtime;				// global param: running time (ms)
extern float g_ratio;				// global param: overall ratio
extern float g_recall;				// global param: recall (%)
extern float g_fraction;			// global param: fraction (%)

// -----------------------------------------------------------------------------
//  uitlity functions
// -----------------------------------------------------------------------------
void create_dir(					// create directory
	char *path);						// input path

// -----------------------------------------------------------------------------
int read_bin_data(					// read data (binary) from disk
	int   n,							// number of data objects
	int   d,							// dimensionality
	bool  is_data,						// is dataset?
	const char *fname,					// address of data
	float *data);						// data (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int    qn,							// number of query objects
	const  char *fname,					// address of truth set
	Result *R);							// ground truth results (return)

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
	const float *data,					// data set
	const float *query,					// query set
	const char  *truth_set);			// address of truth set

// -----------------------------------------------------------------------------
int k_fn_search(					// k-FN search
	int   n, 							// cardinality
	int   d, 							// dimensionality
	const float *data,					// data objects
	const float *query,					// query objects
	MaxK_List *list);					// top-k results
