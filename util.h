#ifndef __UTIL_H
#define __UTIL_H

// -----------------------------------------------------------------------------
//  basic data structures
// -----------------------------------------------------------------------------
class  MaxK_List;

struct Result {						// basic data structure 
	float key_;
	int   id_;
};

// -----------------------------------------------------------------------------
int ResultComp(						// compare function for qsort (ascending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
int ResultCompDesc(					// compare function for qsort (descending)
	const void *e1,						// 1st element
	const void *e2);					// 2nd element

// -----------------------------------------------------------------------------
//  uitlity functions
// -------------------------------------------------------------------------
int create_dir(						// create directory
	char *path);						// input path

// -----------------------------------------------------------------------------
float calc_l2_dist(					// calc L_2 norm (data type is float)
	int dim,							// dimension
	const float* vec1,					// 1st point
	const float* vec2);					// 2nd point

// -----------------------------------------------------------------------------
float calc_recall(					// calc recall (percentage)
	int k,								// top-k value
	const Result *R,					// ground truth results 
	MaxK_List *list);					// results returned by algorithms

// -----------------------------------------------------------------------------
//  functions used for the input/output of data sets and query sets.
// -----------------------------------------------------------------------------
int read_data(						// read data/query set from disk
	int   n,							// number of data/query objects
	int   d,							// dimensionality
	const char *fname,					// address of data/query set
	float **data);						// data/query objects (return)

// -----------------------------------------------------------------------------
int read_ground_truth(				// read ground truth results from disk
	int qn,								// number of query objects
	const char *fname,					// address of truth set
	Result **R);						// ground truth results (return)

#endif
