#ifndef __AFN_H
#define __AFN_H

// -----------------------------------------------------------------------------
int ground_truth(					// find ground truth
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set);				// address of truth set

// -----------------------------------------------------------------------------
int rqalsh_star(					// c-k-AFN search of RQALSH*
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int rqalsh(							// c-k-AFN search of RQALSH
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   beta,							// false positive percentage
	float delta,						// error probability
	float ratio,						// approximation ratio
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder);			// output folder

// -----------------------------------------------------------------------------
int linear_scan(					// k-FN search of linear scan
	int   n,							// number of data  objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const char *data_set,				// address of data  set
	const char *query_set,				// address of query set
	const char *truth_set,				// address of truth set
	const char *output_folder);			// output folder

#endif // __AFN_H
