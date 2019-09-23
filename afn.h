#ifndef __AFN_H
#define __AFN_H

// -----------------------------------------------------------------------------
int linear_scan(					// k-FN search of linear scan
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int ml_rqalsh(						// c-k-AFN search of ML-RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int rqalsh_star(					// c-k-AFN search of RQALSH*
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	int   L,							// number of projection (drusilla)
	int   M,							// number of candidates (drusilla)
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int rqalsh(							// c-k-AFN search of RQALSH
	int   n,							// number of data objects
	int   qn,							// number of query objects
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int drusilla_select(				// c-k-AFN search of Drusilla-Select
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

// -----------------------------------------------------------------------------
int qdafn(							// c-k-AFN search of QDAFN
	int   n,							// number of data objects
	int   qn,							// number of query points
	int   d,							// number of dimensions
	int   L,							// number of projections
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float **data,					// data set
	const float **query,				// query set
	const Result **R, 					// truth set
	const char *out_path);				// output path

#endif // __AFN_H
