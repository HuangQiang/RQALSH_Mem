#ifndef __RQALSH_H
#define __RQALSH_H

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include "def.h"
#include "util.h"
#include "random.h"
#include "pri_queue.h"

struct Result;
class  MaxK_List;

// -----------------------------------------------------------------------------
//  RQALSH: basic data structure for high-dimensional c-k-AFN search
// -----------------------------------------------------------------------------
class RQALSH {
public:
	RQALSH(							// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float ratio);					// approximation ratio

	// -------------------------------------------------------------------------
	RQALSH(							// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// hash table id
		const float *data);				// one data object

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int   top_k,					// top-k value
	    const float *query,				// query object
		MaxK_List *list);				// c-k-AFN results (return)

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int   top_k,					// top-k value
		float R,						// limited search range
		const float *query,				// input query
		std::vector<int> &cand);		// candidate id (return)

	// -------------------------------------------------------------------------
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	float ratio_;					// approximation ratio
	float w_;						// bucket width
	int   m_;						// number of hash tables
	int   l_;						// collision threshold
	const float **data_;			// data objects

	float **a_;						// hash functions
	Result **tables_;				// hash tables

protected:
	// -------------------------------------------------------------------------
	void init();					// init parameters
	
	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius					
		const int *lpos,				// left  position of query in hash table
		const int *rpos,				// right position of query in hash table
		const float *q_val);			// hash value of query
};

#endif // __RQALSH_H
