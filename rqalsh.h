#pragma once

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

// -----------------------------------------------------------------------------
//  RQALSH: basic data structure for high-dimensional c-k-AFN search
// -----------------------------------------------------------------------------
class RQALSH {
public:
	// -------------------------------------------------------------------------
	RQALSH(							// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float ratio,					// approximation ratio
		const int *index,				// index of data objects
		const float *data);				// data objects

	// -------------------------------------------------------------------------
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int   top_k,					// top-k value
		float R,						// limited search range
		const float *query,				// input query
		MaxK_List *list);				// c-k-AFN results (return)

	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		if (proj_a_ != NULL) ret += SIZEFLOAT * m_ * dim_; // proj_a_
		if (tables_ != NULL) ret += sizeof(Result) * m_ * n_pts_; // tables_
		return ret;
	}

protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	float ratio_;					// approximation ratio
	float w_;						// bucket width
	int   m_;						// number of hash tables
	int   l_;						// collision threshold
	const int   *index_;			// index of data objects
	const float *data_;				// data objects

	float  *proj_a_;				// hash functions
	Result *tables_;				// hash tables
	
	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// hash table id
		const float *data);				// one data object

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius					
		const int   *l_pos,				// left  position of query in hash table
		const int   *r_pos,				// right position of query in hash table
		const float *q_val);			// hash value of query
};
