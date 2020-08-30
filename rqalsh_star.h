#pragma once

#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <vector>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "rqalsh.h"

// -----------------------------------------------------------------------------
//  RQALSH_STAR is used for c-k-Approximate Furthest Neighbor (c-k-AFN) search
// -----------------------------------------------------------------------------
class RQALSH_STAR {
public:
	RQALSH_STAR(					// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		int   L,						// number of projection
		int   M,						// number of candidates
		float ratio,					// approximation ratio
		const float *data);				// data objects

	// -------------------------------------------------------------------------
	~RQALSH_STAR();					// destructor

	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int top_k,						// top-k value
		const float *query,				// query object
		MaxK_List *list);				// top-k results (return)
	
	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += SIZEINT * L_ * M_; 			// cand_
		ret += lsh_->get_memory_usage(); 	// lsh_
		return ret;
	}

protected:
	int    n_pts_;					// cardinality
	int    dim_;					// dimensionality
	int    L_;						// number of projections
	int    M_;						// number of candidates for each proj
	const float *data_;				// data objects

	int    *cand_;					// candidate data objects id
	RQALSH *lsh_;					// index of sample data objects

	// -------------------------------------------------------------------------
	void data_dependent_select(		// data dependent selection
		const float *data,				// data objects
		int   *cand);					// candidate id (return)
			
	// -------------------------------------------------------------------------
	void calc_shift_data(			// calculate shift data objects
		int   &max_id,					// data id with max l2-norm (return)
		float &max_norm,				// max l2-norm (return)
		float *norm,					// l2-norm of shift data (return)
		float *shift_data); 			// shift data (return)
};
