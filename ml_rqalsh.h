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
//  ML_RQALSH: Multi-Level RQALSH for c-k-AFN Search 
// -----------------------------------------------------------------------------
class ML_RQALSH {
public:
	ML_RQALSH(						// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float ratio,					// approximation ratio
    	const float *data);				// data objects

	// -------------------------------------------------------------------------
	~ML_RQALSH();					// destructor
	
	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN seach	
		int top_k,		    			// top-k value
		const float *query,				// input query
		MaxK_List *list);				// top-k results (return)

	// -------------------------------------------------------------------------
	int64_t get_memory_usage()		// get memory usage
	{
		int64_t ret = 0;
		ret += sizeof(*this);
		ret += SIZEFLOAT * dim_; 	// centroid_
		ret += SIZEINT * n_pts_;	// sorted_id_
		ret += SIZEFLOAT * radius_.capacity(); // radius_
		for (auto lsh : lsh_) {		// blocks_
			ret += lsh->get_memory_usage();
		}
		return ret;
	}

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	float ratio_;					// approximation ratio
	const float *data_;				// data objects

	int   *sorted_id_;				// sorted data id after ml-partition
	float *centroid_;				// centroid of data objects
	std::vector<float> radius_;		// radius
	std::vector<RQALSH*> lsh_;		// blocks
};
