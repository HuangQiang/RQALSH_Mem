#ifndef __RQALSH_STAR_H
#define __RQALSH_STAR_H

class RQALSH;
class MaxK_List;

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
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~RQALSH_STAR();					// destructor

	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int top_k,						// top-k value
		const float *query,				// query object
		MaxK_List *list);				// top-k results (return)

protected:
	int    n_pts_;					// cardinality
	int    dim_;					// dimensionality
	int    L_;						// number of projections
	int    M_;						// number of candidates for each proj
	float  appr_ratio_;				// approximation ratio

	int    n_cand_;					// number of candidates
	int    *cand_id_;				// candidate data objects id
	float  **cand_data_;			// candidate data objects
	RQALSH *lsh_;					// index of sample data objects

	// -------------------------------------------------------------------------
	void bulkload(					// bulkloading
		const float **data);			// objects after moving to centroid

	// -------------------------------------------------------------------------
	int data_dependent_select(		// data dependent selection
		const float **shift_data);		// shift data objects
};

#endif // __RQALSH_STAR_H
