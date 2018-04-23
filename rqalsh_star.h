#ifndef __RQALSH_STAR_H
#define __RQALSH_STAR_H

class RQALSH;
class MaxK_List;

// -----------------------------------------------------------------------------
//  RQALSH_Star is used to solve the problem of c-k-Approximate Furthest 
//  Neighbor (c-k-AFN) search
// -----------------------------------------------------------------------------
class RQALSH_Star {
public:
	RQALSH_Star();
	~RQALSH_Star();

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// cardinality
		int   d,						// dimensionality
		int   L,						// number of projection
		int   M,						// number of candidates
		int   beta,						// false positive percentage
		float delta,					// error probability
		float ratio,					// approximation ratio
		const float **data);			// data objects

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
	int    beta_;					// false positive percentage
	float  delta_;					// error probability
	float  appr_ratio_;				// approximation ratio

	int    sample_size_;			// number of sample data objects
	int    *sample_id_;				// sample data objects id
	float  **sample_data_;			// sample data objects
	RQALSH *lsh_;					// index of sample data objects

	// -------------------------------------------------------------------------
	int bulkload(					// bulkloading
		const float **data);			// objects after moving to centroid

	// -------------------------------------------------------------------------
	int calc_shift_data(			// calc shift data
		const float **data,  			// data objects
		float *shift_data);  			// shift data objects (return)

	// -------------------------------------------------------------------------
	int data_dependent_select(		// data dependent selection
		const float *shift_data);		// shift data objects

	// -------------------------------------------------------------------------
	void display();			        // display parameters
};

#endif // __RQALSH_STAR_H
