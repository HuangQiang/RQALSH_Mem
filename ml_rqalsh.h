#ifndef __ML_RQALSH_H
#define __ML_RQALSH_H

class  RQALSH;
class  MinK_List;
struct Result;

// -----------------------------------------------------------------------------
//  Blocks: an block which stores hash tables for some of data objects
// -----------------------------------------------------------------------------
struct Blocks {
	int    n_pts_;					// number of data objects
	float  radius_;					// radius of this block
	RQALSH *lsh_;					// index of data objects in this blocks
	vector<int> index_;				// data object id

	Blocks() { n_pts_ = -1; radius_ = -1.0f; lsh_ = NULL; }
	~Blocks() { if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; } }
};

// -----------------------------------------------------------------------------
//  ML_RQALSH: Multi-Level RQALSH for c-k-AFN Search 
// -----------------------------------------------------------------------------
class ML_RQALSH {
public:
	ML_RQALSH(						// constructor
		int   n,						// cardinality
		int   d,						// dimensionality
		float ratio,					// approximation ratio
    	const float **data);			// data objects

	// -------------------------------------------------------------------------
	~ML_RQALSH();					// destructor
	
	// -------------------------------------------------------------------------
	void display();			        // display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN seach	
		int top_k,		    			// top-k value
		const float *query,				// input query
		MaxK_List *list);				// top-k results (return)

protected:
	int   n_pts_;					// number of data objects
	int   dim_;						// dimensionality
	float appr_ratio_;				// approximation ratio
	const float **data_;					// data objects (modified)

	int   num_blocks_;				// number of blocks
	float *centroid_;				// centroid of data objects
	float **shift_data_;			// shift data (move to centroid)
	vector<Blocks*> blocks_;		// blocks
    
    // -------------------------------------------------------------------------
	void bulkload();				// build hash tables
};

#endif // __ML_RQALSH_H