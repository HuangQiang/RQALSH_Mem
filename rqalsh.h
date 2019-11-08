#ifndef __RQALSH_H
#define __RQALSH_H

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
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search
		int   top_k,					// top-k value
	    const float *query,				// query object
		MaxK_List *list);				// c-k-AFN results (return)

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search (for RQALSH*)
		int   top_k,					// top-k value
	    const float *query,				// query object
		const int *object_id,			// object id mapping
		MaxK_List *list);				// c-k-AFN results (return)

	// -------------------------------------------------------------------------
	int kfn(						// c-k-AFN search (for ML_RQALSH)
		int   top_k,					// top-k value
		float R,						// limited search range
		const float *query,				// input query
		const int *object_id,			// objects id mapping
		MaxK_List *list);				// c-k-AFN results (return)

protected:
	int   n_pts_;					// cardinality
	int   dim_;						// dimensionality
	float beta_;					// false positive percentage
	float delta_;					// error probability
	float appr_ratio_;				// approximation ratio
	const float **data_;			// data objects

	float  w_;						// bucket width
	float  p1_;						// positive probability
	float  p2_;						// negative probability
	float  alpha_;					// collision threshold percentage
	int    m_;						// number of hash tables
	int    l_;						// collision threshold
	float  *a_array_;				// hash functions
	Result **tables_;				// hash tables	

	int    *freq_;					// frequency of data objects
	int    *lpos_;					// left  position of hash table
	int    *rpos_;					// right position of hash table
	bool   *checked_;				// whether the data objects are checked
	bool   *bucket_flag_;			// flag of bucket width
	bool   *range_flag_;			// flag of search range
	float  *q_val_;					// hash value of query				

	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)
	
	// -------------------------------------------------------------------------
	int bulkload();					// build hash tables	

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int   tid,						// hash table id
		const float *data);				// one data object

	// -------------------------------------------------------------------------
	float find_radius();			// find proper radius
};

#endif // __RQALSH_H
