#ifndef __RQALSH_H
#define __RQALSH_H

class  MaxK_List;
struct Result;

// -----------------------------------------------------------------------------
//  RQALSH: an LSH scheme for for high-dimensional c-k-AFN search
// -----------------------------------------------------------------------------
class RQALSH {
public:
	RQALSH();						// default constructor
	~RQALSH();						// destructor

	// -------------------------------------------------------------------------
	int build(						// build index
		int   n,						// cardinality
		int   d,						// dimensionality
		int   beta,						// false positive percentage
		float delta,					// error probability
		float ratio,					// approximation ratio
		const float **data);			// data objects

	// -------------------------------------------------------------------------
	int kfn(						// k-FN search
		int top_k,						// top-k value
	    const float *query,				// query object
		MaxK_List *list);				// k-FN results (return)

	// -------------------------------------------------------------------------
	int kfn(						// k-FN search
		int top_k,						// top-k value
	    const float *query,				// query object
		const int *object_id,			// object id mapping
		MaxK_List *list);				// k-FN results (return)

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
	bool   *flag_;					// flag of bucket width
	float  *q_val_;					// hash value of query				

	// -------------------------------------------------------------------------
	void calc_params();				// calc parameters of qalsh

	// -------------------------------------------------------------------------
	float calc_l2_prob(				// calc <p1> and <p2> for L_{2.0} distance
		float x);						// x = w / (2.0 * r)

	// -------------------------------------------------------------------------
	void gen_hash_func();			// generate hash functions
	
	// -------------------------------------------------------------------------
	int bulkload();					// build hash tables	

	// -------------------------------------------------------------------------
	float calc_hash_value(			// calc hash value
		int table_id,					// hash table id
		const float *data);				// one data object

	// -------------------------------------------------------------------------
	void display();					// display parameters

	// -------------------------------------------------------------------------
	float find_radius(				// find proper radius
		const float *q_val,				// hash value of query
		const int *lpos,				// left position of hash table
		const int *rpos);				// right position of hash table
};

#endif // __RQALSH_H
