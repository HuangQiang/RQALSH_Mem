#include "ml_rqalsh.h"

// -----------------------------------------------------------------------------
ML_RQALSH::ML_RQALSH(				// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float *data)					// data objects
	: n_pts_(n), dim_(d), ratio_(ratio), data_(data)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data obejcts
	// -------------------------------------------------------------------------
	centroid_ = new float[d];
	for (int i = 0; i < d; ++i) centroid_[i] = 0.0f;

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < d; ++j) {
			centroid_[j] += data_[i*d+j];
		}
	}
	for (int i = 0; i < d; ++i) centroid_[i] /= n;

	// -------------------------------------------------------------------------
	//  reorder data objects by their l2-dist to centroid (descending order)
	// -------------------------------------------------------------------------
	Result *arr = new Result[n];
	for (int i = 0; i < n; ++i) {
		arr[i].id_  = i;
		arr[i].key_ = calc_l2_dist(d, &data_[i*d], centroid_);
	}
	qsort(arr, n, sizeof(Result), ResultCompDesc);

	sorted_id_ = new int[n];
	for (int i = 0; i < n; ++i) sorted_id_[i] = arr[i].id_;

	// -------------------------------------------------------------------------
	//  multi-level partition
	// -------------------------------------------------------------------------
	int start = 0;	
	while (start < n) {
		//  get index for each block
		int   idx    = start;
		int   cnt    = 0;
		float radius = arr[start].key_;
		float min_r  = LAMBDA * radius;

		while (idx < n && arr[idx].key_> min_r) {
			++idx;
			if (++cnt >= MAX_BLOCK_NUM) break; 
		}
		// build rqalsh for each block
		const int *index = (const int*) sorted_id_ + start;
		RQALSH *lsh = new RQALSH(cnt, d, ratio, index, data);
		
		// update info
		lsh_.push_back(lsh);
		radius_.push_back(radius);
		start += cnt;
	}
	assert(start == n);
	delete[] arr;
}

// -----------------------------------------------------------------------------
ML_RQALSH::~ML_RQALSH()				// destructor
{
	for (auto lsh : lsh_) { delete lsh; lsh = NULL; }
	lsh_.clear();    lsh_.shrink_to_fit();
	radius_.clear(); radius_.shrink_to_fit();

	delete[] sorted_id_; sorted_id_ = NULL;
	delete[] centroid_;  centroid_  = NULL;
}

// -----------------------------------------------------------------------------
void ML_RQALSH::display()			// display parameters
{
	printf("Parameters of ML_RQALSH:\n");
	printf("    n       = %d\n",   n_pts_);
	printf("    d       = %d\n",   dim_);
	printf("    c       = %.1f\n", ratio_);
	printf("    #blocks = %d\n\n", (int) lsh_.size());
}

// -----------------------------------------------------------------------------
int ML_RQALSH::kfn(					// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k results (return)
{
	float dist2ctr = calc_l2_dist(dim_, centroid_, query);
	float radius = MINREAL;

	int cnt = 0;
	for (int i = 0; i < (int) lsh_.size(); ++i) {
		// early stop pruning
		float ub = radius_[i] + dist2ctr;
		if (radius > ub / ratio_) break;

		// k-FN search by rqalsh on each block
		cnt += lsh_[i]->kfn(top_k, radius, query, list);
		radius = list->min_key();
	}
	return cnt;
}
