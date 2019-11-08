#include <algorithm>
#include <cassert>

#include "def.h"
#include "util.h"
#include "pri_queue.h"
#include "rqalsh.h"
#include "ml_rqalsh.h"


// -----------------------------------------------------------------------------
ML_RQALSH::ML_RQALSH(				// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_       = n;
	dim_         = d;
	appr_ratio_  = ratio;
	data_        = data;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
ML_RQALSH::~ML_RQALSH()				// destructor
{
	if (!blocks_.empty()) {
		for (int i = 0; i < num_blocks_; ++i) {
			delete blocks_[i]; blocks_[i] = NULL;
		}
		blocks_.clear(); blocks_.shrink_to_fit();
	}

	delete[] centroid_; centroid_ = NULL;
	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data_[i]; shift_data_[i] = NULL;
	}
	delete[] shift_data_; shift_data_ = NULL;
}

// -----------------------------------------------------------------------------
void ML_RQALSH::bulkload()			// build hash tables
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data obejcts
	// -------------------------------------------------------------------------
	centroid_ = new float[dim_];
	for (int i = 0; i < dim_; ++i) centroid_[i] = 0.0f;

	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid_[j] += data_[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid_[i] /= (float) n_pts_;

	// -------------------------------------------------------------------------
	//  calculate the l2-dist after moving all data objects to their centroid
	// -------------------------------------------------------------------------
	Result *arr = new Result[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		arr[i].id_ = i;
		arr[i].key_ = calc_l2_dist(dim_, data_[i], centroid_);
	}

	// -------------------------------------------------------------------------
	//  reorder the shifted data objects by their l2-dist (descending order)
	// -------------------------------------------------------------------------
	qsort(arr, n_pts_, sizeof(Result), ResultCompDesc);

	shift_data_ = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		int id = arr[i].id_;

		shift_data_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			shift_data_[i][j] = data_[id][j] - centroid_[j];
		}
	}

	// -------------------------------------------------------------------------
	//  divide datasets into blocks and build hash tables for each block
	// -------------------------------------------------------------------------
	int start = 0;	
	while (start < n_pts_) {
		// ---------------------------------------------------------------------
		//  divide one block
		// ---------------------------------------------------------------------
		Blocks *block = new Blocks();

		int   idx       = start;
		int   count     = 0;
		float radius    = arr[start].key_;
		float threshold = LAMBDA * radius;

		while (idx<n_pts_ && arr[idx].key_>threshold && count<MAX_BLOCK_NUM) {
			block->index_.push_back(arr[idx].id_);
			++count;
			++idx;
		}

		// ---------------------------------------------------------------------
		//  update info
		// ---------------------------------------------------------------------
		block->n_pts_  = count;
		block->radius_ = radius;

		if (count > N_THRESHOLD) {
			block->lsh_ = new RQALSH(count, dim_, appr_ratio_, 
				(const float **) shift_data_ + start);
		}
		blocks_.push_back(block);
		start += count;
	}
	// -------------------------------------------------------------------------
	//  get the number of blocks
	// -------------------------------------------------------------------------
	num_blocks_ = (int) blocks_.size();
	assert(start == n_pts_);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] arr; arr = NULL;
}

// -----------------------------------------------------------------------------
void ML_RQALSH::display()			// display parameters
{
	printf("Parameters of ML_RQALSH:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    c          = %.1f\n", appr_ratio_);
	printf("    num_blocks = %d\n\n", num_blocks_);
}

// -----------------------------------------------------------------------------
int ML_RQALSH::kfn(					// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k results (return)
{
	// -------------------------------------------------------------------------
	//  calculate the Euclidean arr of query
	// -------------------------------------------------------------------------
	float dist2ctr_q = 0.0f;		// l2-dist to centroid for query
	float *ctr_q = new float[dim_];	// centroid to query
	for (int i = 0; i < dim_; ++i) {
		ctr_q[i] = query[i] - centroid_[i];
		dist2ctr_q += ctr_q[i] * ctr_q[i];
	}
	dist2ctr_q = sqrt(dist2ctr_q);

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1;
	float radius = MINREAL;

	for (int i = 0; i < num_blocks_; ++i) {
		// ---------------------------------------------------------------------
		//  pruning
		// ---------------------------------------------------------------------
		float ub = blocks_[i]->radius_ + dist2ctr_q;
		if (radius > ub / appr_ratio_) break;

		// ---------------------------------------------------------------------
		//  c-k-AFN search on each block
		// ---------------------------------------------------------------------
		int n = blocks_[i]->n_pts_;
		if (n > N_THRESHOLD) {
			// -----------------------------------------------------------------
			//  larger than <candidates>, use RQALSH for c-k-AFN search
			// -----------------------------------------------------------------	
			blocks_[i]->lsh_->kfn(top_k, radius, (const float *) ctr_q, 
				blocks_[i]->index_.data(), list);
		}
		else {
			// -----------------------------------------------------------------
			//  otherwise, use linear scan for c-k-AFN search
			// -----------------------------------------------------------------
			for (int j = 0; j < n; ++j) {
				int   id   = blocks_[i]->index_[j];
				float dist = calc_l2_dist(dim_, query, data_[id]);
				list->insert(dist, id + 1);
			}
		}
		radius = list->min_key();
	}

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] ctr_q; ctr_q = NULL;

	return 0;
}
