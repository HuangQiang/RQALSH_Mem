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
	n_pts_  = n;
	dim_    = d;
	ratio_  = ratio;
	data_   = data;

	// -------------------------------------------------------------------------
	//  calculate the centroid of data obejcts
	// -------------------------------------------------------------------------
	g_memory += SIZEFLOAT * dim_;
	centroid_ = new float[dim_];
	for (int i = 0; i < dim_; ++i) centroid_[i] = 0.0f;

	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid_[j] += data_[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid_[i] /= n_pts_;

	// -------------------------------------------------------------------------
	//  reorder data objects by their l2-dist to centroid (descending order)
	// -------------------------------------------------------------------------
	Result *arr = new Result[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		arr[i].id_ = i;
		arr[i].key_ = calc_l2_dist(dim_, data_[i], centroid_);
	}
	qsort(arr, n_pts_, sizeof(Result), ResultCompDesc);

	g_memory += SIZEINT * n_pts_;
	new_order_id_ = new int[n_pts_];
	for (int i = 0; i < n_pts_; ++i) new_order_id_[i] = arr[i].id_;

	// -------------------------------------------------------------------------
	//  multi-level partition
	// -------------------------------------------------------------------------
	int start = 0;	
	while (start < n_pts_) {
		// ---------------------------------------------------------------------
		//  get new_order_id for one block
		// ---------------------------------------------------------------------
		int   idx = start, cnt = 0;
		float radius = arr[start].key_;
		float min_r  = LAMBDA * radius;

		while (idx < n_pts_ && arr[idx].key_> min_r) {
			++idx;
			if (++cnt >= MAX_BLOCK_NUM) break; 
		}

		// ---------------------------------------------------------------------
		//  init each block (and build rqalsh for it)
		// ---------------------------------------------------------------------
		Block *block = new Block();
		block->n_pts_  = cnt;
		block->radius_ = radius;
		block->index_  = new_order_id_ + start;

		if (cnt > N_THRESHOLD) {
			block->lsh_ = new RQALSH(cnt, dim_, ratio_);

			int m = block->lsh_->m_;
			for (int i = 0; i < cnt; ++i) {
				int id = block->index_[i];
				for (int j = 0; j < m; ++j) {
					float val = block->lsh_->calc_hash_value(j, data[id]);
					block->lsh_->tables_[j][i].id_  = i;
					block->lsh_->tables_[j][i].key_ = val;
				}
			}
			for (int i = 0; i < m; ++i) {
				qsort(block->lsh_->tables_[i], cnt, sizeof(Result), ResultComp);
			}
		}
		blocks_.push_back(block);
		start += cnt;
	}
	num_blocks_ = (int) blocks_.size();
	assert(start == n_pts_);

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] arr; arr = NULL;
}

// -----------------------------------------------------------------------------
ML_RQALSH::~ML_RQALSH()				// destructor
{
	delete[] new_order_id_; new_order_id_ = NULL; g_memory -= SIZEINT * n_pts_;
	delete[] centroid_; centroid_ = NULL; g_memory -= SIZEFLOAT * dim_;

	if (!blocks_.empty()) {
		for (int i = 0; i < num_blocks_; ++i) {
			delete blocks_[i]; blocks_[i] = NULL;
		}
		blocks_.clear(); blocks_.shrink_to_fit();
	}
}

// -----------------------------------------------------------------------------
void ML_RQALSH::display()			// display parameters
{
	printf("Parameters of ML_RQALSH:\n");
	printf("    n          = %d\n",   n_pts_);
	printf("    d          = %d\n",   dim_);
	printf("    c          = %.1f\n", ratio_);
	printf("    num_blocks = %d\n",   num_blocks_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int ML_RQALSH::kfn(					// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// input query
	MaxK_List *list)					// top-k results (return)
{
	float dist2ctr = calc_l2_dist(dim_, centroid_, query);
	float radius = MINREAL;
	std::vector<int> cand;

	for (int i = 0; i < num_blocks_; ++i) {
		Block *block = blocks_[i];

		// ---------------------------------------------------------------------
		//  pruning
		// ---------------------------------------------------------------------
		float ub = block->radius_ + dist2ctr;
		if (radius > ub / ratio_) break;

		// ---------------------------------------------------------------------
		//  c-k-AFN search on each block
		// ---------------------------------------------------------------------
		int n = block->n_pts_;
		if (n > N_THRESHOLD) {
			// find candidates by rqalsh
			cand.clear();
			block->lsh_->kfn(top_k, radius, query, cand);

			// check candidates
			for (size_t j = 0; j < cand.size(); ++j) {
				int   id   = block->index_[cand[j]];
				float dist = calc_l2_dist(dim_, data_[id], query);
				list->insert(dist, id + 1);
			}	
		}
		else {
			// check all candidates
			for (int j = 0; j < n; ++j) {
				int   id   = block->index_[j];
				float dist = calc_l2_dist(dim_, data_[id], query);
				list->insert(dist, id + 1);
			}
		}
		radius = list->min_key();
	}
	return 0;
}
