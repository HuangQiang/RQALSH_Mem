#include "headers.h"

// -----------------------------------------------------------------------------
RQALSH_STAR::RQALSH_STAR(			// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	int   L,							// number of proj
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	n_pts_      = n;
	dim_        = d;
	L_          = L;
	M_          = M;
	appr_ratio_ = ratio;
	n_cand_     = L_ * M_;
	lsh_        = NULL;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload(data);
}

// -----------------------------------------------------------------------------
RQALSH_STAR::~RQALSH_STAR()			// destructor
{
	delete[] cand_id_; cand_id_ = NULL;
	for (int i = 0; i < n_cand_; ++i) {
		delete[] cand_data_[i]; cand_data_[i] = NULL;
	}
	delete[] cand_data_; cand_data_ = NULL;

	if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::bulkload(			// bulkloading for each block
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  calculate shift data
	// -------------------------------------------------------------------------
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) shift_data[i] = new float[dim_];
	
	vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid[i] /= (float) n_pts_;

	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			shift_data[i][j] = data[i][j] - centroid[j];
		}
	}

	// -------------------------------------------------------------------------
	//  get sample data from data dependent selection
	// -------------------------------------------------------------------------
	cand_id_ = new int[n_cand_];
	data_dependent_select((const float **) shift_data);

	cand_data_ = new float*[n_cand_];
	for (int i = 0; i < n_cand_; ++i) {
		int id = cand_id_[i];
		
		cand_data_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			cand_data_[i][j] = data[id][j];
		}
	}

	// -------------------------------------------------------------------------
	//  build hash tables for objects from drusilla select using RQALSH
	// -------------------------------------------------------------------------
	if (n_cand_ > N_THRESHOLD) {
		lsh_ = new RQALSH(n_cand_, dim_, appr_ratio_, (const float **) cand_data_);
	}
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::data_dependent_select( // drusilla select
	const float **shift_data)			// shift data
{
	// -------------------------------------------------------------------------
	//  calc the norm of data objects and find the data object with max norm
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = -1.0f;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, shift_data[i], shift_data[i]));
		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id   = i;
		}
	}

	float  *proj  = new float[dim_];
	Result *score = new Result[n_pts_];
	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			proj[j] = shift_data[max_id][j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] >= 0.0f) {
				float offset = calc_inner_product(dim_, shift_data[j], proj);

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j][k] - offset * proj[k]);
				}
				score[j].id_  = j;
				score[j].key_ = offset * offset - distortion;
			}
			else {
				score[j].id_  = j;
				score[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < M_; ++j) {
			int id = score[j].id_;

			cand_id_[i * M_ + j] = id;
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > max_norm) {
				max_norm = norm[j];
				max_id = j;
			}
		}
	}
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::display()			// display parameters
{
	printf("Parameters of RQALSH*:\n");
	printf("    n    = %d\n",   n_pts_);
	printf("    d    = %d\n",   dim_);
	printf("    L    = %d\n",   L_);
	printf("    M    = %d\n",   M_);
	printf("    c    = %.1f\n", appr_ratio_);
	printf("    size = %d\n\n", n_cand_);
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::kfn(				// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	if (n_cand_ > N_THRESHOLD) {
		lsh_->kfn(top_k, query, (const int*) cand_id_, list);
	} else {
		for (int i = 0; i < n_cand_; ++i) {
			int   id   = cand_id_[i];
			float dist = calc_l2_dist(dim_, cand_data_[i], query);
			list->insert(dist, id + 1);
		}
	}
	return 0;
}
