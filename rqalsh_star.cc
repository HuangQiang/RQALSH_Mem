#include "headers.h"

// -----------------------------------------------------------------------------
RQALSH_Star::RQALSH_Star()			// constructor
{
	n_pts_       = -1;
	dim_         = -1;
	L_           = -1;
	M_           = -1;
	beta_        = -1;
	delta_       = -1.0f;
	appr_ratio_  = -1.0f;
	sample_size_ = -1;
	sample_id_   = NULL;
	sample_data_ = NULL;
	lsh_         = NULL;
}

// -----------------------------------------------------------------------------
RQALSH_Star::~RQALSH_Star()			// destructor
{
	if (sample_id_ != NULL) {
		delete[] sample_id_; sample_id_ = NULL;
	}
	if (sample_data_ != NULL) {
		for (int i = 0; i < sample_size_; ++i) {
			delete[] sample_data_[i]; sample_data_[i] = NULL;
		}
		delete[] sample_data_; sample_data_ = NULL;
	}
	if (lsh_ != NULL) {
		delete lsh_; lsh_ = NULL;
	}
}

// -----------------------------------------------------------------------------
int RQALSH_Star::build(				// build index	
	int   n,							// cardinality
	int   d,							// dimensionality
	int   L,							// number of projection
	int   M,							// number of candidates
	int   beta,							// false positive percentage
	float delta,						// error probability
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
	beta_       = beta;
	delta_      = delta;
	appr_ratio_ = ratio;

	// -------------------------------------------------------------------------
	//  build hash tables (bulkloading)
	// -------------------------------------------------------------------------
	bulkload(data);
	// display();
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::bulkload(			// bulkloading for each block
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  calculate shift data
	// -------------------------------------------------------------------------
	float *shift_data = new float[n_pts_ * dim_];
	calc_shift_data(data, shift_data);

	// -------------------------------------------------------------------------
	//  get sample data from data dependent selection
	// -------------------------------------------------------------------------
	sample_size_ = L_ * M_;
	sample_id_   = new int[sample_size_];
	data_dependent_select(shift_data);

	sample_data_ = new float*[sample_size_];
	for (int i = 0; i < sample_size_; ++i) {
		int id = sample_id_[i];
		
		sample_data_[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			sample_data_[i][j] = data[id][j];
		}
	}

	// -------------------------------------------------------------------------
	//  build hash tables for objects from drusilla select using RQALSH
	// -------------------------------------------------------------------------
	if (sample_size_ > CANDIDATES) {
		lsh_ = new RQALSH();
		lsh_->build(sample_size_, dim_, beta_, delta_, appr_ratio_, 
			(const float **) sample_data_);
	}
	
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] shift_data; shift_data = NULL;
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::calc_shift_data( 	// calc shift data
	const float **data,					// data objects
	float *shift_data)					// shift data objects (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) {
		centroid[i] /= (float) n_pts_;
	}

	// -------------------------------------------------------------------------
	//  make a copy of data objects which move to the centroid of data objects
	// -------------------------------------------------------------------------
	for (int i = 0; i < n_pts_; ++i) {
		int base = i * dim_;
		for (int j = 0; j < dim_; ++j) {
			shift_data[base + j] = data[i][j] - centroid[j];
		}
	}
	
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH_Star::data_dependent_select( // drusilla select
	const float *shift_data)			// shift data
{
	// -------------------------------------------------------------------------
	//  calc the norm of data objects and find the data object with max norm
	// -------------------------------------------------------------------------
	int   max_id   = -1;
	float max_norm = -1.0f;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		int base = i * dim_;
		for (int j = 0; j < dim_; ++j) {
			float x = shift_data[base + j];
			norm[i] += x * x;
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id   = i;
		}
	}

	vector<bool>  close_angle(n_pts_);
	vector<float> projection(dim_);
	Result *score_pair = new Result[n_pts_];

	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			projection[j] = shift_data[max_id * dim_ + j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			int base = j * dim_;

			if (norm[j] >= 0.0f) {
				float offset = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					offset += (shift_data[base + k] * projection[k]);
				}

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					float x = shift_data[base + k] - offset * projection[k];
					distortion += x * x;
				}

				score_pair[j].id_ = j;
				score_pair[j].key_ = offset * offset - distortion;
			}
			else {
				score_pair[j].id_ = j;
				score_pair[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this projection
		// ---------------------------------------------------------------------
		qsort(score_pair, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < M_; ++j) {
			int id = score_pair[j].id_;

			sample_id_[i * M_ + j] = id;
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

	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] score_pair; score_pair = NULL;

	return 0;
}

// -----------------------------------------------------------------------------
void RQALSH_Star::display()			// display parameters
{
	printf("Parameters of RQALSH_Star:\n");
	printf("    n           = %d\n", n_pts_);
	printf("    d           = %d\n", dim_);
	printf("    L           = %d\n", L_);
	printf("    M           = %d\n", M_);
	printf("    beta        = %d\n", beta_);
	printf("    delta       = %f\n", delta_);
	printf("    c           = %f\n", appr_ratio_);
	printf("    sample_size = %d\n", sample_size_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH_Star::kfn(				// c-k-AFN search
	int top_k,							// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	// -------------------------------------------------------------------------
	//  use index to speed up c-k-AFN search
	// -------------------------------------------------------------------------
	int candidates = CANDIDATES + top_k - 1;
	if (sample_size_ > candidates) {
		return lsh_->kfn(top_k, query, (const int*) sample_id_, list);
	}

	// -------------------------------------------------------------------------
	//  if the number of samples is small enough, linear scan directly
	// -------------------------------------------------------------------------
	for (int i = 0; i < sample_size_; ++i) {
		int   id   = sample_id_[i];
		float dist = calc_l2_dist(dim_, sample_data_[i], query);

		list->insert(dist, id + 1);
	}
	
	return 0;
}
