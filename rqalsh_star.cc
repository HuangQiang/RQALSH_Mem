#include "rqalsh_star.h"

// -----------------------------------------------------------------------------
RQALSH_STAR::RQALSH_STAR(			// constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	int   L,							// number of proj
	int   M,							// number of candidates
	float ratio,						// approximation ratio
	const float *data)					// data objects
	: n_pts_(n), dim_(d), L_(L), M_(M), data_(data), lsh_(NULL)
{
	// get candidates from data dependent selection
	int n_cand = L * M;
	cand_ = new int[n_cand];
	data_dependent_select(data, cand_);

	//  build rqalsh if necessary
	lsh_ = new RQALSH(n_cand, d, ratio, (const int*) cand_, data);
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::data_dependent_select(// data dependent selection
	const float *data,					// data objects
	int   *cand)						// candidate id (return)
{
	// -------------------------------------------------------------------------
	//  calc the shift data
	// -------------------------------------------------------------------------
	int   max_id      = -1;
	float max_norm    = -1.0f;
	float *norm       = new float[n_pts_];
	float *shift_data = new float[n_pts_ * dim_];

	calc_shift_data(max_id, max_norm, norm, shift_data);

	// -------------------------------------------------------------------------
	//  data dependent selection
	// -------------------------------------------------------------------------
	float  *proj  = new float[dim_];
	Result *score = new Result[n_pts_];
	for (int i = 0; i < L_; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < dim_; ++j) {
			proj[j] = shift_data[max_id*dim_+j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] >= 0.0f) {
				float offset = calc_inner_product(dim_, &shift_data[j*dim_], proj);

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j*dim_+k] - offset * proj[k]);
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

			cand[i * M_ + j] = id;
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding object
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > max_norm) { max_norm = norm[j]; max_id = j; }
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;
	delete[] proj;
	delete[] score;
	delete[] shift_data;
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::calc_shift_data(	// calculate shift data objects
	int   &max_id,						// data id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float *shift_data) 					// shift data (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	std::vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data_[i*dim_ + j];
		}
	}
	for (int i = 0; i < dim_; ++i) centroid[i] /= n_pts_;

	// -------------------------------------------------------------------------
	//  calc shift data and their l2-norm and find max l2-norm and its id
	// -------------------------------------------------------------------------
	max_id   = -1;
	max_norm = MINREAL;

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = 0.0f;
		for (int j = 0; j < dim_; ++j) {
			float tmp = data_[i*dim_+j] - centroid[j];
			
			shift_data[i*dim_+j] = tmp;
			norm[i] += tmp * tmp;
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
	centroid.clear(); centroid.shrink_to_fit();
}

// -----------------------------------------------------------------------------
RQALSH_STAR::~RQALSH_STAR()			// destructor
{
	delete   lsh_;  lsh_  = NULL;
	delete[] cand_; cand_ = NULL; 
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::display()			// display parameters
{
	printf("Parameters of RQALSH*:\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    L = %d\n",   L_);
	printf("    M = %d\n\n", M_);

	lsh_->display();
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::kfn(				// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	return lsh_->kfn(top_k, MINREAL, query, list);
}
