#include "rqalsh_star.h"

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
	n_pts_ = n;
	dim_   = d;
	L_     = L;
	M_     = M;
	data_  = data;
	lsh_   = NULL;

	// -------------------------------------------------------------------------
	//  get sample data from data dependent selection
	// -------------------------------------------------------------------------
	int n_cand = L_ * M_;
	g_memory += SIZEINT * n_cand;
	cand_ = new int[n_cand];
	data_dependent_select(data, cand_);

	// -------------------------------------------------------------------------
	//  build rqalsh if necessary
	// -------------------------------------------------------------------------
	if (n_cand > N_THRESHOLD) {
		lsh_ = new RQALSH(n_cand, dim_, ratio);

		int m = lsh_->m_;
		for (int i = 0; i < n_cand; ++i) {
			int id = cand_[i];
			for (int j = 0; j < m; ++j) {
				lsh_->tables_[j][i].id_  = i;
				lsh_->tables_[j][i].key_ = lsh_->calc_hash_value(j, data[id]);
			}
		}
		for (int i = 0; i < m; ++i) {
			qsort(lsh_->tables_[i], n_cand, sizeof(Result), ResultComp);
		}
	}
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::data_dependent_select(// data dependent selection
	const float **data,					// data objects
	int   *cand)						// candidate id (return)
{
	// -------------------------------------------------------------------------
	//  calc the shift data
	// -------------------------------------------------------------------------
	int   max_id = -1;
	float max_norm = -1.0f;
	float *norm = new float[n_pts_];
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) shift_data[i] = new float[dim_];

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
	delete[] norm;  norm  = NULL;
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::calc_shift_data(	// calculate shift data objects
	int   &max_id,						// data id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float **shift_data) 				// shift data (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	std::vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data_[i][j];
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
			float tmp = data_[i][j] - centroid[j];
			shift_data[i][j] = tmp;
			norm[i] += SQR(tmp);
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
	centroid.clear(); centroid.shrink_to_fit();
}

// -----------------------------------------------------------------------------
RQALSH_STAR::~RQALSH_STAR()			// destructor
{
	delete[] cand_; cand_ = NULL; 
	g_memory -= SIZEINT * L_ * M_;

	if (lsh_ != NULL) { delete lsh_; lsh_ = NULL; }
}

// -----------------------------------------------------------------------------
void RQALSH_STAR::display()			// display parameters
{
	printf("Parameters of RQALSH*:\n");
	printf("    n    = %d\n",   n_pts_);
	printf("    d    = %d\n",   dim_);
	printf("    L    = %d\n",   L_);
	printf("    M    = %d\n",   M_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH_STAR::kfn(				// c-k-AFN search
	int   top_k,						// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	int n_cand = L_ * M_;
	if (n_cand > N_THRESHOLD) {
		// find candidates by rqalsh
		std::vector<int> cand_list;
		lsh_->kfn(top_k, MINREAL, query, cand_list);

		// check candidates
		for (size_t i = 0; i < cand_list.size(); ++i) {
			int   id   = cand_[cand_list[i]];
			float dist = calc_l2_dist(dim_, data_[id], query);
			list->insert(dist, id + 1);
		}
	}
	else {
		// check all candidates
		for (int i = 0; i < n_cand; ++i) {
			int   id   = cand_[i];
			float dist = calc_l2_dist(dim_, data_[id], query);
			list->insert(dist, id + 1);
		}
	}
	return 0;
}
