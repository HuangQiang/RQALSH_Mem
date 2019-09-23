#include "headers.h"

// -----------------------------------------------------------------------------
Drusilla_Select::Drusilla_Select(	// default constructor
	int   n,							// number of data idects
	int   d,							// number of dimensions
	int   l,							// number of projections
	int   m,							// number of candidates on each proj
	const float **data)					// data idects
{
	// -------------------------------------------------------------------------
	//  init the parameters
	// -------------------------------------------------------------------------
	n_pts_ = n;
	dim_   = d;
	l_     = l;
	m_     = m;
	data_  = data;

	// -------------------------------------------------------------------------
	//  build hash tables
	// -------------------------------------------------------------------------
	bulkload();
}

// -----------------------------------------------------------------------------
Drusilla_Select::~Drusilla_Select()
{
	delete[] cand_; cand_ = NULL;
}

// -----------------------------------------------------------------------------
void Drusilla_Select::bulkload()		// build hash tables
{
	// -------------------------------------------------------------------------
	//  calculate centroid
	// -------------------------------------------------------------------------
	vector<float> centroid(dim_, 0.0f);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data_[i][j];
		}
	}
	for (int i = 0; i < dim_; ++i) {
		centroid[i] /= (float) n_pts_;
	}

	// -------------------------------------------------------------------------
	//  calc the shift data
	// -------------------------------------------------------------------------
	float **shift_data = new float*[n_pts_];
	for (int i = 0; i < n_pts_; ++i) {
		shift_data[i] = new float[dim_];
		for (int j = 0; j < dim_; ++j) {
			shift_data[i][j] = data_[i][j] - centroid[j];
		}
	}

	// -------------------------------------------------------------------------
	//  find the idect with maximum Euclidean norm
	// -------------------------------------------------------------------------
	int   max_id = -1;
	float max_norm = -1.0f;
	vector<float> norm(n_pts_, 0.0f);

	for (int i = 0; i < n_pts_; ++i) {
		norm[i] = sqrt(calc_inner_product(dim_, shift_data[i], shift_data[i]));
		if (norm[i] > max_norm) {
			max_norm = norm[i];
			max_id = i;
		}
	}

	vector<bool>  close_angle(n_pts_, false);
	Result *score = new Result[n_pts_];
	float  *proj  = new float[dim_];

	cand_ = new int[l_ * m_];
	for (int i = 0; i < l_; ++i) {
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
			score[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = calc_inner_product(dim_, shift_data[j], proj);

				float distortion = 0.0F;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j][k] - offset * proj[k]);
				}
				distortion = sqrt(distortion);

				score[j].key_ = fabs(offset) - fabs(distortion);
				if (atan(distortion / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < FLOATZERO) {
				score[j].key_ = MINREAL + 1.0f;
			}
			else {
				score[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the idects that are well-represented by this proj
		// ---------------------------------------------------------------------
		qsort(score, n_pts_, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < m_; ++j) {
			int id = score[j].id_;
			cand_[i * m_ + j] = id;
			
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding idect
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n_pts_; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) {
				norm[j] = 0.0f;
			}
			if (norm[j] > max_norm) {
				max_norm = norm[j];
				max_id = j;
			}
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] proj;  proj  = NULL;
	delete[] score; score = NULL;

	for (int i = 0; i < n_pts_; ++i) {
		delete[] shift_data[i]; shift_data[i] = NULL;
	}
	delete[] shift_data; shift_data = NULL;
}

// -----------------------------------------------------------------------------
void Drusilla_Select::display()		// display the parameters
{
	printf("Parameters of Drusilla-Select (SISAP2016 paper):\n");
	printf("    n     = %d\n",   n_pts_);
	printf("    d     = %d\n",   dim_);
	printf("    l     = %d\n",   l_);
	printf("    m     = %d\n\n", m_);
}

// -----------------------------------------------------------------------------
int Drusilla_Select::kfn(			// c-k-AFN search
	const float *query,					// query object
	MaxK_List   *list)					// top-k results (return)
{
	int size = l_ * m_;
	for (int i = 0; i < size; ++i) {
		int id = cand_[i];
		float dist = calc_l2_dist(dim_, query, data_[id]);

		list->insert(dist, id + 1);
	}
	return 0;
}