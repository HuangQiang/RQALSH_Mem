#include "drusilla_select.h"

// -----------------------------------------------------------------------------
Drusilla_Select::Drusilla_Select(	// default constructor
	int   n,							// number of data idects
	int   d,							// number of dimensions
	int   l,							// number of projections
	int   m,							// number of candidates on each proj
	const float *data)					// data idects
	: n_pts_(n), dim_(d), l_(l), m_(m), data_(data)
{
	// -------------------------------------------------------------------------
	//  calc the shift data
	// -------------------------------------------------------------------------
	int   max_id      = -1;
	float max_norm    = -1.0f;
	float *norm       = new float[n];
	float *shift_data = new float[n*d];

	calc_shift_data(max_id, max_norm, norm, shift_data);

	// -------------------------------------------------------------------------
	//  drusilla-select
	// -------------------------------------------------------------------------
	float  *proj  = new float[d];
	Result *score = new Result[n];
	bool   *close_angle = new bool[n];

	cand_ = new int[l * m];
	for (int i = 0; i < l; ++i) {
		// ---------------------------------------------------------------------
		//  select the projection vector with largest norm and normalize it
		// ---------------------------------------------------------------------
		for (int j = 0; j < d; ++j) {
			proj[j] = shift_data[max_id*d+j] / norm[max_id];
		}

		// ---------------------------------------------------------------------
		//  calculate offsets and distortions
		// ---------------------------------------------------------------------
		for (int j = 0; j < n; ++j) {
			score[j].id_ = j;
			close_angle[j] = false;

			if (norm[j] > 0.0f) {
				float offset = calc_inner_product(d, &shift_data[j*d], proj);

				float distortion = 0.0f;
				for (int k = 0; k < dim_; ++k) {
					distortion += SQR(shift_data[j*d+k] - offset * proj[k]);
				}
				distortion = sqrt(distortion);

				score[j].key_ = fabs(offset) - fabs(distortion);
				if (atan(distortion / fabs(offset)) < ANGLE) {
					close_angle[j] = true;
				}
			}
			else if (fabs(norm[j]) < CHECK_ERROR) {
				score[j].key_ = MINREAL + 1.0f;
			}
			else {
				score[j].key_ = MINREAL;
			}
		}

		// ---------------------------------------------------------------------
		//  collect the objects that are well-represented by this proj
		// ---------------------------------------------------------------------
		qsort(score, n, sizeof(Result), ResultCompDesc);
		for (int j = 0; j < m; ++j) {
			int id = score[j].id_;

			cand_[i*m+j] = id;
			norm[id] = -1.0f;
		}

		// ---------------------------------------------------------------------
		//  find the next largest norm and the corresponding idect
		// ---------------------------------------------------------------------
		max_id = -1;
		max_norm = -1.0f;
		for (int j = 0; j < n; ++j) {
			if (norm[j] > 0.0f && close_angle[j]) { norm[j] = 0.0f; }
			if (norm[j] > max_norm) { max_norm = norm[j]; max_id = j; }
		}
	}
	// -------------------------------------------------------------------------
	//  release space
	// -------------------------------------------------------------------------
	delete[] norm;
	delete[] close_angle;
	delete[] proj;
	delete[] score;
	delete[] shift_data;
}

// -----------------------------------------------------------------------------
void Drusilla_Select::calc_shift_data( // calculate shift data objects
	int   &max_id,						// data id with max l2-norm (return)
	float &max_norm,					// max l2-norm (return)
	float *norm,						// l2-norm of shift data (return)
	float *shift_data) 					// shift data (return)
{
	// -------------------------------------------------------------------------
	//  calculate the centroid of data objects
	// -------------------------------------------------------------------------
	float *centroid = new float[dim_];
	memset(centroid, 0.0f, dim_ * SIZEFLOAT);
	for (int i = 0; i < n_pts_; ++i) {
		for (int j = 0; j < dim_; ++j) {
			centroid[j] += data_[i*dim_+j];
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
			norm[i] += SQR(tmp);
		}
		norm[i] = sqrt(norm[i]);

		if (norm[i] > max_norm) { max_norm = norm[i]; max_id = i; }
	}
	delete[] centroid;
}

// -----------------------------------------------------------------------------
Drusilla_Select::~Drusilla_Select()
{
	delete[] cand_; cand_ = NULL;
}

// -----------------------------------------------------------------------------
void Drusilla_Select::display()		// display the parameters
{
	printf("Parameters of Drusilla-Select (SISAP2016 paper):\n");
	printf("    n = %d\n",   n_pts_);
	printf("    d = %d\n",   dim_);
	printf("    l = %d\n",   l_);
	printf("    m = %d\n\n", m_);
}

// -----------------------------------------------------------------------------
int Drusilla_Select::kfn(			// c-k-AFN search
	const float *query,					// query object
	MaxK_List   *list)					// top-k results (return)
{
	int size = l_ * m_;
	for (int i = 0; i < size; ++i) {
		int id = cand_[i];
		float dist = calc_l2_dist(dim_, query, &data_[id*dim_]);
		list->insert(dist, id + 1);
	}
	return size;
}
