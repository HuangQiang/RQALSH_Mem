#include "headers.h"

// -----------------------------------------------------------------------------
RQALSH::RQALSH(						// default constructor
	int   n,							// cardinality
	int   d,							// dimensionality
	float ratio,						// approximation ratio
	const float **data)					// data objects
{
	// -------------------------------------------------------------------------
	//  init the input parameters
	// -------------------------------------------------------------------------
	assert(n >= CANDIDATES);
	
	n_pts_      = n;
	dim_        = d;
	beta_       = (float) CANDIDATES / (float) n;
	delta_      = 0.49f;
	appr_ratio_ = ratio;
	data_       = data;
	
	// -------------------------------------------------------------------------
	//  init <w_> <p1_> and <p2_> (auto tuning-w)
	// -------------------------------------------------------------------------
	w_  = sqrt((8.0f * log(appr_ratio_)) / (appr_ratio_ * appr_ratio_ - 1.0f));
	
	p1_ = calc_l2_prob(w_ / 2.0f);
	p2_ = calc_l2_prob(w_ * appr_ratio_ / 2.0f);

	// -------------------------------------------------------------------------
	//  init <alpha_> <m_> and <l_>
	// -------------------------------------------------------------------------
	float para1 = sqrt(log(2.0f / beta_));
	float para2 = sqrt(log(1.0f / delta_));
	float para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);

	float eta = para1 / para2;
	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);

	m_ = (int) ceil((para1 + para2) * (para1 + para2) / para3);
	l_ = (int) ceil(alpha_ * m_);

	// -------------------------------------------------------------------------
	//  calc parameters and generate hash functions
	// -------------------------------------------------------------------------
	int size = m_ * dim_;
	a_array_ = new float[size];
	for (int i = 0; i < size; i++) {
		a_array_[i] = gaussian(0.0f, 1.0f);
	}
	
	// -------------------------------------------------------------------------
	//  bulkloading
	// -------------------------------------------------------------------------
	freq_        = new int[n_pts_];
	lpos_        = new int[m_];
	rpos_        = new int[m_];
	checked_     = new bool[n_pts_];
	bucket_flag_ = new bool[m_];
	range_flag_  = new bool[m_];
	q_val_       = new float[m_];
	tables_      = new Result*[m_];

	for (int i = 0; i < m_; ++i) {
		tables_[i] = new Result[n_pts_];
		for (int j = 0; j < n_pts_; ++j) {
			tables_[i][j].id_ = j;
			tables_[i][j].key_ = calc_hash_value(i, data_[j]);
		}
		qsort(tables_[i], n_pts_, sizeof(Result), ResultComp);
	}
}

// -----------------------------------------------------------------------------
RQALSH::~RQALSH()					// destructor
{
	delete[] a_array_;     a_array_     = NULL;
	delete[] freq_;        freq_        = NULL;
	delete[] lpos_;        lpos_        = NULL;
	delete[] rpos_;        rpos_        = NULL;
	delete[] checked_;     checked_     = NULL;
	delete[] bucket_flag_; bucket_flag_ = NULL;
	delete[] range_flag_;  range_flag_  = NULL;
	delete[] q_val_;       q_val_       = NULL;

	for (int i = 0; i < m_; ++i) {
		delete[] tables_[i]; tables_[i] = NULL;
	}
	delete[] tables_; tables_ = NULL;
}

// -------------------------------------------------------------------------
inline float RQALSH::calc_l2_prob(	// calc <p1> and <p2> for L2 distance
	float x)							// x = w / (2.0 * r)
{
	return 1.0f - new_gaussian_prob(x);
}

// -------------------------------------------------------------------------
inline float RQALSH::calc_hash_value( // calc hash value
	int   tid,							// hash table id
	const float *data)					// one data object
{
	int   base = tid * dim_;
	float ret  = 0.0f;
	for (int i = 0; i < dim_; ++i) {
		ret += a_array_[base + i] * data[i];
	}
	return ret;
}

// -------------------------------------------------------------------------
void RQALSH::display()				// display the parameters
{
	printf("Parameters of RQALSH:\n");
	printf("    n          = %d\n", n_pts_);
	printf("    d          = %d\n", dim_);
	printf("    ratio      = %f\n", appr_ratio_);
	printf("    w          = %f\n", w_);
	printf("    p1         = %f\n", p1_);
	printf("    p2         = %f\n", p2_);
	printf("    alpha      = %f\n", alpha_);
	printf("    beta       = %f\n", beta_);
	printf("    delta      = %f\n", delta_);
	printf("    m          = %d\n", m_);
	printf("    l          = %d\n", l_);
	printf("\n");
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search
	int top_k,							// top-k value
	const float *query,					// query object
	MaxK_List *list)					// k-FN results (return)
{
	int    candidates = CANDIDATES + top_k - 1; // candidate size
	float  kdist      = MINREAL;
	Result *table     = NULL;

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		lpos_[i]  = 0;  
		rpos_[i]  = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   dist_cnt = 0;			// counter of Eucldiean distance compuation
	float radius   = find_radius(); // search radius
	float bucket   = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		int cnt = -1, lpos = -1, rpos = -1, id = -1;
		while (num_bucket < m_) {
			float ldist = -1.0f;	// left  proj dist with query
			float rdist = -1.0f;	// right proj dist with query
			float q_val = -1.0f;	// hash value of 
			float dist  = -1.0f;	// l2 dist

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				q_val = q_val_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				lpos = lpos_[j]; rpos = rpos_[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos < rpos) {
						ldist = fabs(q_val - table[lpos].key_);
					}
					else break;
					if (ldist < bucket) break;

					id = table[lpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, id + 1);

						if (++dist_cnt >= candidates) break;
					}
					++lpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				lpos_[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos < rpos) {
						rdist = fabs(q_val - table[rpos].key_);
					}
					else break;
					if (rdist < bucket) break;

					id = table[rpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, id + 1);

						if (++dist_cnt >= candidates) break;
					}
					--rpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				rpos_[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < bucket && rdist < bucket)) {
					bucket_flag_[j] = false;
					++num_bucket;
				}
				if (num_bucket >= m_) break;
			}
			if (num_bucket >= m_ || dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist > radius / appr_ratio_ && dist_cnt >= top_k) break;
		if (dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search
	int top_k,							// top-k value
	const float *query,					// query object
	const int *object_id,				// object id mapping
	MaxK_List *list)					// k-FN results (return)
{
	int    candidates = CANDIDATES + top_k - 1; // candidate size
	float  kdist      = MINREAL;
	Result *table     = NULL;

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		lpos_[i]  = 0;  
		rpos_[i]  = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   dist_cnt = 0;			// counter of Eucldiean distance compuation
	float radius   = find_radius(); // search radius
	float bucket   = radius * w_ / 2.0f; // bucket width

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		int cnt = -1, lpos = -1, rpos = -1, id = -1;
		while (num_bucket < m_) {
			float ldist = -1.0f;	// left  proj dist with query
			float rdist = -1.0f;	// right proj dist with query
			float q_val = -1.0f;	// hash value of 
			float dist  = -1.0f;	// l2 dist

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				q_val = q_val_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				lpos = lpos_[j]; rpos = rpos_[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos < rpos) {
						ldist = fabs(q_val - table[lpos].key_);
					}
					else break;
					if (ldist < bucket) break;

					id = table[lpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, object_id[id] + 1);

						if (++dist_cnt >= candidates) break;
					}
					++lpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				lpos_[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos < rpos) {
						rdist = fabs(q_val - table[rpos].key_);
					}
					else break;
					if (rdist < bucket) break;

					id = table[rpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						dist = calc_l2_dist(dim_, data_[id], query);
						kdist = list->insert(dist, object_id[id] + 1);

						if (++dist_cnt >= candidates) break;
					}
					--rpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				rpos_[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < bucket && rdist < bucket)) {
					bucket_flag_[j] = false;
					++num_bucket;
				}
				if (num_bucket >= m_) break;
			}
			if (num_bucket >= m_ || dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (kdist > radius / appr_ratio_ && dist_cnt >= top_k) break;
		if (dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}
	return 0;
}

// -----------------------------------------------------------------------------
int RQALSH::kfn(					// k-FN search (for ML_RQALSH)
	int   top_k,						// top-k value
	float R,							// limited search range
	const float *query,					// input query
	const int *object_id,				// objects id mapping
	MaxK_List *list)					// k-FN results (return)
{
	int    candidates = CANDIDATES + top_k - 1; // candidate size
	float  kdist      = MINREAL;
	Result *table     = NULL;

	// -------------------------------------------------------------------------
	//  init parameters
	// -------------------------------------------------------------------------
	memset(freq_, 0, n_pts_ * SIZEFLOAT);
	memset(checked_, false, n_pts_ * SIZEBOOL);
	memset(range_flag_, true, m_ * SIZEBOOL);

	for (int i = 0; i < m_; ++i) {
		q_val_[i] = calc_hash_value(i, query);
		lpos_[i]  = 0;  
		rpos_[i]  = n_pts_ - 1;
	}

	// -------------------------------------------------------------------------
	//  c-k-AFN search
	// -------------------------------------------------------------------------
	int   dist_cnt  = 0;			// counter of Eucldiean distance compuation
	int   num_range = 0; 			// used for search range

	float radius    = find_radius(); // search radius
	float bucket    = radius * w_ / 2.0f; // bucket width
	float range     = -1.0f;			// search range width

	if (R < FLOATZERO) range = 0.0f;
	else range = R * w_ / 2.0f;

	while (true) {
		// ---------------------------------------------------------------------
		//  step 1: initialization
		// ---------------------------------------------------------------------
		int num_bucket = 0;
		memset(bucket_flag_, true, m_ * SIZEBOOL);

		// ---------------------------------------------------------------------
		//  step 2: (R,c)-FN search
		// ---------------------------------------------------------------------
		int cnt = -1, lpos = -1, rpos = -1, id = -1;
		while (num_bucket < m_ && num_range < m_) {
			float ldist = -1.0f;	// left  proj dist with query
			float rdist = -1.0f;	// right proj dist with query
			float q_val = -1.0f;	// hash value of 
			float dist  = -1.0f;	// l2 dist

			for (int j = 0; j < m_; ++j) {
				if (!bucket_flag_[j]) continue;

				table = tables_[j];
				q_val = q_val_[j];
				// -------------------------------------------------------------
				//  step 2.1: scan left part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				lpos = lpos_[j]; rpos = rpos_[j];
				while (cnt < SCAN_SIZE) {
					ldist = MAXREAL;
					if (lpos < rpos) {
						ldist = fabs(q_val - table[lpos].key_);
					}
					else break;
					if (ldist < bucket) break;

					id = table[lpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						// assert(id >= 0 && id < n_pts_);
						dist = calc_l2_dist(dim_, data_[id], query);
						list->insert(dist, object_id[id] + 1);

						if (++dist_cnt >= candidates) break;
					}
					++lpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				lpos_[j] = lpos;

				// -------------------------------------------------------------
				//  step 2.2: scan right part of hash table
				// -------------------------------------------------------------
				cnt = 0;
				while (cnt < SCAN_SIZE) {
					rdist = MAXREAL;
					if (lpos < rpos) {
						rdist = fabs(q_val - table[rpos].key_);
					}
					else break;
					if (rdist < bucket) break;

					id = table[rpos].id_;
					if (++freq_[id] >= l_ && !checked_[id]) {
						checked_[id] = true;
						// assert(id >= 0 && id < n_pts_);
						dist = calc_l2_dist(dim_, data_[id], query);
						list->insert(dist, object_id[id] + 1);

						if (++dist_cnt >= candidates) break;
					}
					--rpos; ++cnt;
				}
				if (dist_cnt >= candidates) break;
				rpos_[j] = rpos;

				// -------------------------------------------------------------
				//  step 2.3: check whether this bucket is finished scanned
				// -------------------------------------------------------------
				if (lpos >= rpos || (ldist < bucket && rdist < bucket)) {
					bucket_flag_[j] = false;
					++num_bucket;
					
					if ((lpos >= rpos || (ldist < range && rdist < range)) && range_flag_[j]) {
						range_flag_[j] = false;
						++num_range;
					}
				}
				if (num_bucket >= m_ || num_range >= m_) break;
			}
			if (num_bucket >= m_ || num_range >= m_) break;
			if (dist_cnt >= candidates) break;
		}

		// ---------------------------------------------------------------------
		//  step 3: stop conditions 1 & 2
		// ---------------------------------------------------------------------
		if (num_range >= m_ || dist_cnt >= candidates) break;

		// ---------------------------------------------------------------------
		//  step 4: update radius
		// ---------------------------------------------------------------------
		radius = radius / appr_ratio_;
		bucket = radius * w_ / 2.0f;
	}

	return 0;
}

// -----------------------------------------------------------------------------
float RQALSH::find_radius()			// find proper radius
{
	// -------------------------------------------------------------------------
	//  find an array of projected distance which is closest to the query in
	//  each of <m> hash tables 
	// -------------------------------------------------------------------------
	vector<float> list;
	for (int i = 0; i < m_; ++i) {
		if (lpos_[i] < rpos_[i]) {
			list.push_back(fabs(tables_[i][lpos_[i]].key_ - q_val_[i]));
			list.push_back(fabs(tables_[i][rpos_[i]].key_ - q_val_[i]));
		}
	}

	// -------------------------------------------------------------------------
	//  sort the array in ascending order 
	// -------------------------------------------------------------------------
	sort(list.begin(), list.end());

	// -------------------------------------------------------------------------
	//  find the median distance and return the new radius
	// -------------------------------------------------------------------------
	int num  = (int) list.size();
	float dist = -1.0f;
	if (num % 2 == 0) dist = (list[num / 2 - 1] + list[num / 2]) / 2.0f;
	else dist = list[num / 2];

	int kappa = (int) ceil(log(2.0f * dist / w_) / log(appr_ratio_));
	return pow(appr_ratio_, kappa);
}
